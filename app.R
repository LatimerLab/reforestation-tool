library(shiny)
library(ggplot2)
library(raster)
library(sf)
library(dplyr)
library(rgdal)
library(lme4)
library(fasterize)
library(viridis)
library(shinyjs)
library(V8)
library(tidyr)
library(scales)

options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize = 100*1024^2)

#### Globals ####

var_lims = read.csv("data/var_lims.csv", header=TRUE)
color_pal <- c("below acceptable density range (even with planting)" = "#fde725", "within acceptable density range if planted" = "#7ad151", "within acceptable density range regardless of planting" = "#22a884",  "below acceptable density range when unplanted; above it when planted" = "#2a788e", "within acceptable density range if unplanted" = "#414487", "above acceptable density range (even without planting)" = "#440154")
mod = readRDS("data/model.rds")
dat = readRDS("data/data.rds")

maps_loaded = reactiveVal(FALSE) # Holds status of whether maps for display are fully computed (so can disable upload box)
severity_uploaded = reactiveVal(FALSE)

df_plot_export = reactiveVal(NULL) # for exporting raster maps in file download handler

experimental_enabled = reactiveVal(FALSE) # stores whether the experimental tools are enabled

jsResetCode <- "shinyjs.reset = function() {history.go(0)}" # Define the js method that resets the page

albers = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "


### A blank ggplot that says "Calculating..."

d = data.frame(x = c(0,0,100,100),
               y = c(0,100,0,100),
               z = c(1,1,1,1))
blank_plot = ggplot(d,aes(x=x,y=y)) +
  theme_void() +
  geom_blank() +
  theme(plot.background = element_rect(fill="grey90",color=NA)) +
  annotate("text",x=0,y=97,label="Calculating...",hjust=0)



#### Fire perim and sev ####

# perim = st_read("data/power_perim.geojson") %>% st_transform(3310)
# sevtmp= st_read("data/rim_sev.gpkg") %>% st_transform(3310)

#### Get predictors and predictions for a given fire ####

prep_mapping_vars = function(perim,sev,input,min_high_sev_manual) {
  
  cat("starting to prep mapping vars")
  
  if(is.null(sev)) return(NULL)
  
  if(input$resolution_input == 1) {
    env = brick("data/env_raster_stack.tif")
  } else if (input$resolution_input == 2) {
    env = brick("data/env_raster_stack_coarse.tif")
  } else {
    env = brick("data/env_raster_stack_extracoarse.tif")
  }
  
  env = crop(env,perim %>% st_transform(projection(env)))
  env = mask(env,perim %>% st_transform(projection(env)))

  
  # get a fire perimeter shapefile and severity raster
  # temp sev = raster("/home/derek/Desktop/ca4005312066920190904_20181007_20191010_rdnbr_ba7.tif")
  #temp perim = st_read("/home/derek/Downloads/ca4005312066920190904_20181007_20191010_ravg_data/bdy.gpkg") %>% st_transform(3310)
  
  # # Get a (smoothed) perimeter from the severity shapefile
  # perim = st_union(sev) %>% st_as_sf() %>% st_simplify(dTolerance = 30) %>% st_buffer(0) %>% smooth(method="ksmooth") %>% st_buffer(0)
  # 
  
  ### Compute seed distance ###
  
  ## get all the non-high-sev area
  sev_nonhigh = sev
  sev_nonhigh = projectRaster(sev_nonhigh,env,method="ngb")
  
  if(min_high_sev_manual == -1) {
    sev_nonhigh[sev_nonhigh >= as.numeric(input$min_high_sev)] = NA
  } else {
    sev_nonhigh[sev_nonhigh >= min_high_sev_manual] = NA
  }
  
  ## comp distance to non-high-sev
  seed_dist = distance(sev_nonhigh)
  

  ## crop to the fire footprint
  seed_dist = crop(seed_dist,perim %>% st_transform(projection(seed_dist)))
  seed_dist = mask(seed_dist,perim %>% st_transform(projection(seed_dist)))
  
  ## Make a seed dist (0 m) mask layer
  non_high_sev_mask = seed_dist
  non_high_sev_mask[seed_dist == 0] = 1
  non_high_sev_mask[seed_dist != 0] = 0
  
  
  
  #### Load and assemble env predictor data for focal area ####
  env = stack(env,seed_dist,non_high_sev_mask)
  
  env_df = as.data.frame(env,xy=TRUE)
  names(env_df) = c("x","y","tpi","ppt","tmean","shrub","elev","eveg","seed_dist","non_high_sev_mask")
  
  env_df = env_df[!is.na(env_df$ppt),]
  
  ## manually set seed dist based on user selection
  if(input$seed_dist_assumption == "30") env_df$seed_dist = 30
  if(input$seed_dist_assumption == "70") env_df$seed_dist = 70
  if(input$seed_dist_assumption == "250") env_df$seed_dist = 250
  # if the final option, "severity", is selected, use the default of estimating from severity raster (distance to non-high severity)
  
  
  ## override shrub cover based on input selection
  if(input$shrub_assumption == "low") env_df$shrub = 3000
  if(input$shrub_assumption == "moderate") env_df$shrub = 6000
  if(input$shrub_assumption == "high") env_df$shrub = 9000
  # otherwise, if it is "predicted", use the default value
  
  env_df = env_df %>%
    rename(normal_annual_precip = ppt,
           tpi2000 = "tpi",
           Shrubs = "shrub") %>%
    ## undo the raster scaling that was done to be able to save layers as int
    mutate(tmean = tmean/100,
           tpi2000 = tpi2000/10,
           Shrubs = Shrubs/100,
           eveg = eveg/100) %>%
    # apply transformation needed by stat model
    ###!!! consider different cap on seed wall max distance
    mutate(seed_dist = ifelse(seed_dist == 0,15,seed_dist)) %>% # correct for a limitation of using remotely sensed data: no plot center is exactly 0 m from a tree. Use 15 since we are focused on high-severity areas, so the closest a tree could be is half the width of the pixel. Also conveniently 15 m was the closest a tree was in our plot dataset.
    mutate(seed_dist = ifelse(seed_dist >= 250,250,seed_dist)) %>% # cap it at 200 since our field data only go that far and we know it tends to level off by then
    mutate(log10SeedWallConifer = log10(seed_dist))
  
  ## make an extrapolation column
  env_df = env_df %>%
    dplyr::mutate(extrap = !(tpi2000 %>% between(var_lims$tpi_min,var_lims$tpi_max)) |
             !(normal_annual_precip %>% between(var_lims$ppt_min,var_lims$ppt_max)) |
             !(tmean %>% between(var_lims$tmean_min,var_lims$tmean_max)) |
               !(elev %>% between(var_lims$elev_min,var_lims$elev_max))   )
  
  
  
  benefit_range = get_benefit_range(env_df)
  
  
  
  mapping_vars = list(env_df = env_df,
                      perim = perim,
                      sev = sev,
                      benefit_range = benefit_range)
  
  cat("mapping vars prepped")
  
  return(mapping_vars)
  
}






# 
# summary(env_df$tmean)
# summary(env_df$tpi2000)
# summary(env_df$normal_annual_precip)
# 
# 
# 

# 
# dat = dat %>%
#   select(tpi2000,facts.planting.first.year,Shrubs,fsplanted,tmin,normal_annual_precip,neglog5SeedWallConifer,ShrubHt) %>%
#   mutate(fsplanted = ifelse(fsplanted == "planted",1,0))
# 
# summary(dat$tmean)
# summary(dat$tpi2000)
# summary(dat$normal_annual_precip)


#### Function to load fire files ####
load_sev = function(input,output) {
  
  # perim = st_read(input$perim_file) %>% st_transform(3310)
  # sev = st_read(input$sev_file) %>% st_transform(3310)
  
  if(is.null(input$sev_file)) {
    sev = NULL
  } else  {
    sev = raster(input$sev_file$datapath)
    severity_uploaded(TRUE)
    cat("Sev raster loaded")
    

    
  }
  
  return(sev)
  
}


#### Function to load fire files ####
load_perim = function(input) {
  
  # perim = st_read(input$perim_file) %>% st_transform(3310)
  # sev = st_read(input$sev_file) %>% st_transform(3310)
  
  if(is.null(input$perim_file)) {
    perim = NULL
  } else  {
    perim = st_read(input$perim_file$datapath) %>% st_transform
    cat("Perim loaded")
  }
  
  return(perim)
  
}


#### Function to load fire files ####
test_perim_uploaded = function(input) {
  
  # perim = st_read(input$perim_file) %>% st_transform(3310)
  # sev = st_read(input$sev_file) %>% st_transform(3310)
  
  return(!is.null(input$perim_file))
  

  
  # str(input$perim_file)
  
}


custom_debug = function(input) {
  
  # perim = st_read(input$perim_file) %>% st_transform(3310)
  # sev = st_read(input$sev_file) %>% st_transform(3310)
  
  str(input$sev_file)
  str(file_uploaded)
  
  reset('sev_file')
  
}


#### Function for computing range of potential planting benefit values on a fire (for scaling relative planting benefit display) ####
get_benefit_range = function(env_df) {

  env_df = env_df %>%
    filter(!is.na(tmean))
  
  ## replicate for each possible value of planting year
  env_df = tidyr::expand_grid(env_df,facts.planting.first.year = c(1,2,3))
  
  # calc only for high sev, YPMC, non-extrapolation
  env_df = env_df[which(env_df$non_high_sev_mask != 1),]
  env_df = env_df[which(env_df$eveg > 0.5),]
  env_df = env_df[which(env_df$extrap != TRUE),]
  
  
  env_df_noplant = env_df %>%
    mutate(fsplanted = "unplanted")
  
  ## for predictions of a planting scenario
  env_df_plant = env_df %>%
    mutate(fsplanted = "planted")
  
  ## put tables together for prediction
  env_df = bind_rows(env_df_noplant,
                     env_df_plant)
  
  # predict seedlings per ha (incl undoing the response variable trasformation)
  pred = predict(mod,env_df,re.form=NA) %>% exp() - 24.99
  
  pred_df = env_df
  
  pred_df$pred = pred / 2.47 # convert seedlings/ha to seedlings/acre
  
  # any model predictions below 0 should be 0
  pred_df[pred_df$pred < 0.1,"pred"] = 0.1
  
  ## put the two data frames side by side
  pred_df_plant = pred_df %>%
    filter(fsplanted == "planted") %>%
    rename("pred_plant" = "pred")
  
  pred_df_noplant = pred_df %>%
    filter(fsplanted == "unplanted") %>%
    rename("pred_noplant" = "pred")
  
  pred_df = bind_cols(pred_df_plant,pred_df_noplant %>% dplyr::select(pred_noplant))
  
  df_plot = pred_df %>%
    ### Compute planting benefit
    mutate(planting_benefit = pred_plant - pred_noplant)

  ## get min and max values
  min_benefit = df_plot$planting_benefit %>% quantile(0.025)
  max_benefit = df_plot$planting_benefit %>% quantile(0.975)
  
  benefit_range = c(min_benefit,max_benefit)

  return(benefit_range)  
  
}




#### Function for making maps based on inputs ####
# This function relies on some globals from above
make_maps = function(mapping_vars, plant_year, density_low, density_high, map_masking) {

  if(is.null(mapping_vars)) {
    
    ret = list("main" = blank_plot, "planting_benefit" = blank_plot, "density_unplanted" = blank_plot, "density_planted" = blank_plot,
               "cover_shrub" = blank_plot,
               "seed_distance" = blank_plot,
               tmean = blank_plot,
               precip = blank_plot,
               tpi = blank_plot)
    
    
  } else {
    
    perim = mapping_vars$perim
    
    env_df = mapping_vars$env_df %>%
      mutate(#neglog5SeedWallConifer = input$seedwall,
        facts.planting.first.year = as.numeric(plant_year)) %>%
      filter(!is.na(tmean)) %>%  ### exclude NA grid cells!
      filter(!is.na(tpi2000))
    # fsplanted = input$planted)
    #  Shrubs = input$shrub_cover,  
    #ShrubHt = input$shrub_height,
    #LitDuff = input$lit_duff)
    
    # env_df = env_df %>%
    #   mutate(normal_annual_precip = 997.3,
    #          tpi2000 = 14.6,
    #          tmean = 2.06)
    
    # env_df = env_df %>%
    #   mutate(Shrubs = asin(sqrt(Shrubs/100)))
    
    ## for predictions of a no-planting scenario
    env_df_noplant = env_df %>%
      mutate(fsplanted = "unplanted")
    
    ## for predictions of a planting scenario
    env_df_plant = env_df %>%
      mutate(fsplanted = "planted")
    
    ## put tables together for prediction
    env_df = bind_rows(env_df_noplant,
                       env_df_plant)
    
    # predict seedlings per ha (incl undoing the response variable trasformation)
    pred = predict(mod,env_df,re.form=NA) %>% exp() - 24.99

    #browser()
    
    ####!!!! get uncertainty following http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html # lme4 section

    if("mask_uncertain" %in% map_masking) {
      
      newdat = env_df %>%
      mutate(ln.dens.planted = 0) %>%
      mutate(fsplanted = as.factor(fsplanted))
      
      mm = model.matrix(terms(mod),newdat)
      tcross = tcrossprod(vcov(mod),mm)
      
      nobs = nrow(newdat)
      block_size = 2000
      nblocks = ceiling(nobs/block_size)
      
      stderr = NULL
      
      for(i in 1:nblocks) {
        start = 1 + (i-1)*block_size
        end = i*block_size
        if(end > nobs) end = nobs
        
        mm_sub = mm[start:end,]
        tcross_sub = tcross[,start:end]
        stderr_sub <- diag(mm_sub %*% tcross_sub) %>% sqrt()
        
        stderr = c(stderr,stderr_sub)
        
      }
    } else {
      stderr = rep(0,nrow(env_df))
    }

    
    pred_df = env_df
    
    pred_df$pred = pred / 2.47 # convert seedlings/ha to seedlings/acre
    pred_df$stderr = stderr
    
    
    # any model predictions below 0 should be 0
    pred_df[which(pred_df$pred < 0.1),"pred"] = 0.1
    
    
    ## put the two data frames side by side
    pred_df_plant = pred_df %>%
      filter(fsplanted == "planted") %>%
      rename("pred_plant" = "pred",
             "stderr_plant" = "stderr")
    
    pred_df_noplant = pred_df %>%
      filter(fsplanted == "unplanted") %>%
      rename("pred_noplant" = "pred",
             "stderr_noplant" = "stderr")
    
    pred_df = bind_cols(pred_df_plant,pred_df_noplant %>% dplyr::select(pred_noplant,stderr_noplant))
    

    ## classify densities
    pred_df = pred_df %>%
      mutate(noplant_class = cut(pred_noplant,breaks=c(-Inf,density_low,density_high,Inf),labels=c("low","good","high"))) %>%
      mutate(plant_class = cut(pred_plant,breaks=c(-Inf,density_low,density_high,Inf),labels=c("low","good","high"))) %>%
      filter(!is.na(pred_noplant))
    
    
    
    pred_df$overall = NA
    ##too low even with planting
    pred_df[which(pred_df$plant == "low"),"overall"] = "below acceptable density range (even with planting)"
    ##too high even without planting
    pred_df[which(pred_df$noplant == "high"),"overall"] = "above acceptable density range (even without planting)"
    ## good with planting only
    pred_df[which(pred_df$noplant != "good" & pred_df$plant == "good"),"overall"] = "within acceptable density range if planted"
    ## good with natural only
    pred_df[which(pred_df$noplant == "good" & pred_df$plant != "good"), "overall"] = "within acceptable density range if unplanted"
    ## good regardless of planting
    pred_df[which(pred_df$noplant == "good" & pred_df$plant == "good"), "overall"] = "within acceptable density range regardless of planting"
    ## too low when unplanted; too high when planted
    pred_df[which(pred_df$noplant == "low" & pred_df$plant == "high"), "overall"] = "below acceptable density range when unplanted; above it when planted"
    
    # make a column with these values as numeric codes
    pred_df = pred_df %>%
      mutate(overall_code = recode(overall,
                                   "below acceptable density range (even with planting)" = 1,
                                   "within acceptable density range if planted" = 2,
                                   "within acceptable density range regardless of planting" = 3,
                                   "below acceptable density range when unplanted; abive it when planted" = 4,
                                   "within acceptable density range if unplanted" = 5,
                                   "above acceptable density range (even without planting)" = 6))
    
    df_plot = pred_df %>%
      mutate(overall = factor(overall, levels = rev(names(color_pal)))) %>%
      mutate(pred_noplant_trunc = ifelse(pred_noplant > 400, 400, pred_noplant)) %>%
      mutate(pred_plant_trunc = ifelse(pred_plant > 400, 400, pred_plant)) %>%
    ### Compute planting benefit
      mutate(planting_benefit = pred_plant - pred_noplant)
    
  
    ## Drop rows of masked-out values
    if("mask_non_high_sev" %in% map_masking) {
      df_plot = df_plot[which(df_plot$non_high_sev_mask != 1),]
    }
    
    if("mask_uncertain" %in% map_masking) {
      df_plot = df_plot[which(df_plot$stderr_noplant < 0.60),]
    }
    
    if("mask_extrapolation" %in% map_masking) {
      df_plot = df_plot[which(df_plot$extrap != TRUE),]
    }
    
    if("mask_non_ypmc" %in% map_masking) {
      df_plot = df_plot[which(df_plot$eveg > 0.5),]
    }
    
    ### get the range of possible benefit for this fire and scale the planting benefit column
    benefit_range = mapping_vars$benefit_range
    
    ## set the max to a reasonable cap (400)
    cap = 300
    if(benefit_range[2] > cap) benefit_range[2] = cap
    if(benefit_range[1] > cap) benefit_range[1] = cap-10
    
    
    df_plot = df_plot %>%
      mutate(planting_benefit = ifelse(planting_benefit > cap,cap,planting_benefit)) %>%
      mutate(planting_benefit_rel = rescale(planting_benefit,from=benefit_range)) %>%
      mutate(planting_benefit_rel = ifelse(planting_benefit_rel > 1,1,planting_benefit_rel),
             planting_benefit_rel = ifelse(planting_benefit_rel < 0,0,planting_benefit_rel))


    df_plot_export(df_plot)
    cat("Saved plot data table")
    
    
    perim = perim %>% st_transform(3310)
    
    main_map = ggplot(data=df_plot,aes(x=x,y=y,fill=overall)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_manual(values = color_pal) +
      theme(legend.position="bottom", legend.title=element_blank(), plot.title = element_text(hjust = 0.5)) +
      labs(title = "Predicted seedling density outcomes") +
      guides(fill = guide_legend(nrow = 6))
    
    
    
    planting_benefit = ggplot(df_plot,aes(x=x,y=y,fill=planting_benefit_rel)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(direction = -1,
                         option = "inferno",
                         begin = 0.2,
                         end = 0.9,
                         limits=c(0,1),
                         labels = c("low","high"),
                         breaks = c(0.1,0.9),
                         guide = guide_colourbar(ticks=FALSE)) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Predicted planting benefit", fill = NULL)
    
    
    
    density_unplanted = ggplot(df_plot,aes(x=x,y=y,fill=pred_noplant_trunc)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(breaks=c(0,200,400),
                           labels=c(0,200,"400+"),
                           limits=c(0,400),
                           direction = -1,
                         option = "inferno",
                         begin = 0.2,
                         end = 0.9) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Predicted seedling density: Natural", fill = "Seedlings / acre") #+
      #guides(fill = guide_legend(reverse=FALSE))
      #guides(fill = guide_legend(nrow = 2))
    
    density_planted = ggplot(df_plot,aes(x=x,y=y,fill=pred_plant_trunc)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(breaks=c(0,200,400),
                         labels=c(0,200,"400+"),
                         limits=c(0,400),
                         direction = -1,
                         option = "inferno",
                         begin = 0.2,
                         end = 0.9) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Predicted seedling density: Natural + planted", fill = "Seedlings / acre") #+
    
    cover_shrub = ggplot(df_plot,aes(x=x,y=y,fill=Shrubs)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(direction = -1) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Shrub cover", fill = "Percent cover")
    
    seed_distance = ggplot(df_plot,aes(x=x,y=y,fill=seed_dist)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis() +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Seed tree distance", fill = "Distance (m)")
    
    tmean = ggplot(df_plot,aes(x=x,y=y,fill=tmean)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(option = "inferno") +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Mean annual temperature", fill = "°C")
    
    precip = ggplot(df_plot,aes(x=x,y=y,fill=normal_annual_precip)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(option = "viridis",
                         direction = -1) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Annual average precipitation", fill = "mm")
    
    tpi = ggplot(df_plot,aes(x=x,y=y,fill=tpi2000)) +
      geom_raster() +
      geom_sf(data=perim,color="black",fill = NA, inherit.aes = FALSE) +
      theme_void(20) +
      scale_fill_viridis(option = "inferno",
                         direction = 1) +
      theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.direction = "vertical") +
      labs(title = "Topographic position index", fill = "TPI")
    
    ret = list("main" = main_map, "planting_benefit" = planting_benefit, "density_unplanted" = density_unplanted, "density_planted" = density_planted,
               "cover_shrub" = cover_shrub,
               "seed_distance" = seed_distance,
               tmean = tmean,
               precip = precip,
               tpi = tpi)
    
    maps_loaded(TRUE)
    
  }
  
  return(ret)
  
}




# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  useShinyjs(),                                           # Include shinyjs in the UI
  extendShinyjs(text = jsResetCode, functions = c("reset")),                      # Add the js code to the page
  
  # App title ----
  titlePanel("Post-fire reforestation success estimation tool"),
  h4('aka "PRESET"'),
  h6('v 0.1 (beta)'),
  h6("Developed by: ",a(href="http://www.changingforests.com","Derek Young",.noWS="after"),", ",a(href="https://twitter.com/qmsorenson?lang=en","Quinn Sorenson",.noWS="after"),", ",a(href="https://www.plantsciences.ucdavis.edu/people/andrew-latimer","Andrew Latimer",.noWS="after")),
  h6(a(href="https://latimer.ucdavis.edu/","Latimer Lab",.noWS="after")," and ",a(href="http://www.changingforests.com","Young Lab",.noWS="after"),", Department of Plant Sciences, UC Davis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      ### Upload fire perim and severity
      #fileInput("perim_file", "Fire perimeter"),
      
      conditionalPanel(
        condition="output.map_loaded ===false",
        
        radioButtons(inputId = "resolution_input",
                     label = "1. Choose spatial resolution:",
                     choices = list("30 m pixels (slowest)" = 1, "60 m pixels" = 2, "120 m pixels (fastest)" = 3),
                     selected = 3),
        
        tags$br(),
        
        fileInput("sev_file", "2. Upload severity raster"),
        
        conditionalPanel( condition = "output.severity_uploaded === true",
              textInput("min_high_sev",label="2b. Minimum value to consider high severity"),
              htmlOutput("severity_text")
 
                          ),
        
        fileInput("perim_file", "3. Upload fire perimeter (or focal area perimeter)"),
        
        HTML('&emsp;'),tags$b("Or,"),"use 2020 Creek Fire demo data: ",
        
        actionButton("upload_button","Go")
      ),
      
      conditionalPanel(
        condition="output.map_loaded === true",
        "Fire perimeter/severity map loaded.",
        actionButton("reset_button","Clear")
      ),

      
      # checkboxInput(inputId = "show_map_controls",
      #               label = "Show map controls",
      #               value = FALSE),
      
      
      
      ### Map display controls
      conditionalPanel(
        condition="output.map_loaded === true",
      
        
        # Input: Slider for the number of bins ----
        # sliderInput(inputId = "seedwall",
        #             label = "Seed wall:",
        #             min = -1.2,
        #             max = -0.6,
        #             value = -0.9),

        # sliderInput(inputId = "max_density",
        #             label = "Maximum seedling density:",
        #             min = 0,
        #             max = 300,
        #             value = 50),
        radioButtons(inputId = "planted_year",
                    label = "Planting year:",
                    choices = list("1 year post-fire" = 1, "2 years post-fire" = 2, "3 years post-fire" = 3),
                    selected = 2),
        radioButtons(inputId = "shrub_assumption",
                     label = "Expected shrub cover:",
                     choices = list("Predicted (varies across space)" = "predicted", "Low (30%) everywhere" = "low", "Moderate (60%) everywhere" = "moderate", "High (90%) everywhere" = "high"),
                     selected = "predicted"),
        checkboxGroupInput(inputId = "map_masking",
                           label = "Map filtering:",
                           choices = list("Show high-severity area only" = "mask_non_high_sev",
                                          "Show yellow pine & mixed-conifer only" = "mask_non_ypmc",
                                          #"Hide model extrapolation areas" = "mask_extrap",
                                          "Hide low model confidence areas (slow)" = "mask_uncertain"
                           ),
                           selected = c("mask_non_high_sev","mask_non_ypmc")),
                           

        conditionalPanel(condition = "output.experimental_enabled === true",
                         sliderInput(inputId = "density_range",
                                     label = "Acceptable seedling density range (seedlings/acre):",
                                     min = 0,
                                     max = 600,
                                     value = c(50,250)),
                         
                         radioButtons(inputId = "seed_dist_assumption",
                                      label = "Seed source distance:",
                                      choices = list("Near (100 ft)" = "30",
                                                     "Moderate (230 ft)" = "70",
                                                     "Far (820 ft)" = "820",
                                                     "Derived from severity map (dist. to non-high sev.)" = "severity"),
                                      selected = 70),
        ),
        conditionalPanel(condition = "output.experimental_enabled === false",
                          checkboxGroupInput(inputId = "map_selection_basic",
                                             label = "Layers to display:",
                                             choices = list("Planting benefit" = "planting_benefit",
                                                            "Shrub cover" = "cover_shrub",
                                                            "Mean temperature" = "tmean",
                                                            "Annual precipitation" = "precip",
                                                            "Topographic position index" = "tpi"
                                             ),
                                             selected = "planting_benefit")
        ),
        
        conditionalPanel(condition = "output.experimental_enabled === true",
                         checkboxGroupInput(inputId = "map_selection_advanced",
                                            label = "Layers to display:",
                                            choices = list("Predicted outcomes" = "main",
                                                           "Planting benefit" = "planting_benefit",
                                                           "Natural seedling density" = "density_unplanted",
                                                           "Planted + natural seedling density" = "density_planted",
                                                           "Modeled shrub cover" = "cover_shrub",
                                                           "Seed tree distance" = "seed_distance",
                                                           "Mean temperature" = "tmean",
                                                           "Annual precipitation" = "precip",
                                                           "Topographic position index" = "tpi"
                                            ),
                                            selected = "main")
                         
                         ),
     
        tags$br(),
        
        conditionalPanel(condition = "output.experimental_enabled === true",
          selectInput("dataset_advanced", "Download data:",
                      choices = c("Predicted outcomes", "Planting benefit", "Natural seedling density", "Planted + natural seedling density"))
        ),
        conditionalPanel(condition = "output.experimental_enabled === false",
                         selectInput("dataset_basic", "Download data:",
                                     choices = c("Planting benefit"))

        ),
        downloadButton("downloadData", "Download"),
        
        

        tags$br(),
        tags$br(),
        
        conditionalPanel(condition = "output.experimental_enabled === false",
          actionButton("experimental_enabled","Enable experimental/beta tools")
        ),
        conditionalPanel(condition = "output.experimental_enabled === true",
                         tags$b("Experimental/beta tools enabled.")
        ),
        
        "Important caveats"
        
        
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(

      # Conditional panel to hide all the maps if user has not yet uploaded a perimeter shapefile
      conditionalPanel(
          condition = "output.perim_uploaded === true",
          
        
          
          ## Main map
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('main') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('main') === true)",      
            plotOutput(outputId = "distPlot",height="600px"),
            tags$br()
          ),
          
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('planting_benefit') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('planting_benefit') === true)",      
            plotOutput(outputId = "plantingBenefitPlot",height="600px"),
            tags$br()
          ),
          
          ## Unplanted density
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('density_unplanted') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('density_unplanted') === true)",      
            plotOutput(outputId = "densityUnplantedPlot",height="600px"),
            tags$br(" ")
          ),
          
          ## Planted density
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('density_planted') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('density_planted') === true)",      
            plotOutput(outputId = "densityPlantedPlot",height="600px"),
            tags$br()
          ),
          
          ## Shrub cover
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('cover_shrub') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('cover_shrub') === true)",      
            plotOutput(outputId = "coverShrubPlot",height="600px"),
            tags$br()
          ),
          
          ## Seed distance
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('seed_distance') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('seed_distance') === true)",      
            plotOutput(outputId = "seedDistancePlot",height="600px"),
            tags$br()
          ),
          
          ## Tmin
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('tmean') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('tmean') === true)",      
            plotOutput(outputId = "tmeanPlot",height="600px"),
            tags$br()
          ),
          
          ## Precip
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('precip') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('precip') === true)",      
            plotOutput(outputId = "precipPlot",height="600px"),
            tags$br()
          ),
          
          ## TPI
          conditionalPanel(
            condition = "(!output.experimental_enabled && input.map_selection_basic.includes('tpi') === true) || (output.experimental_enabled && input.map_selection_advanced.includes('tpi') === true)",      
            plotOutput(outputId = "tpiPlot",height="600px"),
            tags$br()
          )
      )


    )
  )
)





# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  sev = reactiveVal(NULL) # Holds fire seveirty raster once loaded
  perim = reactiveVal(NULL) # Holds fire perimeter sf object once loaded
  min_high_sev_manual = reactiveVal(-1) # Used to manually set the minimum high sev value (for, e.g., when the demo layer is loaded)
  
  output$map_loaded = reactive({ maps_loaded() })
  outputOptions(output, 'map_loaded', suspendWhenHidden=FALSE)
  
  output$severity_uploaded = reactive({ severity_uploaded() })
  outputOptions(output, 'severity_uploaded', suspendWhenHidden=FALSE)
  
  output$experimental_enabled = reactive({ experimental_enabled() })
  outputOptions(output, 'experimental_enabled', suspendWhenHidden=FALSE)
  
  # output$max_sev_value = reactive({ max_sev_value() })
  # outputOptions(output, 'max_sev_value', suspendWhenHidden=FALSE)
  
  max_sev_value = reactive({ 
    
        if(!is.null(sev())) {
          return(max(values(sev())))
        } else {
          return(0)
        }
    })
    
  
  output$severity_text = renderUI({ 
    
    if(!is.null(sev())) {
      
      ## get severity range
      sev_values = values(sev())
      min_sev = min(sev_values, na.rm = TRUE)
      max_sev = max(sev_values, na.rm = TRUE)
      range_text = paste0(min_sev, " to ", max_sev)
    } else {
      range_text = "[no severity uploaded yet]"
    }
    
    HTML('&emsp;', paste0( "Range of severity values in uploaded raster: ",range_text,tags$br(),tags$br()))
    
    
    })
  



  observeEvent(input$sev_file,
               
               { sev(load_sev(input,output))
                 cat("Max sev value is:",max_sev_value())
                 updateTextInput(session, "min_high_sev", value=max_sev_value())
                 }       
               
               )
  
  observeEvent(input$perim_file,
               
               { perim(load_perim(input)) }       
               
  )
  
  observeEvent(input$reset_button,
               
               {
                    maps_loaded(FALSE)
                 severity_uploaded(FALSE)
                    js$reset()
               }
              )
  
  observeEvent(input$experimental_enabled,
               
               {
                 
                 cat("clicket")
                 experimental_enabled(TRUE)
               }
  )
  
  output$perim_uploaded = reactive({ return(!is.null(perim())) })
  outputOptions(output, 'perim_uploaded', suspendWhenHidden=FALSE)

  observeEvent(input$upload_button,
                       #{custom_debug(input)}
               {

                 min_high_sev_manual(300)
                 
                 demofile = NULL
                 demofile$sev_file$datapath = "data/demo-data/creek_dnbr_utm.tif"
                 sev(load_sev(demofile))
                 demofile = NULL
                 demofile$perim_file$datapath = "data/demo-data/creek_perim.kml"
                 perim(load_perim(demofile))
                 
               }
                    )
  
  mapping_vars = reactive({ prep_mapping_vars(perim(),sev(),input,min_high_sev_manual()) })

  maps = reactive({ make_maps(mapping_vars(), input$planted_year, input$density_range[1], input$density_range[2], input$map_masking) })
  

  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    # x    <- faithful$waiting
    # bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 

    ## Make model predictions
    
    maps_list = maps()
    plot(maps_list$main)
    
    
    
    
    # hist(x, breaks = bins, col = "#75AADB", border = "white",
    #      xlab = "Waiting time to next eruption (in mins)",
    #      main = "Histogram of waiting times")
    
  })
  
  output$plantingBenefitPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$planting_benefit)
  })
  
  output$densityUnplantedPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$density_unplanted)
  })
  
  output$densityPlantedPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$density_planted)
  })
  
  output$coverShrubPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$cover_shrub)
  })
  
  output$seedDistancePlot = renderPlot({
    maps_list = maps()
    plot(maps_list$seed_distance)
  })
  
  output$tmeanPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$tmean)
  })
  
  output$precipPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$precip)
  })
  
  output$tpiPlot = renderPlot({
    maps_list = maps()
    plot(maps_list$tpi)
  })
  

  dataset_colname <- reactive({
    switch(ifelse(experimental_enabled(),input$dataset_advanced,input$dataset_basic),
           "Predicted outcomes" = "overall_code",
           "Planting benefit" = "planting_benefit_rel",
           "Natural seedling density" = "pred_noplant",
           "Planted + natural seedling density" = "pred_plant")
  })
  
  dataset_filename <- reactive({
    switch(ifelse(experimental_enabled(),input$dataset_advanced,input$dataset_basic),
           "Predicted outcomes" = "predicted_outcomes",
           "Planting benefit" = "planting_benefit",
           "Natural seedling density" = "density_natural",
           "Planted + natural seedling density" = "density_planted_plus_natural")
  })

  
  
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste0(dataset_filename(), ".tif")
    },
    content = function(file) {
      
      xyz = df_plot_export() %>%
        dplyr::select(x,y,!!!dataset_colname())

      rast = rasterFromXYZ(xyz, crs = albers)
      
      writeRaster(rast, file)
    }
  )
  

  
  
}


shinyApp(ui = ui, server = server)
