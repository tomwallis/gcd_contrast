## Import data into a nice R format, save in "output" directory:
# TSAW.

# full_set gives you all the predictors we've calculated, for all observers, including missing data.

import_data <- function(){
  
  # if you want the full set, also include NAs:
  include_NAs = TRUE
  
  # set up data directory:
  data_directory <- paste(getwd(),'/data/',sep="")
  
  # file to output:
  output_file <- paste0(getwd(),'/output/csf_data_full.RData')
  
  # import basic features--------------------------------------------
  fileName <- paste(data_directory,'csf_basic_data',sep="")
  
  dat <- read.table(fileName,col.names=c("unique_id","subject","target_spatial_band","alpha","stim_pos","resp"))
  
  # change some variable levels into labelled factors:
  subjects <- c("MACD","PJB","TSAW","LAL","AJ","LB","CPT","AM")
  
  dat$subject <- as.factor(dat$subject)
  levels(dat$subject) <- subjects
  
  # Create a column for "correct":
  dat$correct <- as.numeric(dat$stim_pos==dat$resp)
  
  # stim_pos and resp as factors:
  dat$stim_pos_factor <- factor(dat$stim_pos)
  levels(dat$stim_pos_factor) <- c('south','west','east','north')
  dat$resp_factor <- factor(dat$resp)
  levels(dat$resp_factor) <- c('south','west','east','north')
  
  initialValidTrials <- nrow(dat)
  
  # band energy (all bands)--------------------------------------------
  fileName <- paste(data_directory,'csf_bandenergy_data_allbands',sep="")
  
  # from michael's data file:
  # uniqueID
  # ground truth position (0, 1, 2, 3 = 2, 4, 6, 8 = south, west, east, north)
  # response
  # spatial level where contrast was increased (target band, 0 -- 5).
  # 6x2x4 columns: for the four potential target locations
  # (south, west, east, north), mean pixel energy in target 
  # patch + double-sized surround . 
  # Layout: (band 0, 4xcentre energies, 4x surround energies), (band 1, 4x centre, 4x surround), ...
  colNames <- c("unique_id","position","response","spatialLevel")
  
  # set up column labels for all band energy:
  positions <- c('S', 'W', 'E', 'N')
  for(band in 1 : 6){
    for(type in 1 : 2){
      for(position in 1 : 4){
        if (type == 1){
          label <- paste0('c_centre_band_',band-1,'_pos_',positions[position])  
        } else {
          label <- paste0('c_surround_band_',band-1,'_pos_',positions[position])  
        }
        colNames <- c(colNames,label)
      }      
    }
  }
  
  # dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
  
  thisDat <- read.table(fileName,col.names=colNames)
  
  # drop position and response columsn:
  thisDat$position <- NULL
  thisDat$response <- NULL
  thisDat$spatialLevel <- NULL
  
  # # which cases don't match?
  # full_dat <- merge(dat,thisDat,by="unique_id",all.x=TRUE)
  # exclude_dat <- merge(dat,thisDat,by="unique_id")
  # matches <- match(full_dat$unique_id,exclude_dat$unique_id,nomatch=0)
  # non_matches <- full_dat$unique_id[matches==0]
  # write.csv(non_matches,file=paste(getwd(),'/data/non_matches.csv',sep=""))
  
  # merge data frames based on the unique_id column:
  ifelse(include_NAs,{
    dat <- merge(dat,thisDat,by="unique_id",all.x=TRUE)
  },{
    dat <- merge(dat,thisDat,by="unique_id")
  })
  
  rm(thisDat)
  
  # cumulatively add the rest of the features to the basic data frame --------------------------------------------
  addFeatures <- function(fileName,colNames,dat,include_NAs=FALSE){
    # read data file:
    thisDat <- read.table(fileName,col.names=colNames)
    
    # merge data frames based on the unique_id column:
    ifelse(include_NAs,{
      dat <- merge(dat,thisDat,by="unique_id",all.x=TRUE)
    },{
      dat <- merge(dat,thisDat,by="unique_id")
    })
    return(dat)
  }
  
  #--------------------------------------------
  # eye movements:
  fileName <- paste(data_directory,'csf_em_data',sep="")
  
  # from michael's data file:
  # ID
  # cumulative EM during trial (deg)
  # distance (deg) eye from stim on- to offset
  # time (us) from stim onset to previous saccade onset
  # time (us) from stim onset to previous saccade offset
  # amplitude (deg) previous saccade
  # direction previous saccade (0: right; PI/2: down; -PI/2: up; PI: left))
  # time (us) from stim onset to next saccade
  # amplitude (deg) next saccade
  # direction next saccade
  # boolean flag peri-saccadic trial (1=saccade ended after stim onset; 2=saccade started before offset)
  colNames <- c("unique_id","em_cumDist","em_dist_StimOnToOff",
                "em_prevSacc_time_StimOnToOn","em_prevSacc_time_StimOnToOff","em_prevSacc_amp","em_prevSacc_dir",
                "em_nextSacc_time_StimOnToOn","em_nextSacc_amp","em_nextSacc_dir",
                "em_perisacc_boolean",
                "em_time_stimOfftoOn", "em_time_stimOfftoOff")
  dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
  
  #--------------------------------------------
  # brightness:
  fileName <- paste(data_directory,'csf_brightness_data',sep="")
  
  # from michael's data file:
  # ID
  # mean feature intensity target patch south
  # mean feature intensity target patch west
  # mean feature intensity target patch east
  # mean feature intensity target patch north
  colNames <- c("unique_id","lum_S","lum_W","lum_E","lum_N")
  dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
  
  #--------------------------------------------
  # contrast:
  fileName <- paste(data_directory,'csf_rms_data',sep="")
  
  # from michael's data file:
  # ID
  # mean feature intensity target patch south
  # mean feature intensity target patch west
  # mean feature intensity target patch east
  # mean feature intensity target patch north
  colNames <- c("unique_id","rms_S","rms_W","rms_E","rms_N")
  dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
  
  #--------------------------------------------
  # 2D invariants
  # for H and K, at 6 spatial scales.
  
  invariant <- c('H','K')
  scales <- c('00','10','20','30','40','50')
  
  for (i in 1:length(invariant)){
    for (j in 1:length(scales)){
      fileName <- paste(data_directory,'csf_invariants2D_',invariant[i],'_',scales[j],'_data',sep="")
      colNames <- c("unique_id",
                    paste(invariant[i],'_2D_',scales[j],'_S',sep=""),
                    paste(invariant[i],'_2D_',scales[j],'_W',sep=""),
                    paste(invariant[i],'_2D_',scales[j],'_E',sep=""),
                    paste(invariant[i],'_2D_',scales[j],'_N',sep=""))
      dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
    }
  }
  
  #--------------------------------------------
  # 3D invariants
  # for H, S and K, at 6 spatial scales.
  
  invariant <- c('H','S','K')
  spatial_scale <- c('0','1','2','3','4','5')
  temporal_scale <- c('0','1','2')
  
  for (i in 1:length(invariant)){
    for (j in 1:length(spatial_scale)){
      for (k in 1:length(temporal_scale)){
        fileName <- paste(data_directory,'csf_invariants3D_',invariant[i],'_',spatial_scale[j],temporal_scale[k],'_data',sep="")
        colNames <- c("unique_id",
                      paste(invariant[i],'_3D_',spatial_scale[j],temporal_scale[k],'_S',sep=""),
                      paste(invariant[i],'_3D_',spatial_scale[j],temporal_scale[k],'_W',sep=""),
                      paste(invariant[i],'_3D_',spatial_scale[j],temporal_scale[k],'_E',sep=""),
                      paste(invariant[i],'_3D_',spatial_scale[j],temporal_scale[k],'_N',sep=""))
        dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
      }
    }
  }
  
  #--------------------------------------------
  # template matching:
  fileName <- paste(data_directory,'csf_template_data',sep="")
  
  # from michael's data file:
  # ID
  # position
  # response
  # mean feature intensity target patch south
  # mean feature intensity target patch west
  # mean feature intensity target patch east
  # mean feature intensity target patch north
  colNames <- c("unique_id","position","response","template_S","template_W","template_E","template_N")
  # dat <- addFeatures(fileName,colNames,dat,include_NAs=include_NAs)
  
  thisDat <- read.table(fileName,col.names=colNames)
  
  # drop position and response columsn:
  thisDat$position <- NULL
  thisDat$response <- NULL
  
  # merge data frames based on the unique_id column:
  ifelse(include_NAs,{
    dat <- merge(dat,thisDat,by="unique_id",all.x=TRUE)
  },{
    dat <- merge(dat,thisDat,by="unique_id")
  })
  
  # Variable transforms--------------------------------------------
  # create column for single-stream features at target location:
  features <- c('lum','rms','template')
  positions <- c(2,4,6,8)
  position_names <- c('S','W','E','N')
  for (f in 1 : length(features)){
    text <- paste('dat$',features[f],'_target <- NA',sep="") 
    eval(parse(text=text))
    for (pos in 1 : 4){
      text <- paste('dat$',features[f],'_target[dat$stim_pos==',positions[pos],'] <- dat$',features[f],'_',position_names[pos],'[dat$stim_pos==',positions[pos],']',sep="")
      eval(parse(text=text))
    }
  }
  
  # create column for 2D invariants at target location, across spatial scales:
  features <- c('H','K')
  scales <- c('00','10','20','30','40','50')
  positions <- c(2,4,6,8)
  position_names <- c('S','W','E','N')
  for (f in 1 : length(features)){
    for (s in 1 : length(scales)){
      text <- paste('dat$',features[f],'_2D_',scales[s],'_target <- NA',sep="") 
      eval(parse(text=text))
      for (pos in 1 : 4){
        text <- paste('dat$',features[f],'_2D_',scales[s],'_target[dat$stim_pos==',positions[pos],'] <- dat$',features[f],'_2D_',scales[s],'_',position_names[pos],'[dat$stim_pos==',positions[pos],']',sep="")
        eval(parse(text=text))
      }    
    }
  }
  
  
  # create column for 3D invariants at target location, across spatial scales:
  features <- c('H','S','K')
  spatial_scale <- c('0','1','2','3','4','5')
  temporal_scale <- c('0','1','2')
  positions <- c(2,4,6,8)
  position_names <- c('S','W','E','N')
  for (f in 1 : length(features)){
    for (i in 1 : length(spatial_scale)){
      for (j in 1: length(temporal_scale)){
        text <- paste('dat$',features[f],'_3D_',spatial_scale[i],temporal_scale[j],'_target <- NA',sep="") 
        eval(parse(text=text))
        for (pos in 1 : 4){
          text <- paste('dat$',features[f],'_3D_',spatial_scale[i],temporal_scale[j],'_target[dat$stim_pos==',positions[pos],'] <- dat$',features[f],'_3D_',spatial_scale[i],temporal_scale[j],'_',position_names[pos],'[dat$stim_pos==',positions[pos],']',sep="")
          eval(parse(text=text))
        }          
      }
    }
  }
  
  # create column for band energies at target location:
  features <- c('c_centre','c_surround')
  positions <- c(2,4,6,8)
  position_names <- c('S','W','E','N')
  bands <- 0:5
  for (f in 1 : length(features)){
    # target band:
    text <- paste0('dat$',features[f],'_target <- NA') 
    eval(parse(text=text))
    
    for (band in 1 : length(bands)){
      # for each spatial band, including target band:
      text <- paste0('dat$',features[f],'_band_',bands[band],'_target <- NA') 
      eval(parse(text=text))
      
      for (pos in 1 : 4){
        # target band:
        text <- paste0('dat$',features[f],'_target[dat$stim_pos==',positions[pos],' & dat$target_spatial_band==',bands[band],'] <- dat$',features[f],'_band_',bands[band],'_pos_',position_names[pos],'[dat$stim_pos==',positions[pos],' & dat$target_spatial_band==',bands[band],']')
        eval(parse(text=text))
        
        # for each spatial band, including target band:
        text <- paste0('dat$',features[f],'_band_',bands[band],'_target[dat$stim_pos==',positions[pos],'] <- dat$',features[f],'_band_',bands[band],'_pos_',position_names[pos],'[dat$stim_pos==',positions[pos],']')
        eval(parse(text=text))
      }
    }      
  }
  
  
  # change boolean perisaccadic variable to a factor:
  dat$em_perisacc_boolean <- factor(dat$em_perisacc_boolean)
  levels(dat$em_perisacc_boolean) <- c("noSaccade","saccEndAfterStimOnset","saccStartBeforeStimOffset")
  
  # Convert times to milliseconds:
  dat$em_prevSacc_time_StimOnToOn <- dat$em_prevSacc_time_StimOnToOn/1000
  dat$em_prevSacc_time_StimOnToOff <- dat$em_prevSacc_time_StimOnToOff/1000
  dat$em_nextSacc_time_StimOnToOn <- dat$em_nextSacc_time_StimOnToOn/1000
  dat$em_time_stimOfftoOn <- dat$em_time_stimOfftoOn/1000
  dat$em_time_stimOfftoOff <- dat$em_time_stimOfftoOff/1000
  
  # signed angle between stimulus position and saccade direction:
  
  # rotated saccade directions to be relative to target direction (with target direction always north -- -pi/2):
  directions <- c(pi/2,pi,0,-pi/2) # for targets 2, 4, 6 and 8.
  
  targetPositions <- c(2,4,6,8)
  dat$em_prevSacc_dir_relative <- NA
  dat$em_nextSacc_dir_relative <- NA
  
  directions <- directions 
  
  for (targPos in 1 : 4){
    # signed angle = atan2(by,bx) - atan2(ay,ax). In the way I do this here, a positive difference = clockwise.
    prevSacc <- dat$em_prevSacc_dir[dat$stim_pos==targetPositions[targPos]]
    nextSacc <- dat$em_nextSacc_dir[dat$stim_pos==targetPositions[targPos]]
    targAngle <- directions[targPos]
    
    diffs_prev <- atan2(sin(prevSacc),cos(prevSacc)) - atan2(sin(targAngle),cos(targAngle))
    diffs_next <- atan2(sin(nextSacc),cos(nextSacc)) - atan2(sin(targAngle),cos(targAngle))
    
    # add / subtract 2*pi if signed value is larger than pi (don't care about angle so much as wraparound):
    diffs_prev[diffs_prev>pi & is.na(diffs_prev)!=1] <- diffs_prev[diffs_prev>pi & is.na(diffs_prev)!=1] - 2*pi
    diffs_prev[diffs_prev<(-pi) & is.na(diffs_prev)!=1] <- diffs_prev[diffs_prev<(-pi) & is.na(diffs_prev)!=1] + 2*pi
    diffs_next[diffs_next>pi & is.na(diffs_next)!=1] <- diffs_next[diffs_next>pi & is.na(diffs_next)!=1] - 2*pi
    diffs_next[diffs_next<(-pi) & is.na(diffs_next)!=1] <- diffs_next[diffs_next<(-pi) & is.na(diffs_next)!=1] + 2*pi
    
    # stick back into data frame:
    dat$em_prevSacc_dir_relative[dat$stim_pos==targetPositions[targPos]] <- diffs_prev
    dat$em_nextSacc_dir_relative[dat$stim_pos==targetPositions[targPos]] <- diffs_next
  }
  
  # # test code for this:
  # testAngle <- -pi/2
  # saccs <- seq(-pi,pi,by=pi/4)
  # diffs <-  atan2(sin(saccs),cos(saccs)) - atan2(sin(testAngle),cos(testAngle))
  # diffs/pi
  # diffs[diffs>pi] <- diffs[diffs>pi] - 2*pi
  # diffs[diffs<(-pi)] <- diffs[diffs<(-pi)] + 2*pi
  
  
  ############################################################################################################################
  # Other cleanup stuff -------------------------------------------------------------------------------------------------
  ############################################################################################################################
  
  # calculate the increment contrast value by multiplying alpha and pedestal then taking the difference:
  dat$increment <- (dat$alpha*dat$c_centre_target) - dat$c_centre_target
  
  # express sf as the middle frequency of the band:
  sf <- dat$target_spatial_band
  sf[dat$target_spatial_band==0] <- 18 # 12-24 band
  sf[dat$target_spatial_band==1] <- 9 # 6-12
  sf[dat$target_spatial_band==2] <- 4.5 # 3-6
  sf[dat$target_spatial_band==3] <- 2.25 # 1.5-3
  sf[dat$target_spatial_band==4] <- 1.125 # .75-1.5
  sf[dat$target_spatial_band==5] <- 0.5625 # .375-.75
  
  dat$sf <- sf
  
  dat$sf_factor <- factor(dat$sf)
  levels(dat$sf_factor) <- c('.375-.75','.75-1.5','1.5-3','3-6','6-12','12-24')
  
  dat$alpha_factor <- factor(dat$alpha)
  
  
  # save to file.
  save(dat,file=output_file)
  
  # Print summary --------------------------------------------
  library(plyr)
  sum_dat <- ddply(dat,.(subject),summarise,N=length(correct))
  print(sum_dat)
  
  # save reduced subset to file.
  key_data_frame(dat)
  
}



# function to reduce the variables to key ones used in the contrast paper:
key_data_frame <- function(dat){
  dat <- subset(dat, select = c(unique_id : resp_factor, c_centre_target : c_surround_band_5_target, increment : alpha_factor))
  # exclude observers with few trials-----------------------------
  
  # exclude TW, AJ and LB (few trials)
  dat <- dat[dat$subject!='AJ',]
  dat <- dat[dat$subject!='LB',]
  dat <- dat[dat$subject!='TSAW',]
  dat$subject <- factor(dat$subject) # get rid of empty levels.
  
  # rename subjects to anonymous codes:
  levels(dat$subject) <- c('S1','S2','S3','S4','S5') # corresponds to MACD, PJB, LAL, CPT and AM.
  
  
  # drop missing cases (should only be a few trials)
  dat <- na.omit(dat)
  
  output_file <- paste0(getwd(),'/output/csf_data_reduced.RData')
  
  save(dat,file=output_file)
  
  # Print summary --------------------------------------------
  library(plyr)
  sum_dat <- ddply(dat,.(subject),summarise,N=length(correct))
  print(sum_dat)
  
  return(dat)
}
