# Install and load required packages if they are not already installed
R_packages <- c("ggplot2", "chromote", "dplyr", "lubridate", "DT", "tidyr", "plotly", "shiny", "kableExtra", "reticulate", "webshot2", "gdata", "bslib", "pagedown", "shinyMatrix", "shinybusy")
for(i in 1:length(R_packages)){
  is.present <- R_packages[i] %in% data.frame(installed.packages())$Package
  if(isFALSE(is.present)){
    install.packages(R_packages[i])
  }
  lapply(R_packages[i], library, character.only = TRUE)
}

# Create venv if it doesn't exist
if(isFALSE(virtualenv_exists("./python_venv"))){
  install_python(version = "3.10.11")
  virtualenv_create("./python_venv")
}

# install python packages into venv if they aren't already
py_packages <- c("pydicom==2.4.0", "rdsr_navigator==0.5.0", "pandas==2.2.3")
for(i in 1:length(py_packages)){
  is.present <- py_packages[i] %in% py_list_packages("./python_venv")$requirement
  if(isFALSE(is.present)){
    virtualenv_install("./python_venv", packages = py_packages[i])
  }
}

options(chromote.headless = "new")
# options(shiny.error = browser)

use_virtualenv("./python_venv", required = TRUE)
pydicom <- import("pydicom")
nav <- import("rdsr_navigator")
pandas <- import("pandas")

load("environment.RData")
reticulate::source_python("read_RDSR_irradiation_events.py")
reticulate::source_python("read_RDSR_acc_dose_data.py")

do_math <- function(RDSR,
                    Manufacturer,
                    Height_Of_System,
                    Table_Lateral_home_position,
                    Table_Longitudinal_home_position,
                    Table_height_home_position,
                    dh,
                    HR,
                    TL,
                    table_width,
                    DAPerrorA,
                    DAPerrorB,
                    hs_length,
                    hs_width) {
  
  d <- drop_na(RDSR, Dose_RP) # Remove rows containing NA values
  d <- filter(d, Dose_RP > 0, Dose_Area_Product > 0) # remove rows with 0 values in the Dose_RP and Dose_Area_Product columns as these doses are essentially negligible and field size is unable to be calculated
  d$Irradiation_Event <- seq(1:nrow(d))
  
  planes <- sort(unique(d$Acquisition_Plane_in_Irradiation_Event))
  
  # Define the home position of the table
  d$Table_Lateral_home_position <- Table_Lateral_home_position*1000 #mm
  d$Table_Longitudinal_home_position <- Table_Longitudinal_home_position*1000 #mm
  d$Table_height_home_position <- Table_height_home_position*1000 #mm
  d$Distance_Source_to_Detector <- ifelse(is.na(d$Distance_Source_to_Detector), round(mean(drop_na(d, Distance_Source_to_Detector)$Distance_Source_to_Detector), 0), d$Distance_Source_to_Detector) # if Distance_Source_to_Detector contains NA values, replace with the average of the non-NA values (usually Philips issue).
  d$Distance_Source_to_Detector <- ifelse(d$Distance_Source_to_Detector == 0, round(mean(drop_na(d, Distance_Source_to_Detector)$Distance_Source_to_Detector), 0), d$Distance_Source_to_Detector) # if Distance_Source_to_Detector contains values of 0, replace with the average of the non-0 values (usually Philips issue).
  d$Distance_Source_to_Detector <- d$Distance_Source_to_Detector/1000 # convert mm to m
  
  d$KAP <- d$Dose_Area_Product # Gy.m2
  d$K_air <- d$Dose_RP # Gy
  d$FS_RP <- d$KAP/d$K_air  # calculate field size at reference point m2
  d$FS_RP <- ifelse(is.na(d$FS_RP), 0.01, d$FS_RP)
  d$FS_RP <- ifelse(is.infinite(d$FS_RP), 0.01, d$FS_RP)
  d$FS_detector <- sqrt(d$FS_RP*((d$Distance_Source_to_Detector/((d$Distance_Source_to_Isocenter-150)/1000))^2)) # calculate the field size at the detector (m2)
  d$KVP <- ifelse(is.na(d$KVP), round(mean(drop_na(d, KVP)$KVP), 0), d$KVP) # if kVp is NA, assign it the average value of all the non-NA values
  d$Acquisition_Plane_in_Irradiation_Event <- ifelse(is.na(d$Acquisition_Plane_in_Irradiation_Event), "Single Plane", d$Acquisition_Plane_in_Irradiation_Event) # if plane is NA, assign it "Single Plane"
  planes <- sort(unique(d$Acquisition_Plane_in_Irradiation_Event))
  d$cal <- 0
  
  for(i in 1:length(planes)){
    plane <- planes[i]
    if (i == 1) {DAPerror <- DAPerrorA}
    else {DAPerror <- DAPerrorB}
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP < 60, DAPerror[1],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 60 & d$KVP <= 70, DAPerror[8],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 70 & d$KVP <= 80, DAPerror[15],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 80 & d$KVP <= 90, DAPerror[22],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 90 & d$KVP <= 100, DAPerror[29],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 100 & d$KVP <= 110, DAPerror[36],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 110 & d$KVP <= 120, DAPerror[43],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector >= 0.35 & d$KVP >= 120, DAPerror[50], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP < 60, DAPerror[2],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 60 & d$KVP < 70, DAPerror[9],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 70 & d$KVP < 80, DAPerror[16],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 80 & d$KVP < 90, DAPerror[23],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 90 & d$KVP < 100, DAPerror[30],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 100 & d$KVP < 110, DAPerror[37],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 110 & d$KVP < 120, DAPerror[44],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.35 & d$FS_detector >= 0.3 & d$KVP >= 120, DAPerror[51], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP < 60, DAPerror[3],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 60 & d$KVP < 70, DAPerror[10],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 70 & d$KVP < 80, DAPerror[17],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 80 & d$KVP < 90, DAPerror[24],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 90 & d$KVP < 100, DAPerror[31],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 100 & d$KVP < 110, DAPerror[38],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 110 & d$KVP < 120, DAPerror[45],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.3 & d$FS_detector >= 0.25 & d$KVP >= 120, DAPerror[52], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP < 60, DAPerror[4],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 60 & d$KVP < 70, DAPerror[11],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 70 & d$KVP < 80, DAPerror[18],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 80 & d$KVP < 90, DAPerror[25],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 90 & d$KVP < 100, DAPerror[32],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 100 & d$KVP < 110, DAPerror[39],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 110 & d$KVP < 120, DAPerror[46],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.25 & d$FS_detector >= 0.2 & d$KVP >= 120, DAPerror[53], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP < 60, DAPerror[5],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 60 & d$KVP < 70, DAPerror[12],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 70 & d$KVP < 80, DAPerror[19],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 80 & d$KVP < 90, DAPerror[26],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 90 & d$KVP < 100, DAPerror[33],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 100 & d$KVP < 110, DAPerror[40],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 110 & d$KVP < 120, DAPerror[47],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.2 & d$FS_detector >= 0.15 & d$KVP >= 120, DAPerror[54], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP < 60, DAPerror[6],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 60 & d$KVP < 70, DAPerror[13],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 70 & d$KVP < 80, DAPerror[20],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 80 & d$KVP < 90, DAPerror[27],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 90 & d$KVP < 100, DAPerror[34],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 100 & d$KVP < 110, DAPerror[41],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 110 & d$KVP < 120, DAPerror[48],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.15 & d$FS_detector >= 0.1 & d$KVP >= 120, DAPerror[55], d$cal))))))))
    
    d$cal <- ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP < 60, DAPerror[7],
                    ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 60 & d$KVP < 70, DAPerror[14],
                           ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 70 & d$KVP < 80, DAPerror[21],
                                  ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 80 & d$KVP < 90, DAPerror[28],
                                         ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 90 & d$KVP < 100, DAPerror[35],
                                                ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 100 & d$KVP < 110, DAPerror[42],
                                                       ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 110 & d$KVP < 120, DAPerror[49],
                                                              ifelse(d$Acquisition_Plane_in_Irradiation_Event == plane & d$FS_detector < 0.1 & d$KVP >= 120, DAPerror[56], d$cal))))))))
  }
  d$cal <- as.numeric(d$cal) # not sure why this is not numeric already. It causes an error later in the code for being a character class
  d$`X-Ray_Filter_Thickness_Minimum` <- ifelse(is.na(d$`X-Ray_Filter_Thickness_Minimum`), 0, d$`X-Ray_Filter_Thickness_Minimum`) # reassign NA values as 0
  
  # Calculate the nearest BSF to the factors in each irradiation event
  for(i in 1:nrow(d)) {
    # The order in which the data is filtered is based on the order of highest to lowest impact on the BSF
    # sort the data for the nearest field size
    BSF_FS_nearest <- BSF_data$field_size_m2[which.min(abs(BSF_data$field_size_m2-d$FS_RP[i]))]
    BSF_data_filtered <- subset(BSF_data, field_size_m2 == BSF_FS_nearest)
    # sort the data for the nearest kV
    BSF_kV_nearest <- BSF_data_filtered$kVp[which.min(abs(BSF_data_filtered$kVp-d$KVP[i]))]
    BSF_data_filtered <- subset(BSF_data_filtered, kVp == BSF_kV_nearest)
    # sort the data for the nearest Cu filtration
    BSF_Cu_nearest <- BSF_data_filtered$mm_Cu[which.min(abs(BSF_data_filtered$mm_Cu-d$`X-Ray_Filter_Thickness_Minimum`[i]))]
    BSF_data_filtered <- subset(BSF_data_filtered, mm_Cu == BSF_Cu_nearest)
    # sort the data for the nearest inherent filtration
    # BSF_IF_nearest <- BSF_data_filtered$mm_Al_inherent_filtration[which.min(abs(BSF_data_filtered$mm_Al_inherent_filtration-inherent_filtration))]
    # BSF_data_filtered <- subset(BSF_data_filtered, mm_Al_inherent_filtration == BSF_IF_nearest)
    BSF <- BSF_data_filtered$BSF
    f <- BSF_data_filtered$f
    d$BSF[i] <- BSF
    d$f[i] <- f
    TTF_kV_nearest <- TTF_data$kVp[which.min(abs(TTF_data$kVp-d$KVP[i]))]
    TTF_Cu_nearest <- TTF_data$mm_Cu[which.min(abs(TTF_data$mm_Cu-d$`X-Ray_Filter_Thickness_Minimum`[i]))]
    TTF_FS_nearest <- TTF_data$field_size_m2[which.min(abs(TTF_data$field_size_m2-d$FS_RP[i]))]
    d$TTF[i] <- subset(TTF_data, kVp == TTF_kV_nearest & mm_Cu == TTF_Cu_nearest & field_size_m2 == TTF_FS_nearest)$TTF
  }
  
  # Define variables outlined in Table 1
  d$a1 <- -d$Positioner_Primary_Angle*(pi/180) # radians
  d$a2 <- d$Positioner_Secondary_Angle*(pi/180) # radians
  d$Distance_Source_to_Isocenter <- ifelse(is.na(d$Distance_Source_to_Isocenter), round(mean(drop_na(d, Distance_Source_to_Isocenter)$Distance_Source_to_Isocenter), 0), d$Distance_Source_to_Isocenter) # if Distance_Source_to_Isocenter contains NA values, replace with the average of the non-NA values (usually Philips issue).
  d$Distance_Source_to_Isocenter <- ifelse(d$Distance_Source_to_Isocenter == 0, round(mean(drop_na(d, Distance_Source_to_Isocenter)$Distance_Source_to_Isocenter), 0), d$Distance_Source_to_Isocenter) # if Distance_Source_to_Isocenter contains values of 0, replace with the average of the non-0 values (usually Philips issue).
  d$radius <- d$Distance_Source_to_Isocenter/1000 # Dist Source to the C-arm Isocenter
  d$TLongInc <- (d$Table_Longitudinal_Position+d$Table_Longitudinal_home_position)/1000 # Table longitudinal increment
  d$TLatInc <- (d$Table_Lateral_Position+d$Table_Lateral_home_position)/1000 # Table lateral increment
  d$THeightInc <- ifelse(Manufacturer == "Philips" | Manufacturer == "Philips Medical Systems", (Height_Of_System-d$Table_Height_Position+d$Table_height_home_position)/1000,
                         (d$Table_Height_Position+d$Table_height_home_position)/1000) # Table vertical displacement (m)
  
  d$dh <- dh # Horiz dist from patient origin to the head of table
  d$TL <- TL # Table Length
  d$table_width <- table_width
  d$HR <- HR # Hit rejection plane
  
  # Define variables outlined in Table 2
  d$IstoTable <- d$TL/2 # IStoTable is the Horizontal distance from isocenter toward patient head to the edge of table—head first position
  d$distance <- abs(d$IstoTable) # Distance from isocenter to top head of the patient
  
  # GEOMETRIC CALCULATIONS
  # Isocenter
  d$x_Is <- 0 # eq 1 (modified)
  d$y_Is <- 0 # eq 2
  d$z_Is <- 0
  
  # X-ray source
  d$x_spot <- -d$radius*sin(d$a2) # eq 4 (modified)
  d$y_spot <- d$radius*sin(d$a1)*cos(d$a2) # eq 5 (modified)
  d$z_spot <- d$radius*cos(d$a1)*cos(d$a2) # eq 6
  
  # Patient Entrance Reference Point (PERP)
  d$x_RF <- -0.15*sin(d$a2) # eq 7 (modified)
  d$y_RF <- 0.15*sin(d$a1)*cos(d$a2) # eq 8
  d$z_RF <- 0.15*cos(d$a1)*cos(d$a2) # eq 9
  # The length 0.15m in eq 7,8,9 is the distance from PERP to the X-ray source
  
  # eq 14
  d$v1_x <- cos(d$a2)
  d$v1_y <- sin(d$a1)*sin(d$a2)
  d$v1_z <- cos(d$a1)*sin(d$a2)
  
  # eq 15
  d$v2_x <- 0
  d$v2_y <- cos(d$a1)
  d$v2_z <- -sin(d$a1)
  
  # side length of irradiated area
  d$l <- sqrt(d$FS_RP) # eq 16
  
  # Detector Coords
  d$x_det <- (-d$x_spot+d$x_RF)/(d$radius-0.15)*d$Distance_Source_to_Detector+d$x_spot # eq 7 (modified)
  d$y_det <- (-d$y_spot+d$y_RF)/(d$radius-0.15)*d$Distance_Source_to_Detector+d$y_spot # eq 7 (modified)
  d$z_det <- (-d$z_spot+d$z_RF)/(d$radius-0.15)*d$Distance_Source_to_Detector+d$z_spot # eq 7 (modified)
  
  # side length at detector
  d$l_d <- d$l/(d$radius-0.15)*d$Distance_Source_to_Detector
  
  # corners at detector
  # corner 1, eq 17
  d$x_dc1 <- d$x_det-d$l_d/2*d$v1_x+d$l_d/2*d$v2_x
  d$y_dc1 <- d$y_det-d$l_d/2*d$v1_y+d$l_d/2*d$v2_y
  d$z_dc1 <- d$z_det-d$l_d/2*d$v1_z+d$l_d/2*d$v2_z
  
  # corner 2, eq 18
  d$x_dc2 <- d$x_det+d$l_d/2*d$v1_x-d$l_d/2*d$v2_x
  d$y_dc2 <- d$y_det+d$l_d/2*d$v1_y-d$l_d/2*d$v2_y
  d$z_dc2 <- d$z_det+d$l_d/2*d$v1_z-d$l_d/2*d$v2_z
  
  # corner 3, eq 19
  d$x_dc3 <- d$x_det-d$l_d/2*d$v1_x-d$l_d/2*d$v2_x
  d$y_dc3 <- d$y_det-d$l_d/2*d$v1_y-d$l_d/2*d$v2_y
  d$z_dc3 <- d$z_det-d$l_d/2*d$v1_z-d$l_d/2*d$v2_z
  
  # corner 4, eq 20
  d$x_dc4 <- d$x_det+d$l_d/2*d$v1_x+d$l_d/2*d$v2_x
  d$y_dc4 <- d$y_det+d$l_d/2*d$v1_y+d$l_d/2*d$v2_y
  d$z_dc4 <- d$z_det+d$l_d/2*d$v1_z+d$l_d/2*d$v2_z
  
  # Hit rejection Coords
  d$x_HR <- (-d$x_spot+d$x_RF)/(d$radius-0.15)*(d$Distance_Source_to_Detector-HR)+d$x_spot # eq 7 (modified)
  d$y_HR <- (-d$y_spot+d$y_RF)/(d$radius-0.15)*(d$Distance_Source_to_Detector-HR)+d$y_spot # eq 7 (modified)
  d$z_HR <- (-d$z_spot+d$z_RF)/(d$radius-0.15)*(d$Distance_Source_to_Detector-HR)+d$z_spot # eq 7 (modified)
  
  # side length at hit rejection point
  d$l_HR <- d$l/(d$radius-0.15)*(d$Distance_Source_to_Detector-HR)
  
  # corners at hit rejection point
  # corner 1, eq 17
  d$x_HR1 <- d$x_HR-d$l_HR/2*d$v1_x+d$l_HR/2*d$v2_x
  d$y_HR1 <- d$y_HR-d$l_HR/2*d$v1_y+d$l_HR/2*d$v2_y
  d$z_HR1 <- d$z_HR-d$l_HR/2*d$v1_z+d$l_HR/2*d$v2_z
  
  # corner 2, eq 18
  d$x_HR2 <- d$x_HR+d$l_HR/2*d$v1_x-d$l_HR/2*d$v2_x
  d$y_HR2 <- d$y_HR+d$l_HR/2*d$v1_y-d$l_HR/2*d$v2_y
  d$z_HR2 <- d$z_HR+d$l_HR/2*d$v1_z-d$l_HR/2*d$v2_z
  
  # corner 3, eq 19
  d$x_HR3 <- d$x_HR-d$l_HR/2*d$v1_x-d$l_HR/2*d$v2_x
  d$y_HR3 <- d$y_HR-d$l_HR/2*d$v1_y-d$l_HR/2*d$v2_y
  d$z_HR3 <- d$z_HR-d$l_HR/2*d$v1_z-d$l_HR/2*d$v2_z
  
  # corner 4, eq 20
  d$x_HR4 <- d$x_HR+d$l_HR/2*d$v1_x+d$l_HR/2*d$v2_x
  d$y_HR4 <- d$y_HR+d$l_HR/2*d$v1_y+d$l_HR/2*d$v2_y
  d$z_HR4 <- d$z_HR+d$l_HR/2*d$v1_z+d$l_HR/2*d$v2_z
  
  # Beam pyramid (this pyramid is different from that in the paper as it uses the source and detector to define the pyramid to define the corners, not the coordinates at the PERP.)
  # Face 1 (c1, c4)
  d$A1 <- d$y_spot*(d$z_dc1-d$z_dc4)+d$y_dc1*(d$z_dc4-d$z_spot)+d$y_dc4*(d$z_spot-d$z_dc1) # eq 21
  d$B1 <- d$z_spot*(d$x_dc1-d$x_dc4)+d$z_dc1*(d$x_dc4-d$x_spot)+d$z_dc4*(d$x_spot-d$x_dc1) # eq 22
  d$C1 <- d$x_spot*(d$y_dc1-d$y_dc4)+d$x_dc1*(d$y_dc4-d$y_spot)+d$x_dc4*(d$y_spot-d$y_dc1) # eq 23
  d$D1 <- -d$x_spot*(d$y_dc1*d$z_dc4-d$y_dc4*d$z_dc1)-d$x_dc1*(d$y_dc4*d$z_spot-d$y_spot*d$z_dc4)-d$x_dc4*(d$y_spot*d$z_dc1-d$y_dc1*d$z_spot) # eq 24
  
  # Face 2 (c4, c2)
  d$A2 <- d$y_spot*(d$z_dc4-d$z_dc2)+d$y_dc4*(d$z_dc2-d$z_spot)+d$y_dc2*(d$z_spot-d$z_dc4) # eq 21
  d$B2 <- d$z_spot*(d$x_dc4-d$x_dc2)+d$z_dc4*(d$x_dc2-d$x_spot)+d$z_dc2*(d$x_spot-d$x_dc4) # eq 22
  d$C2 <- d$x_spot*(d$y_dc4-d$y_dc2)+d$x_dc4*(d$y_dc2-d$y_spot)+d$x_dc2*(d$y_spot-d$y_dc4) # eq 23
  d$D2 <- -d$x_spot*(d$y_dc4*d$z_dc2-d$y_dc2*d$z_dc4)-d$x_dc4*(d$y_dc2*d$z_spot-d$y_spot*d$z_dc2)-d$x_dc2*(d$y_spot*d$z_dc4-d$y_dc4*d$z_spot) # eq 24
  
  # Face 3 (c2, c3)
  d$A3 <- d$y_spot*(d$z_dc2-d$z_dc3)+d$y_dc2*(d$z_dc3-d$z_spot)+d$y_dc3*(d$z_spot-d$z_dc2) # eq 21
  d$B3 <- d$z_spot*(d$x_dc2-d$x_dc3)+d$z_dc2*(d$x_dc3-d$x_spot)+d$z_dc3*(d$x_spot-d$x_dc2) # eq 22
  d$C3 <- d$x_spot*(d$y_dc2-d$y_dc3)+d$x_dc2*(d$y_dc3-d$y_spot)+d$x_dc3*(d$y_spot-d$y_dc2) # eq 23
  d$D3 <- -d$x_spot*(d$y_dc2*d$z_dc3-d$y_dc3*d$z_dc2)-d$x_dc2*(d$y_dc3*d$z_spot-d$y_spot*d$z_dc3)-d$x_dc3*(d$y_spot*d$z_dc2-d$y_dc2*d$z_spot) # eq 24
  
  # Face 4 (c3, c1)
  d$A4 <- d$y_spot*(d$z_dc3-d$z_dc1)+d$y_dc3*(d$z_dc1-d$z_spot)+d$y_dc1*(d$z_spot-d$z_dc3) # eq 21
  d$B4 <- d$z_spot*(d$x_dc3-d$x_dc1)+d$z_dc3*(d$x_dc1-d$x_spot)+d$z_dc1*(d$x_spot-d$x_dc3) # eq 22
  d$C4 <- d$x_spot*(d$y_dc3-d$y_dc1)+d$x_dc3*(d$y_dc1-d$y_spot)+d$x_dc1*(d$y_spot-d$y_dc3) # eq 23
  d$D4 <- -d$x_spot*(d$y_dc3*d$z_dc1-d$y_dc1*d$z_dc3)-d$x_dc3*(d$y_dc1*d$z_spot-d$y_spot*d$z_dc1)-d$x_dc1*(d$y_spot*d$z_dc3-d$y_dc3*d$z_spot) # eq 24
  
  # Define body support coordinates
  d$x_table1 <- d$TLatInc+0.5-hs_length
  d$y_table1 <- d$table_width/2+d$TLongInc
  
  d$x_table4 <- d$TLatInc+0.5-hs_length
  d$y_table4 <- -d$table_width/2+d$TLongInc
  
  d$x_table3 <- d$TLatInc+0.5-d$TL-hs_length
  d$y_table3 <- d$table_width/2+d$TLongInc
  
  d$x_table2 <- d$TLatInc+0.5-d$TL-hs_length
  d$y_table2 <- -d$table_width/2+d$TLongInc
  
  d$z_table <-  d$THeightInc
  
  # Body support pyramid
  # Face 1 (c1, c4)
  d$table_A1 <- d$y_spot*(d$z_table-d$z_table)+d$y_table1*(d$z_table-d$z_spot)+d$y_table4*(d$z_spot-d$z_table) # eq 21
  d$table_B1 <- d$z_spot*(d$x_table1-d$x_table4)+d$z_table*(d$x_table4-d$x_spot)+d$z_table*(d$x_spot-d$x_table1) # eq 22
  d$table_C1 <- d$x_spot*(d$y_table1-d$y_table4)+d$x_table1*(d$y_table4-d$y_spot)+d$x_table4*(d$y_spot-d$y_table1) # eq 23
  d$table_D1 <- -d$x_spot*(d$y_table1*d$z_table-d$y_table4*d$z_table)-d$x_table1*(d$y_table4*d$z_spot-d$y_spot*d$z_table)-d$x_table4*(d$y_spot*d$z_table-d$y_table1*d$z_spot) # eq 24
  
  # Face 2 (c4, c2)
  d$table_A2 <- d$y_spot*(d$z_table-d$z_table)+d$y_table4*(d$z_table-d$z_spot)+d$y_table2*(d$z_spot-d$z_table) # eq 21
  d$table_B2 <- d$z_spot*(d$x_table4-d$x_table2)+d$z_table*(d$x_table2-d$x_spot)+d$z_table*(d$x_spot-d$x_table4) # eq 22
  d$table_C2 <- d$x_spot*(d$y_table4-d$y_table2)+d$x_table4*(d$y_table2-d$y_spot)+d$x_table2*(d$y_spot-d$y_table4) # eq 23
  d$table_D2 <- -d$x_spot*(d$y_table4*d$z_table-d$y_table2*d$z_table)-d$x_table4*(d$y_table2*d$z_spot-d$y_spot*d$z_table)-d$x_table2*(d$y_spot*d$z_table-d$y_table4*d$z_spot) # eq 24
  
  # Face 3 (c2, c3)
  d$table_A3 <- d$y_spot*(d$z_table-d$z_table)+d$y_table2*(d$z_table-d$z_spot)+d$y_table3*(d$z_spot-d$z_table) # eq 21
  d$table_B3 <- d$z_spot*(d$x_table2-d$x_table3)+d$z_table*(d$x_table3-d$x_spot)+d$z_table*(d$x_spot-d$x_table2) # eq 22
  d$table_C3 <- d$x_spot*(d$y_table2-d$y_table3)+d$x_table2*(d$y_table3-d$y_spot)+d$x_table3*(d$y_spot-d$y_table2) # eq 23
  d$table_D3 <- -d$x_spot*(d$y_table2*d$z_table-d$y_table3*d$z_table)-d$x_table2*(d$y_table3*d$z_spot-d$y_spot*d$z_table)-d$x_table3*(d$y_spot*d$z_table-d$y_table2*d$z_spot) # eq 24
  
  # Face 4 (c3, c1)
  d$table_A4 <- d$y_spot*(d$z_table-d$z_table)+d$y_table3*(d$z_table-d$z_spot)+d$y_table1*(d$z_spot-d$z_table) # eq 21
  d$table_B4 <- d$z_spot*(d$x_table3-d$x_table1)+d$z_table*(d$x_table1-d$x_spot)+d$z_table*(d$x_spot-d$x_table3) # eq 22
  d$table_C4 <- d$x_spot*(d$y_table3-d$y_table1)+d$x_table3*(d$y_table1-d$y_spot)+d$x_table1*(d$y_spot-d$y_table3) # eq 23
  d$table_D4 <- -d$x_spot*(d$y_table3*d$z_table-d$y_table1*d$z_table)-d$x_table3*(d$y_table1*d$z_spot-d$y_spot*d$z_table)-d$x_table1*(d$y_spot*d$z_table-d$y_table3*d$z_spot) # eq 24
  
  # Define head support coordinates
  d$x_hs1 <- d$TLatInc+0.5
  d$y_hs1 <- hs_width/2+d$TLongInc
  
  d$x_hs4 <- d$TLatInc+0.5
  d$y_hs4 <- -hs_width/2+d$TLongInc
  
  d$x_hs3 <- d$TLatInc+0.5-hs_length
  d$y_hs3 <- hs_width/2+d$TLongInc
  
  d$x_hs2 <- d$TLatInc+0.5-hs_length
  d$y_hs2 <- -hs_width/2+d$TLongInc
  
  # Head support pyramid
  # Face 1 (c1, c4)
  d$hs_A1 <- d$y_spot*(d$z_table-d$z_table)+d$y_hs1*(d$z_table-d$z_spot)+d$y_hs4*(d$z_spot-d$z_table) # eq 21
  d$hs_B1 <- d$z_spot*(d$x_hs1-d$x_hs4)+d$z_table*(d$x_hs4-d$x_spot)+d$z_table*(d$x_spot-d$x_hs1) # eq 22
  d$hs_C1 <- d$x_spot*(d$y_hs1-d$y_hs4)+d$x_hs1*(d$y_hs4-d$y_spot)+d$x_hs4*(d$y_spot-d$y_hs1) # eq 23
  d$hs_D1 <- -d$x_spot*(d$y_hs1*d$z_table-d$y_hs4*d$z_table)-d$x_hs1*(d$y_hs4*d$z_spot-d$y_spot*d$z_table)-d$x_hs4*(d$y_spot*d$z_table-d$y_hs1*d$z_spot) # eq 24
  
  # Face 2 (c4, c2)
  d$hs_A2 <- d$y_spot*(d$z_table-d$z_table)+d$y_hs4*(d$z_table-d$z_spot)+d$y_hs2*(d$z_spot-d$z_table) # eq 21
  d$hs_B2 <- d$z_spot*(d$x_hs4-d$x_hs2)+d$z_table*(d$x_hs2-d$x_spot)+d$z_table*(d$x_spot-d$x_hs4) # eq 22
  d$hs_C2 <- d$x_spot*(d$y_hs4-d$y_hs2)+d$x_hs4*(d$y_hs2-d$y_spot)+d$x_hs2*(d$y_spot-d$y_hs4) # eq 23
  d$hs_D2 <- -d$x_spot*(d$y_hs4*d$z_table-d$y_hs2*d$z_table)-d$x_hs4*(d$y_hs2*d$z_spot-d$y_spot*d$z_table)-d$x_hs2*(d$y_spot*d$z_table-d$y_hs4*d$z_spot) # eq 24
  
  # Face 3 (c2, c3)
  d$hs_A3 <- d$y_spot*(d$z_table-d$z_table)+d$y_hs2*(d$z_table-d$z_spot)+d$y_hs3*(d$z_spot-d$z_table) # eq 21
  d$hs_B3 <- d$z_spot*(d$x_hs2-d$x_hs3)+d$z_table*(d$x_hs3-d$x_spot)+d$z_table*(d$x_spot-d$x_hs2) # eq 22
  d$hs_C3 <- d$x_spot*(d$y_hs2-d$y_hs3)+d$x_hs2*(d$y_hs3-d$y_spot)+d$x_hs3*(d$y_spot-d$y_hs2) # eq 23
  d$hs_D3 <- -d$x_spot*(d$y_hs2*d$z_table-d$y_hs3*d$z_table)-d$x_hs2*(d$y_hs3*d$z_spot-d$y_spot*d$z_table)-d$x_hs3*(d$y_spot*d$z_table-d$y_hs2*d$z_spot) # eq 24
  
  # Face 4 (c3, c1)
  d$hs_A4 <- d$y_spot*(d$z_table-d$z_table)+d$y_hs3*(d$z_table-d$z_spot)+d$y_hs1*(d$z_spot-d$z_table) # eq 21
  d$hs_B4 <- d$z_spot*(d$x_hs3-d$x_hs1)+d$z_table*(d$x_hs1-d$x_spot)+d$z_table*(d$x_spot-d$x_hs3) # eq 22
  d$hs_C4 <- d$x_spot*(d$y_hs3-d$y_hs1)+d$x_hs3*(d$y_hs1-d$y_spot)+d$x_hs1*(d$y_spot-d$y_hs3) # eq 23
  d$hs_D4 <- -d$x_spot*(d$y_hs3*d$z_table-d$y_hs1*d$z_table)-d$x_hs3*(d$y_hs1*d$z_spot-d$y_spot*d$z_table)-d$x_hs1*(d$y_spot*d$z_table-d$y_hs3*d$z_spot) # eq 24
  
  return(d)
}

calculate_IE <- function(d, patient_phantom) {
  
  for (i in 1:nrow(d)) {
    patient_phantom$z_flat <- patient_phantom$z + 0.03
    patient_phantom$z_flat <- ifelse(patient_phantom$z_flat > 0, 0, patient_phantom$z_flat)
    
    patient_phantom$x_adj <- patient_phantom$x+d$TLatInc[i]-d$dh[i] # Table lateral increment
    patient_phantom$y_adj <- patient_phantom$y+d$TLongInc[i] # Table longitudinal increment
    patient_phantom$z_adj <- patient_phantom$z_flat+d$THeightInc[i] # Table vertical displacement (m)
    
    patient_phantom$eff_angle <- (pi/2) - atan(abs(d$z_spot[i]-patient_phantom$z_adj)/sqrt(abs(d$x_spot[i]-patient_phantom$x_adj)^2 + abs(d$y_spot[i]-patient_phantom$y_adj)^2)) # effective angle which x-rays would travel through the patient support
    patient_phantom$f_theta <- exp(log(d$TTF[i])*((1/cos(patient_phantom$eff_angle))-1)) # table transmission due to oblique angle
    
    patient_phantom$Dist <- sqrt((patient_phantom$x_adj-d$x_spot[i])^2+(patient_phantom$y_adj-d$y_spot[i])^2+(patient_phantom$z_adj-d$z_spot[i])^2) # eq 10
    patient_phantom$`is point irradiated?` <- ifelse(d$A1[i]*patient_phantom$x_adj+d$B1[i]*patient_phantom$y_adj+d$C1[i]*patient_phantom$z_adj+d$D1[i] > 0 & # Face 1
                                                       d$A2[i]*patient_phantom$x_adj+d$B2[i]*patient_phantom$y_adj+d$C2[i]*patient_phantom$z_adj+d$D2[i] > 0 & # Face 2
                                                       d$A3[i]*patient_phantom$x_adj+d$B3[i]*patient_phantom$y_adj+d$C3[i]*patient_phantom$z_adj+d$D3[i] > 0 & # Face 3
                                                       d$A4[i]*patient_phantom$x_adj+d$B4[i]*patient_phantom$y_adj+d$C4[i]*patient_phantom$z_adj+d$D4[i] > 0 & # Face 4
                                                       patient_phantom$Dist < d$Distance_Source_to_Detector[i]-d$HR[i], # Hit rejection plane
                                                     T, F)
    
    
    patient_phantom$`table attenuation` <- ifelse(d$table_A1[i]*patient_phantom$x_adj+d$table_B1[i]*patient_phantom$y_adj+d$table_C1[i]*patient_phantom$z_adj+d$table_D1[i] > 0 & # body support face 1
                                                    d$table_A2[i]*patient_phantom$x_adj+d$table_B2[i]*patient_phantom$y_adj+d$table_C2[i]*patient_phantom$z_adj+d$table_D2[i] > 0 & # body support face 2
                                                    d$table_A3[i]*patient_phantom$x_adj+d$table_B3[i]*patient_phantom$y_adj+d$table_C3[i]*patient_phantom$z_adj+d$table_D3[i] > 0 & # body support face 3
                                                    d$table_A4[i]*patient_phantom$x_adj+d$table_B4[i]*patient_phantom$y_adj+d$table_C4[i]*patient_phantom$z_adj+d$table_D4[i] > 0 | # body support face 4
                                                    d$hs_A1[i]*patient_phantom$x_adj+d$hs_B1[i]*patient_phantom$y_adj+d$hs_C1[i]*patient_phantom$z_adj+d$hs_D1[i] > 0 & # head support face 1
                                                    d$hs_A2[i]*patient_phantom$x_adj+d$hs_B2[i]*patient_phantom$y_adj+d$hs_C2[i]*patient_phantom$z_adj+d$hs_D2[i] > 0 & # head support face 2
                                                    d$hs_A3[i]*patient_phantom$x_adj+d$hs_B3[i]*patient_phantom$y_adj+d$hs_C3[i]*patient_phantom$z_adj+d$hs_D3[i] > 0 & # head support face 3
                                                    d$hs_A4[i]*patient_phantom$x_adj+d$hs_B4[i]*patient_phantom$y_adj+d$hs_C4[i]*patient_phantom$z_adj+d$hs_D4[i] > 0, # head support face 4
                                                  TRUE, FALSE)
    
    patient_phantom$TTF <- ifelse(patient_phantom$`table attenuation`, d$TTF[i], 1)
    patient_phantom$f_theta <- ifelse(patient_phantom$`table attenuation`, patient_phantom$f_theta, 1)
    patient_phantom$K_skin <- ifelse(patient_phantom$`is point irradiated?`, d$K_air[i]*((d$radius[i]-0.15)^2)/(patient_phantom$Dist^2), 0) # eq 11
    patient_phantom$irradiation_event_skin_dose <- patient_phantom$K_skin*((100-d$cal[i])/100)*d$f[i]*d$BSF[i]*patient_phantom$TTF*patient_phantom$f_theta# eq 13
    patient_phantom$irradiation_event_skin_dose <- ifelse(is.na(patient_phantom$irradiation_event_skin_dose), 0, patient_phantom$irradiation_event_skin_dose)
    patient_phantom$cumulative_skin_dose <- patient_phantom$cumulative_skin_dose+patient_phantom$irradiation_event_skin_dose
  }
  return(patient_phantom)
}

ui <- fluidPage(
  add_busy_spinner(spin = "circle"),
  titlePanel("skin dose calculator"),
  modalDialog(title = "Warning",
              HTML("Calculations performed by this software are an estimation.<br><br>"),
              size = "l",
              easyClose = FALSE,
              modalButton("Acknowledge"),
              footer = NULL),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      strong("Import system configuration"),
      fileInput("import_config", label  = NULL, accept = c(".rds")),
      fileInput("RDSR_dcm", "Import RDSR (.dcm)"),
      radioButtons("input_phantom", label = NULL, choices = c("New skin dose estimation",
                                                              "Append to existing skin dose estimation")),
      conditionalPanel(condition = "input.input_phantom == 'New skin dose estimation'",
                       selectInput("patient_phantom", "Select patient phantom", choices = c("", "ICRP145 Adult Male", "ICRP145 Adult Female", "130 cm, 30 kg", "140 cm, 40 kg", "150 cm, 60 kg","160 cm, 70 kg", "170 cm, 80 kg", "180 cm, 90 kg", "190 cm, 100 kg", "200 cm, 110 kg", "cow"), selected = "")),
      conditionalPanel(condition = "input.input_phantom == 'Append to existing skin dose estimation'",
                       fileInput("patient_phantom", "Upload previous skin dose estimation", accept = c(".rds"))),
      actionButton("calculate_PSD", "Calculate"),
      downloadButton("report", "Download Report"),
      downloadButton("phantom", "Download Data"),
      HTML("<br><br>"),
      splitLayout(cellWidths = c(75, 400),
                  numericInput("HR", NULL, value = 0.15, step = 0.05),
                  strong("Distance from hit rejection plane to detector (m)")),
      splitLayout(cellWidths = c(75, 400),
                  numericInput("dh", NULL, value = 0.15, step = 0.05),
                  strong("Distance from patient head to end of table (m)"))),
    mainPanel(
      width = 9,
      tabsetPanel(type = "tabs",
                  tabPanel("Cumulative Skin Dose",
                           plotlyOutput("cumulative", height = "800px")),
                  tabPanel("Irradiation Event Visualisation",
                           splitLayout(cellWidths = c(120, 80, 100),
                                       strong(HTML("Irradiation Event")),
                                       numericInput("IE", NULL, value = 1, min = 1, width = '75px'),
                                       textOutput("no_of_IE")),
                           plotlyOutput("visualisation", height = "700px")),
                  tabPanel("Procedure Details",
                           splitLayout(cellWidths = c("80%", "20%"),
                                       DT::dataTableOutput("summary"),
                                       plotlyOutput("summary_figure", height = "700px"))),
                  tabPanel("System Configuration",
                           HTML("<br>"),
                           downloadButton("export_config", "Export system config"),
                           HTML("<br><br>"),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("Table_Lateral_home_position", NULL, value = 0, step = 0.05),
                                       strong("Table lateral (head to feet) home position (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("Table_Longitudinal_home_position", NULL, value = 0, step = 0.05),
                                       strong("Table longitudinal (left to right) home position (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("Table_height_home_position", NULL, value = 0, step = 0.05),
                                       strong("Table height home position (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("hs_length", NULL, value = 0.3, min = 0, step = 0.05),
                                       strong("Head support length (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("hs_width", NULL, value = 0.2, min = 0, step = 0.05),
                                       strong("Head support width (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("TL", NULL, value = 2.5, min = 0, step = 0.05),
                                       strong("Body support length (m)")),
                           splitLayout(cellWidths = c(75, 400),
                                       numericInput("table_width", NULL, value = 0.5, min = 0, step = 0.05),
                                       strong("Body support width (m)")),
                           layout_columns(
                             card(card_header(strong("Plane A/Single Plane")),
                                  matrixInput("DAPerrorA",
                                              value = matrix(0, 7, 8, dimnames = list(c("FS≥ 35x35cm", "30x30≤ FS< 35x35cm", "25x25≤ FS< 30x30cm", "20x20≤ FS< 25x25cm", "15x15≤ FS< 20x20cm", "10x10≤ FS< 15x15cm", "FS< 10x10cm"),
                                                                                      c("kV<60", "60≤kV<70", "70≤kV<80", "80≤kV<90", "90≤kV <100", "100≤kV <110", "110≤kV <120", "kV≥120"))),
                                              rows = list(names = TRUE),
                                              cols = list(names = TRUE),
                                              cells = list(editableCells = TRUE))),
                             card(card_header(strong("Plane B (if applicable)")),
                                  matrixInput("DAPerrorB",
                                              value = matrix(0, 7, 8, dimnames = list(c("FS≥ 35x35cm", "30x30≤ FS< 35x35cm", "25x25≤ FS< 30x30cm", "20x20≤ FS< 25x25cm", "15x15≤ FS< 20x20cm", "10x10≤ FS< 15x15cm", "FS< 10x10cm"),
                                                                                      c("kV<60", "60≤kV<70", "70≤kV<80", "80≤kV<90", "90≤kV <100", "100≤kV <110", "110≤kV <120", "kV≥120"))),
                                              rows = list(names = TRUE),
                                              cols = list(names = TRUE),
                                              cells = list(editableCells = TRUE))),
                             card(plotlyOutput("DAPerrorA")),
                             card(plotlyOutput("DAPerrorB")))),
                  tabPanel("Mathematics/Physics/Validation",
                           HTML("Use the arrow keys to navigate the below presentation. This works best in Mozilla Firefox.<br>"),
                           htmlOutput("presentation"))))))


server <- function(input, output, session) {
  options(shiny.maxRequestSize=20*1024^2)
  
  output$presentation <- renderUI({
    tags$iframe(seamless="seamless", src = "presentation.html", width="100%", height=850, frameborder = "0")
  })
  
  RDSR <- eventReactive(input$RDSR_dcm,
                        read_RDSR_irradiation_events(input$RDSR_dcm$datapath))
  
  study_details <- eventReactive(input$RDSR_dcm,
                                 read_RDSR_acc_dose_data(input$RDSR_dcm$datapath))
  
  patient_phantom_df <- reactive({ifelse(input$input_phantom == "New skin dose estimation",
                                         get(input$patient_phantom),
                                         readRDS(input$patient_phantom$datapath))})
  
  system_config <- eventReactive(input$import_config,
                                 readRDS(input$import_config$datapath))
  
  observeEvent(system_config(), {
    updateNumericInput(session, "TL", value = system_config()[["TL"]])
    updateNumericInput(session, "table_width", value = system_config()[["table_width"]])
    updateNumericInput(session, "Table_Lateral_home_position", value = system_config()[["Table_Lateral_home_position"]])
    updateNumericInput(session, "Table_Longitudinal_home_position", value = system_config()[["Table_Longitudinal_home_position"]])
    updateNumericInput(session, "Table_height_home_position", value = system_config()[["Table_height_home_position"]])
    updateMatrixInput(session, "DAPerrorA", value = system_config()[["DAPerrorA"]])
    updateMatrixInput(session, "DAPerrorB", value = system_config()[["DAPerrorB"]])
    updateNumericInput(session, "hs_length", value = system_config()[["hs_length"]])
    updateNumericInput(session, "hs_width", value = system_config()[["hs_width"]])
  })
  
  previous_summary <- reactive(drop_na(select(patient_phantom_df()[[1]], 
                                              Study_Date,
                                              Study_Description,
                                              Patient_Name,
                                              Patient_ID,
                                              Gender,
                                              Birth_Date,
                                              Performing_Physician,
                                              Referring_Physician,
                                              Manufacturer,
                                              Model,
                                              Serial_Number,
                                              Fluoro_Dose_RP_Total,
                                              Fluoro_Dose_Area_Product_Total,
                                              Total_Fluoro_Time,
                                              Acquisition_Dose_RP_Total,
                                              Acquisition_Dose_Area_Product_Total,
                                              Total_Acquisition_Time,
                                              Height_Of_System), Patient_ID))
  
  output$no_of_IE <- reactive(paste("of", nrow(d()))) # number of irradiation events
  
  summary <- reactive({
    if (nrow(study_details()) > 1) { # if there biPlane systems return a dataframe with duplicate rows of study details because of the 2 planes. This block of code consolidates this into a single row.
      study_details_consolidated <- data.frame("Study_Date" = 1)
      study_details_consolidated$Study_Date <- study_details()$Study_Date[1]
      study_details_consolidated$Study_Description <- study_details()$Study_Description[1]
      study_details_consolidated$Patient_Name <- study_details()$Patient_Name[1]
      study_details_consolidated$Patient_ID <- study_details()$Patient_ID[1]
      study_details_consolidated$Gender <- study_details()$Gender[1]
      study_details_consolidated$Birth_Date <- study_details()$Birth_Date[1]
      study_details_consolidated$Performing_Physician <- study_details()$Performing_Physician[1]
      study_details_consolidated$Referring_Physician <- study_details()$Referring_Physician[1]
      study_details_consolidated$Manufacturer <- study_details()$Manufacturer[1]
      study_details_consolidated$Model <- study_details()$Model[1]
      study_details_consolidated$Serial_Number <- study_details()$Serial_Number[1]
      study_details_consolidated$Fluoro_Dose_RP_Total <- sum(study_details()$Fluoro_Dose_RP_Total)
      study_details_consolidated$Fluoro_Dose_Area_Product_Total <- sum(study_details()$Fluoro_Dose_Area_Product_Total)
      study_details_consolidated$Total_Fluoro_Time <- sum(study_details()$Total_Fluoro_Time)
      study_details_consolidated$Acquisition_Dose_RP_Total <- sum(study_details()$Acquisition_Dose_RP_Total)
      study_details_consolidated$Acquisition_Dose_Area_Product_Total <- sum(study_details()$Acquisition_Dose_Area_Product_Total)
      study_details_consolidated$Total_Acquisition_Time <- sum(study_details()$Total_Acquisition_Time)
      study_details_consolidated$Height_Of_System <- study_details()$Height_Of_System[1]
      summary <- rbind(study_details_consolidated, previous_summary())
    } else {
      summary <- rbind(study_details(), previous_summary())
    } 
    # summary$Study_Date <- format(ymd(summary$Study_Date), "%d/%m/%Y")
    # summary$Birth_Date <- format(ymd(summary$Birth_Date), "%d/%m/%Y")
    return(summary)
  })
  
  summary_display <- eventReactive(d(), {
    summary_display <- summary()
    summary_display$Study_Date <- format(ymd(summary_display$Study_Date), "%d/%m/%Y")
    summary_display$Birth_Date <- format(ymd(summary_display$Birth_Date), "%d/%m/%Y")
    summary_display <- t(summary_display)
    rownames(summary_display) <- c("Study Date",
                                   "Study Description",
                                   "Patient Name",
                                   "Patient ID",
                                   "Gender",
                                   "Birth Date",
                                   "Performing Physician",
                                   "Referring Physician",
                                   "Manufacturer",
                                   "Model",
                                   "Serial Number",
                                   "Total Fluoro Dose (Gy)",
                                   "Total Fluoro DAP (Gy.cm2)",
                                   "Total Fluoro Time (s)",
                                   "Total Acquisition Dose (Gy)",
                                   "Total Acquisition DAP (Gy.cm2)",
                                   "Total Acquisition Time (s)",
                                   "Height of System (private Philips tag)")
    column_names <- NA
    for(i in 1:ncol(summary_display)) {
      column_names[i] <- paste("Procedure", ncol(summary_display)+1-i)
    }
    colnames(summary_display) <- column_names
    return(summary_display)
  })
  
  output$summary = DT::renderDataTable(summary_display(), options = list(pageLength = 18, lengthChange = FALSE, bFilter = 0, bInfo=0, paging=FALSE))
  
  patient_phantom <- reactive(select(patient_phantom_df()[[1]], x,	y,	z,	i,	j,	k,	cumulative_skin_dose))
  
  d <- eventReactive(input$calculate_PSD, do_math(RDSR(),
                                                  study_details()$Manufacturer[1],
                                                  study_details()$Height_Of_System[1],
                                                  input$Table_Lateral_home_position,
                                                  input$Table_Longitudinal_home_position,
                                                  input$Table_height_home_position,
                                                  input$dh,
                                                  input$HR,
                                                  input$TL,
                                                  input$table_width,
                                                  input$DAPerrorA,
                                                  input$DAPerrorB,
                                                  input$hs_length,
                                                  input$hs_width))
  
  cumulative_patient_phantom <- eventReactive(d(), calculate_IE(d(), patient_phantom()))
  
  output$cumulative <- renderPlotly(plot_ly(cumulative_patient_phantom(),
                                            name = "Patient phantom",
                                            x = ~x,
                                            y = ~y,
                                            z = ~z,
                                            i = ~i,
                                            j = ~j,
                                            k = ~k,
                                            type = "mesh3d",
                                            intensity = ~cumulative_skin_dose,
                                            showscale = TRUE,
                                            colorbar = list(title = "Gy", len = 1),
                                            text = cumulative_patient_phantom()$cumulative_skin_dose,
                                            hovertemplate = 'Skin Dose: %{text:.2f} Gy') %>% 
                                      add_annotations(xref = "paper",
                                                      yref = "paper",
                                                      zref = "paper",
                                                      x = 1,
                                                      y = 1,
                                                      text = paste("Peak Skin Dose:", signif(max(cumulative_patient_phantom()$cumulative_skin_dose), 3), "Gy"),showarrow = F) %>% 
                                      layout(paper_bgcolor='transparent',
                                             scene = list(xaxis = list(title = "Lateral (x)", range = c(-1.5, 1), tickvals = c(-1.5 -1, -0.5, 0, 0.5, 1)),
                                                          yaxis = list(title = "Longitudinal (y)", range =c(1, -1), tickvals = c(1, 0.5, 0, -0.5, -1)),
                                                          zaxis = list(title = "Height (z)", range =c(1, -1), tickvals = c(1, 0.5, 0, -0.5, -1)),
                                                          camera = list(eye = list(x = 0.1, y = 0, z = -2)),
                                                          aspectratio = list(x =2.5, y = 2, z = 2))))
  
  IE_patient_phantom <- reactive(calculate_IE(d()[input$IE ,], patient_phantom()))
  
  output$DAPerrorA <- renderPlotly(plot_ly(x = c("kV<60", "60≤kV<70", "70≤kV<80", "80≤kV<90", "90≤kV<100", "100≤kV<110", "110≤kV<120", "kV≥120"),
                                           y = c("FS≥35x35cm", "30x30≤FS<35x35cm", "25x25≤FS<30x30cm", "20x20≤FS<25x25cm", "15x15≤FS<20x20cm", "10x10≤FS<15x15cm", "FS<10x10cm"),
                                           z = input$DAPerrorA,
                                           type = 'heatmap') %>%
                                     layout(paper_bgcolor='transparent',
                                            yaxis = list(autorange = "reversed")))
  
  output$DAPerrorB <- renderPlotly(plot_ly(x = c("kV<60", "60≤kV<70", "70≤kV<80", "80≤kV<90", "90≤kV<100", "100≤kV<110", "110≤kV<120", "kV≥120"),
                                           y = c("FS≥35x35cm", "30x30≤FS<35x35cm", "25x25≤FS<30x30cm", "20x20≤FS<25x25cm", "15x15≤FS<20x20cm", "10x10≤FS<15x15cm", "FS<10x10cm"),
                                           z = input$DAPerrorB,
                                           type = 'heatmap') %>%
                                     layout(paper_bgcolor='transparent',
                                            yaxis = list(autorange = "reversed")))
  
  output$visualisation <-renderPlotly(plot_ly(data = d()) %>%
                                        add_trace(name = "Beam line 1",
                                                  x = c(~x_spot[input$IE], ~x_dc1[input$IE]),
                                                  y = c(~y_spot[input$IE], ~y_dc1[input$IE]),
                                                  z = c(~z_spot[input$IE], ~z_dc1[input$IE]),
                                                  type = "scatter3d",
                                                  mode = "lines",
                                                  line = list(color = "red", width = 3),
                                                  hoverinfo = "none") %>% 
                                        add_trace(name = "Beam line 2",
                                                  x = c(~x_spot[input$IE], ~x_dc2[input$IE]),
                                                  y = c(~y_spot[input$IE], ~y_dc2[input$IE]),
                                                  z = c(~z_spot[input$IE], ~z_dc2[input$IE]),
                                                  type = "scatter3d",
                                                  mode = "lines",
                                                  line = list(color = "red", width = 3),
                                                  hoverinfo = "none") %>% 
                                        add_trace(name = "Beam line 3",
                                                  x = c(~x_spot[input$IE], ~x_dc3[input$IE]),
                                                  y = c(~y_spot[input$IE], ~y_dc3[input$IE]),
                                                  z = c(~z_spot[input$IE], ~z_dc3[input$IE]),
                                                  type = "scatter3d",
                                                  mode = "lines",
                                                  line = list(color = "red", width = 3),
                                                  hoverinfo = "none") %>% 
                                        add_trace(name = "Beam line 4",
                                                  x = c(~x_spot[input$IE], ~x_dc4[input$IE]),
                                                  y = c(~y_spot[input$IE], ~y_dc4[input$IE]),
                                                  z = c(~z_spot[input$IE], ~z_dc4[input$IE]),
                                                  type = "scatter3d",
                                                  mode = "lines",
                                                  line = list(color = "red", width = 3),
                                                  hoverinfo = "none") %>% 
                                        add_trace(name = "Patient Support (body)",
                                                  x = c(~x_table1[input$IE], ~x_table2[input$IE], ~x_table3[input$IE], ~x_table4[input$IE]),
                                                  y = c(~y_table1[input$IE], ~y_table2[input$IE], ~y_table3[input$IE], ~y_table4[input$IE]),
                                                  z = c(~z_table[input$IE], ~z_table[input$IE], ~z_table[input$IE], ~z_table[input$IE]),
                                                  type = "mesh3d",
                                                  opacity = 0.5,
                                                  flatshading = T,
                                                  facecolor = rep("black", 12),
                                                  hovertext = ~paste("Patient support (body)"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Patient Support (head)",
                                                  x = c(~x_hs1[input$IE], ~x_hs2[input$IE], ~x_hs3[input$IE], ~x_hs4[input$IE]),
                                                  y = c(~y_hs1[input$IE], ~y_hs2[input$IE], ~y_hs3[input$IE], ~y_hs4[input$IE]),
                                                  z = c(~z_table[input$IE], ~z_table[input$IE], ~z_table[input$IE], ~z_table[input$IE]),
                                                  type = "mesh3d",
                                                  opacity = 0.5,
                                                  flatshading = T,
                                                  facecolor = rep("black", 12),
                                                  hovertext = ~paste("Patient support (head)"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Focal Spot",
                                                  x = ~x_spot[input$IE],
                                                  y = ~y_spot[input$IE],
                                                  z = ~z_spot[input$IE],
                                                  type = "scatter3d",
                                                  mode = "markers",
                                                  marker = list(color = "red"),
                                                  hovertext = ~paste("Focal spot"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Isocentre",
                                                  x = ~x_Is[input$IE],
                                                  y = ~y_Is[input$IE],
                                                  z = ~z_Is[input$IE],
                                                  type = "scatter3d",
                                                  mode = "markers",
                                                  marker = list(color = "magenta"),
                                                  hovertext = ~paste("Isocentre"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Exclusion plane",
                                                  x = c(~x_HR1[input$IE], ~x_HR2[input$IE], ~x_HR3[input$IE], ~x_HR4[input$IE]),
                                                  y = c(~y_HR1[input$IE], ~y_HR2[input$IE], ~y_HR3[input$IE], ~y_HR4[input$IE]),
                                                  z = c(~z_HR1[input$IE], ~z_HR2[input$IE], ~z_HR3[input$IE], ~z_HR4[input$IE]),
                                                  type = "mesh3d",
                                                  facecolor = rep("blue", 12),
                                                  flatshading = T,
                                                  hovertext = ~paste("Exclusion plane"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Reference Point",
                                                  x = ~x_RF[input$IE],
                                                  y = ~y_RF[input$IE],
                                                  z = ~z_RF[input$IE],
                                                  type = "scatter3d",
                                                  mode = "markers",
                                                  marker = list(color = "orange"),
                                                  hovertext = ~paste("Reference Point"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Detector",
                                                  x = c(~x_dc1[input$IE], ~x_dc2[input$IE], ~x_dc3[input$IE], ~x_dc4[input$IE]),
                                                  y = c(~y_dc1[input$IE], ~y_dc2[input$IE], ~y_dc3[input$IE], ~y_dc4[input$IE]),
                                                  z = c(~z_dc1[input$IE], ~z_dc2[input$IE], ~z_dc3[input$IE], ~z_dc4[input$IE]),
                                                  type = "mesh3d",
                                                  facecolor = rep("green", 12),
                                                  flatshading = T,
                                                  hovertext = ~paste("Detector"),
                                                  hoverinfo = "text") %>% 
                                        add_trace(name = "Patient phantom",
                                                  x = ~ IE_patient_phantom()$x_adj,
                                                  y = ~ IE_patient_phantom()$y_adj,
                                                  z = ~ IE_patient_phantom()$z_adj,
                                                  i = ~ IE_patient_phantom()$i,
                                                  j = ~ IE_patient_phantom()$j,
                                                  k = ~ IE_patient_phantom()$k,
                                                  type = "mesh3d",
                                                  opacity = 0.5,
                                                  intensity = ~IE_patient_phantom()$irradiation_event_skin_dose*1000,
                                                  showscale = TRUE,
                                                  colorbar = list(title = "mGy", yanchor = "top", len = 1, yref = "container", y = 1),
                                                  hovertext = ~paste(" Skin irradiated:", IE_patient_phantom()$`is point irradiated?`, "<br>",
                                                                     "Table attenuation:", IE_patient_phantom()$`table attenuation`, "<br>",
                                                                     "<i>f</i> =", signif(d()$f[input$IE], 3), "<br>",
                                                                     "<i>BSF</i> =", signif(d()$BSF[input$IE],3 ), "<br>",
                                                                     "<i>TTF</i> =", signif(IE_patient_phantom()$TTF, 3), "<br>",
                                                                     "Angle of incident X-rays through table =", signif(IE_patient_phantom()$eff_angle * (180/pi), 3), "&deg;", "<br>",
                                                                     "<i>f<sub>&#952;</sub></i> =", signif(IE_patient_phantom()$f_theta, 3), "<br>",
                                                                     "Source to skin distance =", signif(IE_patient_phantom()$Dist, 3), "m", "<br>",
                                                                     "Kerma at source to skin distance (<i>K<sub>skin</sub></i>) =", signif(IE_patient_phantom()$K_skin, 3)*1000, "mGy", "<br>",
                                                                     "Correction factor (<i>CF</i>) =", signif(((100-d()$cal[input$IE])/100), 3), "<br>",
                                                                     "Skin dose = <i>f &#215; BSF &#215; TTF &#215; f<sub>&#952;</sub> &#215; K<sub>skin</sub> &#215; CF <i>=", signif(IE_patient_phantom()$irradiation_event_skin_dose, 3)*1000, "mGy", "<br>"),
                                                  hoverinfo = "text") %>% 
                                        add_annotations(xref = "paper",
                                                        yref = "paper",
                                                        zref = "paper",
                                                        x = 1,
                                                        y = 0,
                                                        text = paste0("Peak Skin Dose For Irradiation Event: ", signif(max(IE_patient_phantom()$irradiation_event_skin_dose)*1000, 3), " mGy", "<br>",
                                                                      "Dose at RP: ", signif(d()$Dose_RP[input$IE]*1000, 3), " mGy", "<br>",
                                                                      "Dose Area Product: ", signif(d()$Dose_Area_Product[input$IE]*10000, 3), " Gy.cm<sup>2</sup>", "<br>",
                                                                      "Tube Voltage: ", signif(d()$KVP[input$IE], 3),  " kV", "<br>",
                                                                      "Field size at RP: ", signif(d()$FS_RP[input$IE], 3), " m<sup>2</sup>", "<br>",
                                                                      "Exposure: ", d()$Exposure[input$IE]/1000, " mAs", "<br>",
                                                                      "Average Pulse Width: ", d()$Pulse_Width[input$IE], " ms", "<br>",
                                                                      "Focal Spot Size: ", d()$Focal_Spot_Size[input$IE], "<br>",
                                                                      "Pulse Rate: ", d()$Pulse_Rate[input$IE], " per s", "<br>",
                                                                      "X-Ray Filter Thickness: ", d()$`X-Ray_Filter_Thickness_Minimum`[input$IE], " mm", "<br>",
                                                                      "X-Ray Filter Material: ", d()$`X-Ray_Filter_Material`[input$IE], "<br>",
                                                                      "Irradiation Event Type: ", d()$Irradiation_Event_Type[input$IE], "<br>",
                                                                      "Acquisition Protocol: ", d()$Acquisition_Protocol[input$IE], "<br>",
                                                                      "Positioner Primary Angle: ", d()$Positioner_Primary_Angle[input$IE], "&deg;", "<br>",
                                                                      "Positioner Secondary Angle: ", d()$Positioner_Secondary_Angle[input$IE], "&deg;", "<br>",
                                                                      "Table Lateral Position: ", d()$Table_Lateral_Position[input$IE], " mm", "<br>",
                                                                      "Table Longitudinal Position: ", d()$Table_Longitudinal_Position[input$IE], " mm", "<br>",
                                                                      "Table Height Position: ", d()$Table_Height_Position[input$IE], " mm", "<br>",
                                                                      "Distance Source to Detector: ", d()$Distance_Source_to_Detector[input$IE]*1000, " mm", "<br>",
                                                                      "Distance Source to Isocenter: ", d()$radius[input$IE]*1000, " mm", "<br>",
                                                                      "Acquisition Plane: ", d()$Acquisition_Plane_in_Irradiation_Event[input$IE], "<br>",
                                                                      "Irradiation Event UID: ", d()$Irradiation_Event_UID[input$IE], "<br>"),
                                                        align = "right",
                                                        showarrow = F) %>% 
                                        layout(paper_bgcolor='transparent',
                                               showlegend = F, 
                                               scene = list(xaxis = list(title = "Lateral (x)", 
                                                                         range = c(-2, 2), 
                                                                         tickvals = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)),
                                                            yaxis = list(title = "Longitudinal (y)", 
                                                                         range = c(1, -1), 
                                                                         tickvals = c(1, 0.5, 0, -0.5, -1)),
                                                            zaxis = list(title = "Height (z)", 
                                                                         range = c(1, -1), 
                                                                         tickvals = c(1, 0.5, 0, -0.5, -1)),
                                                            camera = list(eye = list(x = -0.75, y = 1, z = 0.75),
                                                                          center = list(x = -0.25, y = 0, z = 0)),
                                                            aspectratio = list(x = 2, y = 1, z = 1))))
  
  summary_figure_data <- reactive(data.frame("Value" = c(sum(summary()$Fluoro_Dose_RP_Total), 
                                                         sum(summary()$Acquisition_Dose_RP_Total),
                                                         sum(summary()$Fluoro_Dose_Area_Product_Total), 
                                                         sum(summary()$Acquisition_Dose_Area_Product_Total), 
                                                         sum(summary()$Total_Fluoro_Time), 
                                                         sum(summary()$Total_Acquisition_Time)),
                                             "Metric" = c("Dose at Reference Point (Gy)","Dose at Reference Point (Gy)",
                                                          "Dose Area Product (Gy.cm2)", "Dose Area Product (Gy.cm2)",
                                                          "Screening Time (s)", "Screening Time (s)"),
                                             "Irradiation_Type" = c("Fluoroscopy", "Acquisition",
                                                                    "Fluoroscopy", "Acquisition",
                                                                    "Fluoroscopy", "Acquisition")))
  
  output$summary_figure <- renderPlotly(ggplot(summary_figure_data(), aes(x = Irradiation_Type, y = Value, fill = Irradiation_Type, group = Metric)) +
                                          geom_bar(stat = 'identity') + 
                                          facet_wrap(~Metric, ncol = 1, scales = 'free') +
                                          theme(legend.position="none",
                                                axis.title.y = element_blank(),
                                                axis.title.x = element_blank()))
  
  output$phantom <- downloadHandler(
    filename = "skin_dose_data.rds",
    content = function(file) {
      patient_phantom <- list(cbindX(cumulative_patient_phantom(), summary()))
      saveRDS(patient_phantom, file = file)})
  
  output$export_config <- downloadHandler(
    filename = "system_config.rds",
    content = function(file) {
      system_config_export <- list("DAPerrorA" = input$DAPerrorA,
                                   "DAPerrorB" = input$DAPerrorB,
                                   "Table_Lateral_home_position" = input$Table_Lateral_home_position,
                                   "Table_Longitudinal_home_position" = input$Table_Longitudinal_home_position,
                                   "Table_height_home_position" = input$Table_height_home_position,
                                   "TL" = input$TL,
                                   "table_width" = input$table_width,
                                   "hs_length" = input$hs_length,
                                   "hs_width" = input$hs_width)
      saveRDS(system_config_export, file = file)})
  
  output$report <- downloadHandler(filename = "Peak Skin Dose Estimate Report.docx",
                                   content = function(file) {
                                     rmarkdown::render(
                                       "report.Rmd",
                                       params = list(
                                         Patient_ID = unique(summary()$Patient_ID),
                                         DOB = unique(summary()$`Birth_Date`),
                                         Patient_Name = unique(summary()$Patient_Name),
                                         Gender = unique(summary()$Gender),
                                         phantom = cumulative_patient_phantom(),
                                         summary = summary(),
                                         skin_reactions = skin_reactions),
                                       output_file = file,
                                       envir = new.env(parent = globalenv()))})
}
shinyApp(ui, server)
