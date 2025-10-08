library(plotly)
library(tidyverse)
phantoms <- list()
cow <- read_csv("data/cow.csv")
cow$x <- cow$x-2.762274+1
cow$z <- cow$z-0.9000000
phantoms$cow <- cow

library(readobj)
f <- "data/MRCP_AM_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_AM.148"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_AM.148"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.5
positions$z <- min(positions$z)+positions$z
positions[41003:82000,] <- NA
`ICRP145 Adult Male` <- cbind(positions, indicies)
phantoms$`ICRP145 Adult Male` <- `ICRP145 Adult Male`

f <- "data/MRCP_AF_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_AF.150"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_AF.150"]][["positions"]]))

indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.5
positions$z <- min(positions$z)+positions$z
positions[38486:76966,] <- NA

`ICRP145 Adult Female` <- cbind(positions, indicies)
phantoms$`ICRP145 Adult Female` <- `ICRP145 Adult Female`

f <- "data/MRCP_00F_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_00F.166"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_00F.166"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.2
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 00F` <- cbind(positions, indicies)
phantoms$`ICRP156 00F` <- `ICRP156 00F`

f <- "data/MRCP_00M_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_00M.165"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_00M.165"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.2
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 00M` <- cbind(positions, indicies)
phantoms$`ICRP156 00M` <- `ICRP156 00M`

f <- "data/MRCP_01F_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_01F.164"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_01F.164"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.25
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 01F` <- cbind(positions, indicies)
phantoms$`ICRP156 01F` <- `ICRP156 01F`

f <- "data/MRCP_01M_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_01M.163"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_01M.163"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.25
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 01M` <- cbind(positions, indicies)
phantoms$`ICRP156 01M` <- `ICRP156 01M`

f <- "data/MRCP_05F_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_05F.162"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_05F.162"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.3
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 05F` <- cbind(positions, indicies)
phantoms$`ICRP156 05F` <- `ICRP156 05F`

f <- "data/MRCP_05M_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_05M.161"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_05M.161"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.3
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 05M` <- cbind(positions, indicies)
phantoms$`ICRP156 05M` <- `ICRP156 05M`

f <- "data/MRCP_10F_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_10F.163"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_10F.163"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.35
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 10F` <- cbind(positions, indicies)
phantoms$`ICRP156 10F` <- `ICRP156 10F`

f <- "data/MRCP_10M_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_10M.161"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_10M.161"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.35
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 10M` <- cbind(positions, indicies)
phantoms$`ICRP156 10M` <- `ICRP156 10M`

f <- "data/MRCP_15F_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_15F.162"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_15F.162"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.4
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 15F` <- cbind(positions, indicies)
phantoms$`ICRP156 15F` <- `ICRP156 15F`

f <- "data/MRCP_15M_skin_surface.obj"
obj <- read.obj(f, materialspath = NULL, convert.rgl = FALSE, triangulate = TRUE)
indicies <- data.frame()
indicies <- data.frame(t(obj[["shapes"]][["MRCP_15M.162"]][["indices"]]))
positions <- data.frame(t(obj[["shapes"]][["MRCP_15M.162"]][["positions"]]))
indicies <- data.frame(i = indicies$X1,
                       j = indicies$X2,
                       k = indicies$X3)
positions <- data.frame(x = positions$X3/100,
                        y = positions$X1/100,
                        z = positions$X2/100)
positions$x <- min(positions$x)+positions$x+0.4
positions$z <- min(positions$z)+positions$z
positions[(nrow(positions)+1):nrow(indicies),] <- NA

`ICRP156 15M` <- cbind(positions, indicies)
phantoms$`ICRP156 15M` <- `ICRP156 15M`

for(i in 1:length(phantoms)){
  phantoms[[i]]$cumulative_skin_dose <- 0
  phantoms[[i]]$Study_Date <- NA
  phantoms[[i]]$Study_Description <- NA
  phantoms[[i]]$Patient_Name <- NA
  phantoms[[i]]$Patient_ID <- NA
  phantoms[[i]]$Gender <- NA
  phantoms[[i]]$Birth_Date <- NA
  phantoms[[i]]$Performing_Physician <- NA
  phantoms[[i]]$Referring_Physician <- NA
  phantoms[[i]]$Manufacturer <- NA
  phantoms[[i]]$Model <- NA
  phantoms[[i]]$Serial_Number <- NA
  phantoms[[i]]$Fluoro_Dose_RP_Total <- NA
  phantoms[[i]]$Fluoro_Dose_Area_Product_Total <- NA
  phantoms[[i]]$Total_Fluoro_Time <- NA
  phantoms[[i]]$Acquisition_Dose_RP_Total <- NA
  phantoms[[i]]$Acquisition_Dose_Area_Product_Total <- NA
  phantoms[[i]]$Total_Acquisition_Time <- NA
  phantoms[[i]]$Height_Of_System <- NA
}

`cow` <- list(phantom = phantoms$`cow`)
`ICRP145 Adult Male` <- list(phantom = phantoms$`ICRP145 Adult Male`)
`ICRP145 Adult Female` <- list(phantom = phantoms$`ICRP145 Adult Female`)
`ICRP156 00F` <- list(phantom = phantoms$`ICRP156 00F`)
`ICRP156 00M` <- list(phantom = phantoms$`ICRP156 00M`)
`ICRP156 01F` <- list(phantom = phantoms$`ICRP156 01F`)
`ICRP156 01M` <- list(phantom = phantoms$`ICRP156 01M`)
`ICRP156 05F` <- list(phantom = phantoms$`ICRP156 05F`)
`ICRP156 05M` <- list(phantom = phantoms$`ICRP156 05M`)
`ICRP156 10F` <- list(phantom = phantoms$`ICRP156 10F`)
`ICRP156 10M` <- list(phantom = phantoms$`ICRP156 10M`)
`ICRP156 15F` <- list(phantom = phantoms$`ICRP156 15F`)
`ICRP156 15M` <- list(phantom = phantoms$`ICRP156 15M`)

plot_ly(phantoms[[2]],
        x = ~x,
        y = ~y,
        z = ~z,
        i = ~i,
        j = ~j,
        k = ~k,
        type = "mesh3d") %>%
  layout(scene = list(xaxis = list(title = "Lateral (x)", range = c(-2, 2)),
                      yaxis = list(title = "Longitudinal (y)", range =c(1, -1)),
                      zaxis = list(title = "Height (z)", range = c(1, -1)),
                      aspectratio = list(x =2, y = 1, z = 1)))

rm(phantom, phantoms, i, shift_y, shift_z)

skin_reactions <- read_csv("data/skin_reactions.csv")
BSF_data <- read_csv("data/BSF and f.csv")

# Calculate BSF for an average inherent filtration
BSF_data_average <- data.frame()
field_sizes <- unique(BSF_data$field_size_m2)
for(i in 1:length(field_sizes)) {
  field_size_data <- subset(BSF_data, field_size_m2 == field_sizes[i])
  kVps <- unique(field_size_data$kVp)
  for(j in 1:length(kVps)) {
    kVp_data <- subset(field_size_data, kVp == kVps[j])
    mm_Cus <- unique(BSF_data$mm_Cu)
    for(k in 1:length(mm_Cus)) {
      filtered_data <- subset(kVp_data, mm_Cu == mm_Cus[k])

      newrow <- data.frame(field_size_m2 = c(unique(BSF_data$field_size_m2)[i]),
                           kVp = c(unique(BSF_data$kVp)[j]),
                           mm_Cu = c(unique(BSF_data$mm_Cu)[k]),
                           BSF = mean(filtered_data$BSF),
                           f = mean(filtered_data$f))
      BSF_data_average <- rbind(BSF_data_average, newrow)
    }
  }
}
BSF_data <- BSF_data_average
# BSF_data <- rbind(BSF_data_average, BSF_data)
TTF_data <- read_csv("data/TTF.csv")
rm(newrow, indicies, positions, obj, filtered_data, BSF_data_average, field_size_data, kVp_data, f, i, j, k, kVps, mm_Cus, field_sizes)
save.image("environment.RData")
rm(list = ls())
