
library(dplyr)
library(ggplot2)

## read in data from Harmony 
#single cell data
nuclei_table <- read.delim("Objects_Population - Nuclei_within_gastruloids.txt", header=TRUE, sep="\t")
#gastruloid (colony) data
gastruloid_table <- read.delim("Objects_Population - gastruloids.txt", header=TRUE, sep="\t")

#Adjust WellID column
nuclei_table_2 <- nuclei_table
nuclei_table_2[1] <- lapply(nuclei_table_2[1], chartr, old = '234567', new = 'BCDEFG')
nuclei_table_2[1] <- paste(nuclei_table_2[,1], nuclei_table_2[,2], sep = '')
nuclei_table_3 <- nuclei_table_2[,-2] 
colnames(nuclei_table_3)[which(names(nuclei_table_3) == "Row")] <- "Well.ID"

gastruloid_table_2 <- gastruloid_table
gastruloid_table_2[1] <- lapply(gastruloid_table_2[1], chartr, old = '234567', new = 'BCDEFG')
gastruloid_table_2[1] <- paste(gastruloid_table_2[,1], gastruloid_table_2[,2], sep = '')
gastruloid_table_3 <- gastruloid_table_2[,-2]
colnames(gastruloid_table_3)[which(names(gastruloid_table_3) == "Row")] <- "Well.ID"

#select columns of interest
gastruloid_keep <- c("Well.ID", "Object.No", "Position.X..µm.", "Position.Y..µm.")
nuclei_keep <- c("Well.ID", "Object.No", "Nuclei_within_gastruloids...ROI.No", "Nuclei_within_gastruloids...Intensity.Nucleus.Alexa.488..global..Mean", "Nuclei_within_gastruloids...Intensity.Nucleus.DAPI..global..Mean",
                 "Position.X..µm.", "Position.Y..µm.")

#remove redundant columns
gastruloid_table_4 <- gastruloid_table_3[gastruloid_keep] 
nuclei_table_4 <- nuclei_table_3[nuclei_keep] 

#rename columns - gastruloid table
colnames(gastruloid_table_4)[which(names(gastruloid_table_4) == "Object.No")] <- "Gastruloid.No"
colnames(gastruloid_table_4)[which(names(gastruloid_table_4) == "Position.X..µm.")] <- "Gastruloid.X"
colnames(gastruloid_table_4)[which(names(gastruloid_table_4) == "Position.Y..µm.")] <- "Gastruloid.Y"

#rename columns - nuclei table
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Object.No")] <- "Nuclei.No"
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Nuclei_within_gastruloids...ROI.No")] <- "Gastruloid.No"
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Nuclei_within_gastruloids...Intensity.Nucleus.Alexa.488..global..Mean")] <- "Alexa.488.Mean"
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Nuclei_within_gastruloids...Intensity.Nucleus.DAPI..global..Mean")] <- "DAPI.Mean"
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Position.X..µm.")] <- "Nuclei.X"
colnames(nuclei_table_4)[which(names(nuclei_table_4) == "Position.Y..µm.")] <- "Nuclei.Y"

#merge gastruloid and nuclei table
merged_file1 <- merge(x = gastruloid_table_4, y = nuclei_table_4, by.x = c("Well.ID", "Gastruloid.No"), by.y =  c("Well.ID","Gastruloid.No"))

#calculate distance Nuclei-Gastruloids
merged_file1[,"Nuclei.Distance"] <- sqrt((merged_file1[,"Nuclei.X"] - merged_file1[,"Gastruloid.X"])^2 + (merged_file1[,"Nuclei.Y"] - merged_file1[,"Gastruloid.Y"])^2) 

#calculate bins
x <- max(merged_file1$Nuclei.Distance)
merged_file1$bin <- ifelse(merged_file1$Nuclei.Distance < 25, 0, 
                           ifelse(merged_file1$Nuclei.Distance < 50, 25, 
                                ifelse(merged_file1$Nuclei.Distance < 75, 50, 
                                     ifelse(merged_file1$Nuclei.Distance < 100, 75, 
                                          ifelse(merged_file1$Nuclei.Distance < 125, 100, 
                                               ifelse(merged_file1$Nuclei.Distance < 150, 125, 
                                                    ifelse(merged_file1$Nuclei.Distance < 175, 150, 
                                                         ifelse(merged_file1$Nuclei.Distance < 200, 175, 
                                                              ifelse(merged_file1$Nuclei.Distance < 225, 200,
                                                                   ifelse(merged_file1$Nuclei.Distance < 250, 225,
                                                                        ifelse(merged_file1$Nuclei.Distance < 275, 250,
                                                                             ifelse(merged_file1$Nuclei.Distance < 300, 275,
                                                                                  ifelse(merged_file1$Nuclei.Distance < 325, 300,
                                                                                       ifelse(merged_file1$Nuclei.Distance < 350, 325,
                                                                                            ifelse(merged_file1$Nuclei.Distance < 375, 350,
                                                                                                 ifelse(merged_file1$Nuclei.Distance < 400, 375,
                                                                                                      ifelse(merged_file1$Nuclei.Distance < 425, 400,
                                                                                                           ifelse(merged_file1$Nuclei.Distance < 450, 425,
                                                                                                                ifelse(merged_file1$Nuclei.Distance < 475, 450,
                                                                                                                     ifelse(merged_file1$Nuclei.Distance < 500, 475,
                                                                                                                          ifelse(merged_file1$Nuclei.Distance < 525, 500,
                                                                                                                               ifelse(merged_file1$Nuclei.Distance < 550, 525,
                                                                                                                                    ifelse(merged_file1$Nuclei.Distance <= x, 550)))))))))))))))))))))))
#remove outer bin
merged_file <- subset(merged_file1, !bin %in% c(525,550))

#add column to dataframe to indicate protein marker (SOX2, BRA or SOX17), line number (i.e. position in plate) (1, 2 or 3) and biological replicate
SOX2_1_n1 <- c("B2", "B3", "B4")
SOX17_1_n1 <- c("B5", "B6", "B7")
BRA_1_n1 <- c("B8", "B9", "B10")

SOX2_2_n1 <- c("D2", "D3", "D4")
SOX17_2_n1 <- c("D5", "D6", "D7")
BRA_2_n1 <- c("D8", "D9", "D10")

SOX2_3_n1 <- c("F2", "F3", "F4")
SOX17_3_n1 <- c("F5", "F6", "F7")
BRA_3_n1 <- c("F8", "F9", "F10")

merged_file$marker_line <- ifelse(merged_file$Well.ID %in% SOX2_1_n1, "SOX2_1_n1", 
                                  ifelse(merged_file$Well.ID %in% BRA_1_n1, "BRA_1_n1",
                                         ifelse(merged_file$Well.ID %in% SOX17_1_n1, "SOX17_1_n1",
                                                ifelse(merged_file$Well.ID %in% SOX2_2_n1, "SOX2_2_n1",
                                                       ifelse(merged_file$Well.ID %in% BRA_2_n1, "BRA_2_n1",
                                                              ifelse(merged_file$Well.ID %in% SOX17_2_n1, "SOX17_2_n1",
                                                                     ifelse(merged_file$Well.ID %in% SOX2_3_n1, "SOX2_3_1",
                                                                            ifelse(merged_file$Well.ID %in% BRA_3_n1, "BRA_3_n1",
                                                                                   ifelse(merged_file$Well.ID %in% SOX17_3_n1, "SOX17_3_n1", "negative_control_n1")))))))))


## HISTOGRAMS OF ALEXA488 INTENSITY (SINGLE CELL DATA)
histograms2 <- "insert_folder_location"

histogram_plots <- function(merged_file, na.rm = TRUE, ...){
  
  Well.ID_list <- unique(merged_file$Well.ID)
  
  for (i in seq_along(Well.ID_list)) { 
    
    histo <- 
      ggplot(subset(merged_file, merged_file$Well.ID==Well.ID_list[i]),
             aes(Alexa.488.Mean)) + 
      geom_histogram(color="black", fill="white", bins=1000) + 
      ggtitle(paste(Well.ID_list[i]), "") +
      theme_bw()
    
    ggsave(histo, file=paste(histograms2, Well.ID_list[i], ".png", sep=''), width = 50, height = 20, units = "cm")

  }
}

histogram_plots(merged_file)

#split merged_file by well 
split_merged_file <- split(merged_file, merged_file$Well.ID)
wells <- names(split_merged_file)
for(i in seq_along(wells)){
  assign(wells[i], split_merged_file[[i]])
}

#remove background for each well based on histogram plots
#LINE 1 
#SOX2 
B2_expression <- B2 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 7500, 0)) 

B3_expression <- B3 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 8000, 0)) 

B4_expression <- B4 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 10000, 0)) 

SOX2_1_collated <- rbind(B2_expression, B3_expression, B4_expression)

#SOX17
B5_expression <- B5 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 32000, 0)) 

B6_expression <- B6 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 32000, 0)) 

B7_expression <- B7 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 22000, 0)) 

SOX17_1_collated <- rbind(B5_expression, B6_expression, B7_expression)

#BRA
B8_expression <- B8 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 30000, 0)) 

B9_expression <- B9 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 30000, 0)) 

B10_expression <- B10 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 25000, 0)) 

BRA_1_collated <- rbind(B8_expression, B9_expression, B10_expression)

#LINE 2 
#SOX2 
D2_expression <- D2 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 8000, 0)) 

D3_expression <- D3 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 8000, 0)) 

D4_expression <- D4 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 7500, 0)) 

SOX2_2_collated <- rbind(D2_expression, D3_expression, D4_expression)

#SOX17
D5_expression <- D5 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 8000, 0)) 

D6_expression <- D6 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 8000, 0)) 

D7_expression <- D7 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 12000, 0)) 

SOX17_2_collated <- rbind(D5_expression, D6_expression, D7_expression)

#BRA
D8_expression <- D8 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 30000, 0)) 

D9_expression <- D9 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 35000, 0)) 

D10_expression <- D10 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 30000, 0)) 

BRA_2_collated <- rbind(D8_expression, D9_expression, D10_expression)

#LINE 3 
#SOX2 
F2_expression <- F2 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 10000, 0)) 

F3_expression <- F3 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 10000, 0)) 

F4_expression <- F4 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 10000, 0)) 

SOX2_3_collated <- rbind(F2_expression, F3_expression, F4_expression)

#SOX17
F5_expression <- F5 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 6000, 0)) 

F6_expression <- F6 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 5000, 0)) 

F7_expression <- F7 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 6000, 0)) 

SOX17_3_collated <- rbind(F5_expression, F6_expression, F7_expression)

#BRA
F8_expression <- F8 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 10000, 0)) 

F9_expression <- F9 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 15000, 0)) 

F10_expression <- F10 %>%
  mutate(background_removed = replace(Alexa.488.Mean, Alexa.488.Mean < 12000, 0)) 

BRA_3_collated <- rbind(F8_expression, F9_expression, F10_expression)
write.table(BRA_3_collated, file = "BRA_3_data_n1.txt", sep = "\t")

#combine SOX2, BRA + SOX17 datasets
all_data <- rbind(SOX2_1_collated, BRA_1_collated, SOX17_1_collated, SOX2_2_collated, BRA_2_collated, SOX17_2_collated, SOX2_3_collated, BRA_3_collated, SOX17_3_collated)

#normalise Alexa 488 mean to DAPI mean 
all_data[,"Normalised.Alexa488"] <- all_data[,"background_removed"]/all_data[,"DAPI.Mean"]

#find mean per bin per gastruloid per well 
mean_data <- all_data %>%
  group_by(Well.ID, bin, Gastruloid.No, marker_line) %>%
  summarise(Alexa488_normalised = mean(Normalised.Alexa488))

#normalise to max Alexa 488 intensity within each well
normalised_data <- mean_data %>%
  group_by(Well.ID) %>%
  mutate(mean_intensity = (Alexa488_normalised/(max(Alexa488_normalised))))

#find mean per bin per well (averaged across gastruloids)
well_means <- normalised_data %>%
  group_by(Well.ID, bin, marker_line) %>%
  summarise(well_mean = mean(mean_intensity),
            sd_intensity = sd(mean_intensity))

#find mean per bin across triplicate wells 
marker_means <- normalised_data %>%
  group_by(marker_line, bin) %>%
  summarise(marker_mean = mean(mean_intensity),
            sd_intensity = sd(mean_intensity))

write.table(marker_means, file = "marker_means_n1.txt", sep = "\t")

## Plots of individual colonies (grouped by well)
gastruloids <- "insert_folder_location"

graphs <- function(normalised_data, na.rm = TRUE, ...){
  
  Well.ID_list <- unique(normalised_data$Well.ID)
  
  for (i in seq_along(Well.ID_list)) { 
    
    plot <- 
      ggplot(subset(normalised_data, normalised_data$Well.ID==Well.ID_list[i]),
             aes(bin, mean_intensity, group=Gastruloid.No, color=Gastruloid.No)) + 
      geom_line(size=1.5) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"), axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold")) + 
      labs(x = "Distance from colony centre (μM)", y = "Fluorescence intensity (a.u.)")
    ggtitle(paste(Well.ID_list[i]), "")
    
    ggsave(plot, file=paste(gastruloids, Well.ID_list[i], ".pdf", sep=''), dpi=300)

  }
}

graphs(normalised_data)


## Plots of triplicate wells 
well_plots <- "insert_folder_location"

well_graphs <- function(well_means, na.rm = TRUE, ...){
  
  marker_list <- unique(well_means$marker_line)
  
  for (i in seq_along(marker_list)) { 
    
    wellplot <- 
      ggplot(subset(well_means, well_means$marker_line==marker_list[i]),
             aes(bin, well_mean, group=Well.ID, colour=Well.ID)) + 
      geom_line(size=1.5) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold")) + 
      labs(x = "Distance from colony centre (μM)", y = "Fluorescence intensity (a.u.)") +
      ggtitle(paste(marker_list[i]), "")

    ggsave(wellplot, file=paste(well_plots, marker_list[i], ".pdf", sep=''), dpi=300)

  }
}

well_graphs(well_means)


#calculate % of cells expressing protein of interest within each colony 
all_data_expression <- all_data %>%
  group_by(Well.ID, Gastruloid.No, marker_line) %>%
  summarise(n_cells = n(),
            n_gt0 = sum(background_removed > 0),
            p_gt0 = (n_gt0 / n_cells) * 100)

write.table(all_data_expression, file = "all_data_expression_n1.txt", sep = "\t")



