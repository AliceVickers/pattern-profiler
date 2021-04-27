
library(dplyr) 
library(ggplot2)

## read in files 
n1_data <- read.delim("marker_means_n1.txt", header=TRUE, sep="\t")
n2_data  <- read.delim("marker_means_n2.txt", header=TRUE, sep="\t")
n3_data  <- read.delim("marker_means_n3.txt", header=TRUE, sep="\t")

#merge tables
all_data <- rbind(n1_data, n2_data, n3_data)

#add column to group by protein marker 
#cell line 1
SOX2_1 <- c("SOX2_1_n1", "SOX2_1_n2", "SOX2_1_n3")
BRA_1 <- c("BRA_1_n1", "BRA_1_n2", "BRA_1_n3")
SOX17_1 <- c("SOX17_1_n1", "SOX17_1_n2", "SOX17_1_n3")

#cell line 2
SOX2_2 <- c("SOX2_2_n1", "SOX2_2_n2", "SOX2_2_n3")
BRA_2 <- c("BRA_2_n1", "BRA_2_n2", "BRA_2_n3")
SOX17_2 <- c("SOX17_2_n1", "SOX17_2_n2", "SOX17_2_n3")

#cell line 3
SOX2_3 <- c("SOX2_3_n1", "SOX2_3_n2", "SOX2_3_n3")
BRA_3 <- c("BRA_3_n1", "BRA_3_n2", "BRA_3_n3")
SOX17_3 <- c("SOX17_3_n1", "SOX17_3_n2", "SOX17_3_n3")

all_data$group <- ifelse(all_data$marker_line %in% SOX2_1, "SOX2_1", 
                         ifelse(all_data$marker_line %in% BRA_1, "BRA_1",
                                ifelse(all_data$marker_line %in% SOX17_1, "SOX17_1",
                                       ifelse(all_data$marker_line %in% SOX2_2, "SOX2_2",
                                               ifelse(all_data$marker_line %in% BRA_2, "BRA_2",
                                                      ifelse(all_data$marker_line %in% SOX17_2, "SOX17_2",
                                                             ifelse(all_data$marker_line %in% SOX2_3, "SOX2_3",
                                                                     ifelse(all_data$marker_line %in% BRA_3, "BRA_3",
                                                                             ifelse(all_data$marker_line %in% SOX17_3, "SOX17_3", "NA")))))))))

#find mean of biological replicates 
averaged_reps <- all_data %>%
  group_by(group, bin) %>%
  summarise(mean_reps = mean(marker_mean),
            sd_reps = sd(marker_mean)) 

cell_line_1 <- c("SOX2_1", "BRA_1", "SOX17_1")
cell_line_2 <- c("SOX2_2", "BRA_2", "SOX17_2")
cell_line_3 <- c("SOX2_3", "BRA_3", "SOX17_3")

averaged_reps$cell_line <- ifelse(averaged_reps$group %in% cell_line_1, "cell_line_1",
                                  ifelse(averaged_reps$group %in% cell_line_2, "cell_line_2",
                                         ifelse(averaged_reps$group %in% cell_line_3, "cell_line_3", "NA")))

##Plot spatial patterns of each marker on one graph (grouped by cell line)
all_markers <- "insert_folder_location"

graphs <- function(averaged_reps, na.rm = TRUE, ...){
  
  group_list <- unique(averaged_reps$cell_line)
  
  for (i in seq_along(group_list)) { 
    
    plot <- 
      ggplot(subset(averaged_reps, averaged_reps$cell_line==group_list[i]),
             aes(bin, mean_reps, group=group, color=group)) + 
      geom_line(size=2) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
      scale_x_continuous(expand = c(0, 0)) +
      geom_errorbar(aes(ymin=mean_reps-sd_reps, ymax=mean_reps+sd_reps), width=.2, position=position_dodge(0.05)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"), axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold"), legend.title = element_blank()) + 
      labs(x = "Distance from colony centre (μM)", y = "Fluorescence intensity (a.u.)") +
      ggtitle(paste(group_list[i]), "")

    ggsave(plot, file=paste(all_markers, group_list[i], ".pdf", sep=''), dpi=300)

  }
}

graphs(averaged_reps)


## % positive plots

## read in files 
n1_datafr <- read.delim("all_data_expression_n1.txt", header=TRUE, sep="\t")
n2_datafr  <- read.delim("all_data_expression_n1.txt", header=TRUE, sep="\t")
n3_datafr  <- read.delim("all_data_expression_n1.txt", header=TRUE, sep="\t")

#merge tables
datafr <- rbind(n1_datafr, n2_datafr, n3_datafr)

#add column to group by protein marker 
#cell line 1
SOX2_1 <- c("SOX2_1_n1", "SOX2_1_n2", "SOX2_1_n3")
BRA_1 <- c("BRA_1_n1", "BRA_1_n2", "BRA_1_n3")
SOX17_1 <- c("SOX17_1_n1", "SOX17_1_n2", "SOX17_1_n3")

#cell line 2
SOX2_2 <- c("SOX2_2_n1", "SOX2_2_n2", "SOX2_2_n3")
BRA_2 <- c("BRA_2_n1", "BRA_2_n2", "BRA_2_n3")
SOX17_2 <- c("SOX17_2_n1", "SOX17_2_n2", "SOX17_2_n3")

#cell line 3
SOX2_3 <- c("SOX2_3_n1", "SOX2_3_n2", "SOX2_3_n3")
BRA_3 <- c("BRA_3_n1", "BRA_3_n2", "BRA_3_n3")
SOX17_3 <- c("SOX17_3_n1", "SOX17_3_n2", "SOX17_3_n3")

datafr$group <- ifelse(datafr$marker_line %in% SOX2_1, "SOX2_1",
                       ifelse(datafr$marker_line %in% BRA_1, "BRA_1",
                              ifelse(datafr$marker_line %in% SOX17_1, "SOX17_1",
                                     ifelse(datafr$marker_line %in% SOX2_2, "SOX2_2",
                                            ifelse(datafr$marker_line %in% BRA_2, "BRA_2",
                                                   ifelse(datafr$marker_line %in% SOX17_2, "SOX17_2",
                                                          ifelse(datafr$marker_line %in% SOX2_3, "SOX2_3",
                                                                 ifelse(datafr$marker_line %in% BRA_3, "BRA_3",
                                                                        ifelse(datafr$marker_line %in% SOX17_3, "SOX17_3", "NA")))))))))

##add column for cell names 
cell_line_1 <- c("SOX2_1", "BRA_1", "SOX17_1")
cell_line_2 <- c("SOX2_2", "BRA_2", "SOX17_2")
cell_line_3 <- c("SOX2_3", "BRA_3", "SOX17_3")
datafr$cell_line <- ifelse(datafr$group %in% cell_line_1, "cell_line_1",
                           ifelse(datafr$group %in% cell_line_2, "cell_line_2",
                                  ifelse(datafr$group %in% cell_line_3, "cell_line_3", "NA")))

#split dataframe by cell line 
split_datafr <- split(datafr, datafr$cell_line)
cell_lines <- names(split_datafr)
for(i in seq_along(cell_lines)){
  assign(cell_lines[i], split_datafr[[i]])
}

#calculate mean and SD for each protein 
cell_line_1_data <- cell_line_1 %>%
  group_by(group) %>%
  summarise(mean = mean(p_gt0),
            sd = sd(p_gt0))
cell_line_1 <- merge(cell_line_1, cell_line_1_data, by="group")
write.table(cell_line_1, file = "cell_line_1_data.txt", sep = "\t")

#plot protein expression data 
cell_line_1$group2 <- factor(cell_line_1$group, levels = c("SOX2_1", "BRA_1", "SOX17_1"))
cell_line_1_plot <- ggplot(cell_line_1, aes(group2, p_gt0)) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 1, fill = "grey39") +
  ylim(0,100) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, color="black", size=1.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), axis.title.y = element_text(size=16, face="bold"), axis.text=element_text(size=14, face="bold", colour = "black")) +
  labs(x = NULL, y = "Percentage positive (%)") +
  scale_x_discrete(labels=c("SOX2","BRA","SOX17")) +
  ggtitle("cell_line_1")
ggsave(cell_line_1_plot, file="cell_line_1_plot.tiff", width = 15, height = 15, units = "cm")

#calculate mean and SD for each protein 
cell_line_2_data <- cell_line_2 %>%
  group_by(group) %>%
  summarise(mean = mean(p_gt0),
            sd = sd(p_gt0))
cell_line_2 <- merge(cell_line_2, cell_line_2_data, by="group")
write.table(cell_line_2, file = "cell_line_2_data.txt", sep = "\t")

#plot protein expression data 
cell_line_2$group2 <- factor(cell_line_2$group, levels = c("SOX2_1", "BRA_1", "SOX17_1"))
cell_line_2_plot <- ggplot(cell_line_2, aes(group2, p_gt0)) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 1, fill = "grey39") +
  ylim(0,100) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, color="black", size=1.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), axis.title.y = element_text(size=16, face="bold"), axis.text=element_text(size=14, face="bold", colour = "black")) +
  labs(x = NULL, y = "Percentage positive (%)") +
  scale_x_discrete(labels=c("SOX2","BRA","SOX17")) +
  ggtitle("cell_line_2")
ggsave(cell_line_2_plot, file="cell_line_2_plot.tiff", width = 15, height = 15, units = "cm")


#calculate mean and SD for each protein 
cell_line_3_data <- cell_line_3 %>%
  group_by(group) %>%
  summarise(mean = mean(p_gt0),
            sd = sd(p_gt0))
cell_line_3 <- merge(cell_line_3, cell_line_3_data, by="group")
write.table(cell_line_3, file = "cell_line_3_data.txt", sep = "\t")

#plot protein expression data 
cell_line_3$group2 <- factor(cell_line_3$group, levels = c("SOX2_1", "BRA_1", "SOX17_1"))
cell_line_3_plot <- ggplot(cell_line_3, aes(group2, p_gt0)) + 
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 1, fill = "grey39") +
  ylim(0,100) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, color="black", size=1.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), axis.title.y = element_text(size=16, face="bold"), axis.text=element_text(size=14, face="bold", colour = "black")) +
  labs(x = NULL, y = "Percentage positive (%)") +
  scale_x_discrete(labels=c("SOX2","BRA","SOX17")) +
  ggtitle("cell_line_3")
ggsave(cell_line_3_plot, file="cell_line_3_plot.tiff", width = 15, height = 15, units = "cm")






