setwd("/Users/christina/Desktop/gal")
install.packages("naniar")
library(tidyverse)
library(naniar)
library('dplyr')
Z_score_analysis <- read_csv("Corrected_Z_score.csv")

# get rid of blanks
Z_score_analysis<-Z_score_analysis[!(Z_score_analysis$`Sys. Name`=='Blank'),] 
# turn 0 to NA in WT; turn NA to false
Z_score_analysis[,3:8][Z_score_analysis[,3:8]=='#N/A']<-NA  
Z_score_analysis <- Z_score_analysis[complete.cases(Z_score_analysis[, 3:8]), ]
Z_score_analysis$`Exp 1_Z` <- as.numeric(as.character(Z_score_analysis$`Exp 1_Z`))
Z_score_analysis$`Exp 2_Z` <- as.numeric(as.character(Z_score_analysis$`Exp 2_Z`))
Z_score_analysis$`Exp 3_Z` <- as.numeric(as.character(Z_score_analysis$`Exp 3_Z`))
Z_score_analysis$`ctrl 1_Z` <- as.numeric(as.character(Z_score_analysis$`ctrl 1_Z`))
Z_score_analysis$`ctrl 2_Z` <- as.numeric(as.character(Z_score_analysis$`ctrl 2_Z`))
Z_score_analysis$`ctrl 3_Z` <- as.numeric(as.character(Z_score_analysis$`ctrl 3_Z`))


#t-test for WT vs. GG control or SM control

p.value_GG<-c()
p.value_GG <- numeric(length(rownames(Z_score_analysis)))

for (i in 1:length(rownames(Z_score_analysis)))
{
  if (sum(!is.na(Z_score_analysis[i,3:8]))>=2) {
    
    GGp.value<-t.test(Z_score_analysis[i,6:8],Z_score_analysis[i,3:5])$p.value
    
    p.value_GG[i]<-GGp.value
    
    
  }
  
  
}


#FDR calculation

FDR<-unlist(p.adjust(p.value_GG,method="fdr"))

analysis_tab<-cbind(Z_score_analysis,p.value_GG,FDR)

#q value calculation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")
n
library(qvalue)
qval_GG<-qvalue(p.value_GG)
q_GG<-qval_GG$qvalues


localFDR<-qval_GG$lfdr


analysis_tab<-cbind(Z_score_analysis,p.value_GG,FDR, q_GG, localFDR)
names(analysis_tab)[names(analysis_tab) == "Z difference"] <- "Z_difference"

#Calculate Z score difference
analysis_tab$Exp_mean_Z <- rowMeans(analysis_tab[, c("Exp 1_Z", "Exp 2_Z", "Exp 3_Z")], na.rm = TRUE)
analysis_tab$Ctrl_mean_Z <- rowMeans(analysis_tab[, c("ctrl 1_Z", "ctrl 2_Z", "ctrl 3_Z")], na.rm = TRUE)
analysis_tab$Z_difference <- analysis_tab$Exp_mean_Z - analysis_tab$Ctrl_mean_Z
analysis_tab$abs_Z_difference <- ifelse(analysis_tab$Z_difference < 0, -analysis_tab$Z_difference, analysis_tab$Z_difference)

analysis_tab <- analysis_tab %>%
  mutate(
    positive_fold_change = ifelse(Z_difference >= 0, abs_Z_difference, NA),
    negative_fold_change = ifelse(Z_difference < 0, -abs_Z_difference, NA)
  )

ordered_analysis <- analysis_tab %>%
  arrange(desc(Z_difference))
#Fix the wrong gene names
#Pdr3YBL005W/AIM33YML087C
row_indices <- which(ordered_analysis$'Sys. Name' == 'YIL077C')
print(row_indices)
row_index <- 4
ordered_analysis[row_index, c('Sys. Name', 'Name')] <- c('YBL005W', 'PDR3')
row_index <- 17
ordered_analysis[row_index, c('Sys. Name', 'Name')] <- c('YML087C', 'AIM33')
row_index <- 19
ordered_analysis[row_index, c('Name')] <- c('RCI37')

#Get rid of empty sys name
ordered_analysis<-ordered_analysis[!(ordered_analysis$`Sys. Name`=='Empty'),] 

# Create a basic volcano plot -------------------------------------------------
vol_plot <- analysis_tab %>%
  ggplot(aes(x = Z_difference,
             y = -log10(q_GG))) + 
  geom_point() 

vol_plot # Visualise ggplot output

#Mark points
# Highlight the Top 20 hits
# Identify the top 20 points
top_20_points <- ordered_analysis[order(ordered_analysis$Z_difference, decreasing = TRUE), ][1:20, ]

# Highlight the top 20 points
highlighted_points <- geom_point(data = top_20_points, aes(x = Z_difference, y = -log10(q_GG)), color = "palegreen", size = 2)
# Highlight the mito points
mitohits <- subset(top_20_points, Name %in% c("MGE1", "RCI37", "COX9", "PET100", "JID1", "INA22", "TIM21","TIM54", "AIM33"))
mge1 <- subset(mitohits, Name %in% c("MGE1"))

# Add subplot layer to the main volcano plot -----------------------------------
#install ggrepel
library(ggrepel)
#load ggplot2

Vol_Plot_2 <- 
  ggplot(data = ordered_analysis, # Original data  
       aes(Z_difference,
           y = -log10(q_GG))) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "black") +
  geom_point(colour = "grey", alpha = 0.5) +
  geom_point(data = top_20_points, # New layer containing data subset mitohits30      
             size = 2,
             shape = 21,
             fill = "olivedrab",
             colour = "olivedrab") +
  geom_point(data = mitohits, # New layer containing data subset mitohits30      
             size = 2,
             shape = 21,
             fill = "maroon",
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_label_repel(size = 2, data = mitohits,# Add labels last so they appear as the top layer  
                   aes(label = Name),
                   force = 1,
                   nudge_y = 0.25) +
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                     limits = c(-5, 10)) +
  labs(x = "mean Z-score difference", y = "-log(q-value)")
Vol_Plot_2
dev.off()
tiff("vocanoplot-1C.tiff", units="in", width=5, height=5, res=300)
library(dplyr)
ordered_analysis <- ordered_analysis %>%
  mutate(logq_GG = -log(p.value_GG))

install.packages("openxlsx")
library(openxlsx)
write.xlsx(ordered_analysis, "ordered_z_score.xlsx")


 