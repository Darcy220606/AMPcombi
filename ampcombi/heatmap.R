# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(d3heatmap)

# Load data 
amp_summary <- readr::read_csv("AMPcombi_summary.csv",show_col_types = FALSE) %>%
  unique()
#colnames(amp_summary) <- gsub("\\.", " ", colnames(amp_summary))

# get unique activities found by diamond alignment in column "Activity"
#splitlst <- strsplit(unique(amp_summary$Activity), ", ", fixed=T)
#activities <- splitlst[!is.na(splitlst)]
#activities <- unique(unlist(splitlst, recursive=FALSE))
#activities <- activities[!(activities %in% c("Not found", "Unknown")) & !is.na(activities)]

samples = c(unique(amp_summary[["name"]]))

# Select data (name, contig_id and activity) and delete rows with activity NA
data <- amp_summary[,c("name","contig_id","Activity")]
data$Activity[data$Activity == "Not found"] <- 'Activity unknown'
data$Activity[data$Activity =="Unknown"] <- 'Activity unknown'
data$Activity[is.na(data$Activity)] <- 'Not found in reference db'

#c("Not found", "Unknown") "activity unknown"
#NA "not found in reference db"
activities <- unique(data[["Activity"]])

#student[!is.na(student$science),]
df_sample1 = data[data$name %in% samples[1],]
df_sample2 = data[data$name %in% samples[2],]
#df_act2 <-  subset(df_sample2, grepl(paste(activities[2], collapse = "|"), Activity))
#df_a2 <- data[data$Activity %in% activities[2]]

# Create a matrix with sample_names and activities containing list of corresponding contigs
count_mat <- matrix(nrow=length(activities),ncol=length(samples))
node_mat <- matrix(nrow=length(activities),ncol=length(samples))
colnames(count_mat); colnames(node_mat) = samples
rownames(count_mat), rownames(node_mat) = activities

for(i in 1:length(samples)){
  df_sub = data[data$name == samples[i],]
  for(j in 1:length(activities)){
    df_sub2 <- df_sub[df_sub$Activity == activities[j],]
    count <- nrow(df_sub2)
    contig_names <- df_sub$contig_id
    nodes <- paste(contig_names, collapse=",")
    node_mat[j,i] <- nodes
    count_mat[j,i] <- count
  }
}

rownames(mat) <- mat[,1]
mat <- mat %>% dplyr::select(-Country, -Group, -Continent)
mat <- as.matrix(mat)

# Heatmap
#d3heatmap(mtr, scale="column", dendrogram = "none", width="800px", height="80Opx", colors = "Blues")

heatmap(count_mat)
        
library(heatmaply)
p <- heatmaply(count_mat, 
               dendrogram = "none",
               #xlab = "", ylab = "", 
               #main = "",
               scale = "column",
               #margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               #hide_colorbar = TRUE,
               branches_lwd = 0.1,
               #label_names = c("Country", "Feature:", "Value"),
               fontsize_row = 5, fontsize_col = 5,
               labCol = colnames(count_mat),
               labRow = rownames(count_mat),
               heatmap_layers = theme(axis.line=element_blank())
)
