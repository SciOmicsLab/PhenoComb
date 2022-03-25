
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(flowCore)
library(plotly)


folder <- "../PhenoCombAnalysis/input"

sample_data <- read.csv(file.path(folder,"sample_data.csv"))

flow_data <- as.data.frame(read.FCS(file.path(folder,"concat_1.fcs"), truncate_max_range = FALSE)@exprs)


aux_sample <- as.data.frame(table(flow_data$Sample_ID))

colnames(aux_sample) <- c("Sample_ID","cells")

aux_sample <- aux_sample[order(-aux_sample$cells),]

fig <- plot_ly(aux_sample,
         x = 1:nrow(aux_sample),
         y = ~cells,
         #error_y = ~list(array = preprocess_runtime_sd,color = '#111111'),
         type = 'bar',showlegend = F)

fig <- fig %>% layout(xaxis = list(title = 'Sample'), yaxis = list(type = "log"))

fig


filtered_samples <- as.numeric(as.character(aux_sample[aux_sample$cells >= 5000,"Sample_ID"]))

filtered_samples_data <- sample_data[sample_data$Sample_ID %in% filtered_samples,]

write.csv(filtered_samples_data,file.path(folder,"ncell_filtered_sample_data.csv"),row.names = F)
