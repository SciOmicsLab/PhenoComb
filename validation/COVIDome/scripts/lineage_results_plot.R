
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)
library(ggtext)
library(gridExtra)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  t_test <- t.test(eval(as.name(measurevar)) ~ eval(as.name(groupvars[1])), data=data, conf.level=conf.interval)$conf.int
  
  datac$low_mean_diff_ci <- t_test[1]
  
  datac$high_mean_diff_ci <- t_test[2]
  
  return(datac)
}



sample_file <- "../PhenoCombAnalysis/input/Sample_Info.csv"

compiled_results_file <- "../results/compiled_lineage_results.csv"

sample_data <- read.csv(sample_file)
sample_data$Sample_ID <- as.character(sample_data$Sample_ID)

compiled_results <- read.csv(compiled_results_file,check.names = FALSE)


plot_phenotype <- function(sample_data, phenotype_data, result_row){
  
  plot_data <- data.frame(covid_status = sample_data$COVID_status, phenotype_frequency = unlist(phenotype_data[result_row,c(sample_data$Sample_ID)]))
  
  data_summary <- summarySE(plot_data, measurevar="phenotype_frequency", groupvars=c("covid_status"))
  
  p <- ggplot(plot_data, aes(x = covid_status, y = phenotype_frequency, color = covid_status)) +
    
    scale_color_manual(values=c("#69b3a2", "purple")) +
    
    geom_jitter(position=position_jitter(0.1), alpha=0.7) +
    
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .3, size = 0.9, color = "black") +
    
    geom_errorbar(data = data_summary, aes(ymin=phenotype_frequency-ci, ymax=phenotype_frequency+ci), width=.15, size = .6, color = "black") +
    
    
    
    
    labs(title = compiled_results[result_row,"Phenotype"], x = "COVID Status", y = "Phenotype Frequency") +
    
    theme_light() +
    
    theme(legend.position="none",
          plot.title = element_markdown(hjust = 0.5,face="bold",size=ifelse(result_row<=2,6,8)),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 8)
    )
  
  
  
  p <- p + annotate("text", x = 1.5, y = 0.9*layer_scales(p)$y$range$range[2],
                    label = paste('95% CI Mean Diff\n',format(round(data_summary[1,"low_mean_diff_ci"], 3), nsmall = 3)," ~ ",format(round(data_summary[1,"high_mean_diff_ci"], 3), nsmall = 3)),
                    size = 2)
  
  
  p
  
}


ggsave("../results/phenotype_plots.pdf",
       grid.arrange(grobs = lapply(1:nrow(compiled_results), function(i) plot_phenotype(sample_data, compiled_results, i)),                             
                    nrow = ,
                    padding	= 1),
       width = 9, height = 8, units = "in"
       )


