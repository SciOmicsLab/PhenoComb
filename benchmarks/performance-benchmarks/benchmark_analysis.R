
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(plotly)
library(dplyr)



get_run_time <- function(file_path, type = "full"){
  
  txt <- readLines(file_path)
  
  if(type == "full"){
    
    start <- txt[1]
    end <- txt[length(txt)]  
    
  }
  
  if(type == "pre_process"){
    
    start <- txt[1]
    
    for(i in 1:length(txt)){
      
      if(grepl( "Getting unique phenotypes", txt[i], fixed = TRUE)){
        
        end <- txt[i]
        break
        
      }
      
    }
    
  }
  
  
  if(type == "combinatorics_run"){
    
    for(i in 1:length(txt)){
      
      if(grepl( "Getting unique phenotypes", txt[i], fixed = TRUE)){
        
        start <- txt[i]
        
      }
      
    }
    
    end <- txt[length(txt)] 
    
  }
  
  start <- regmatches(start, gregexpr("(?<=\\[).*?(?=\\])",start, perl=T))[[1]]
  start <- as.POSIXct(start,format="%Y/%m/%d %H:%M:%OS")
  end <- regmatches(end, gregexpr("(?<=\\[).*?(?=\\])",end, perl=T))[[1]]
  end <- as.POSIXct(end,format="%Y/%m/%d %H:%M:%OS")
  
  return(as.numeric(difftime(end, start, units = "sec")))
}





get_phenotypes_generated <- function(file_path){
  
  txt <- readLines(file_path)
  
  phenotypes_generated <- 0
  
  for(i in 1:length(txt)){
    
    if(grepl( "Writing", txt[i], fixed = TRUE)){
      
      if(grepl( "done", txt[i], fixed = TRUE)) next  #Fixing one other log line that thas "writing" word. Can be done more elegantly.
      
      phenotypes_generated <- phenotypes_generated + as.integer(sub(".*?Writing .*?(\\d+).*", "\\1", txt[i]))
      
    }
    
  }
  
  rm(txt)
  
  return(phenotypes_generated)
  
}


get_chunk_runtime <- function(file_path){
  
  txt <- readLines(file_path)
  
  chunk_runtime <- NULL
  
  for(i in 1:length(txt)){
    
    if(grepl( "Computing marker combinations", txt[i], fixed = TRUE)){
      
      start <- txt[i]
      
    }
    
    if(grepl( "Marker combinations computed", txt[i], fixed = TRUE)){
      
      end <- txt[i]
      
      start <- regmatches(start, gregexpr("(?<=\\[).*?(?=\\])",start, perl=T))[[1]]
      start <- as.POSIXct(start,format="%Y/%m/%d %H:%M:%OS")
      end <- regmatches(end, gregexpr("(?<=\\[).*?(?=\\])",end, perl=T))[[1]]
      end <- as.POSIXct(end,format="%Y/%m/%d %H:%M:%OS")
      
      chunk_runtime <- c(chunk_runtime, as.numeric(difftime(end, start, units = "sec")))
      
    }
    
  }
  
  return(chunk_runtime)
  
}

get_chunk_phenotypes_generated <- function(file_path){
  
  txt <- readLines(file_path)
  
  phenotypes_generated <- NULL
  
  for(i in 1:length(txt)){
    
    if(grepl( "Writing", txt[i], fixed = TRUE)){
      
      phenotypes_generated <- c(phenotypes_generated, as.integer(sub(".*?Writing .*?(\\d+).*", "\\1", txt[i])))
      
    }
    
  }
  
  return(phenotypes_generated)
    
}

get_chunk_info <- function(results_folder,experiment_id,replicates){
  
  file_name <- file.path(results_folder, paste("experiment_",experiment_id,"_replicate_1.log",sep = ""))
  
  chunk_results <- get_chunk_phenotypes_generated(file_name)
  chunk_results <- data.frame(chunk = 1:length(chunk_results), phenotypes_generated = chunk_results)
  
  chunk_results <- merge(chunk_results,data.frame(Replicate = 1:replicates))
  
  chunk_runtimes <- NULL
  for(i in 1:replicates){
    file_name <- file.path(results_folder, paste("experiment_",experiment_id,"_replicate_",i,".log",sep = ""))
    chunk_runtimes <- c(chunk_runtimes,get_chunk_runtime(file_name))
  }
  chunk_results$runtime <- chunk_runtimes
  
  return(chunk_results)
}


read_results <- function(experiments,results_folder,replicates){
  
  results <- merge(experiments,data.frame(Replicate = 1:replicates))
  
  
  for(i in 1:nrow(results)){
    
    file_name <- paste("experiment_",results[i,"experiment_ID"],"_replicate_",results[i,"Replicate"],".log",sep = "")
    file_path <- file.path(results_folder,file_name)
    
    if(file.exists(file_path)){
      
      results[i,"preprocess_runtime"] <- get_run_time(file_path,type = "pre_process")
      results[i,"combinatorics_runtime"] <- get_run_time(file_path,type = "combinatorics_run")
      results[i,"runtime"] <- get_run_time(file_path)
      results[i,"phenotypes_generated"] <- get_phenotypes_generated(file_path)  
      
    }
    
  }
  
  return(results)
  
}



plot_nphenotypes <- function(){
  
  
  
  phen_data <- expand.grid(4:28,2:3)
  colnames(phen_data) <- c("markers","states")
  phen_data$phenotypes <- (phen_data$states + 1)^(phen_data$markers)
  
  
  
  fig <- plot_ly(phen_data, x = ~markers, y = ~phenotypes, color = ~factor(states), type = 'scatter', mode = 'markers+lines',
                 colors = c('#29a329','#006450'),
                 marker = list(),
                 showlegend = TRUE)
  
  
  
  
  fig <- fig %>% layout(yaxis = list(title = "Number of Phenotypes", type = "log"),
                        xaxis = list(
                          title = "Number of Markers",
                          tickmode = "array",
                          nticks = 7,
                          tickvals = c(4,10,15,20,25,30,35),
                          range = c(3, 35)
                          ),
                        legend = list(x = 0.08, y = 0.98,title=list(text='States'))
  )
  
  return(fig)
  
}




plot_nmarker_experiment_runtime <- function(results){
  
  runtime_data <- summarise(group_by(filter(results,experiment == "n_markers"), markers),
                            mean(combinatorics_runtime)/3600,
                            max(combinatorics_runtime)/3600,
                            min(combinatorics_runtime)/3600
  )
  
  
  
  colnames(runtime_data) <- c("x",
                              "mean",
                              "max",
                              "min"
  )
  
  
  fig <- plot_ly(runtime_data, x = ~x, y = ~max, type = 'scatter', mode = 'lines',
                 line = list(color = 'rgba(1,1,1,0)'),
                 showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~min, type = 'scatter', mode = 'lines',
                           
                           fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'rgba(1,1,1,0)'),
                           
                           showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~mean, type = 'scatter', mode = 'lines + markers',
                           line = list(color = 'rgba(0,100,80,1)'),
                           marker = list(color = 'rgba(41, 163, 41,1)'),
                           showlegend = FALSE)
  
  
  fig <- fig %>% layout(yaxis = list(title = "Runtime (hours)", type = "log"),
                        xaxis = list(
                          title = "Number of Markers",
                          tickmode = "array",
                          nticks = 7,
                          tickvals = c(4,10,15,20,25,30,35),
                          range = c(3, 35)
                          )
  )
  
  return(fig)
  
}











plot_max_phen_len_experiment_runtime <- function(results){
  

  results[results$experiment == "n_markers" & results$markers == 15, "experiment"] <- "max_phen_len"
  
  runtime_data <- summarise(group_by(filter(results,experiment == "max_phen_len"), max_phen_len),
                            mean(combinatorics_runtime),
                            max(combinatorics_runtime),
                            min(combinatorics_runtime)
  )
  
  
  
  colnames(runtime_data) <- c("x",
                              "mean",
                              "max",
                              "min"
  )
  
  
  fig <- plot_ly(runtime_data, x = ~x, y = ~max, type = 'scatter', mode = 'lines',
                 line = list(color = 'rgba(1,1,1,0)'),
                 showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~min, type = 'scatter', mode = 'lines',
                           
                           fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'rgba(1,1,1,0)'),
                           
                           showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~mean, type = 'scatter', mode = 'lines + markers',
                           line = list(color = 'rgba(0,100,80,1)'),
                           marker = list(color = 'rgba(41, 163, 41,1)'),
                           showlegend = FALSE)
  
  
  fig <- fig %>% layout(yaxis = list(title = "Runtime (seconds)", type = "log",
                                     range=c(1,2.4)),
                        xaxis = list(
                          title = "Maximum Phenotype Length"
                           )
  )
  
  return(fig)
  
}






plot_n_cells_experiment_runtime <- function(results){
  
  
  runtime_data <- summarise(group_by(filter(results,experiment == "n_cells"), cells),
                            mean(combinatorics_runtime),
                            max(combinatorics_runtime),
                            min(combinatorics_runtime)
  )
  
  
  
  colnames(runtime_data) <- c("x",
                              "mean",
                              "max",
                              "min"
  )
  
  runtime_data <- runtime_data[runtime_data$x >= 10000,]
  
  
  fig <- plot_ly(runtime_data, x = ~x, y = ~max, type = 'scatter', mode = 'lines',
                 line = list(color = 'rgba(1,1,1,0)'),
                 showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~min, type = 'scatter', mode = 'lines',
                           fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'rgba(1,1,1,0)'),
                           showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~mean, type = 'scatter', mode = 'lines + markers',
                           line = list(color = 'rgba(0,100,80,1)'),
                           marker = list(color = 'rgba(41, 163, 41,1)'),
                           showlegend = FALSE)
  
  
  fig <- fig %>% layout(yaxis = list(title = "Runtime (seconds)",range = c(0,200)),
                        xaxis = list(
                          title = "Number of Cells",
                          tickmode = "array",
                          nticks = nrow(runtime_data),
                          tickvals = runtime_data$x
                        )
  )
  
  return(fig)
}







plot_n_samples_experiment_runtime <- function(results){
  
  
  runtime_data <- summarise(group_by(filter(results,experiment == "n_samples"), samples),
                            mean(combinatorics_runtime),
                            max(combinatorics_runtime),
                            min(combinatorics_runtime)
  )
  
  
  
  colnames(runtime_data) <- c("x",
                              "mean",
                              "max",
                              "min"
  )
  
  
  
  fig <- plot_ly(runtime_data, x = ~x, y = ~max, type = 'scatter', mode = 'lines',
                 line = list(color = 'rgba(1,1,1,0)'),
                 showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~min, type = 'scatter', mode = 'lines',
                           
                           fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'rgba(1,1,1,0)'),
                           
                           showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~mean, type = 'scatter', mode = 'lines + markers',
                           line = list(color = 'rgba(0,100,80,1)'),
                           marker = list(color = 'rgba(41, 163, 41,1)'),
                           showlegend = FALSE)
  
  
  fig <- fig %>% layout(yaxis = list(title = ""),
                        xaxis = list(
                          title = "Number of Samples",
                          tickmode = "array",
                          nticks = 5,
                          tickvals = c(20,50,80,100,200)
                        )
  )
  
  return(fig)
  
}





plot_n_threads_experiment_runtime <- function(results){
  
  
  runtime_data <- summarise(group_by(filter(results,experiment == "n_threads"), threads),
                            mean(combinatorics_runtime),
                            max(combinatorics_runtime),
                            min(combinatorics_runtime)
  )
  
  
  
  colnames(runtime_data) <- c("x",
                              "mean",
                              "max",
                              "min"
  )
  
  
  
  fig <- plot_ly(runtime_data, x = ~x, y = ~max, type = 'scatter', mode = 'lines',
                 line = list(color = 'rgba(1,1,1,0)'),
                 showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~min, type = 'scatter', mode = 'lines',
                           
                           fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'rgba(1,1,1,0)'),
                           
                           showlegend = FALSE)
  
  fig <- fig %>% add_trace(y = ~mean, type = 'scatter', mode = 'lines + markers',
                           line = list(color = 'rgba(0,100,80,1)'),
                           marker = list(color = 'rgba(41, 163, 41,1)'),
                           showlegend = FALSE)
  
  
  fig <- fig %>% layout(
                        xaxis = list(
                          title = "Number of Threads",
                          tickmode = "array",
                          nticks = 6,
                          tickvals = c(1,4,8,16,32,50)
                        )
  )
  
  return(fig)
  
}


add_point_to_figure <- function(fig,x,y,ann, color='purple'){
  
  dt <- data.frame(x = x,y = y, ann = ann)
  
  t <- list(
    
    family = "PT Sans Narrow",
    
    size = 7,
    
    color = toRGB("grey50"))
  
  fig <- fig %>% add_markers(data = dt, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
                           marker = list(color = color),
                           showlegend = FALSE,
                           inherit = FALSE)
  
  fig <- fig %>% add_text(data = dt, x = ~x, y = ~y, text = ~ann, textfont = t, textposition = "top right",showlegend = FALSE,
                          inherit = FALSE)
  
  return(fig)
  
}



experiments <- read.csv("experiments.csv")

replicates <- 3

results_folder <- "outputs"

results <- read_results(experiments,results_folder,replicates)


real_datasets_result_folder <- "real_datasets_results"

real_datasets <- c("combinatorial_phenotypes_parent_CD45+CD14+.log",
                    "combinatorial_phenotypes_parent_CD45+CD3-CD19-CD56-CD14-.log",
                    "combinatorial_phenotypes_parent_CD45+CD19+.log",
                    "combinatorial_phenotypes_parent_CD45+CD56+.log",
                    "combinatorial_phenotypes_parent_CD45+CD3+.log",
                    "HIV_combinatorial_phenotypes.log")


dataset_names <- c("CD45+CD14+",
                   "CD45+CD3-CD19-CD56-CD14-",
                   "CD45+CD19+",
                   "cCD45+CD56+",
                   "CD45+CD3+",
                   "HIV")

real_datasets_n_markers <- c(16,18,22,16,26,12)

real_datasets_runtimes <- unlist(lapply(real_datasets, function(i) get_run_time(file.path(real_datasets_result_folder,i),type = "combinatorics_run")/3600))

real_datasets_phenotypes <- unlist(lapply(real_datasets, function(i) get_phenotypes_generated(file.path(real_datasets_result_folder,i))))



n_phenotypes <- plot_nphenotypes()

n_phenotypes <- add_point_to_figure(n_phenotypes,real_datasets_n_markers,real_datasets_phenotypes,dataset_names)


n_marker <- plot_nmarker_experiment_runtime(results)

n_marker <- add_point_to_figure(n_marker,real_datasets_n_markers,real_datasets_runtimes,dataset_names)

# adding flowType in orange for comparison
n_marker <- add_point_to_figure(n_marker,  flowtype_results$markers, flowtype_results$combinatorics_runtime, "", color='orange')

row1 <- subplot(
        n_phenotypes,
        n_marker,
        plot_max_phen_len_experiment_runtime(results),
        nrows = 1,
        widths = c(0.3,0.3,0.3),
        margin = 0.05,
        titleX = TRUE, titleY = TRUE
        )


row2 <- subplot(
        plot_n_cells_experiment_runtime(results),
        plot_n_samples_experiment_runtime(results),
        plot_n_threads_experiment_runtime(results),
        nrows = 1,
        widths = c(0.3,0.3,0.3),
        margin = 0.05,
        titleX = TRUE, titleY = TRUE)

subplot(row1,row2,nrows = 2,margin = 0.08, titleX = TRUE, titleY = TRUE)

save_image(subplot(row1,row2,nrows = 2,margin = 0.08, titleX = TRUE, titleY = TRUE), "benchmark.pdf",width = 800, height = 500)
