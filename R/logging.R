# Start of log
# Fix filecreate TRUE - Fixed 6/17/2 Ann S.

start_log <- function(log_path = NULL, append = FALSE, verbose = FALSE){
  
  options(PhenoComb.verbose = verbose)
  options(PhenoComb.log.file = log_path)
  
  if(!is.null(log_path)){
    
    if(!file.exists(log_path) | !append){
      
      file.create(log_path)
      
    }
    
    if(!append){
      
      print_log("Starting log...")
      
    }
    
  }
  
}


print_log <- function(..., new_line = TRUE){
  
  if(getOption("PhenoComb.verbose") | !is.null(getOption("PhenoComb.log.file"))){
    
    time <- paste(c("[", format(Sys.time(), "%y/%m/%d %X"),"] "), collapse = "")
    
    txt <- paste(c(time, ...), collapse = "")
    
    if(getOption("PhenoComb.verbose")) message(txt)
    
    if(!is.null(getOption("PhenoComb.log.file"))){
      
      if(new_line){
        txt <- paste(c(txt, "\n"), collapse = "")
      }
      
      cat(txt, file = getOption("PhenoComb.log.file"), append = TRUE)
      
    }
    
  }
  
}

stop_log <- function(){
  
  print_log("End logging.")
  
  options(PhenoComb.verbose = FALSE)
  options(PhenoComb.log.file = NULL)
  
}