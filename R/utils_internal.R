### Function to colapse into a single data.frame data stored in multiple lists
### Writen by Gustavo Paterno (2017)
recombine <- function(list, slot1, slot2 = NULL){
  ### One level indexing list[[slot1]]:
  if (is.null(slot2)){
    nam<-c(1:length(list))
    x <- lapply(list, function(x) x[[slot1]])
    
    if(is(x[[1]],"list"))
      stop("Please also provide a value for slot2")
    if(is.na(match(class(x[[1]]), c("data.frame", "matrix"))))
      stop("Check if slot indexing returns a valid data.frame or numeric value withn your list")
    
    n.rows<- nrow(x[[1]])
    N <- rep(nam, each = n.rows)
    res <- data.frame(iteration = N, do.call("rbind", x))
    return(res)
  }
  
  ### Two levels indixing list[[slot1]][[slot2]]
  if(is.null(slot2) == FALSE){ 
    nam<-c(1:length(list))
    x <- lapply(list, function(x) x[[slot1]][[slot2]])
    
    if(class(x[[1]]) == "data.frame" | class(x[[1]]) == "matrix" | class(x[[1]]) == "numeric" ){
      n.rows<- nrow(x[[1]])
      
      ### If there is only a single value (e.g. AIC, optpar)
      if(is.null(n.rows)) {
        n.rows <- length(x[[1]])
        est.name <- names(list[[1]][[slot1]])[[slot2]]
        N <- rep(nam, each = n.rows)
        res <- data.frame(iteration = N, do.call("rbind", x))
        colnames(res)[2] <- est.name
        return(res)
      }
      
      else
        N <- rep(nam, each = n.rows)
      res <- data.frame(iteration = N, do.call("rbind", x))
      return(res)
    }
    else stop("Check if slot indexing returns a valid data.frame or numeric value withn your list") 
  }
}