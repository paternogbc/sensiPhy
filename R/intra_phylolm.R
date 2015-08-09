intra_phylolm <- function(formula,data,phy,nintra=1,
                          vari.resp=NA,vari.pred=NA,minmax=FALSE,
                          method="uniform",model="lambda",track=TRUE){

  #Error check
  if (!inherits(phy, "phylo"))
    stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
  
  if (method=="normal" && minmax==T)
    stop("Cannot generate normal distribution from min and max values!")
  
  #Matching tree and phylogeny using utils.R
  datphy<-match_dataphy(formula,data,phy)
  data<-datphy[[1]]
  phy<-datphy[[2]]
  
  resp<-all.vars(formula)[1]
  pred<-all.vars(formula)[2]
  
  if(exists("vari.resp") && sum(is.na(data[,vari.resp]))!=0){
    data[is.na(data[,vari.resp]),] <- 0}
    
    if(exists("vari.pred") && sum(is.na(data[,vari.pred]))!=0){
      data[is.na(data[,vari.pred]),] <- 0}



  #Function to pick a random value in the interval
  if (method=="normal") funr <- function(a, b) {rnorm(1,a,b)}
  else  funr <- function(a,b) {runif(1,a-b,a+b)}
  

  #Create the results data.frame
  intra.model.estimates<-data.frame("n.intra"=numeric(),"intercept"=numeric(),"se.intercept"=numeric(),
                                   "pval.intercept"=numeric(),"slope"=numeric(),"se.slope"=numeric(),
                                   "pval.slope"=numeric(),"aic"=numeric(),"optpar"=numeric())
  

  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  
  for (i in 1:nintra) {

      ##Set response and predictor variables
      #vari.resp is not provided or is not numeric, do not pick random value
      if(!inherits(data[,resp], c("numeric","integer")) || !exists("vari.resp")) {data$respV<-data[,resp]}
      
      #choose a random value in min/max if vari.resp is provided and minmax=T
      if(exists("vari.resp") && minmax==T)
      {data$respV<-apply(data[,c(resp,vari.resp)],1,function(x)runif(1,x[1],x[2]))}
      
      #choose a random value in [mean-se,mean+se] if vari.resp is provided and minmax=F
      if(exists("vari.resp") && minmax==F)
      {data$respV<-apply(data[,c(resp,vari.resp)],1,function(x)funr(x[1],x[2]))}
      
    
      #vari.pred is not provided or is not numeric, do not pick random value
      if(!inherits(data[,pred], c("numeric","integer")) || !exists("vari.pred")) {data$predV<-data[,pred]}
      
      #choose a random value in min/max if vari.pred is provided and minmax=T
      if(exists("vari.pred") && minmax==T)
      {data$predV<-apply(data[,c(pred,vari.pred)],1,function(x)runif(1,x[1],x[2]))}
      
      #choose a random value in [mean-se,mean+se] if vari.pred is provided
      if(exists("vari.pred") && is.null(dim(vari.pred)))
      {data$predV<-apply(data[,c(pred,vari.pred)],1,function(x)funr(x[1],x[2]))}

      #model
      mod = try(phylolm::phylolm(respV~predV, data=data, model=model,phy=phy),TRUE)

      if(isTRUE(class(mod)=="try-error")) {
        error <- i
        names(error) <- rownames(c.data$data)[i]
        errors <- c(errors,error)
        next }
      
      
      else{
        intercept            <- phylolm::summary.phylolm(mod)$coefficients[[1,1]]
        se.intercept         <- phylolm::summary.phylolm(mod)$coefficients[[1,2]]
        slope                <- phylolm::summary.phylolm(mod)$coefficients[[2,1]]
        se.slope             <- phylolm::summary.phylolm(mod)$coefficients[[2,2]]
        pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
        pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
        aic.mod              <- mod$aic
        n                    <- mod$n
        d                    <- mod$d
        
        if (model == "BM"){
          optpar <- NA
        }
        if (model != "BM"){
          optpar               <- mod$optpar
        }
        
        if(track==TRUE) print(paste("intra: ",i,sep=""))
        
        #write in a table
        intra.model.estimates[counter,1] <- i
        intra.model.estimates[counter,2] <- intercept
        intra.model.estimates[counter,3] <- se.intercept
        intra.model.estimates[counter,4] <- pval.intercept
        intra.model.estimates[counter,5] <- slope
        intra.model.estimates[counter,6] <- se.slope
        intra.model.estimates[counter,7] <- pval.slope
        intra.model.estimates[counter,8] <- aic.mod
        intra.model.estimates[counter,9]<-  optpar
        
        
        counter=counter+1
        
      }
     }

  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval<-stats::aggregate(.~n.intra, data=intra.model.estimates, mean)
  
  statresults<-data.frame(min=apply(intra.model.estimates,2,min),
                          max=apply(intra.model.estimates,2,max),
                          mean=apply(intra.model.estimates,2,mean),
                          sd_intra=apply(mean_by_randomval,2,sd))[-(1:2),]
  
  
  output <- list(analysis.type="intra_phylolm",formula=formula,
                 model_results=intra.model.estimates,N.obs=n,
                 stats=statresults)
  
  return(output)
}

