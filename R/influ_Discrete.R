
#' @export

influ_Discrete <- function(data,phy,model="ARD",
                           transform = "none",bounds = list(),
                           cutoff=2,track=TRUE,...){
  
        #Error check
        if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
        if(length(levels(data))>2) stop("discrete data can have maximal two levels")
        if(class(phy)!="phylo") stop("phy must be class 'phylo'")
        if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
        if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
        else
          
        #Matching tree and phylogeny using utils.R
        full.data<-data
        phy<-phy
          
        #Calculates the full model, extracts model parameters
        N                   <- nrow(full.data)
        mod.0               <- geiger::fitDiscrete(phy = phy,dat = full_data,
                                                   model = model,transform = transform,
                                                   bounds = bounds,ncores = NULL,...)
        q12.0               <- mod.0$opt$q12
        q21.0               <- mod.0$opt$q12
        aicc.0              <- mod.0$opt$aicc
        if (transform == "none"){
          optpar.0 <- NA
        }
        if (transform == "EB"){
          optpar.0               <- mod.0$opt$a
        }
        if (transform == "lambda"){
          optpar.0               <- mod.0$opt$lambda
        }
        if (transform == "kappa"){
          optpar.0               <- mod.0$opt$kappa
        }
        if (transform == "delta"){
          optpar.0               <- mod.0$opt$delta
        }
        
        #Creates empty data frame to store model outputs
        sensi.estimates<-data.frame("species" =numeric(),
                                    "q12"=numeric(),"DIFq12"= numeric(),"q12.perc"= numeric(),
                                    "q21"=numeric(),"DIFq21"= numeric(),"q21.perc"= numeric(),
                                    "aicc"=numeric(),"optpar"=numeric()) 

        #Loops over all species, and removes each one individually
        counter <- 1
        errors <- NULL
        if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = N, style = 3)
        for (i in 1:N){
                
                crop.data <- full.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,rownames(crop.data)))
                
                mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.data,
                                              model = model,transform = transform,
                                              bounds = bounds,ncores = NULL,...),TRUE)
                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(full.data$data)[i]
                        errors <- c(errors,error)
                        next }
                else {  sp                   <- phy$tip.label[i]
                        intercept            <- mod$coefficients[[1]]
                        estimate                <- mod$coefficients[[2]]
                        DIFintercept          <- intercept - intercept.0
                        DIFestimate              <- estimate - estimate.0
                        intercept.perc       <- round((abs(DIFintercept/intercept.0))*100,digits=1)
                        estimate.perc           <- round((abs(DIFestimate/estimate.0))*100,digits=1)
                        pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
                        pval.estimate           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
                        aic.mod              <- mod$aic
                        if (model == "BM" | model == "trend"){
                            optpar <- NA
                        }
                        if (model != "BM" & model != "trend" ){
                            optpar               <- mod$optpar
                        }

                        if(track==TRUE) utils::setTxtProgressBar(pb, i)
                        # Stores values for each simulation
                        estim.simu <- data.frame(sp, intercept, DIFintercept, intercept.perc,
                                                 pval.intercept, estimate, DIFestimate, estimate.perc,
                                                 pval.estimate, aic.mod, optpar,
                                                 stringsAsFactors = F)
                        sensi.estimates[counter, ]  <- estim.simu
                        counter=counter+1
                }
        }
        if(track==TRUE) on.exit(close(pb))
        #Calculates Standardized DFbeta and DIFintercept
        sDIFintercept <- sensi.estimates$DIFintercept/
                stats::sd(sensi.estimates$DIFintercept)
        sDIFestimate     <- sensi.estimates$DIFestimate/
                stats::sd(sensi.estimates$DIFestimate)

        sensi.estimates$sDIFestimate     <- sDIFestimate
        sensi.estimates$sDIFintercept <- sDIFintercept

        #Creates a list with full model estimates:
        param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                       aic=phylolm::summary.phylolm(mod.0)$aic,
                       optpar=mod.0$optpar)

       #Identifies influencital species (sDF > cutoff) and orders by influence
       reorder.on.estimate         <-sensi.estimates[order(abs(
               sensi.estimates$sDIFestimate),decreasing=T),c("species","sDIFestimate")]
       influ.sp.estimate           <-as.character(reorder.on.estimate$species[abs(
               reorder.on.estimate$sDIFestimate)>cutoff])
       reorder.on.intercept     <-sensi.estimates[order(abs(
               sensi.estimates$sDIFintercept),decreasing=T),c("species","sDIFintercept")]
       influ.sp.intercept       <-as.character(reorder.on.intercept$species[abs(
               reorder.on.intercept$sDIFintercept)>cutoff])

        #Generates output:
        res <- list(call = match.call(),
                    cutoff=cutoff,
                    formula=formula,
                    full.model.estimates=param0,
                    influential.species= list(influ.sp.estimate=influ.sp.estimate,
                                              influ.sp.intercept=influ.sp.intercept),
                    sensi.estimates=sensi.estimates,
                    data=full.data,errors=errors)
        class(res) <- "sensiInflu"
        ### Warnings:
        if (length(res$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                res$errors <- "No errors found."
        }

        return(res)

}

