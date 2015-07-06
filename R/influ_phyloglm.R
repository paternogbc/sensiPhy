

influ_phyloglm <- function(formula,data,phy,btol=50,cutoff=2,...){
        if(class(formula)!="formula") stop("formula must be class 'formula'")
        if(class(data)!="data.frame") stop("data must be class 'data.frame'")
        if(class(phy)!="phylo") stop("phy must be class 'phylo'")
        else

        #Calculates the full model, extracts model parameters
        full.data       <- data
        N               <- nrow(full.data)
        mod.0           <- phylolm::phyloglm(formula, data=full.data,
                                   phy=phy,method="logistic_MPLE",btol=btol,...)
        intercept.0      <- mod.0$coefficients[[1]]
        slope.0          <- mod.0$coefficients[[2]]
        pval.intercept.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]]
        pval.slope.0     <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]]
        optpar.0         <- mod.0$alpha
        if(isTRUE(mod.0$convergence!=0)) stop("Full model failed to converge,
                                        consider changing btol. See ?phyloglm")
        else

        #Creates empy data frame to store model outputs
        influ.model.estimates<-
                data.frame("species" =numeric(), "intercept"=numeric(),
                           "DFintercept"=numeric(),"intercept.perc"=numeric(),
                           "pval.intercept"=numeric(),"slope"=numeric(),
                           "DFslope"=numeric(),"slope.perc"=numeric(),
                           "pval.slope"=numeric(),"AIC"=numeric(),
                           "optpar" = numeric())

        #Loops over all species, and removes each one individually
        counter <- 1
        errors <- NULL
        for (i in 1:N){
                crop.data <- full.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])
                mod=try(phylolm::phyloglm(formula, data=crop.data,phy=crop.phy,
                                          method="logistic_MPLE",btol=btol,...),
                        TRUE)
                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(full.data$data)[i]
                        errors <- c(errors,error)
                        next }
                else {  sp                   <- phy$tip.label[i]
                        intercept            <- mod$coefficients[[1]]
                        slope                <- mod$coefficients[[2]]
                        DFintercept          <- intercept - intercept.0
                        DFslope              <- slope - slope.0
                        intercept.perc       <- round((abs(DFintercept/intercept.0))*100,digits=1)
                        slope.perc           <- round((abs(DFslope/slope.0))*100,digits=1)
                        pval.intercept       <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
                        pval.slope           <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
                        aic.mod              <- mod$aic
                        optpar               <- mod$alpha
                        print(paste(i," / ",N,sep=""))

                        #Stores values for eacht simulation
                        influ.model.estimates[counter,1]  <- sp
                        influ.model.estimates[counter,2]  <- intercept
                        influ.model.estimates[counter,3]  <- DFintercept
                        influ.model.estimates[counter,4]  <- intercept.perc
                        influ.model.estimates[counter,5]  <- pval.intercept
                        influ.model.estimates[counter,6]  <- slope
                        influ.model.estimates[counter,7]  <- DFslope
                        influ.model.estimates[counter,8]  <- slope.perc
                        influ.model.estimates[counter,9]  <- pval.slope
                        influ.model.estimates[counter,10] <- aic.mod
                        influ.model.estimates[counter,11] <- optpar
                        counter=counter+1
                }
        }

        ### Calculating Standardized DFbeta and DFintercept
        sDFintercept <- influ.model.estimates$DFintercept/
                sd(influ.model.estimates$DFintercept)
        sDFslope     <- influ.model.estimates$DFslope/
                sd(influ.model.estimates$DFslope)

        influ.model.estimates$sDFslope     <- sDFslope
        influ.model.estimates$sDFintercept <- sDFintercept

        #Creates a list with full model estimates:
        param0 <- list(coef=phylolm::summary.phyloglm(mod.0)$coefficients,
                       aic=phylolm::summary.phyloglm(mod.0)$aic,
                       optpar=phylolm::summary.phyloglm(mod.0)$alpha)

        #Identifies influencital species (sDF > cutoff) and orders by influence
        reorder.on.slope         <-influ.model.estimates[order(abs(
                influ.model.estimates$sDFslope),decreasing=T),c("species","sDFslope")]
        influ.sp.slope           <-as.character(reorder.on.slope$species[abs(
                reorder.on.slope$sDFslope)>cutoff])
        reorder.on.intercept     <-influ.model.estimates[order(abs(
                influ.model.estimates$sDFintercept),decreasing=T),c("species","sDFintercept")]
        influ.sp.intercept       <-as.character(reorder.on.intercept$species[abs(
                reorder.on.intercept$sDFintercept)>cutoff])

        #Generates output:
        res <- list(analysis.type="influ_phyloglm",
                    cutoff=cutoff,
                    formula=formula,
                    full.model.estimates=param0,
                    influential.species= list(influ.sp.slope=influ.sp.slope,
                                              influ.sp.intercept=influ.sp.intercept),
                    influ.model.estimates=influ.model.estimates,
                    data=full.data,errors=errors)

        ### Warnings:
        if (length(res$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                message("No errors found. All single deletions were performed and stored successfully. Please, check outpu$influ.model.estimates.")
                res$errors <- "No errors found."
        }

        return(res)

}



