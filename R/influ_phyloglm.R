#' Leave-one-out-deletion analysis for gls phylogenetic regression. (GW: Make helpfile still!)

influ_phyloglm <- function(formula,data,phy,btol=50,cutoff=2,...)
{
        # Basic error checking:
        if(class(formula)!="formula") stop("Please formula must be class
                                           'forumla'")
        if(class(data)!="data.frame") stop("Please data must be class
                                           'data.frame'")
        if(class(phy)!="phylo") stop("Please phy must be class
                                           'phylo'")
        if (sum(rownames(data) != phy$tip.label) > 0) stop("Species must be at the same order
                                                      in data and phy")
        else

        # FULL MODEL calculations:
        full.data <- data
        N         <- nrow(full.data)
        mod.0 <- phylolm::phyloglm(formula, data=full.data,
                                   phy=phy,method="logistic_MPLE",btol=btol,...)
        intercept.0      <- mod.0$coefficients[[1]]             # Intercept (full model)
        slope.0          <- mod.0$coefficients[[2]]             # Slope (full model)
        pval.intercept.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]] # p.value (intercept)
        pval.slope.0     <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]] # p.value (slope)
        optpar.0         <- mod.0$alpha
        if(isTRUE(mod.0$convergence!=0)) stop("Full model failed to converge, consider changing btol. See ?phyloglm")
        else
                #Create the influ.model.estimates data.frame
                influ.model.estimates<-data.frame("species" =numeric(), "intercept"=numeric(),
                                                  "DFintercept"=numeric(),"intercept.perc"=numeric(),"pval.intercept"=numeric(),
                                                  "slope"=numeric(),"DFslope"=numeric(),"slope.perc"=numeric(),
                                                  "pval.slope"=numeric(),"AIC"=numeric(),
                                                  "optpar" = numeric())

        #Loop:
        counter <- 1
        errors <- NULL

        for (i in 1:N){
                crop.data <- full.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])

                mod=try(phylolm::phyloglm(formula, data=crop.data,
                                          phy=crop.phy,method="logistic_MPLE",btol=btol,...),TRUE)
                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(full.data$data)[i]
                        errors <- c(errors,error)
                        next }
                else {
                        ### Extracting model estimates:

                        sp                   <- phy$tip.label[i]      # species removed
                        intercept            <- mod$coefficients[[1]] # Intercept (crop model)
                        slope                <- mod$coefficients[[2]] # Beta (crop model)
                        DFintercept          <- intercept - intercept.0 # DF intercept
                        DFslope              <- slope - slope.0 # DF beta
                        intercept.perc       <- round((abs(DFintercept/intercept.0))*100,digits=1)  # Percentage of intercept change
                        slope.perc           <- round((abs(DFslope/slope.0))*100,digits=1)  # Percentage of beta change
                        pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]] # p.value (intercept)
                        pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]] # p.value
                        aic.mod              <- mod$aic # Model AIC
                        optpar               <- mod$alpha
                        print(i)

                        ### Storing values for each simulation
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

        sDFbetas <- DFbetas/sd(DFbetas)
        sDFintercepts <- DFintercepts/sd(DFintercepts)
        # Dataframe with results:
        estimates <- data.frame(species,betas,DFbetas,sDFbetas,intercepts,DFintercepts,sDFintercepts,
                                p.values)

        param0 <- data.frame(intercept.0,beta.0)

        ### Statistically Influential species for Beta (sDFbetas > 2)
        sb.ord <- which(abs(estimates$sDFbetas) > 2)
        influ.sp.b <- as.character(estimates$species[sb.ord])

        ### Statistically Influential species for intercept (sDFintercepts > 2)
        si.ord <- which(abs(estimates$sDFintercepts) > 2)
        influ.sp.i <- as.character(estimates$species[si.ord])
        influ.sp.i <- as.character(estimates[order(estimates$DFintercepts,decreasing=T)[1:5],]$species)
        #Should we add these to the input too? Only first 5?

        output <- list(errors=errors,formula=formula,
                       model_estimates=param0,
                       influential_species= influ.sp.b,
                       results=estimates,data=c.data)

        if (length(output$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                print("No errors found. All single deletions were performed and stored successfully")
                output$errors <- "No errors found."
        }

        return(output)

}
