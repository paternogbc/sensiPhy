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

        # Loop:

        for (i in 1:nrow(c.data)){
                exclude <- c(1:nrow(c.data))[-i]
                crop.data <- c.data[exclude,]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])
                mod=try(phylolm::phyloglm(formula, data=crop.data,
                                          phy=crop.phy,method="logistic_MPLE",btol=btol,...),TRUE)

                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(c.data$data)[i]
                        errors <- c(errors,error)
                        next }

                else {

                        ### Calculating model estimates:
                        intercept <-    mod$coefficients[[1]]       # Intercept (crop model)
                        beta <-    mod$coefficients[[2]]            # Beta (crop model)
                        alpha <-    mod$alpha                #Alpha (phylogenetic correlation parameter)
                        pval <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
                        DFbeta <- beta - beta.0
                        DFint  <- intercept - intercept.0
                        sp <- phy$tip.label[i]
                        pval.intercept <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
                        alpha<-mod$alpha

                        ### Storing values for each simulation
                        betas <- c(betas,beta)
                        intercepts <- c(intercepts,intercept)
                        DFbetas <- c(DFbetas,DFbeta)
                        DFintercepts <- c(DFintercepts,DFint)
                        species <- c(species,sp)
                        p.values <- c( p.values,pval)
                        print(i)
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
