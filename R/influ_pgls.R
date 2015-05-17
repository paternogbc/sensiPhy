#' Leave-one-out-deletion analysis for gls phylogenetic regression.
#'
#' \code{influ_pgls} performs leave-one-out-deletion analysis for
#' \code{\link{gls}} phylogenetic linear regression. It removes one species at a time, fits
#' a gls model without the species and store the results of the model
#' estimates. The function repeats this procedure for all species in the dataset
#' and store all the results in a data.frame.
#' @aliases influ_pgls
#' @param formula A model formula
#' @param data A data frame containing variables that can be attributed to the taxa
#' at the tips of a phylogeny
#' @param phy A phylogeny (class 'phylo') to be matched to the data above
#' @details This functions only works for simple linear regression \eqn{y = bx +
#'   a}. Future implementation will deal with more complex models. This functions only
#'   works with Pagel correlation (see \code{corPagel}). If any error occur during simulation,
#'   check 'output$errors' for details.
#' @return The function \code{influ_gls} returns a list with the following
#'   components:
#'
#' @return \code{formula} The model formula
#' @return \code{model_estimates} Full model estimates
#' @return \code{influential_species} Most influential species for beta
#' @return \code{results} A data frame with all simulation estimates. DFbeta and
#'   DFintercept represent absolute difference between full model and simulation.
#'   sDFintercept and sDFbeta represent the standardized differences.
#' @return \code{data} Original dataset
#' @return \code{errors} Species that showed erros during pgls fit
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{samp_pgls}}
#' @export

influ_pgls <- function(formula,data,phy)
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
        c.data <- data
        N <- nrow(c.data)
        cor.0 <- ape::corPagel(1,phy=phy,fixed=F)

        mod.0 <- nlme::gls(formula, data=c.data,method="ML",correlation=cor.0)
        sumMod <- as.data.frame(summary(mod.0)$tTable)

        a.0 <-           sumMod[1,1]            # Intercept (full model)
        b.0 <-           sumMod[2,1]            # Beta (full model)
        p.val.a0 <-    anova(mod.0)[1,3]        # p.value (intercept)
        p.val.b0 <-    anova(mod.0)[2,3]        # p.value (beta)
        lambda.0 <-      as.numeric(mod.0$model)# Estimated lambda

        #Create the results data.frame
        results<-data.frame("species" =numeric(), "intercept"=numeric(),
                            "DFintercept"=numeric(),"beta"=numeric(),
                            "DFbeta"=numeric(),"pval.intercept"=numeric(),
                            "pval.beta"=numeric(),"AIC"=numeric(),
                            "Lambda"=numeric())
        # Loop:
        counter <- 1
        errors <- NULL

        for (i in 1:nrow(c.data)){
                exclude <- c(1:nrow(c.data))[-i]
                crop.data <- c.data[exclude,]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])
                crop.cor <- ape::corPagel(1,phy=crop.phy,fixed=F)
                mod=try(nlme::gls(formula, data=crop.data,method="ML",correlation=crop.cor),TRUE)

                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(c.data$data)[i]
                        errors <- c(errors,error)
                        next }

                else {
                        ### Calculating model estimates:
                        sumMod.crop <- as.data.frame(summary(mod)$tTable)

                        sp <- phy$tip.label[i]         # species removed
                        a <-    sumMod.crop[1,1]       # Intercept (crop model)
                        b <-    sumMod.crop[2,1]       # Beta (crop model)
                        DFa <- a - a.0                 # DF intercept
                        DFb <- b - b.0                 # DF beta
                        a.change <- round((abs(DFa/a.0))*100,digits=1)  # Percentage of intercept change
                        b.change <- round((abs(DFb/b.0))*100,digits=1)  # Percentage of beta change
                        pval.a <-    anova(mod)[1,3]        # p.value (intercept)
                        pval.b <-    anova(mod)[2,3]        # p.value (beta)
                        aic.mod <-   AIC(mod)            # Model AIC
                        lambda <-    as.numeric(mod$model)# Estimated lambda
                        print(i)

                        ### Storing values for each simulation
                        #write in a table
                        results[counter,1]<- sp
                        results[counter,2]<- a
                        results[counter,3]<- DFa
                        results[counter,4]<- b
                        results[counter,5]<- DFb
                        results[counter,6]<- pval.a
                        results[counter,7]<- pval.b
                        results[counter,8]<- aic.mod
                        results[counter,9]<- lambda

                        counter=counter+1
                }
        }

        ### Calculating Standardized DFbeta and DFintercept
        sDFintercept <- results$DFintercept/sd(results$DFintercept)
        sDFbeta <- results$DFbeta/sd(results$DFbeta)


        results$sDFbeta <- sDFbeta;
        results$sDFintercept <- sDFintercept

        ### Original model estimates:
        param0 <- data.frame(intercept=a.0,beta=b.0)

        ### Statistically Influential species for Beta (sDFbetas > 2)
        sb.ord <- which(abs(results$sDFbeta) > 2)
        influ.sp.b <- as.character(results$species[sb.ord])

        ### Output:
        res <- list(output="influ_pgls",
             formula=formula,
             model_estimates=param0,
             influential_species= influ.sp.b,
             results=results,data=c.data,errors=errors)
        ### Warnings:
        if (length(res$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                message("No erros found. All single deletions were performed and stored successfully. Please, check outpu$results.")
                res$errors <- "No erros found."
        }

        return(res)

}
