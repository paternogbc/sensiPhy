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
#' @seealso \code{\link[caper]{pgls}}, \code{\link{samp_pgls}}
#' @export

influ_pgls <- function(formula,data,phy,model="lambda",...)
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
        if ((model == "trend") & (is.ultrametric(phy)))
                stop("the trend is unidentifiable for ultrametric trees.")
        else

        # FULL MODEL calculations:
        c.data <- data
        N <- nrow(c.data)
        mod.0 <- phylolm::phylolm(formula, data=c.data,model=model,phy=phy)

        a.0 <- mod.0$coefficients[[1]]             # Intercept (full model)
        b.0 <- mod.0$coefficients[[2]]             # Beta (full model)
        p.val.a0 <-    phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] # p.value (intercept)
        p.val.b0 <-    phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] # p.value (slope)
        optpar.0 <- mod.0$optpar

        #Create the results data.frame
        results<-data.frame("species" =numeric(), "intercept"=numeric(),
                            "DFintercept"=numeric(),"beta"=numeric(),
                            "DFbeta"=numeric(),"pval.intercept"=numeric(),
                            "pval.beta"=numeric(),"AIC"=numeric(),
                             "optpar" = numeric())

        # Loop:
        counter <- 1
        errors <- NULL

        for (i in 1:N){
                crop.data <- c.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])

                mod=try(phylolm::phylolm(formula, data=crop.data,model=model,phy=crop.phy),TRUE)

                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(c.data$data)[i]
                        errors <- c(errors,error)
                        next }

                else {
                        ### Calculating model estimates:

                        sp <- phy$tip.label[i]         # species removed
                        a <-    mod$coefficients[[1]]          # Intercept (crop model)
                        b <-    mod$coefficients[[2]]       # Beta (crop model)
                        DFa <- a - a.0                 # DF intercept
                        DFb <- b - b.0                 # DF beta
                        a.change <- round((abs(DFa/a.0))*100,digits=1)  # Percentage of intercept change
                        b.change <- round((abs(DFb/b.0))*100,digits=1)  # Percentage of beta change
                        pval.a <-    phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] # p.value (intercept)
                        pval.b <-    phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] # p.value
                        aic.mod <-   mod$aic            # Model AIC
                        optpar <-    mod$optpar# Estimated lambda
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
                        results[counter,9]<- optpar

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
                message("No errors found. All single deletions were performed and stored successfully. Please, check outpu$results.")
                res$errors <- "No errors found."
        }

        return(res)

}
