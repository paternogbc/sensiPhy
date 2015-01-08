#' Leave-one-out-deletion analysis for gls phylogenetic regression.
#'
#' \code{influ_gls} performs leave-one-out-deletion analysis for
#' \code{\link{gls}} phylogenetic linear regression. It removes one species at a time, fits
#' a gls model without the species and store the results of the model
#' estimates. The function repeats this procedure for all species in the dataset
#' and store all the results in a data.frame.
#' @aliases influ_gls
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
#' @return \code{errors} Species that showed erros during pgls fit
#' @return \code{formula} The model formula
#' @return \code{model_estimates} Full model estimates
#' @return \code{beta95_IC} Full model beta 95 confidence interval
#' @return \code{influential_species} Most influential species for beta and
#'   intercept
#' @return \code{results} A data frame with all simulation estimates. DFbeta and
#'   DFintercept represent absolute difference between full model and simulation.
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{samp_pgls}}
#' @examples
#' library(caper);library(ggplot2);library(gridExtra)
#' data(shorebird)
#' comp.data <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
#' # First we need to match tip.labels with rownames in data:
#' sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
#' shorebird.data <- shorebird.data[sp.ord,]
#' # Now we can run the function influ_gls:
#' influ <- influ_gls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)
#' # Estimated parameters:
#' head(influ$results)
#' # Most influential species:
#' influ[[5]]
#' # Check for species with erros erros:
#' influ$errors
#' @export


influ_gls <- function(formula,data,phy)
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
        N <- nrow(c.data$data)
        cor.0 <- ape::corPagel(1,phy=phy,fixed=F)

        mod.0 <- nlme::gls(formula, data=c.data,method="ML",correlation=cor.0)
        sumMod <- as.data.frame(summary(mod.0)$tTable)

        intercept.0 <-    sumMod[1,1]       # Intercept (full model)
        beta.0 <-    sumMod[2,1]            # Beta (full model)
        pval.0 <-    sumMod[2,4]            # p.value (full model)
        sd.beta.0 <- sumMod[2,2]            # Standart Error (full model)
        df.0 <- N-1                         # Degree if Freedon (full model))
        beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
        beta.0.low <- beta.0 - beta.IC      # Low limit of beta CI (full model)
        beta.0.up <-  beta.0 + beta.IC      # Up limit of beta CI (full model)

        # Sampling effort analysis:
        betas <- as.numeric()
        intercepts <- as.numeric()
        DFbetas <- as.numeric()
        DFintercepts <- as.numeric()
        p.values <- as.numeric()
        species <- as.character()
        errors <- as.numeric()
        # Loop:

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
                        sum.Mod <- as.data.frame(summary(mod)$tTable)
                        intercept <-    sum.Mod[1,1]       # Intercept (full model)
                        beta <-    sum.Mod[2,1]            # Beta (full model)
                        pval <-    sum.Mod[2,4]            # p.value (full model)
                        DFbeta <- beta - beta.0
                        DFint  <- intercept - intercept.0
                        sp <- phy$tip.label[i]

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
        # Dataframe with results:
        estimates <- data.frame(species,betas,DFbetas,intercepts,DFintercepts,
                                p.values)
        param0 <- data.frame(intercept.0,beta.0)
        influ.sp.b <- as.character(estimates[order(estimates$DFbetas,decreasing=T)[1:5],]$species)
        influ.sp.i <- as.character(estimates[order(estimates$DFintercepts,decreasing=T)[1:5],]$species)
        beta_IC <- data.frame(beta.low=beta.0.low,beta.up=beta.0.up)
        output <- list(errors=errors,formula=formula,
                       model_estimates=param0,
                       beta95_IC=beta_IC,
                       influential_species=rbind(beta=influ.sp.b,
                                                 intercept=influ.sp.i),
                       results=estimates,data=c.data)

        if (length(output$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                print("No erros found. All single deletions were performed and stored successfully")
                output$errors <- "No erros found."
        }

        return(output)

}
