#' Leave-one-out-deletion analysis for pgls regression.
#'
#' \code{influ_pgls} performs leave-one-out-deletion analysis for
#' \code{\link{pgls}} linear regression. It removes one species at a time, fits
#' a pgls model without the species and store the results of the model
#' estimates. The function repeats this procedure for all species in the dataset
#' and store all the results in a data.frame.
#' @aliases influ_pgls
#' @param formula A model formula
#' @param data A 'comparative.data' object containing the covariance matrix and
#'   data to be used in the model.
#' @param lambda A value for the lambda transformation. If NULL, lambda= "ML
#' @details This functions only works for simple linear regression \eqn{y = bx +
#'   a}. Future implementation will deal with more complex models.
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{[caper]pgls}}
#' @return The function \code{influpgls} returns a list with the following
#'   components:
#' @return \code{errors} Species removal that showed erros during pgls fit
#' @return \code{formula} The model formula
#' @return \code{model_estimates} Full model estimates
#' @return \code{beta95_IC} Full model beta 95 confidence interval
#' @return \code{influential_species} Most influential species for beta and
#'   intercept
#' @return \code{results} A data frame with all simulation estimates. DFbeta and
#'   DFintercept represent absolute difference between full model and simulation
#'   estimate, beta and intercept respectively.
#' @export


influ_pgls <- function(formula,data,lambda="ML")
{
        # Basic error checking:
        if(class(formula)!="formula") stop("Please formula must be class
                                           'forumla'")
        if(class(data)!="comparative.data") stop("data object must be of class
        'comparative.data. See ?comparative.data, package (caper) for details")
        else

        # FULL MODEL calculations:
        c.data <- data
        N <- nrow(c.data$data)             # Sample size
        mod.0 <- pgls(formula, data=c.data,lambda=lambda)
        sumMod <- summary(mod.0)
        intercept.0 <-    sumMod[[c(5,1)]]  # Intercept (full model)
        beta.0 <-    sumMod[[c(5,2)]]       # Beta (full model)
        pval.0 <-    sumMod[[c(5,8)]]       # p.value (full model)
        sd.beta.0 <- sumMod[[c(5,4)]]       # Standart Error (full model)
        df.0 <- sumMod[[2]][2]              # Degree if Freedon (full model))
        beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
        beta.0.low <- beta.0 - beta.IC      # Low limit of beta CI (full model)
        beta.0.up <- beta.0 + beta.IC       # Up limit of beta CI (full model)

        # Sampling effort analysis:
        betas <- as.numeric()
        intercepts <- as.numeric()
        DFbetas <- as.numeric()
        DFintercepts <- as.numeric()
        p.values <- as.numeric()
        species <- as.character()
        errors <- as.numeric()
        # Loop:

        for (i in 1:nrow(c.data$data)){
                exclude <- c(1:nrow(c.data$data))[-i]
                crop.data <- c.data[exclude,]

                mod=try(pgls(formula, data=crop.data,lambda),TRUE)
                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(c.data$data)[i]
                        errors <- c(errors,error)
                        next }

                else {
                        ### Calculating model estimates:
                        sum.Mod <- summary(mod)
                        beta <-    sum.Mod[[c(5,2)]]      # Beta
                        intercept <-    sum.Mod[[c(5,1)]] # Intercept
                        pval <-    sum.Mod[[c(5,8)]]      # p.value
                        DFbeta <- beta - beta.0
                        DFint  <- intercept - intercept.0
                        sp <- c.data$phy$tip.label[i]
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
        results=estimates,data=c.data$data)

        if (length(output$errors) >0){ print("Some species deletion
                presented errors,please check: output$errors to see
                which species deletion showed error")}
          else {
                print("No erros found. All single deletions were performed
                      and stored successfully")
                output$errors <- c("No erros found")}

        return(output)

}

