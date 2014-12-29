#' Estimates time of samp_pgls and influ_pgls
#'
#' \code{sensi_time} Estimates simulation time to run \code{samp_pgls} and
#' \code{influ_pgls}
#' @aliases sensi_time
#' @export

sensi_time <- function(formula,data,times=20,breaks=seq(.1,.7,.1),lambda="ML"){

          ### Basic error checking:
          if(class(formula)!="formula") stop("Please formula must be class 'forumla'")
          if(class(data)!="comparative.data") stop("data data must be of class 'comparative.data'. See function `comparative.data`.")
          if(length(breaks)<2) stop("please include more then one break (eg. breaks=c(.3,.5)")
          else

          c.data <- data
          N <- nrow(c.data$data)             # Sample size
          tm0 <- system.time(mod.0 <- pgls(formula, data=c.data,lambda=lambda))
          limit <- sort(round((breaks)*nrow(c.data$data),digits=0))
          tmCs <- as.numeric()
          for (i in limit) {
                    crop.data <- c.data[-(1:i),]
                    tmC <- system.time(mod.0 <- pgls(formula, data=crop.data,lambda=lambda))
                    tmCs <- c(tmCs,tmC[3])
          }

          time.s <- sum((tmCs*times)) + tm0[3]
          time.m <- round(time.s/60,digits=2)

          names(time.m) <- "Time estimated (min)"
          time.influ <- round((tm0[3]*N)/60,digits=2)
          output <- data.frame(sampling=time.m,influence=time.influ)
          return(output)
}


