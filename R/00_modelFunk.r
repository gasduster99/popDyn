#' A list of functions for internal use that formats a generic function for internal use
#'
#' @param self a ddModel object
#' @param data data to be used with likelihood
dLikes = list(
        #
        LN = function(self, data){
                dnorm(log(data), self$lq+log(self$B), exp(self$lsdo), log=T)
        },
        #
        N = function(self, data){
                dnorm(data, exp(self$lq)*self$B, exp(self$lsdo), log=T)
        }
)

#' A list of functions for internal use that formats a generic function for internal use
#'
#' @param self a ddModel object
#' @param prob probability to be used with likelihood to get quantile 
qLikes = list(
        #
        LN = function(self, prob){
                qlnorm(prob, self$lq+log(self$B), exp(self$lsdo))
        },
        #
        N = function(self, prob){
                qnorm(prob, exp(self$lq)*self$B, exp(self$lsdo))
        }
)
