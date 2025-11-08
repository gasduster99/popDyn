##
#suppressWarnings(suppressMessages( library(GGally, quietly=TRUE) ))
#suppressWarnings(suppressMessages( library(mvtnorm, quietly=TRUE) ))

#
#

#' A function primarily for internal use that make colors transparent.
#'
#' @param  someColor a color given however R handles a color
#' @param  alpha     the alpha channel value that wil be added
#'
#' @return a version of someColor with the alpha channel added
makeTransparent = function(someColor, alpha=100){
        newColor = col2rgb(someColor)
        apply(newColor, 2, function(curcoldata){
                rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3], alpha=alpha, maxColorValue=255)
        })
}

#
#

#' A function to print
#'
#' @param ins  
#' @param outs 
printSelf = function(ins, outs=c()){
        #
        n = 5
        #
        nome = names(self)
	extraNames = c(outs, ".__enclos_env__", "initialize", "printSelf", "clone") #"iterate", "optimize", "clone", "model", "prior", "plotMean", "plotBand", "plotRS", "N0Funk", "save", "load")
        display = nome[!nome%in%extraNames | nome%in%ins]
	display = display[order(nchar(display))]
	if(length(ins)>0){ display=display[display%in%ins] } 
        #
        for(d in display){
                #
                text = sprintf("self$%s", d)
                if( typeof(eval(parse(text=text)))=='list' ){
                        #
                        cat( sprintf('%s\t:\n', d) )
                        print( eval(parse(text=text)) )
                } else{
                        #       
                        if(length(eval(parse( text=text )))>n){
                                cat( sprintf('%s\t:', d), eval(parse( text=text ))[1:n], '...\n' )
                        } else{
                                cat( sprintf('%s\t:', d), eval(parse( text=text )), '\n' )
                        }
                }
        }
}

#' A function to plot the mean of the dynamics
#'
#' @param col   
#' @param alpha 
#' @param lwd   
#' @param add
plotMean = function(col='black', alpha=100, lwd=3, add=F){
        #arguments passed directly to plotting functions

        #       
        if(!add){
                plot(self$time, exp(self$lq+log(self$N)),
                        type='l',
                        lwd=lwd,
                        col=col
                )
        } else{
                lines(self$time, exp(self$lq+log(self$N)),
                        lwd=lwd,
                        col=col
                )
        }
}

#' A function to plot uncertainty bands 
#'
#' @param prob  Size of the uncertainty bands as defined as posterior probability
#' @param col   
#' @param alpha 
plotBand = function(prob=0.95, col='black', alpha=100){
        #arguments passed directly to plotting functions

        #
        left = (1-prob)/2
        right = prob+left
        #
        polygon( c(self$time, rev(self$time)),
                c(
                        private$qLikes[[self$model$observation]](self, left),
                        rev(private$qLikes[[self$model$observation]](self, right))
                ),
                col = makeTransparent(col, alpha=alpha),
                border = NA
        )
}

#' A function to plot repeated sampling posterior-like parameter distributions
#'
#' @param m      the number of samples from the repeated sampling distribution
#' @param sample a boolean to indicate if samples should be returned
#' @param save   a boolean to indicate if samples should be saved 
plotRS = function(m=10^4, sample=F, save=F){
        #m      : how many samples
        #save   : FALSE or a filename

        #
        parNames = colnames(self$rsCov)
        sam = mvtnorm::rmvnorm(m, private$selfToPar(parNames), self$rsCov)
        #
        if(save==F){
                GGally::print_if_interactive(GGally::ggpairs(as.data.frame(sam)))
        } else{
		ggplot2::ggsave(filename=save, plot=GGally::ggpairs(as.data.frame(sam)))
        }
        #
        if(sample){ return(sam) }
}
