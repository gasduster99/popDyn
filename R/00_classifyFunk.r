#' A function for internal use that formats the dNdt function for internal use
#'
#' @param fun     the function to internalize
#' @param numPars number of parameters used by fun
dNdt_classify = function(fun, numPars=2){
        #
        name = as.character(substitute(fun))
        arg = formals(fun)
        bdy = body(fun)[-1]
        mainArgs = names(arg[1:(numPars)])
        xtraArgs = names(arg[(numPars+1):length(arg)])
        for( ea in xtraArgs ){
                #NOTE: i need to identify variable better from coincidental text
                var = sprintf("self$%s", ea)
                bdy = gsub(ea, var, bdy)
        }
        #
        eval(parse( text=sprintf("private$dNdt=function(%s, %s, dNdt_par){}", names(arg[1]), names(arg[2])) ))
        body(private$dNdt) = parse( text=c('{', bdy, '}') )
}

#' A function for internal use that formats the N0 function for internal use
#'
#' @param fun     the function to internalize
N0_classify = function(fun){
        #

        #
        name = as.character(substitute(fun))
	arg = formals(fun)
        bdy = body(fun)[-1]
        xtraArgs = names(arg) 
	for( ea in xtraArgs ){
                #
                var = sprintf("self$%s", ea)
                bdy = gsub(ea, var, bdy)
        }

        #
        eval(parse( text=sprintf("self$%s=function(){}", name) ))
        eval(parse( text=sprintf("body(self$%s)=parse(text=c('{', bdy, '}'))", name) ))
}

#' A function for internal use that formats a generic function for internal use
#'
#' @param fun     the function to internalize
#' @param numPars number of parameters used by fun
classify = function(fun, numPars=1){
        #

        #
        name = as.character(substitute(fun))
        arg = formals(fun)
        bdy = body(fun)[-1]
        mainArgs = names(arg[1:(numPars)])
        xtraArgs = names(arg[(numPars+1):length(arg)])
        for( ea in xtraArgs ){
                #
                var = sprintf("self$%s_%s", name, ea)
                bdy = gsub(ea, var, bdy)
        }

        #
        eval(parse( text=sprintf("self$%s=function(%s){}", name, paste(mainArgs, collapse=',')) ))
        eval(parse( text=sprintf("body(self$%s)=parse(text=c('{', bdy, '}'))", name) ))
}
