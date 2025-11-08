#' A function primarily for internal use 
#'
#' @param  selfToPar 
#'
#' @return parameter values taken from self
selfToPar = function(parNames){
        #check if variable names exist

        #
        parNome = c()
        parValues = c()
        for(pn in parNames){
                eval(parse( text=sprintf("parValues=c(parValues,self$%s)", pn) ))
                l = eval(parse( text=sprintf("length(self$%s)", pn) ))
                parNome = c(parNome, pn, rep("", l-1))
        }
        names(parValues) = parNome
        #
        return(parValues)
}

#NOTE: parValues should be passed with names

#' A function primarily for internal use 
#'
#' @param  parValues values of parameters given as named objects to be registered with self
parToSelf = function(parValues){
        #check is names exist

        #
        parNames = names(parValues)
        #
        i = 1

        for(pn in parNames){
                #update self
                if(pn!=""){
                        nome = pn
                        eval(parse( text=sprintf("self$%s=parValues[pn]", pn) ))
                } else{
                        eval(parse( text=sprintf("self$%s=c(self$%s, parValues[i])", nome, nome) ))
                }
                #
                i = i+1
        }
}



##
#selfToPar = function(parNames){
#        #check if variable names exist
#
#        #
#        parValues = c()
#        for(pn in parNames){
#                eval(parse( text=sprintf("parValues['%s']=self$%s", pn, pn) ))
#        }
#        #
#        return(parValues)
#}
#
##NOTE: parValues should be passed with names
#parToSelf = function(parValues){
#        #check is names exist
#
#        #
#        parNames = names(parValues)
#        #
#        for(pn in parNames){
#                #update self
#                eval(parse( text=sprintf("self$%s=parValues[pn]", pn) ))
#        }
#}
