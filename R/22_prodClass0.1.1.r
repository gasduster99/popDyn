##
#suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
#suppressWarnings(suppressMessages( library(deSolve, quietly=TRUE) ))
##
#source('modelFunk.r')
#source('prodOutputFunk.r')
#source('classifyFunk.r')
#source('parSelfParFunk.r')
#source('prodOptimizeFunk.r')

#
#CLASS
#

#' @title ProdModel
#'
#' @description The main production modeling class definition
#'
#' @export
prodModel = R6::R6Class("ProdModel", lock_objects=FALSE,
	#
	public = list(
		#POP
		#' @field N      a vector of population numbers/biomass
		N  = NA,
		#' @field N0     a scalar initial condition for numbers/biomass
		N0 = NA,
		#' @field N0Funk a function defining the equilibrium vigin numbers as a function of the productivity parameters
		N0Funk = NA,

		#TIME
		#' @field time   a vector of times for which to integrate the dynamics to
		time = NA,

		#MODEL
		#' @field lsdo    a scalar of the log scale residual sd
                lsdo = NA,
		#' @field lq      a scalar of the log proportionality constant between cpue and N
		lq = 0,
		#' @field model   a list of strings indicating the residual model, etc.
		model = list(observation="LN"),
                #prior = list(),
		
		#COMPUTATION
		#' @field ODE_method a string specifing the integration method to be used by dede
		ODE_method = 'rk4',
		#' @field OPT_method a string specifing the optimization method to be used by optim
		OPT_method = 'L-BFGS-B',
		
		#' @description Create a new instance of the ddModel class
                #'
                #' @param N0 a scalar initial condition for numbers
                #' @param time a vector of times for which to integrate the dynamics to
                #' @param dNdt a function that evaluates the delay differential equation that defines numbers/biomass dynamics.
                #' @param N0Funk a function defining the equilibrium vigin numbers as a function of the productivity parameters
                #'
                #' @return An instance of the prodModel class, for which data can be provided, parameters optimized and quantities plotted.
		initialize = function( 	N0   = NA,
					time = NA,
					dNdt = NA,
					N0Funk = NA,
					...
		){	
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			#dNdt
			stopifnot(is.function(dNdt))
			private$dNdt_classify(dNdt)	
			
			#N0
			if( is.function(N0Funk) ){ 	
				private$N0_classify(N0Funk)
				N0 = self$N0Funk()
			}
			
			#preallocate N
			self$N0 = N0
			self$time = time
			self$N  = matrix(NA, nrow=length(time), ncol=1)
			rownames(self$N) = sprintf("TIME %d", time)
		},
		
		#' @description A function to integrate the dynamics equations over the times in 'time'
                #'
                #' @param method an optional string specifing the integration method to be used by dede 
                #'
                #' @return NULL B and N are updated interally.
		iterate = function(method=self$ODE_method){
			#method	: optionally update the ODE method to be handed to the ode function (default 'rk4')
			#
			#value	: nothing returned, but self$N is updated with current values

			#prechecking 
			
			#digest and possibly change ode method	
			self$ODE_method = method	
	
			#N0
                        if( is.function(self$N0Funk) ){ N0=self$N0Funk() }
				
			#last minute allocation and variable updates
			self$N0 = N0
                        self$N = matrix(NA, nrow=length(self$time), ncol=1)
                        rownames(self$N) = sprintf("TIME %d", self$time)
			
			#solve 
        		capture.output( self$N <- deSolve::ode(self$N0, self$time, private$dNdt, parms=NULL, method=method)[,2], file="/dev/null" )
		},
		
		#' @description A function to optimize model parameter given the data provided
                #'
                #' @param data     A vector of data used to fit specified model.
                #' @param parNames A vector of strings matching the names of the parameters to be optimized
                #' @param lower    A vector of lower bounds for searching parameter space.
                #' @param upper    A vector of upper bounds for searching parameter space.
                #' @param method   A string optionally defining the method of optimization to be handed to optim (default 'L-BFGS-B')
                #' @param cov      A logical optionally indicating if hessian should be computed and inverted in optimization process
                #' @param fitQ     A logical indicating whether to find the MLE of log(q):lq via profile MLE, or if False don't change the initialized value of lq.
                #' @param gaBoost  A logical optionally (default F) indicating if a persistent genetic algorithm should be used to assist local optimization. Genetic algorithm repeates until first and second finite difference derivatives are successful and hessian is inverted. Optionally gaBoost may be given as a list containting names values for list(popSize, maxiter, run).
                #' @param persistFor A numeric indicating how many iterations of optimization tryCatch to engange in.
                #' @param control  Additional control parameters to be passed to optim
                #'
                #' @return Optimization objects are returned. Parameters values are updated inplace. rsCov is defined to self.
                optimize = optimize,
		
		#' @description A function to print
                #' @param ins ins 
                #' @param outs outs
                printer   = printSelf,
                #' @description The main function for printing this class
                #' @param ins ins
		printSelf = function(ins=c()){
			#NOTE: figure out a workaround for printing functions, here I just remove them
			self$printer(ins, outs=c(
				"iterate", "optimize", "model",	"prior", "like", 
				"plotMean", "plotBand", "plotRS", "N0Funk", 
				"save", "load", "printer"
			))
		},
		#' @description A function to plot the mean of the dynamics
                #' @param col   A color given however R colors are given
                #' @param alpha An alpha level between 0 and 255
                #' @param lwd   Line width
                #' @param add   A boolean indicating if the plot should be added to the existing device or not
		plotMean = plotMean,
		#' @description A function to plot uncertainty bands
                #' @param prob  Size of the uncertainty bands as defined as posterior probability
                #' @param col   A color given however R colors are given
                #' @param alpha An alpha level between 0 and 255
		plotBand = plotBand,
		#' @description  A function to plot repeated sampling posterior-like parameter distributions
                #' @param m      The number of samples from the repeated sampling distribution
                #' @param sample A boolean to indicate if samples should be returned
                #' @param save   A boolean to indicate if samples should be saved
		plotRS 	 = plotRS,
		
		#' @description    A function to save the ddModel class as an rds file.
                #' @param fileName A string defining the name and path of the file to be saved. It should probably end with the .rds extension.
		save = function(fileName){ saveRDS(self, file=fileName) },
		#' @description    A function to read in a ddModel class from .rds file.
                #' @param fileName A string of the name and path of the file to be loaded.
		load = function(fileName){ readRDS(fileName) },
		
		#' @description A function to calculate the likelihood of the model given the data
                #' @param data  A vector of the given index data.
		like = function(data){ sum(private$dLikes[[self$model$observation]](self, data)) }
	),
	
	#
	private = list(
		##
		#dNdt = NA,
		
		#
		selfToPar = selfToPar,	
		parToSelf = parToSelf, #NOTE: parValues should be passed with names
		
		#
		dNdt_classify = dNdt_classify,
		N0_classify = N0_classify,
		classify = classify,
		
		#
		dLikes = dLikes,
		qLikes = qLikes
	)
)




