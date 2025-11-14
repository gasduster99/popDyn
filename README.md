
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fishKitLearn

To install the latest version of the `fishKitLearn` package you can run
the following in an R shell:

    install.packages("devtools") #if required
    devtools::install_github("gasduster99/fishKitLearn")

## Simple Production Model

### Automatic Schnute `prodModel`

`fishKitLearn` ships with several automatically confingured models. The
`prodModel` class is a platform for specifying generic production
models, as will be described in the next section about manually
configuring models. That said, a Schnute production model is also
automatically configured and provided by default.

This model uses a three parameter Schnute production[^1], natural
mortality, and fishing mortality as follows:

$$
\begin{align*}
&I_t = q B_t e^\epsilon ~~~ \epsilon\sim N(0, \sigma^2)\\
&\frac{dB}{dt} = P_{S}(B;[\alpha, \beta]) -(M+F)B\\
&P_{S}(B;[\alpha, \beta, \gamma]) = \alpha B(1-\beta\gamma B)^{\frac{1}{\gamma}}
\end{align*}
$$

This is a versatile model that allows for the specification of many
common models by fixing $\gamma$. The BH and Logistic (Schaefer) Models
arise when $\gamma$ is fixed to -1 or 1 respectively. The Ricker model
is a limiting case as $\gamma\rightarrow0$. For $\gamma<-1$ a family of
strictly increasing Cushing-like curves arise, culminating in linear
production as $\gamma\rightarrow-\infty$. More Information about this
model and how it estimates reference points see[^2].

This Schnute production model is initialized to the `schnuteProdMod`
default instance of the `prodModel` class. You can see the default state
of the model by running the following. This code prints the model state
to your R shell and plots the mean response value.

    #A default Schnute model configuration is provided in the package
    schnuteProdMod$printSelf()

    #Plot
    schnuteProdMod$plotMean()

To parameterize this model to your species of interest, you must add
your own catch data and time settings. At which point, the model is
ready to optimize whichever parameters you would like from the given
index data.

    #Update with your own data, here is Nimibian Hake for example.
    cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
    catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
    TT = length(catch)
    #
    schnuteProdMod$time = 1:TT
    schnuteProdMod$catch = catch

    #optimize parameters
    opt = schnuteProdMod$optimize(cpue,
                  c('lsdo'    , 'lalpha', 'lbeta'   ),
            lower   = c(log(0.001), log(0.1), log(1000) ),
            upper   = c(log(0.3)  , 0       , log(10000)),
            cov     = T
    )

    #
    schnuteProdMod$printSelf()
    #Update Plot
    schnuteProdMod$plotMean(add=T, col='blue')
    schnuteProdMod$plotBand(col='blue')

At this time it is only possible to fit to a single index. However you
can optimize any parameter that you can see in the `schnuteProdMod`
instance by simply adding its name to the optimization and adding
bounds.

### Example of Manual Schaefer Model Instantiation

To manually specify a production model, you must specify a function for
integrating the ODE that defines your model. Further you should provide
a function that defines the virgin population state in terms of the
model parameters or constants. Below you can see a typical Schaefer
model specified using the `prodModel` class.

    #define Schaefer model ODE 
    #log productivity parameterization improves optimization.
    dNdtSchaefer = function(t, N, lalpha, lbeta, catch){
            #linearly interpolate catches
            ft = floor(t)
            q  = (t-ft)
            Cl = catch[ft]
            Cu = catch[ft+1]
            C = q*Cu + (1-q)*Cl
            if(q==0){ C=Cl }
            #
            R = exp(lalpha)*N*(1-(N/exp(lbeta)))
            out = R - C 
            #
            return( list(out) )
    }

    #initialize prodModel class
    schaeferModel = prodModel$new(
        dNdt=dNdtSchaefer, N0Funk=function(lbeta){1/exp(lbeta)},    #Dynamics 
        time=1:TT, catch=catch, #Constants
        lalpha=-1, lbeta=8, #Productivity Parameters
        lq=log(0.0005), lsdo=-2.1176758 #Nuisance Parameters
    )

    #optimize
    optS = schaeferModel$optimize(cpue,
            c('lsdo', 'lalpha', 'lbeta'),
            lower   = c(log(0.001), log(0.1), log(1000)),
            upper   = c(log(0.3), 0, log(10000)),
            gaBoost = T,
            cov     = T,
            fitQ    = T
    )

    #plot
    plot(schaeferModel$time, cpue)
    fitPT$plotMean(add=T, col="blue")
    fitPT$plotBand(col="blue")

## Delay Differential Models

More complex models can be fit in `fishKitLearn` by using the `ddModel`
class that is designed to work with delay differential models (DDM).
Manual specification of DDMs follow the same general motif as shown
above, except that lag values and lagged derivatives are specified with
the `lagvalue` and `lagderiv` functions as would be used with `dede`
from the `deSolve` package.

### Automatic Schnute-Deriso Delay Differential Model

The Schnute-Deriso DDM as specified by Walters[^3] is automatically
provided here in the `schnuteDDMod` object. In this package the model is
constructed with the three parameter Schnute Stock Recruitment
relationship so as to allow the same versatile access to Ricker, BH and
Logistic recruitment model through $\gamma$. For more information about
this model and how it estimates reference points see
[\[2\]](https://escholarship.org/uc/item/1th4n7kd).
<!--[Grunloh, N. (2024) A Metamodeling Approach for Bias Estimation of Biological Reference Points. (Doctoral dissertation, University of California Santa Cruz).](https://escholarship.org/uc/item/1th4n7kd)].-->

    #A default Schnute model configuration is provided in the package
    schnuteDDModel$printSelf()

    #Plot
    schnuteDDModel$plotMean()

    #Optimize Parameters
    optDD = schnuteDDModel$optimize(cpue,
                  c('lsdo'    , 'lalpha', 'lbeta'   ),
            lower   = c(log(0.001), log(0.1), log(1000) ),
            upper   = c(log(0.3)  , 0       , log(10000)),
            cov     = T
    )

    #Update Plot
    schnuteDDModel$plotMean(add=T, col='blue')
    schnuteDDModel$plotBand(col='blue')

#### Shiny App

A shiny app is provided here to visualize the most common metrics one
may be interested for this model. The shiny app allows the user to
perform real-time sensativities of all metrics with respect to the model
parameters.

    schnuteDDModel$launchShinyApp()

![](https://github.com/gasduster99/fishKitLearn/blob/main/shiny.png?raw=true)

<!--
### Manual BH ddModel configuration
&#10;```
```
-->

## Bibliography

[^1]: Jon Schnute. A General Theory for Analysis of Catch and Effort
    Data. Canadian Journal of Fisheries and Aquatic Sciences,
    42(3):414–429, March 1985.

[^2]: [Grunloh, N. (2024) A Metamodeling Approach for Bias Estimation of
    Biological Reference Points. (Doctoral dissertation, University of
    California Santa Cruz).](https://escholarship.org/uc/item/1th4n7kd)

[^3]: [Walters, The Continuous Time Schnute-Deriso Delaydifference Model
    for Age-Structured Population Dynamics, with Example Application to
    the Peru Anchoveta
    Stock.](https://fisheries-2023.sites.olt.ubc.ca/files/2020/06/1Continuous-time-Schnute-Deriso-model-Final.pdf)
