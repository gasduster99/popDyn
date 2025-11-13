
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fishKitLearn

To install the latest version of the `fishKitLearn` package you can run
the following in an R shell:
`install.packages("devtools") #if required devtools::install_github("gasduster99/fishKitLearn")`

## Simple Production Model

### Automatic Schnute prodModel Provided

`fishKitLearn` ships with several automatically confingured models. The
`prodModel` class is a platform for specifying generic production
models, as will be described in the next section about manually
configuring models. That said, a Schnute production model is also
automatically configured and provided by default.

This model includes three parameter Schnute production[^1], natural
mortality, and fishing mortality as follows.

$$
\begin{align*}
&I_t = q B_t e^\epsilon ~~~ \epsilon\sim N(0, \sigma^2)\\
&\frac{dB}{dt} = P_{S}(B;[\alpha, \beta]) -(M+F)B\\~\\
&P_{S}(B;[\alpha, \beta, \gamma]) = \alpha B(1-\beta\gamma B)^{\frac{1}{\gamma}}
\end{align*}
$$

Logistic (Schaefer) Model `\gamma=1`, Ricker Model `\gamma\rightarrow0`,
Beverton-Holt `\gamma=-1`

\`\` \#A default Schnute model configuration is provided in the package
schnuteProdModel\$printSelf()

\#Plot schnuteProdModel\$plotMean()

\#Update with your own data cpue = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90,
0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53,
0.58, 0.64, 0.66, 0.65, 0.63) catch = c(94, 212, 195, 383, 320, 402,
366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211,
231, 223) TT = length(catch) \# schnuteProdModel$time = 1:TT
schnuteProdModel$catch = catch

\#optimize parameters opt = schnuteProdModel\$optimize(cpue, c(‘lsdo’ ,
‘lalpha’, ‘lbeta’ ), lower = c(log(0.001), log(0.1), log(1000) ), upper
= c(log(0.3) , 0 , log(10000)), cov = T )

\#Update Plot schnuteProdModel$plotMean(add=T, col='blue')
schnuteProdModel$plotBand(col=‘blue’) \`\`

### Example of Manual Schaefer Model Instantiation

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

### Automatic Schnute-Deriso Delay Differential Model [^2]

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

    schnuteDDModel$shiny()

<figure>
<img
src="https://github.com/gasduster99/fishKitLearn/blob/main/shiny.png?raw=true"
alt="Example Shiny App" />
<figcaption aria-hidden="true">Example Shiny App</figcaption>
</figure>

<!--
### Manual BH ddModel configuration
&#10;```
```
-->

## Bibliography

[^1]: [cite](link.pdf)

[^2]: [Walters, The Continuous Time Schnute-Deriso Delaydifference Model
    for Age-Structured Population Dynamics, with Example Application to
    the Peru Anchoveta
    Stock.](https://fisheries-2023.sites.olt.ubc.ca/files/2020/06/1Continuous-time-Schnute-Deriso-model-Final.pdf)
