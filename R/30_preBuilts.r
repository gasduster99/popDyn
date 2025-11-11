#
#PRODUCTION MODEL
#

#PT: exp(lalpha)*N/(gamma-1)*(1-(N/exp(lbeta))^(gamma-1))

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

#SCHNUTE

#
dNdtSchnute = function(t, P, lalpha, lbeta, gamma, M, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        if( gamma==-1 ){ Fmsy = sqrt(exp(lalpha)*M)-M
        } else {
                #
                FUpper = abs(exp(lalpha)-M)
                Fmsy = uniroot.all(function(ff){ a(gamma, ff, M)-exp(lalpha) }, c(0, FUpper))[1]
        }
        #
        R = exp(lalpha)*P*(1-exp(lbeta)*gamma*P)^(1/gamma)
        out = R - (M+Fmsy*C)*P
        #
        return( list(out) )
}

#
PBarSchnute = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M+

#
getBetaSchnute = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
bSchnute = Vectorize(getBetaSchnute, c("alpha", "gamma"))

#
getAlphaSchnute = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
aSchnute = Vectorize(getAlphaSchnute, "gamma")

#
getZetaSchnute = function(gamma, ff, M){
        #fgfm = ff*gamma/(ff+M)
        #fgfm/(1 + fgfm - (M/ff+M)^gamma)

        #
        (1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
}
zSchnute = Vectorize(getZetaSchnute, "gamma")

#
FMsySchnute = function(alpha, gamma, M){
        #
        FUpper = alpha-M
        root = rootSolve::uniroot.all(function(ff){ aSchnute(gamma, ff, M)-alpha }, c(0, FUpper))

        #
        return(root)
}

#EXPORT

#' @description Pre-built Schnute Production Model
#'
#' @field dNdt The ode has been initialized with the schnute model of production
#' @field N0Funk Vigin number/biomass for the schnute model 
#' @field defaultParams All other parameters are specified generically and should be replaced with actual species values. 
#'
#' @export
schnuteProdMod = prodModel$new(
        dNdt=dNdtSchnute, N0Funk=function(lalpha, lbeta, gamma, M){PBarSchnute(exp(lalpha), exp(lbeta), gamma, 0, M)},  #Dynamics
        time=1:80, catch=rep(1, 80), M=0.2,	#Constants
        lalpha=-1, lbeta=8, gamma=-1,	#Productivity Parameters
        lq=log(0.0005), #, lsdo=-2.1176758	#Nuisance Parameters
	alphaGiven=aSchnute, betaGiven=bSchnute,	#Productivity Given Functions
	zetaGiven=zSchnute, NBar=PBarSchnute, FMsy=FMsySchnute	#RP Functions
)
schnuteProdMod$iterate()

#
#DD MODEL
#

#Equilibrium Equations

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = Ryacas::ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = Ryacas::as_r( Ryacas::with_value(Bbar, "F", Ryacas::ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = Ryacas::ysym("F*B")
FBbar = Ryacas::with_value(FBbar, "B", Bbar)
FBbar = Ryacas::with_value(FBbar, "F", Ryacas::ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = Ryacas::as_r(dFBdF)
FDebug = function(FF, M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        eval(dFBdFexpr)
}
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10))$root
}

#beta does not matter for either of getAlphaFmsy or getGammaFmsy
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
        #
        #capture.output(
        out <- tryCatch({
                uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
        }, error=function(err){
                out = NA
        })#, file="/dev/null")
        #
        return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        uniroot(f, c(0, 10), tol=pracma::eps())$root
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
        }
        uniroot(f, c(-10, 10))$root
}
  
#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
        BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

#Differential Equation 

#
der = function(t, Y, lalpha, lbeta, gamma, aS, a0, WW, kappa, catch, B0){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        N = Y[1]
        B = Y[2]
        #
        if( (t-aS)<1){
                Blag = B0
        }else{
                Blag = deSolve::lagvalue(t-aS)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
        FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B
        #
        return( list(out) )
}

#SHINY

#
makeUI = function(lalpha, lbeta, gamma, M, kappa, aS){
	shinydashboard::dashboardPage(
        shinydashboard::dashboardHeader(),
        shinydashboard::dashboardSidebar(
                #Base Model NOTE: need to figure out how to update sliders to the new selected model
                #varSelectInput("who", "Base Model", whoMight, selected=startWho),

                ##Contrast
                #sliderInput("sliderChi", "chi", min=0, max=1, step=0.05, value=1),
                #Growth & Maturity
                shiny::sliderInput("sliderAS","As", min=0.1, max=10., step=0.1, value=aS), #schnuteDDMod$aS),
                shiny::sliderInput("sliderKappa","Kappa", min=0.1, max=5, step=0.1, value=kappa), #schnuteDDMod$kappa),
                #Recruitment
                shiny::sliderInput("sliderAlpha", "Alpha", min=M+0.1, max=exp(lalpha)*3, step=0.05, value=exp(lalpha)), #,min=schnuteDDMod$M+0.1, max=exp(schnuteDDMod$lalpha)*3, step=0.05, value=exp(schnuteDDMod$lalpha)),
                shiny::sliderInput("sliderBeta" , "Beta" , min=0, max=exp(lbeta)*1.5, step=10^(floor(log10(exp(lbeta)))-1), value=exp(lbeta)), #max=exp(schnuteDDMod$lbeta)*1.5, step=10^(floor(log10(exp(schnuteDDMod$lbeta)))-1), value=exp(schnuteDDMod$lbeta)),
                shiny::sliderInput("sliderGamma", "Gamma", min=-2, max=2, step=0.10001, value=gamma) #value=schnuteDDMod$gamma)

        ),
        shinydashboard::dashboardBody(
                shiny::fluidRow(shiny::column(12, shiny::plotOutput('rowOne'))),
                shiny::fluidRow(shiny::column(12, shiny::plotOutput('rowTwo'))),
                shiny::fluidRow(shiny::column(12, shiny::plotOutput('rowThree'))),
                shiny::fluidRow(shiny::column(12, shiny::plotOutput('rowFour'))),
                shiny::fluidRow(shiny::column(12, shiny::plotOutput('rowFive')))
        )
)
}

#
server = function(input, output, session){
        #
        reactiveDat = shiny::reactive({

                ##
                #if(input$who!=dat$whoAmI){
                #       whoNew = sprintf("%s/%s", path, as.character(input$who))
                #       dat = readRDS(whoNew)
                #       dat$whoAmI = whoNew
                #}

                ##Contrast
                #con = input$sliderChi
                #Growth
                schnuteDDMod$aS = input$sliderAS
                schnuteDDMod$kappa = input$sliderKappa
                #Recruitment
                schnuteDDMod$alpha = input$sliderAlpha
                schnuteDDMod$lalpha= log(input$sliderAlpha)
                schnuteDDMod$beta  = input$sliderBeta
                schnuteDDMod$lbeta = log(input$sliderBeta)
                schnuteDDMod$gamma = input$sliderGamma

                #
                ww = vbGrow(schnuteDDMod$aS, schnuteDDMod$kappa, schnuteDDMod$WW, schnuteDDMod$a0) #WW*(1-exp(-kappa*a0))
                schnuteDDMod$FMsy = FMsy(schnuteDDMod$M, schnuteDDMod$kappa, ww, schnuteDDMod$WW, exp(schnuteDDMod$lalpha), exp(schnuteDDMod$lbeta), schnuteDDMod$gamma)
                schnuteDDMod$xi   = rbind(schnuteDDMod$xi, schnuteDDMod$FMsy/schnuteDDMod$M)
                schnuteDDMod$zeta = rbind(schnuteDDMod$zeta, getZeta(schnuteDDMod$FMsy, schnuteDDMod$M, schnuteDDMod$kappa, ww, schnuteDDMod$WW, exp(schnuteDDMod$lalpha), exp(schnuteDDMod$lbeta), schnuteDDMod$gamma))
                #dat$catch = fContrast(con)
                #
                schnuteDDMod$iterate()

                #return
                schnuteDDMod
        })
        #
        output$rowOne = shiny::renderPlot({
                #
                dat = reactiveDat()
                layout(t(1:2))
                dat$plotQuan( function(B){B}, main="Biomass", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass")
                dat$plotQuan( function(N){N}, main="Numbers", ylim=c(0,max(dat$N)), xlab="Time", ylab="Numbers")
        })
        #
        output$rowTwo = shiny::renderPlot({
                #
                dat = reactiveDat()
                #
                layout(t(1:2))
                ##
                #curve(SRR(x, dat), 0, 3*dat$B0, lwd=3, xlab="Biomass", ylab="Recruitment", main="Stock-Recruitment", n=1000)
                #abline(0, dat$M, col='red')
                #
                ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
                BMsy = BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
                g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
                g = Vectorize(g, 'x')
                maxF = uniroot(g, c(dat$FMsy, exp(dat$lalpha)))$root
                #print(maxF)
                FFs = seq(0, maxF, length.out=1000)
                #f = function(x){surplus(x, dat)/surplus(BMsy, dat)*BMsy*dat$FMsy}
                #curve(f(x), 0, dat$B0, lwd=3, xlab="Biomass", ylab="Equilibrium Surplus Biomass", main="Yield Curve", n=1000)
                plot(g(FFs), FFs*g(FFs), type='l', lwd=3, xlab="Biomass", ylab="Equilibrium Surplus Biomass", main="Yield Curve")
                segments(BMsy, 0, BMsy, BMsy*dat$FMsy)
                points(BMsy, BMsy*dat$FMsy, pch=19)
                #abline(0, dat$FMsy)
                #curve(BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma))
                #rug( BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma), lwd=3 )
                #rug( BBar(0, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma), lwd=3 )
                #abline(h=BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)*dat$FMsy)
                #
                curve(vbGrow(x, dat$kappa, dat$WW, dat$a0), 0, 15, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, dat$WW))
                segments(dat$aS, 0, dat$aS, ww)
                segments(0, ww, dat$aS, ww)
                points(dat$aS, ww, pch=19)
        })
        #
        output$rowThree = shiny::renderPlot({
                #
                dat = reactiveDat()
                #layout(cbind(c(1,2), 3))
                #
                #ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #WW*(1-exp(-kappa*a0))
                #dat$FMsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
                #
                #par(mar=c(4,5,2,3))
                #dat$plotQuan( function(catch, FMsy){catch*FMsy}, main="Fishing", ylim=c(0,1.5), xlab="Time", ylab="F")
                #dat$plotQuan( function(B, catch, FMsy){B*catch*FMsy}, ylim=c(0,max(dat$B*dat$catch*dat$FMsy)*1.1), xlab="Time", ylab="Catch")
                #par(mar=c(5,5,5,3))
                #
                layout(t(1:2))
                #
                ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
                #
                curve(SRR(x, dat), 0, 3*dat$B0, lwd=3, xlab="Biomass", ylab="Recruitment #s", main="Stock-Recruitment", n=1000)
                abline(0, dat$M*(dat$M+dat$kappa)/dat$kappa/dat$WW/(1+dat$M*ww/dat$kappa/dat$WW), col='red')
                #abline(v=dat$B0)
                segments(dat$B0, 0, dat$B0, SRR(dat$B0, dat))
                #segments(0, SRR(dat$B0, dat), dat$B0, SRR(dat$B0, dat))
                points(dat$B0, SRR(dat$B0, dat), pch=19)
                #
                dat$plotQuan( function(B, N){B/N}, main="Average Size", ylim=c(0,max(dat$B/dat$N)), xlab="Time", ylab="Biomass Per Individual")
        })
        #
        output$rowFour = shiny::renderPlot({
                #
                dat = reactiveDat()
                #
                layout(t(1:2))
		#NOTE: Rethink how catch is parameterized. Inside both models.
                dat$plotQuan( function(catch, FMsy){catch*FMsy}, main="Fishing", ylim=c(0,max(dat$catch*dat$FMsy)), xlab="Time", ylab="F")
                dat$plotQuan( function(B, catch, FMsy){B*catch*FMsy}, ylim=c(0,max(dat$B*dat$catch*dat$FMsy)*1.1), xlab="Time", ylab="Biomass", main="Catch")
        })
        #
        output$rowFive = shiny::renderPlot({
                #
                dat = reactiveDat()
                #
                howManyRP = length(dat$xi)
                howManyGrey = min(howManyRP, 50)
                greys = rev(81-round(pracma::logseq(80, 1, howManyGrey))) #seq(60, 2, -1)
                nWhite = max(howManyRP-howManyGrey+1, 0)
                #
                plot(dat$xi, dat$zeta, pch=19, col=c(rep("white", nWhite), sprintf("grey%d",greys), "black"), xlab="Fmsy/M", ylab="Bmsy/B0", main="Reference Points")
                points(dat$xi[howManyRP], dat$zeta[howManyRP], pch=19, col='black')
                points(dat$xi[howManyRP], dat$zeta[howManyRP], col='red')
        })
}

#
SRR = function(B, mod){
        #
        alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
        gamma = mod$gamma
        #
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
        #
        return(R)
}
SRR = Vectorize(SRR, 'B')

#EXPORT

#
aS = 0.1                #2 #10 #0.1
a0 = -1         #-0.25 #-0.5 #-1   #-2
M  = 0.2
kappa = 1      #10 #0.1
WW = 1
ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
#
B0 = 10000
#
aMin = M*(M+kappa)/kappa/WW/(1+M*ww/kappa/WW)
#
lalpha=-1; gamma=-1
lbetaBH = log(getBeta(B0, M, kappa, ww, WW, exp(lalpha), gamma))

#' @description Pre-built Schnute Production Model
#'
#' @field derivs The dde has been initialized with the schnute model of production
#' @field N0Funk Vigin Numbers for the schnute model 
#' @field B0Funk Vigin Biomass for the schnute model
#' @field defaultParams All other parameters are specified generically and should be replaced with actual species values. 
#'
#' @export
schnuteDDMod = ddModel$new( derivs=der,
        N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin numbers
                #
                alpha = exp(lalpha)
                beta  = exp(lbeta)
                ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
                #
                BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
                (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M

        },
        B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin biomass
                #
                alpha = exp(lalpha)
                beta = exp(lbeta)
                ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
                #
                BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
        },
        time=1:80, catch=rep(1, 80), aS=aS, a0=a0, M=M, WW=WW, kappa=kappa,	#constants
        lalpha=lalpha, lbeta=lbetaBH, gamma=-1,	#parameters
        lq=log(0.00049),	#nuisance parameters
        #other incidentals to carry along
	UIFunk=makeUI
)
schnuteDDMod$iterate()


#schnuteDDMod$launchShinyApp = function(host="0.0.0.0", port=5050){
#	shiny::runApp(shiny::shinyApp(ui, server), host=host, port=port)
#}

