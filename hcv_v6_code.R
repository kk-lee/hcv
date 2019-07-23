
# load packages####
library(tableone)
library(caret)
library(reshape2)
library(tidyr)
library(dplyr)
library(binom)
library(survival)
library(gridExtra)
library(grid)
library(survMisc)
library(scales)
library(ggplot2)
library(tidyr)
library(metafor)
library(xlsx)
library(xlsxjars)
library(rJava)
library(readxl)
library(stringr)
library(reshape)
library(data.table)
library(plyr)
library(RColorBrewer)
library(rworldmap)
library(classInt)
library(gdata)
library(forestplot)
library(Hmisc)
library(data.table)

# load beta pharm function ####

beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.2.2 (December 2012)
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Belisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@clinepi.mcgill.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        k <- min(last.theta/change)
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}

###====================================
### Part 1 -  estimate pooled rate and risk ratio
###====================================

# load files

dat.ratio <- read.csv("data_raw/hcv.ratio.v5.csv",header = T,na.strings = "")

baseline <- dat.ratio %>% 
  distinct(Study.id, .keep_all=TRUE) %>% 
  select(Study.id, Pub.Year, Author, N, male.percentage, Cohort.country, N.HCVpositive) %>% 
  mutate(n.male=as.integer((male.percentage/100)*N)) %>% 
  arrange(desc(Pub.Year), Author) %>% 
  mutate(n.percentage.male=paste0(n.male, " (", male.percentage, ")"))

events <- dat.ratio %>% 
  filter(Overall.estimate==1) %>% 
  select(Study.id, Pub.Year, Author, n.events) %>% 
  group_by(Study.id) %>% 
  dplyr::summarise(sum.events=sum(n.events))

baseline <- left_join(baseline, events)

write.csv(baseline, "data_derived/table1_n.men.csv")

baseline <- baseline %>% 
  mutate(us.taiwan=case_when(
    Cohort.country %in% c("Taiwan", "USA")~"us.taiwan",
    TRUE~"rest"
  )) %>% 
  group_by(us.taiwan) %>% 
  dplyr::summarise(sum.events=sum(N))


##============================================
## calculate meta-estimates for risk ratio
##============================================

names(dat.ratio) <- tolower(names(dat.ratio))
dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$outcome))
dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))
dat.ratio$author <- paste0(dat.ratio$author," ","et al")
dat.ratio[dat.ratio$study.id==783,"author"]<- "DAD cohort"

dat.ratio$lgrr <- logb(dat.ratio$ratio.ce)
dat.ratio$lgll <- logb(dat.ratio$ratio.ll)
dat.ratio$lgul <- logb(dat.ratio$ratio.ul)
dat.ratio$se <- (dat.ratio$lgul-dat.ratio$lgll)/(2*1.96)

num_iter <- 10000


dat.ratio <- dat.ratio %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
            funs(format(round(as.numeric(.),digits=2), nsmall=2)))

dat.ratio$rr <- paste0(as.character(dat.ratio$ratio.ce)," [",as.factor(dat.ratio$ratio.ll),"-",as.factor(dat.ratio$ratio.ul),"]")


#### Plot by endpoint ####

dat.ratio.strat <- dat.ratio %>% 
  filter(outcome!=0)

rr <- rma(yi = dat.ratio.strat$lgrr,sei = dat.ratio.strat$se, method = "ML", slab=dat.ratio.strat$author) 
rr.exp <- predict(rr, transf = exp, digits=2)

rr.sim <- exp(rnorm(num_iter, as.vector(rr$b), rr$se))
(quantile(rr.sim, probs = c(0.025,0.5, 0.975)))
hist((rr.sim))

pdf("Figures/HCV Strat estimate.v2.pdf", height=10, width=12)

par(mar=c(4,4,1,4), font=1)

forest(rr, xlim=c(-4, 2), at=log(c(0.5,1,2,3,3.5)), atransf=exp,
       ilab=cbind(dat.ratio.strat$pub.year, dat.ratio.strat$n, dat.ratio.strat$n.events, dat.ratio.strat$rr),
       ilab.xpos=c(-2.9,-2.4, -1.9, -1.3),cex=0.75, ylim=c(-1, 65),
       order=rev(order(dat.ratio.strat$outcome,dat.ratio.strat$pub.year)),
       rows=c(3:15,21:33,39:61),
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")


### add text with Q-value, dfs, p-value, and I^2 statistic
text(-4, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (Q = ",
                                           .(formatC(rr$QE, digits=2, format="f")), ", df = ", .(rr$k - rr$p),
                                           "; ", I^2, " = ",.(formatC(rr$I2, digits=1, format="f")), "%"," [95% CI ",.(formatC(confint(rr)[[1]][[7]], digits=1, format="f")),
                                           " - ",
                                           .(formatC(confint(rr)[[1]][[11]], digits=1, format="f")),
                                           "])")))

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(62,34,16), pos=4, c("Cardiovascular events",
                               "Myocardial infarction",
                               "Stroke"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,64.5, "Author(s)",  pos=4)
text(-2.9,64.5,"Year")
text(-2.4,64.5,"Total number \n of participants")
text(-1.9,64.5,"Number of \n events")
text(-1.3,64.5, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


### fit random-effects model in the three subgroups
res.c <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome=="CVD",],method = "ML") 
res.m <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome=="MI",],method = "ML") 
res.s <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome=="Stroke",],method = "ML") 

res.c.exp <- predict(res.c, transf = exp, digits=2)
res.m.exp <- predict(res.m, transf = exp, digits=2)
res.s.exp <- predict(res.s, transf = exp, digits=2)


### add summary polygons for the three subgroups
addpoly(res.c, row=37.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(res.m, row= 19.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(res.s, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")

### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups
text(-4, 37.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(res.c$k - res.c$p),
                                             "; ", I^2, " = ",.(formatC(res.c$I2, digits=1, format="f")), "%), ","Pooled risk ratio = ", .(formatC(res.c.exp$pred, digits=2, format="f")), " [95% CI ",.(formatC(res.c.exp$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(res.c.exp$ci.ub, digits=2, format="f")),
                                             "]")))
text(-4, 19.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(res.m$k - res.m$p),
                                             "; ", I^2, " = ",.(formatC(res.m$I2, digits=1, format="f")), "%), ","Pooled risk ratio = ", .(formatC(res.m.exp$pred, digits=2, format="f")), " [95% CI ",.(formatC(res.m.exp$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(res.m.exp$ci.ub, digits=2, format="f")),
                                             "]")))
text(-4, 1.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(res.s$k - res.s$p),
                                            "; ", I^2, " = ",.(formatC(res.s$I2, digits=1, format="f")), "%), ","Pooled risk ratio = ", .(formatC(res.s.exp$pred, digits=2, format="f")), " [95% CI ",.(formatC(res.s.exp$ci.lb, digits=2, format="f")),
                                            " - ",
                                            .(formatC(res.s.exp$ci.ub, digits=2, format="f")),
                                            "]")))

dev.off()

#### Calculate overall ratio ####
dat.ratio.overall <- dat.ratio %>% filter(overall.estimate==1)

rr <- rma(yi = dat.ratio.overall$lgrr,sei = dat.ratio.overall$se, method = "ML", slab=dat.ratio.overall$author) 
rr.exp <- predict(rr, transf = exp, digits=2)

rr.sim <- exp(rnorm(num_iter, as.vector(rr$b), rr$se))
(quantile(rr.sim, probs = c(0.025,0.5, 0.975)))
hist((rr.sim))


#### Plot for overall ratio ####
### decrease margins so the full space is used
pdf("Figures/HCV Overall estimate.pdf", height=8, width=12)

par(mar=c(4,4,1,4), font=1)

forest(rr, xlim=c(-3, 2), at=log(c(0.5,1,2,3,3.5)), atransf=exp,
       ilab=cbind(dat.ratio.overall$pub.year,dat.ratio.overall$rr),
       ilab.xpos=c(-2,-1.3),cex=0.75, ylim=c(-1, 50),
       order=rev(order(dat.ratio.overall$outcome,dat.ratio.overall$pub.year)),
       rows=c(1:47),
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")


### add text with Q-value, dfs, p-value, and I^2 statistic

text(-3, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr$k - rr$p),
                                           "; ", I^2, " = ",.(formatC(rr$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.exp$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.exp$ci.lb, digits=2, format="f")),
                                           " - ",
                                           .(formatC(rr.exp$ci.ub, digits=2, format="f")),
                                           "]")))

### set font expansion factor (as in forest() above) and use bold italic

### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)


### switch to bold font
par(font=2)

### add column headings to the plot
text(-3,49, "Author(s)",  pos=4)
text(-2,49,"Year")
text(-1.3,49, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)

dev.off()

#### publication bias ####
regtest(rr)
predict(trimfill(rr))
funnel(rr)
funnel(trimfill(rr))

#heterogeneity
confint(rr)


#### subgroup analysis ####

rr.m <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$hospitalisation.mortality=="mortality",],method = "ML")
rr.mort <- predict(rr.m, transf = exp, digit=2)
rr.mort <- as.data.frame(rr.mort)
rr.mort$subgroup <- "Mortality"
rr.mort$i2 <- round(rr.m$I2, digits=1)
rr.mort$estimates <- rr.m$k

rr.h <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$hospitalisation.mortality!="mortality",],method = "ML")
rr.hosp <- predict(rr.h, transf = exp, digit=2)
rr.hosp <- as.data.frame(rr.hosp)
rr.hosp$subgroup <- "Hospitalisation"
rr.hosp$i2 <- round(rr.h$I2, digits=1)
rr.hosp$estimates <- rr.h$k

rr.co <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$hiv.hcv.coinfection.1.0.==1,],method = "ML") 
rr.coinf <- predict(rr.co, transf = exp, digit=2)
rr.coinf <- as.data.frame(rr.coinf)
rr.coinf$subgroup <- "HIV co-infected"
rr.coinf$i2 <- round(rr.co$I2, digits=1)
rr.coinf$estimates <- rr.co$k

rr.mo <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$hiv.hcv.coinfection.1.0.==0,],method = "ML") 
rr.mono <- predict(rr.mo, transf = exp, digit=2)
rr.mono <- as.data.frame(rr.mono)
rr.mono$subgroup <- "HIV mono-infected"
rr.mono$i2 <- round(rr.mo$I2, digits=1)
rr.mono$estimates <- rr.mo$k

rr.pr <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$pub.year<=2013,],method = "ML") 
rr.pre <- predict(rr.pr, transf = exp, digit=2)
rr.pre <- as.data.frame(rr.pre)
rr.pre$subgroup <- "pre-2014"
rr.pre$i2 <- round(rr.pr$I2, digits=1)
rr.pre$estimates <- rr.pr$k

rr.po <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$pub.year>2013,],method = "ML") 
rr.post <- predict(rr.po, transf = exp, digit=2)
rr.post <- as.data.frame(rr.post)
rr.post$subgroup <- "post-2014"
rr.post$i2 <- round(rr.po$I2, digits=1)
rr.post$estimates <- rr.po$k

rr.lb <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$rob==0,],method = "ML") 
rr.lbias <- predict(rr.lb, transf = exp, digit=2)
rr.lbias <- as.data.frame(rr.lbias)
rr.lbias$subgroup <- "low bias"
rr.lbias$i2 <- round(rr.lb$I2, digits=1)
rr.lbias$estimates <- rr.lb$k

rr.hb <- rma(yi = lgrr,sei = se, data = dat.ratio[!dat.ratio$rob==0,],method = "ML") 
rr.hbias <- predict(rr.hb, transf = exp, digit=2)
rr.hbias <- as.data.frame(rr.hbias)
rr.hbias$subgroup <- "high bias"
rr.hbias$i2 <- round(rr.hb$I2, digits=1)
rr.hbias$estimates <- rr.hb$k

dat.ratio.icd <- dat.ratio %>% 
  filter(stringr::str_detect(definition.outcome, 'ICD'))

rr.ic <- rma(yi = lgrr,sei = se, data = dat.ratio.icd,method = "ML") 
rr.icd <- predict(rr.ic, transf = exp, digit=2)
rr.icd <- as.data.frame(rr.icd)
rr.icd$subgroup <- "ICD diagnostic code"
rr.icd$i2 <- round(rr.ic$I2, digits=1)
rr.icd$estimates <- rr.ic$k

dat.ratio.phy <- dat.ratio %>% 
  filter(definition.outcome=="physician diagnosis")

rr.ph <- rma(yi = lgrr,sei = se, data = dat.ratio.phy,method = "ML") 
rr.phy <- predict(rr.ph, transf = exp, digit=2)
rr.phy <- as.data.frame(rr.phy)
rr.phy$subgroup <- "physician diagnosis"
rr.phy$i2 <- round(rr.ph$I2, digits=1)
rr.phy$estimates <- rr.ph$k

subgroups <- rbind (rr.hosp, rr.mort, rr.mono, rr.coinf, rr.pre, rr.post, rr.lbias, rr.hbias, rr.icd, rr.phy)

subgroups <- subgroups[,c(6,8,1:3,7)]
subgroups[,c(3:5)] <- round(subgroups[,c(3:5)], digits=2)

write.csv(subgroups, "data_derived/eTable 1.csv")


## sensitivity analysis for studies with viremia

viremia <- dat.ratio %>% 
  filter(hcv.status=="viraemic")

rr.vi <- rma(yi = lgrr,sei = se, data = viremia, method = "ML") 
rr.vir <- predict(rr.vi, transf = exp, digit=2)
rr.vir <- as.data.frame(rr.vir)

## subgroup analysis by risk ratio and odds ratio
rr.or <- rma(yi = lgrr,sei = se, data = dat.ratio.overall[dat.ratio.overall$risk.measure=="OR",], method = "ML") 
rr.or2 <- predict(rr.or, transf = exp, digit=2)

rr.rr <- rma(yi = lgrr,sei = se, data = dat.ratio.overall[dat.ratio.overall$risk.measure!="OR",], method = "ML") 
rr.rr2 <- predict(rr.rr, transf = exp, digit=2)

## subgroup analysis by risk ratio and odds ratio
dat.ratio.overall <- dat.ratio.overall %>% 
  mutate(us.taiwan=case_when(
    cohort.country %in% c("Taiwan", "USA")~"us.taiwan",
    TRUE~"rest"
  ))

rr.us.taiwan <- rma(yi = lgrr,sei = se, data = dat.ratio.overall[dat.ratio.overall$us.taiwan=="us.taiwan",], method = "ML")
rr.us.taiwan2 <- predict(rr.us.taiwan, transf = exp, digit=2)

rr.rest <- rma(yi = lgrr,sei = se, data = dat.ratio.overall[dat.ratio.overall$us.taiwan=="rest",], method = "ML")
rr.rest2 <- predict(rr.rest, transf = exp, digit=2)


#__________________________________________________________________________


###=========================================================
### Part 2: estimating burden using daly and prevalence data
###=========================================================

##==============================
## Part 2a - sort IHME CVD data 
##==============================

# load country data into R with ISO country code

country.iso <- read.csv("data_raw/country.isocode.csv",header = T)
country.iso$countryid <- 1:247

# combine ISO country code and GBD regions
country.gbd.regions <- read.csv("data_raw/ihme.gbd.regions.csv", header=T)
combined.iso.region <- merge(x = country.gbd.regions, y = country.iso, by.x ="location_name", by.y = "country",all.x=T, all.y=T)
names(combined.iso.region) <- c("country", "gbd.region", "country.iso", "countryid")

# load daly data - country specific but in one file

dat <- read.csv("data_raw/IHME-GBD_2015_DATA.csv",header = T)
dat$location_name <- as.character(dat$location_name)

#match country names to UN and relabel unmatched names 

combined <- merge(x = combined.iso.region,y = dat,by.x = "country",by.y ="location_name" ,all.y=T)
str(combined)
combined[,c("country","gbd.region", "country.iso")] <- lapply(combined[,c("country","gbd.region", "country.iso")], as.character)
names(table(combined$country.iso))

combined$country.iso <- ifelse(combined$country=="Global","global",combined$country.iso)
combined$gbd.region <- ifelse(combined$country=="Global","global",combined$gbd.region)

combined <- subset(combined, !is.na(combined$country.iso))
names(table(combined$country))

combined.countrylevel <- combined


# file for all patients

names(combined) <- c("country", "gbd.region", "country.iso", "countryid", 
                     "measure_id", "measure_name", "location_id", "sex_id", "sex_name", 
                     "age_id", "age_name", "cause_id", "cause_name", "metric_id", 
                     "metric_name", "year_id", "nm_mean", "nm_upper", "nm_lower")

# file for all patients
x <- function(x){rbind(x, c("total",names(table(y$gbd.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.all <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Both"),c("country", "gbd.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.all <- split(combined.all,list(combined.all$gbd.region,combined.all$year_id))
combined.all <- lapply(combined.all,function(x){rbind(x, c("total",names(table(x$gbd.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})

combined.all <- rbind.fill((combined.all))

total.both <- combined.all[combined.all$country=="total",]

# file for male patients
x <- function(x){rbind(x, c("total",names(table(y$gbd.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.male <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Male"),c("country", "gbd.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.male <- split(combined.male,list(combined.male$gbd.region,combined.male$year_id))
combined.male <- lapply(combined.male,function(x){rbind(x, c("total.male",names(table(x$gbd.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})
combined.male <- rbind.fill((combined.male))
total.males <- combined.male[combined.male$country=="total.male",]


# file for female patients
x <- function(x){rbind(x, c("total",names(table(y$gbd.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.female <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Female"),c("country", "gbd.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.female <- split(combined.female,list(combined.female$gbd.region,combined.female$year_id))
combined.female <- lapply(combined.female,function(x){rbind(x, c("total.female",names(table(x$gbd.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})
combined.female <- rbind.fill((combined.female))
total.females <- combined.female[combined.female$country=="total.female",]

# file for all patients, stratified by gender

total.all <- as.data.frame(rbind(total.both,total.males,total.females))
total.all <- total.all[,c("country", "gbd.region", "year_id", "nm_mean", "nm_lower", "nm_upper")]

total.all$compound.key <- paste0(total.all$country,total.all$gbd.region)
total.all$nm_mean <- as.numeric(total.all$nm_mean)
total.all$nm_lower <- as.numeric(total.all$nm_lower)
total.all$nm_upper <- as.numeric(total.all$nm_upper)


# daly.ce - estimating central daly estimate startified by sex, year and region

daly.ce <-total.all[,c(3,4,7)] 

daly.ce <- split(daly.ce, daly.ce$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_mean"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ce <- lapply(daly.ce,fun1)
daly.ce <- do.call(rbind.data.frame, daly.ce)
daly.ce$V3 <- row.names(daly.ce)
daly.ce$V3 <- gsub('[[:digit:]]+', '', daly.ce$V3)
daly.ce$V3 <- str_sub(daly.ce$V3,-nchar(daly.ce$V3),-2)
names(daly.ce) <- c("year","daly.ce","subgroup")

# daly.ll - estimating daly.ll estimate startified by sex, year and region

daly.ll <-total.all[,c(3,5,7)] 

daly.ll <- split(daly.ll, daly.ll$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_lower"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ll <- lapply(daly.ll,fun1)
daly.ll <- do.call(rbind.data.frame, daly.ll)
daly.ll$V3 <- row.names(daly.ll)
daly.ll$V3 <- gsub('[[:digit:]]+', '', daly.ll$V3)
daly.ll$V3 <- str_sub(daly.ll$V3,-nchar(daly.ll$V3),-2)
names(daly.ll) <- c("year","daly.ll","subgroup")

# daly.ul - estimating daly.ul estimate startified by sex, year and region

daly.ul <-total.all[,c(3,6,7)] 

daly.ul <- split(daly.ul, daly.ul$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_upper"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ul <- lapply(daly.ul,fun1)
daly.ul <- do.call(rbind.data.frame, daly.ul)
daly.ul$V3 <- row.names(daly.ul)
daly.ul$V3 <- gsub('[[:digit:]]+', '', daly.ul$V3)
daly.ul$V3 <- str_sub(daly.ul$V3,-nchar(daly.ul$V3),-2)
names(daly.ul) <- c("year","daly.ul","subgroup")

# combined daly

daly <- as.data.frame(cbind(daly.ce,daly.ll,daly.ul))

daly <- daly[,c("year", "daly.ce", "daly.ll", "daly.ul", "subgroup")]


rownames(daly) <- 1:dim(daly)[1]

daly$area <- gsub("[.]",'',daly$subgroup)
daly$area <- gsub("totalmale|totalfemale|total",'',daly$area)

pattern1=c("total")
pattern2=c("totalmale")
pattern3=c("totalfemale")

daly$sex1 <- gsub("[.]",'',daly$subgroup)
daly$sex <- str_extract(daly$sex1,paste(rev(mget(ls(pattern = "pattern\\d+"))), collapse="|"))


daly$sex <- factor(daly$sex,
                   levels=c("total","totalmale","totalfemale"),
                   labels=c("adults","male","female"))

daly$area.id <- factor(daly$area,
                       levels=c("global", "Andean Latin America", "Australasia", "Caribbean", "Central Asia", 
                                "Central Europe", "Central Latin America", "Central Sub-Saharan Africa", 
                                "East Asia", "Eastern Europe", "Eastern Sub-Saharan Africa", 
                                "High-income Asia Pacific", "High-income North America", 
                                "North Africa and Middle East", "Oceania", "South Asia", "Southeast Asia", 
                                "Southern Latin America", "Southern Sub-Saharan Africa", "Tropical Latin America", 
                                "Western Europe", "Western Sub-Saharan Africa"),
                       labels=c("All countries", "Andean Latin America", "Australasia", "Caribbean", "Central Asia", 
                                "Central Europe", "Central Latin America", "Central Sub-Saharan Africa", 
                                "East Asia", "Eastern Europe", "Eastern Sub-Saharan Africa", 
                                "High-income Asia Pacific", "High-income North America", 
                                "North Africa and Middle East", "Oceania", "South Asia", "Southeast Asia", 
                                "Southern Latin America", "Southern Sub-Saharan Africa", "Tropical Latin America", 
                                "Western Europe", "Western Sub-Saharan Africa"))

##=================================
## Part 2b - load prevalence files
##=================================

hcv.prev <- read.csv("data_raw/prev.gbd.region.csv",header = T,stringsAsFactors = F)

##=========================================
## Part 2c - merge prevelance data to daly
##=========================================

final <- merge(x=hcv.prev,y=daly,by.x = c("subgroup","gbd.region","year"), by.y=c("sex","area.id","year"))

## derive alpha and beta using beta pharm function for prevalence data 

# convert percentage into decimals

final[,c("ce","ll","ul")] <- lapply(final[,c("ce","ll","ul")],function(x){x/100})

alpha <- function(x){
  alpha <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[5]),as.numeric(x[6])))[1])
}

beta <- function(x){
  beta <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[5]),as.numeric(x[6])))[2])
}

final$alpha <- as.vector(unlist(apply(final, 1,alpha)))
final$beta <- as.vector(unlist(apply(final, 1,beta)))

# check beta results
final$uniq_id <- seq(1000, 1000+nrow(final)-1, 1)

check_res <- by(final, final$uniq_id, function(x) qbeta(shape1 = x$alpha, shape2= x$beta, 
                                                        p = c(0.025, 0.5, 0.975)))
check_res <- do.call(rbind, check_res)
colnames(check_res) <- c("est", "lci", "uci")
check_res <- cbind(final, check_res)
plot(check_res$ll, check_res$lci)
rm(check_res)
final$uniq_id <- NULL

a <- as.list(NULL)
num_iter <- 10000

## derive central estimate and 95% uncertainty term for CVD daly using log distribution

final[,c("daly.ce","daly.ll","daly.ul")] <- lapply(final[,c("daly.ce","daly.ll","daly.ul")],function(x){log(x,base = exp(1))})


final$nm_se <- (as.numeric(final$daly.ul)-as.numeric(final$daly.ll))/3.92
final$nm_mean <- as.numeric(final$daly.ce)

# check if log distribution appropriate

final2 <- final
final2[, c('ce', 'll', 'ul')] <- lapply(final2[, c('ce', 'll', 'ul')],log) 
final2$upper_diff <- final2$ul - final2$ce
final2$lower_diff <- final2$ce - final2$ll
final2$diff_diff <- final2$upper_diff - final2$lower_diff
hist(final2$diff_diff, breaks = 100)


###=======================================================
### Part 3 = calculate paf and corresponding CVD burden due to HCV 
###=======================================================
newvars <- data.frame(id = 1:nrow(final), 
                      burdence = NA,
                      burdenll = NA,
                      burdenul = NA,
                      pafce= NA,
                      pafll = NA,
                      paful = NA)
newvars$id <- NULL


final <- cbind(final, newvars)
for(i in 1:nrow(final)){
  x <- rbeta(num_iter, shape1 = final$alpha[i], shape2 = final$beta[i])
  y <- rnorm(num_iter,final$nm_mean[i],final$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  
  
  final$burdence[i] <- quantile(burden,0.5)
  final$burdenll[i] <- quantile(burden,0.025)
  final$burdenul[i] <- quantile(burden,0.975)
  final$pafce[i] <- quantile(p,0.5)
  final$pafll[i] <- quantile(p,0.025)
  final$paful[i] <- quantile(p,0.975)
}


# in thousands

final$burdence <-final$burdence/1000
final$burdenll <-final$burdenll/1000
final$burdenul <-final$burdenul/1000

##==============================================
##                  COUNTRY LEVEL
##==============================================

##sort prevalence by country, age, gender ####

prev.males <- read.csv("data_raw/polaris.males.csv",header = T,na.strings = "")
prev.females <- read.csv("data_raw/polaris.females.csv", header=T, na.strings = "")
prev.both<- read.csv("data_raw/polaris.both.csv", header=T, na.strings = "")


country <- colnames(prev.males[, !names(prev.males) %in% c("group", "age.cat")]) 
prev.males <- melt(prev.males, id.vars= c("group", "age.cat"), measure.vars= country)
prev.males <- dcast(prev.males, age.cat+variable~group, value="value")
xm <- filter(prev.males, is.na(prev.males$"Total Population"))
ym <- filter(prev.males, is.na(prev.males$"Viremic Population"))
zm <- cbind(xm,ym)
prev.males <- zm[,c(1,2,4,7)]
prev.males$gender <- rep("male", nrow(prev.males))

country <- colnames(prev.females[, !names(prev.females) %in% c("group", "age.cat")]) 
prev.females <- melt(prev.females, id.vars= c("group", "age.cat"), measure.vars= country)
prev.females <- dcast(prev.females, age.cat+variable~group, value="value")
xf <- filter(prev.females, is.na(prev.females$"Total Population"))
yf <- filter(prev.females, is.na(prev.females$"Viremic Population"))
zf <- cbind(xf,yf)
prev.females <- zf[,c(1,2,4,7)]
prev.females$gender <- rep("female", nrow(prev.females))

country <- colnames(prev.both[, !names(prev.both) %in% c("group", "age.cat")]) 
prev.both <- melt(prev.both, id.vars= c("group", "age.cat"), measure.vars= country)
prev.both <- dcast(prev.both, age.cat+variable~group)
xb <- filter(prev.both, is.na(prev.both$"Total Population"))
yb <- filter(prev.both, is.na(prev.both$"Viremic Population"))
zb <- cbind(xb,yb)
pb <- zb[,c(1,2,4,7)]
wb <- prev.both[complete.cases(prev.both),]
wb <- wb[,c(1,2,4,3)]
prev.both <- rbind(wb,pb)
prev.both$gender <- rep("both", nrow(prev.both))

prev.hcv <- rbind(prev.males, prev.females, prev.both)
prev.hcv <- prev.hcv[,c(2,1,5,3,4)]
colnames(prev.hcv)[colnames(prev.hcv)=="variable"] <- "country"


prev.hcv$age.cat <- revalue(prev.hcv$age.cat, c("  0-4 "="0-4", "  10-14 "="10-14", "  15-19 "="15-19", "  20-24 "="20-24", "  25-29 "="25-29", "  30-34 "="30-34", 
                                                        "  35-39 "="35-39", "  40-44 "="40-44", "  45-49 "="45-49", "  5-9 "="5-9", "  50-54 "="50-54", "  55-59 "="55-59", 
                                                        "  60-64 "="60-64", "  65-69 "="65-69", "  70-74 "="70-74", "  75-79 "="75-79", "  80-84 "="80-84", "  85+ "="85+", 
                                                        " 0-4"="0-4", " 10-14"="10-14", " 15-19"="15-19", " 20-24"="20-24", " 25-29"="25-29", " 30-34"="30-34", " 35-39"="35-39", 
                                                        " 40-44"="40-44", " 45-49"="45-49", " 5-9"="5-9", " 50-54"="50-54", " 55-59"="55-59", " 60-64"="60-64", " 65-69"="65-69", 
                                                        " 70-74"="70-74", " 75-79"="75-79", " 80-84"="80-84", " 85+"="85+", " Total "="Total", "Total"="Total", "Totals"="Total"))
prev.hcv[,4] <- as.numeric(gsub(",","",prev.hcv[,4]))
prev.hcv[,5] <- as.numeric(gsub(",","",prev.hcv[,5]))
prev.hcv[,c(1,2)]<- sapply(prev.hcv[,c(1,2)], as.character)

#combine 80-84 and 85+
prev.hcv.over80 <- filter(prev.hcv, prev.hcv$age.cat=="80-84"| prev.hcv$age.cat=="85+")
prev.hcv.over80 <- split(prev.hcv.over80,list(prev.hcv.over80$country,prev.hcv.over80$gender))
totalover80 <- function(z){rbind(z, c(names(table(z$country)),"80+",names(table(z$gender)),colSums(z[4:5])))}
prev.hcv.over80 <- lapply(prev.hcv.over80, totalover80)
prev.hcv.over80 <- rbind.fill((prev.hcv.over80))
prev.hcv.over80 <- prev.hcv.over80[prev.hcv.over80$age.cat=="80+",]

#country level prevalence
prev.hcv.under80 <- filter(prev.hcv, !(prev.hcv$age.cat=="80-84"|prev.hcv$age.cat=="85+"))
prev.hcv.final <- rbind(prev.hcv.under80,prev.hcv.over80)
prev.hcv.final[,1] <- as.factor(prev.hcv.final[,1])

prev.hcv.final$country <- revalue(prev.hcv.final$country, c("Burkina.Faso"="Burkina Faso", "Central.African.Republic"="Central African Republic", "Czech.Republic"="Czech Republic",
                                                            "Dominican.Republic"="Dominican Republic", "Hong.Kong"="Hong Kong", "New.Zealand"="New Zealand", "Papua.New.Guinea"="Papua New Guinea",
                                                            "Puerto.Rico"="Puerto Rico",  "Saudi.Arabia"="Saudi Arabia", "South.Africa"="South Africa", "UAE"= "United Arab Emirates", "Viet.Nam"="Vietnam",
                                                            "US"="United States", "UK"="United Kingdom", "Samoa"="American Samoa", "Korea"="South Korea", "Gambia"="The Gambia"))   
prev.hcv.final <- filter(prev.hcv.final, !prev.hcv.final$age.cat=="Total")

nodata <- c('Hong Kong', 'Guadeloupe')

prev.hcv.final <- prev.hcv.final[ !grepl(paste(nodata, collapse="|"), prev.hcv.final$country),]

prev.hcv.final <- filter(prev.hcv.final, !(prev.hcv.final$age.cat=="0-4"| prev.hcv.final$age.cat=="5-9"|prev.hcv.final$age.cat=="10-14"|prev.hcv.final$age.cat=="15-19"))

## sort cvd dalys
countries.wprev <- c("Afghanistan", 
                     "Albania", "Algeria", "Argentina", "Australia", "Austria", "Azerbaijan", 
                     "Bahrain", "Belgium", "Brazil", "Bulgaria", "Burkina Faso", "Burundi", 
                     "Cambodia", "Cameroon", "Canada", "Central African Republic", 
                     "Chad", "Chile", "China", "Colombia", "Croatia", "Cuba", "Czech Republic", 
                     "Denmark", "Dominican Republic", "Egypt", "Estonia", "Ethiopia", 
                     "Fiji", "Finland", "France", "Gabon", "The Gambia", "Georgia", "Germany", 
                     "Ghana", "Greece", "Guadeloupe", "Hong Kong", "Hungary", "Iceland", 
                     "India", "Indonesia", "Iran", "Iraq", "Ireland", "Israel", "Italy", 
                     "Japan", "Jordan", "Kazakhstan", "Kenya", "South Korea", "Latvia", 
                     "Lebanon", "Libya", "Lithuania", "Luxembourg", "Madagascar", 
                     "Malaysia", "Malta", "Mexico", "Mongolia", "Morocco", "Netherlands", 
                     "New Zealand", "Nigeria", "Norway", "Oman", "Pakistan", "Panama", 
                     "Papua New Guinea", "Peru", "Philippines", "Poland", "Portugal", 
                     "Puerto Rico", "Qatar", "Romania", "Russia", "American Samoa", "Saudi Arabia", 
                     "Slovakia", "Slovenia", "South Africa", "Spain", "Sweden", "Switzerland", 
                     "Syria", "Taiwan", "Thailand", "Tunisia", "Turkey", "United Arab Emirates", "United Kingdom", 
                     "Ukraine", "United States", "Uzbekistan", "Venezuela", "Vietnam", "Yemen")
countries.wprev <- paste(countries.wprev, collapse="|")  
  
# file for all patients
combined.pc <- filter(combined, combined$year_id==2015)
combined.pc <- combined.pc[(!combined.pc$age_name=="Age-standardized"),c("country", "gbd.region","cause_name","sex_name","age_name","nm_mean", "nm_lower", "nm_upper")]
combined.pc[,c("nm_mean", "nm_lower", "nm_upper")] <- sapply(combined.pc[,c("nm_mean", "nm_lower", "nm_upper")], as.numeric)
combined.pc[,c("cause_name","sex_name","age_name")] <- sapply(combined.pc[,c("cause_name","sex_name","age_name")], as.character)
combined.pc$age_name <- revalue(combined.pc$age_name, c("20 to 24"="20-24", "25 to 29"="25-29", "30 to 34"="30-34", "35 to 39"="35-39", "40 to 44"="40-44", "45 to 49"="45-49", "50 to 54"="50-54", 
                                "55 to 59"="55-59", "60 to 64"="60-64", "65 to 69"="65-69", "70 to 74"="70-74", "75 to 79"="75-79", "80 plus"="80+"))
combined.pc <- split(combined.pc,list(combined.pc$country,combined.pc$sex_name, combined.pc$age_name))
totalcvd <- function(z){rbind(z, c(names(table(z$country)),names(table(z$gbd.region)),"CVD",names(table(z$sex_name)),names(table(z$age_name)),colSums(z[6:8])))}
combined.pc <- lapply(combined.pc, totalcvd)
combined.pc <- rbind.fill((combined.pc))
combined.pc <- combined.pc[combined.pc$cause_name=="CVD",]

combined.pc <- combined.pc[grepl(countries.wprev,combined.pc$country)==T,]
combined.pc$country <- as.factor(combined.pc$country)

combined.pc$sex_name <- revalue(combined.pc$sex_name, c("Male"="male", "Female"="female", "Both"="both"))

combined.pc <- setnames(combined.pc, c("sex_name","age_name", "nm_mean", "nm_lower", "nm_upper"), c("gender", "age.cat", "daly.ce", "daly.ll", "daly.ul"))
combined.pc <- combined.pc[c(2,1,5,4,6,7,8)]

##combine country level prevalence and dalys
countryleveldata <- merge(combined.pc, prev.hcv.final, by=c("country","age.cat","gender"), all=T)

countryleveldata2 <- countryleveldata

## Calculate sex and age strat PAF and burden
countryleveldata$alpha <- as.numeric(countryleveldata$'Viremic Population')
countryleveldata$beta <- as.numeric(countryleveldata$'Total Population')-as.numeric(countryleveldata$'Viremic Population')
countryleveldata <- transform(countryleveldata, daly.ce = as.numeric(daly.ce), 
          daly.ll = as.numeric(daly.ll), daly.ul=as.numeric(daly.ul))

tab4 <- countryleveldata

num_iter <- 10000

countryleveldata.total <- countryleveldata

countryleveldata[,c("daly.ce","daly.ul","daly.ll")] <- lapply(countryleveldata[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})

countryleveldata$nm_se <- (as.numeric(countryleveldata$daly.ul)-as.numeric(countryleveldata$daly.ll))/3.92
countryleveldata$nm_mean <- as.numeric(countryleveldata$daly.ce)

for(i in 1:nrow(countryleveldata)){
  x <- rbeta(num_iter, shape1 = countryleveldata$alpha[i], shape2 = countryleveldata$beta[i])
  y <- rnorm(num_iter,countryleveldata$nm_mean[i],countryleveldata$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  countryleveldata$burdence[i] <- quantile(burden,0.5,na.rm = T)
  countryleveldata$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  countryleveldata$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  countryleveldata$pafce[i] <- quantile(p,0.5,na.rm = T)
  countryleveldata$pafll[i] <- quantile(p,0.025,na.rm = T)
  countryleveldata$paful[i] <- quantile(p,0.975,na.rm = T)
}

## total age and sex stratified Burden by country ####
countryleveldata$'Total.Population' <- as.numeric(countryleveldata$'Total.Population')
countryleveldata$'Viremic.Population' <- as.numeric(countryleveldata$'Viremic.Population')
countryleveldata$country <- as.character(countryleveldata$country)
countryleveldata <- countryleveldata[!countryleveldata$gender=="both",]
countryleveldata <- split(countryleveldata,list(countryleveldata$country))
totalburden <- function(z){rbind(z, c(names(table(z$country)),names(table(z$gbd.region)),rep(NA,2),colSums(z[5:19])))}
countryleveldata <- lapply(countryleveldata, totalburden)
countryleveldata <- rbind.fill((countryleveldata))


countryleveldata$burdence <-as.numeric(countryleveldata$burdence)/1000
countryleveldata$burdenll <-as.numeric(countryleveldata$burdenll)/1000
countryleveldata$burdenul <-as.numeric(countryleveldata$burdenul)/1000

ukraine <- countryleveldata[countryleveldata$country=="Ukraine",]

countryleveldata <- countryleveldata[is.na(countryleveldata$gender),]
countryleveldata <- countryleveldata[,c(1,2,9,14:19)]

## calculate total PAF by country
countryleveldata.total$'Total.Population' <- as.numeric(countryleveldata.total$'Total.Population')
countryleveldata.total$'Viremic.Population' <- as.numeric(countryleveldata.total$'Viremic.Population')
countryleveldata.total$country <- as.character(countryleveldata.total$country)
countryleveldata.total <- countryleveldata.total[!countryleveldata.total$gender=="both",]
countryleveldata.total <- split(countryleveldata.total,list(countryleveldata.total$country))
totalburden.total <- function(z){rbind(z, c(names(table(z$country)),names(table(z$gbd.region)),rep(NA,2),colSums(z[5:11])))}
countryleveldata.total <- lapply(countryleveldata.total, totalburden.total)
countryleveldata.total <- rbind.fill((countryleveldata.total))

countryleveldata.total <- countryleveldata.total[is.na(countryleveldata.total$gender),]
countryleveldata.total[,c(5:11)] <- lapply(countryleveldata.total[,c(5:11)], as.numeric)


countryleveldata.total[,c("daly.ce","daly.ul","daly.ll")] <- lapply(countryleveldata.total[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})

countryleveldata.total$nm_se <- (as.numeric(countryleveldata.total$daly.ul)-as.numeric(countryleveldata.total$daly.ll))/3.92
countryleveldata.total$nm_mean <- as.numeric(countryleveldata.total$daly.ce)

for(i in 1:nrow(countryleveldata.total)){
  x <- rbeta(num_iter, shape1 = countryleveldata.total$alpha[i], shape2 = countryleveldata.total$beta[i])
  y <- rnorm(num_iter,countryleveldata.total$nm_mean[i],countryleveldata.total$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  countryleveldata.total$burdence[i] <- quantile(burden,0.5,na.rm = T)
  countryleveldata.total$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  countryleveldata.total$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  countryleveldata.total$pafce[i] <- quantile(p,0.5,na.rm = T)
  countryleveldata.total$pafll[i] <- quantile(p,0.025,na.rm = T)
  countryleveldata.total$paful[i] <- quantile(p,0.975,na.rm = T)
}

countryleveldata.total[,c("pafce","paful","pafll")] <- lapply(countryleveldata.total[,c("pafce","paful","pafll")],function(x){x*100})

countryleveldata.total[,c("pafce", "pafll", "paful")] <- lapply(countryleveldata.total[,c("pafce", "pafll", "paful")], function(x){round(as.numeric(x),digits=2)})
countryleveldata.total$paf <-paste0(countryleveldata.total$pafce," ","(",countryleveldata.total$pafll,"-",countryleveldata.total$paful,")") 

write.csv(countryleveldata.total, "data_derived/countrylevelpaf.csv")

## Plot cartogram for country level data ####
countryleveldata$burdenC <- (as.numeric(countryleveldata$burdence)/as.numeric(countryleveldata$Total.Population))*100000000
countryleveldata.table <- countryleveldata

countryleveldata$burdenCcat <- cut2(countryleveldata$burdenC,g=10)

countryleveldata <- merge(x = countryleveldata, y = country.iso, by.x ="country", by.y = "country",all.x=T)

countryData<-countryleveldata
data(countryData)
sPDF<-joinCountryData2Map(countryData, joinCode= "ISO3", nameJoinColumn= "country.iso")
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

colfunc<- colorRampPalette(c("#FFF7EC", "#FDD49E", "#fc8d59", "#D7301F", "#7F0000"))
levels(sPDF@data[["cat"]]) <- c('<3.16', '3.17 to 5.80', '5.81 to 7.70', '7.71 to 11.07', '11.08 to 15.79', '15.80 to 21.02', '21.03 to 30.55', '30.56 to 46.84', '46.85 to 84.36', '84.37 to 585.89')
palette<- colfunc(1000)
mapParams <-mapCountryData(sPDF, nameColumnToPlot="burdence" ,numCats= 10, colourPalette = palette, missingCountryCol='dark grey', catMethod="logFixedWidth", addLegend=FALSE, mapTitle=' ')
do.call (addMapLegend, c(mapParams, legendLabels="all", legendWidth=0.5, legendIntervals="data",legendMar= 2))


countryData<-countryleveldata
data(countryData)
countryData
sPDF<-joinCountryData2Map(countryData, joinCode= "ISO3", nameJoinColumn= "country.iso")
mapCountryData(sPDF, nameColumnToPlot= "burdenCcat")
colourPalette<-brewer.pal(10, "OrRd")
levels(sPDF@data[["burdenCcat"]]) <- c('<3.16', '3.17 to 5.80', '5.81 to 7.70', '7.71 to 11.07', '11.08 to 15.79', '15.80 to 21.02', '21.03 to 30.55', '30.56 to 46.84', '46.85 to 84.36', '84.37 to 585.89')
par(mai=c(0,0,1,0),xaxs="i",yaxs="i")
mapParams<-mapCountryData(sPDF, nameColumnToPlot="burdenCcat", mapTitle=' ', colourPalette=colourPalette, missingCountryCol='dark grey')



## Plot cartogram for regional data ####
regionaldata <- merge(x= combined.iso.region, y= final, by="gbd.region")
regionaldata$pafp <- as.numeric(regionaldata$pafce)*100000
regionaldata[regionaldata$gbd.region=="Central Asia", "pafp"] <- 700
countryData<-regionaldata
data(countryData)
sPDF<-joinCountryData2Map(countryData, joinCode= "ISO3", nameJoinColumn= "country.iso")
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

colfunc<- colorRampPalette(c("#FFF7EC", "#FDD49E", "#fc8d59", "#D7301F", "#7F0000"))
palette<- colfunc(1000)
mapParams <-mapCountryData(sPDF, nameColumnToPlot="pafp", numCats=1000, colourPalette = palette, missingCountryCol='dark grey',catMethod=c(0:700), addLegend=FALSE, mapTitle=' ')
do.call (addMapLegend, c(mapParams, legendWidth=0.5, legendMar= 2))


##tables
final[,c("burdence", "burdenll", "burdenul")] <- lapply(final[,c("burdence", "burdenll", "burdenul")], function(x){round(as.numeric(x),digits=1)})
final$burden <-paste0(final$burdence," ","(",final$burdenll,"-",final$burdenul,")") 

final[,c("pafce", "pafll", "paful")] <- lapply(final[,c("pafce", "pafll", "paful")], function(x){round(as.numeric(x)*100,digits=2)})
final$paf <-paste0(final$pafce," ","(",final$pafll,"-",final$paful,")") 

table1 <- final[,c(2,23,24)]
write.csv(table1, "data_derived/regions.csv")

countryleveldata.table

countryleveldata.table$burdenC.ll <- (as.numeric(countryleveldata.table$burdenll)/as.numeric(countryleveldata$Total.Population))*100000000
countryleveldata.table$burdenC.ul <- (as.numeric(countryleveldata.table$burdenul)/as.numeric(countryleveldata$Total.Population))*100000000

countryleveldata.table[,c("burdence", "burdenll", "burdenul")] <- lapply(countryleveldata.table[,c("burdence", "burdenll", "burdenul")], function(x){round(as.numeric(x),digits=1)})
countryleveldata.table$burden <-paste0(countryleveldata.table$burdence," ","(",countryleveldata.table$burdenll,"-",countryleveldata.table$burdenul,")") 

countryleveldata.table[,c("burdenC", "burdenC.ll", "burdenC.ul")] <- lapply(countryleveldata.table[,c("burdenC", "burdenC.ll", "burdenC.ul")], function(x){round(as.numeric(x),digits=1)})
countryleveldata.table$burden.C <-paste0(countryleveldata.table$burdenC," ","(",countryleveldata.table$burdenC.ll,"-",countryleveldata.table$burdenC.ul,")") 

countryleveldata.table[,c("pafce", "pafll", "paful")] <- lapply(countryleveldata.table[,c("pafce", "pafll", "paful")], function(x){round(as.numeric(x),digits=2)})
countryleveldata.table$paf <-paste0(countryleveldata.table$pafce," ","(",countryleveldata.table$pafll,"-",countryleveldata.table$paful,")") 

countryleveldata.table.final <- countryleveldata.table[,c(1,3,13:15)]
write.csv(countryleveldata.table.final, "data_derived/burden.bycountry.csv")



## Calculate PAF by age strata ####
cld2 <- countryleveldata2

#both genders
cld2.both <- subset(cld2, cld2$gender=="both")
cld2.both[,c(5:9)] <- lapply(cld2.both[,c(5:9)], as.numeric)
cld2.both$country <- as.character(cld2.both$country)

cld2.both <- split(cld2.both,list(cld2.both$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,2),colSums(z[5:9])))}
cld2.both <- lapply(cld2.both, totalburden.total)
cld2.both <- rbind.fill((cld2.both))
cld2.both <- cld2.both[is.na(cld2.both$gender),]

cld2.both[,c(5:9)] <- lapply(cld2.both[,c(5:9)], as.numeric)
cld2.both[,c("daly.ce","daly.ul","daly.ll")] <- lapply(cld2.both[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})
cld2.both$nm_se <- (as.numeric(cld2.both$daly.ul)-as.numeric(cld2.both$daly.ll))/3.92
cld2.both$nm_mean <- as.numeric(cld2.both$daly.ce)
cld2.both$alpha <- as.numeric(cld2.both$'Viremic Population')
cld2.both$beta <- as.numeric(cld2.both$'Total Population')-as.numeric(cld2.both$'Viremic Population')

for(i in 1:nrow(cld2.both)){
  x <- rbeta(num_iter, shape1 = cld2.both$alpha[i], shape2 = cld2.both$beta[i])
  y <- rnorm(num_iter,cld2.both$nm_mean[i],cld2.both$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  cld2.both$burdence[i] <- quantile(burden,0.5,na.rm = T)
  cld2.both$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  cld2.both$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  cld2.both$pafce[i] <- quantile(p,0.5,na.rm = T)
  cld2.both$pafll[i] <- quantile(p,0.025,na.rm = T)
  cld2.both$paful[i] <- quantile(p,0.975,na.rm = T)
}


#female
cld2.female <- subset(cld2, cld2$gender=="female")
cld2.female[,c(5:9)] <- lapply(cld2.female[,c(5:9)], as.numeric)
cld2.female$country <- as.character(cld2.female$country)

cld2.female <- split(cld2.female,list(cld2.female$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,2),colSums(z[5:9])))}
cld2.female <- lapply(cld2.female, totalburden.total)
cld2.female <- rbind.fill((cld2.female))
cld2.female <- cld2.female[is.na(cld2.female$gender),]

cld2.female[,c(5:9)] <- lapply(cld2.female[,c(5:9)], as.numeric)
cld2.female[,c("daly.ce","daly.ul","daly.ll")] <- lapply(cld2.female[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})
cld2.female$nm_se <- (as.numeric(cld2.female$daly.ul)-as.numeric(cld2.female$daly.ll))/3.92
cld2.female$nm_mean <- as.numeric(cld2.female$daly.ce)
cld2.female$alpha <- as.numeric(cld2.female$'Viremic Population')
cld2.female$beta <- as.numeric(cld2.female$'Total Population')-as.numeric(cld2.female$'Viremic Population')

for(i in 1:nrow(cld2.female)){
  x <- rbeta(num_iter, shape1 = cld2.female$alpha[i], shape2 = cld2.female$beta[i])
  y <- rnorm(num_iter,cld2.female$nm_mean[i],cld2.female$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  cld2.female$burdence[i] <- quantile(burden,0.5,na.rm = T)
  cld2.female$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  cld2.female$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  cld2.female$pafce[i] <- quantile(p,0.5,na.rm = T)
  cld2.female$pafll[i] <- quantile(p,0.025,na.rm = T)
  cld2.female$paful[i] <- quantile(p,0.975,na.rm = T)
}

#male
cld2.male <- subset(cld2, cld2$gender=="male")
cld2.male[,c(5:9)] <- lapply(cld2.male[,c(5:9)], as.numeric)
cld2.male$country <- as.character(cld2.male$country)

cld2.male <- split(cld2.male,list(cld2.male$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,2),colSums(z[5:9])))}
cld2.male <- lapply(cld2.male, totalburden.total)
cld2.male <- rbind.fill((cld2.male))
cld2.male <- cld2.male[is.na(cld2.male$gender),]

cld2.male[,c(5:9)] <- lapply(cld2.male[,c(5:9)], as.numeric)
cld2.male[,c("daly.ce","daly.ul","daly.ll")] <- lapply(cld2.male[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})
cld2.male$nm_se <- (as.numeric(cld2.male$daly.ul)-as.numeric(cld2.male$daly.ll))/3.92
cld2.male$nm_mean <- as.numeric(cld2.male$daly.ce)
cld2.male$alpha <- as.numeric(cld2.male$'Viremic Population')
cld2.male$beta <- as.numeric(cld2.male$'Total Population')-as.numeric(cld2.male$'Viremic Population')

for(i in 1:nrow(cld2.male)){
  x <- rbeta(num_iter, shape1 = cld2.male$alpha[i], shape2 = cld2.male$beta[i])
  y <- rnorm(num_iter,cld2.male$nm_mean[i],cld2.male$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  cld2.male$burdence[i] <- quantile(burden,0.5,na.rm = T)
  cld2.male$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  cld2.male$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  cld2.male$pafce[i] <- quantile(p,0.5,na.rm = T)
  cld2.male$pafll[i] <- quantile(p,0.025,na.rm = T)
  cld2.male$paful[i] <- quantile(p,0.975,na.rm = T)
}

#graph of PAF by age group ####
cld2.both[,17:19] <- cld2.both[,17:19]*100
cld2.female[,17:19] <- cld2.female[,17:19]*100
cld2.male[,17:19] <- cld2.male[,17:19]*100

cld2.both$age <- seq(from=20, to=80, by=5)

dat <-  cld2.both
p5 <- ggplot(data <- dat) 
p5 <- p5+
  geom_line(aes(y=dat$pafce, x=dat$age), colour = "blue")+
  geom_line(aes(y=dat$pafll, x=dat$age), colour = "red",linetype="dashed")+
  geom_line(aes(y=dat$paful, x=dat$age), colour = "red",linetype="dashed")+
  scale_x_continuous(name="\nAge (years)")+
  scale_y_continuous(name="Population attributable fraction, %\n", limits=c(0,0.8),breaks=seq(0,0.8,0.2),labels=c("0","0.2","0.4","0.6","0.8"))+
  theme(legend.position="none")+
  theme_light()

p5


dat <-  cld2.both
p5 <- ggplot(data <- dat) 
p5 <- p5+
  stat_smooth(aes(y=dat$pafce, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "blue", size=0.5)+
  stat_smooth(aes(y=dat$pafll, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "red",linetype="dashed", size=0.5)+
  stat_smooth(aes(y=dat$paful, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "red",linetype="dashed", size=0.5)+
  scale_x_continuous(name="\nAge (years)", breaks=seq(20,80,10))+
  scale_y_continuous(name="Population attributable fraction, %\n", limits=c(0,0.8),breaks=seq(0,0.8,0.2),labels=c("0","0.2","0.4","0.6","0.8"))+
  theme(legend.position="none")+
  theme_light()

p5


cld2.male$gender <- "male"
cld2.male$age <- seq(from=20, to=80, by=5)
cld2.female$gender <- "female"
cld2.female$age <- seq(from=20, to=80, by=5)
cld2.gender <- rbind(cld2.male, cld2.female)
cld2.gender$gender <- as.factor(cld2.gender$gender)

dat <-  cld2.gender

p6 <- ggplot(data <- dat, aes(y=dat$pafce, x=dat$age, ymin=dat$pafll, ymax=dat$paful))+
  geom_ribbon(aes(fill=dat$gender), alpha = 0.2)+
  geom_line(aes(color=dat$gender))+
  scale_x_continuous(name="\nAge (year)")+
  scale_y_continuous(name="Population attributable fraction, %\n", limits=c(0,0.8),breaks=seq(0,0.8,0.2),labels=c("0","0.2","0.4","0.6","0.8"))+
  theme_light()+theme(legend.position='none')

p6


smooth <- data.frame(Age= seq( min(cld2.male$age), max(cld2.male$age), length.out=100))
smooth$mEst <- predict(loess(pafce ~ age, cld2.male, span=0.5), newdata=smooth$Age)
smooth$mLB  <-  predict(loess(pafll ~ age, cld2.male, span=0.5), newdata=smooth$Age)
smooth$mUB  <-  predict(loess(paful ~ age, cld2.male, span=0.5), newdata=smooth$Age)

smooth$fEst <- predict(loess(pafce ~ age, cld2.female, span=0.5), newdata=smooth$Age)
smooth$fLB  <-  predict(loess(pafll ~ age, cld2.female, span=0.5), newdata=smooth$Age)
smooth$fUB  <-  predict(loess(paful ~ age, cld2.female, span=0.5), newdata=smooth$Age)

ggplot(data=smooth,aes(y=Est,x=Age)) +  geom_line(colour="mediumturquoise") +
  geom_ribbon(aes(ymin=LB,ymax=UB),alpha=0.2, fill="mediumturquoise")+
  geom_line(data=smooth,aes(y=fEst,x=Age, colour="lightcoral"), show.legend = F) +
  geom_ribbon(aes(ymin=fLB,ymax=fUB),alpha=0.2, fill="lightcoral")+
scale_x_continuous(name="\nAge (years)")+
  scale_y_continuous(name="Population attributable fraction, %\n", limits=c(0,0.8),breaks=seq(0,0.8,0.2),labels=c("0","0.2","0.4","0.6","0.8"))+
  theme(legend.position="none")+
  theme_light()


## Burden by gender ####
cld2[,c(5:9)] <- lapply(cld2[,c(5:9)], as.numeric)
cld2[,c("daly.ce","daly.ul","daly.ll")] <- lapply(cld2[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})
cld2$nm_se <- (as.numeric(cld2$daly.ul)-as.numeric(cld2$daly.ll))/3.92
cld2$nm_mean <- as.numeric(cld2$daly.ce)
cld2$alpha <- as.numeric(cld2$'Viremic Population')
cld2$beta <- as.numeric(cld2$'Total Population')-as.numeric(cld2$'Viremic Population')

for(i in 1:nrow(cld2)){
  x <- rbeta(num_iter, shape1 = cld2$alpha[i], shape2 = cld2$beta[i])
  y <- rnorm(num_iter,cld2$nm_mean[i],cld2$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  cld2$burdence[i] <- quantile(burden,0.5,na.rm = T)
  cld2$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  cld2$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  cld2$pafce[i] <- quantile(p,0.5,na.rm = T)
  cld2$pafll[i] <- quantile(p,0.025,na.rm = T)
  cld2$paful[i] <- quantile(p,0.975,na.rm = T)
}


#burden for both
cld2.bd.both <- subset(cld2, cld2$gender=="both")
cld2.bd.both$country <- as.character(cld2.bd.both$country)
cld2.bd.both[,14:19] <- lapply(cld2.bd.both[,14:19], as.numeric)

cld2.bd.both <- split(cld2.bd.both,list(cld2.bd.both$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,11),colSums(z[14:19])))}
cld2.bd.both <- lapply(cld2.bd.both, totalburden.total)
cld2.bd.both <- rbind.fill((cld2.bd.both))
cld2.bd.both <- cld2.bd.both[is.na(cld2.bd.both$gender),]
cld2.bd.both[,14:16] <- lapply(cld2.bd.both[,14:16], function(x){round((as.numeric(x)/1000), digits=2)})
cld2.bd.both$gender <- "both"
cld2.bd.both$age <- seq(from=20, to=80, by=5)

#burden for males
cld2.bd.male <- subset(cld2, cld2$gender=="male")
cld2.bd.male$country <- as.character(cld2.bd.male$country)
cld2.bd.male[,14:19] <- lapply(cld2.bd.male[,14:19], as.numeric)

cld2.bd.male <- split(cld2.bd.male,list(cld2.bd.male$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,11),colSums(z[14:19])))}
cld2.bd.male <- lapply(cld2.bd.male, totalburden.total)
cld2.bd.male <- rbind.fill((cld2.bd.male))
cld2.bd.male <- cld2.bd.male[is.na(cld2.bd.male$gender),]
cld2.bd.male[,14:16] <- lapply(cld2.bd.male[,14:16], function(x){round((as.numeric(x)/1000), digits=2)})
cld2.bd.male$gender <- "male"
cld2.bd.male$age <- seq(from=20, to=80, by=5)

#burden for females
cld2.bd.female <- subset(cld2, cld2$gender=="female")
cld2.bd.female$country <- as.character(cld2.bd.female$country)
cld2.bd.female[,14:19] <- lapply(cld2.bd.female[,14:19], as.numeric)

cld2.bd.female <- split(cld2.bd.female,list(cld2.bd.female$age.cat))
totalburden.total <- function(z){rbind(z, c("total", names(table(z$age.cat)),rep(NA,11),colSums(z[14:19])))}
cld2.bd.female <- lapply(cld2.bd.female, totalburden.total)
cld2.bd.female <- rbind.fill((cld2.bd.female))
cld2.bd.female <- cld2.bd.female[is.na(cld2.bd.female$gender),]
cld2.bd.female[,14:16] <- lapply(cld2.bd.female[,14:16], function(x){round((as.numeric(x)/1000), digits=2)})
cld2.bd.female$gender <- "female"
cld2.bd.female$age <- seq(from=20, to=80, by=5)

cld2.bd.gender <- rbind(cld2.bd.female,cld2.bd.male)

p7 <- ggplot(cld2.bd.gender, aes(x=as.numeric(as.character(age)), y=burdence,order=gender,fill=gender))
p7 <- p7 +
  stat_smooth(geom = 'area', method = 'loess', span = 0.4,
    alpha = 0.3)+
  scale_x_continuous(name="\nAge (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,150),breaks=seq(0,150,25),labels=c("0","25", "50","75","100","125","150"))+
  theme_light()+
  theme(legend.position = c(0.2, 0.83),legend.text = element_text(size=10),legend.title = element_blank())


p7


dat <-  cld2.bd.both
p8 <- ggplot(data <- dat) 
p8 <- p8+
  stat_smooth(aes(y=dat$burdence, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "blue", size=0.5)+
  stat_smooth(aes(y=dat$burdenll, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "red",linetype="dashed", size=0.5)+
  stat_smooth(aes(y=dat$burdenul, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "red",linetype="dashed", size=0.5)+
  scale_x_continuous(name="\nAge (years)", breaks=seq(20,80,10))+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,300),breaks=seq(0,300,50),labels=c("0","50","100","150","200", "250", "300"))+
  theme(legend.position="none")+
  theme_light()

p8

## Forest plots by age strata

cld2.bd.both$burden <-paste0(cld2.bd.both$burdence," ","(",cld2.bd.both$burdenll,"-",cld2.bd.both$burdenul,")") 


#forest plot
cld2.bd.both$age <- as.factor(cld2.bd.both$age)

p9 <- ggplot(data=cld2.bd.both, aes(x=age, y=burdence))+
  geom_errorbar(aes(ymin=burdenll, ymax = burdenul, width=0.3, colour="red"))+
  geom_point(colour="red",size=2)
p9 <- p9+scale_x_discrete(labels=c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"))+
  scale_y_continuous(breaks=seq(0,300,50))


p9 <- p9+theme_classic()+theme(axis.title.x = element_blank(), legend.position="none") + ylab("DALY, thousands (95% CI)\n")+xlab("\nAge Category (years)")
p9 <- p9+theme(axis.text.x = element_text(size=12,hjust=.5,vjust=.8,face="plain", angle=45),axis.text.y = element_text(size=12,hjust=.5,vjust=.8,face="plain"))
p9 <- p9+theme(axis.title.x = element_text(size=14,face="bold"),axis.title.y = element_text(size=14,face="bold"))
p9

## plot by income ####
wb <- read.csv("data_raw/WorldBankIncome.csv")
cld2.income <- merge(x = cld2, y = wb, by.x ="country", by.y = "Economy",all.x=T, all.y=F)
cld2.income <- cld2.income[,c(1:9,14:16,25)]
cld2.income$Income.group <- revalue(cld2.income$Income.group, c("High income"= "H", "Low income"="LM", "Lower middle income"="LM", "Upper middle income"="LM"))
cld2.income$Income.group <- as.character(cld2.income$Income.group)
cld2.income <- subset(cld2.income, cld2.income$gender=="both")
cld2.income$country <- as.character(cld2.income$country)

cld2.income <- split(cld2.income,list(cld2.income$Income.group,cld2.income$age.cat))
cld2.income <- lapply(cld2.income,function(x){rbind(x, c("total",names(table(x$age.cat)),rep(NA,5),colSums(x[8:12]),names(table(x$Income.group))))})
cld2.income <- rbind.fill((cld2.income))
cld2.income <- cld2.income[cld2.income$country=="total",]
cld2.income$burdence <- as.numeric(cld2.income$burdence)/1000
cld2.income$burdenll <- as.numeric(cld2.income$burdenll)/1000
cld2.income$burdenul <- as.numeric(cld2.income$burdenul)/1000

dat <-subset(cld2.income, cld2.income$Income.group=="H")
dat$age <- seq(from=20, to=80, by=5)
dat2 <- subset(cld2.income, cld2.income$Income.group=="LM")
dat2$age <- seq(from=20, to=80, by=5)

p10 <- ggplot(data <- dat) 
p10 <- p10+
  stat_smooth(aes(y=dat$burdence, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "mediumturquoise", size=0.5)+
  stat_smooth(aes(y=dat$burdenll, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "mediumturquoise",linetype="dashed", size=0.5)+
  stat_smooth(aes(y=dat$burdenul, x=dat$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "mediumturquoise",linetype="dashed", size=0.5)+
  stat_smooth(aes(y=dat2$burdence, x=dat2$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "lightcoral", size=0.5)+
  stat_smooth(aes(y=dat2$burdenll, x=dat2$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "lightcoral",linetype="dashed", size=0.5)+
  stat_smooth(aes(y=dat2$burdenul, x=dat2$age), method=lm, formula=y~poly(x,7), se=FALSE, colour = "lightcoral",linetype="dashed", size=0.5)+
  scale_x_continuous(name="\nAge (years)", breaks=seq(20,80,10))+
  scale_y_continuous(name="DALYs in thousands\n")+
  theme_light()


p10

x <- filter(cld2.income, cld2.income$age.cat=="80+")
x$Income.group <- as.factor(x$Income.group)
nrow(filter(x, x$Income.group=="H"))

## eTable Burden of HCV by country ####
tab4$'Total.Population' <- as.numeric(tab4$'Total.Population')
tab4$'Viremic.Population' <- as.numeric(tab4$'Viremic.Population')
tab4$country <- as.character(tab4$country)
tab4 <- tab4[!tab4$gender=="both",]
tab4 <- split(tab4,list(tab4$country))
totalburden <- function(z){rbind(z, c(names(table(z$country)),names(table(z$gbd.region)),rep(NA,2),colSums(z[5:9])))}
tab4 <- lapply(tab4, totalburden)
tab4 <- rbind.fill((tab4))
tab4 <- tab4[is.na(tab4$gender),]

tab4[,5:9] <- lapply(tab4[,5:9], function(x){as.numeric(x)})
tab4$daly.ce <- (tab4$daly.ce/tab4$Total.Population)*100000
tab4$daly.ll <- (tab4$daly.ll/tab4$Total.Population)*100000
tab4$daly.ul <- (tab4$daly.ul/tab4$Total.Population)*100000

tab4[,c("daly.ce", "daly.ll", "daly.ul")] <- lapply(tab4[,c("daly.ce", "daly.ll", "daly.ul")], function(x){round(as.numeric(x),digits=2)})
tab4$CVD <-paste0(tab4$daly.ce," ","(",tab4$daly.ll,"-",tab4$daly.ul,")") 
tab4$HCV <- countryleveldata.table$burden.C
tab4$prevalence <- round(((tab4$Viremic.Population/tab4$Total.Population)*100), digits=2)
tab4 <- tab4[,c(1,9,8,14,12,13)]

write.csv(tab4, "data_derived/eTable.burdenbycountry.csv")


## eTable Burden of HCV by age group ####
tab3 <- cld2.both[,c(2,8,9,17:19)]
tab3[,c("pafce", "pafll", "paful")] <- lapply(tab3[,c("pafce", "pafll", "paful")], function(x){round(as.numeric(x),digits=2)})
tab3$paf <-paste0(tab3$pafce," ","(",tab3$pafll,"-",tab3$paful,")") 
tab3$burden <- cld2.bd.both[,"burden"]
tab3$prevalence <- round((tab3$`Viremic Population`/tab3$`Total Population`)*100, digits=2)
tab3 <- tab3[,c(1,3,2,9,7:8)]
write.csv(tab3, "data_derived/eTable.burdenbyage.csv")


## eTable Burden of HCV by income group ####
dat[,c("burdence", "burdenll", "burdenul")] <- lapply(dat[,c("burdence", "burdenll", "burdenul")], function(x){round(as.numeric(x),digits=2)})
dat$burden <-paste0(dat$burdence," ","(",dat$burdenll,"-",dat$burdenul,")") 
dat2[,c("burdence", "burdenll", "burdenul")] <- lapply(dat2[,c("burdence", "burdenll", "burdenul")], function(x){round(as.numeric(x),digits=2)})
dat2$burden <-paste0(dat2$burdence," ","(",dat2$burdenll,"-",dat2$burdenul,")") 

tab5 <- dat[,c(2,15)]
names(tab5)[2] <- "High-income"
tab5$LMIC <- dat2[,"burden"]
write.csv(tab5, "data_derived/eTable.burdenbyincome.csv")

##Regional total burden
regional_total <- countryleveldata.table %>% 
  select(age.cat, burdence, burdenll, burdenul) %>% 
  group_by(age.cat) %>% 
  summarise_all(sum)

regional_total$burden <-paste0(regional_total$burdence," ","(",regional_total$burdenll,"-",regional_total$burdenul,")") 

sum(regional_total$burdence)

write.csv(regional_total, "data_derived/eTable.burdenbyregion.csv")


## Ukraine calculation
dev.off()

hist((rr.sim), breaks=100, border="blue", main= "Histogram for Risk Ratio", xlab="Risk ratio")

x <- rbeta(10000, 83800, 418200)
x <- x*100
hist(x, breaks=100, border="blue", main= "Histogram for Prevalence", xlab="Prevalence, %")

y <- rnorm(10000,13.605712,0.03534759)
y <- exp(y)
y <- y/1000
hist(y, breaks=100, border="blue", main= "Histogram for all CVD DALYs", xlab="DALYs, thousands")

p <- ((x/100)*(rr.sim-1))/(1+(x/100)*(rr.sim-1))
p <- p*100
hist(p, breaks=100, border="blue", main= "Histogram for PAF", xlab="PAF, %")

burden <- p*y
hist(burden, breaks=100, border="blue", main= "Histogram for CVD DALYs attributable to HCV", xlab="DALYs, thousands")


ukraine[,c("burdence", "burdenll", "burdenul")] <- lapply(ukraine[,c("burdence", "burdenll", "burdenul")], function(x){round(as.numeric(x),digits=2)})
ukraine$burden <-paste0(ukraine$burdence," ","(",ukraine$burdenll,"-",ukraine$burdenul,")") 


ukraine[,c("pafce", "pafll", "paful")] <- lapply(ukraine[,c("pafce", "pafll", "paful")], function(x){round(as.numeric(x)*100,digits=2)})
ukraine$paf <-paste0(ukraine$pafce," ","(",ukraine$pafll,"-",ukraine$paful,")") 


ukraine[,c("daly.ce","daly.ul","daly.ll")] <- lapply(ukraine[,c("daly.ce","daly.ul","daly.ll")],function(x){exp(as.numeric(x))})
ukraine[,c("daly.ce", "daly.ll", "daly.ul")] <- lapply(ukraine[,c("daly.ce", "daly.ll", "daly.ul")], function(x){round(as.numeric(x)/1000,digits=2)})
ukraine$daly <-paste0(ukraine$daly.ce," ","(",ukraine$daly.ll,"-",ukraine$daly.ul,")") 

ukraine$Viremic.Population <- as.numeric(ukraine$Viremic.Population)
ukraine$Total.Population <- as.numeric(ukraine$Total.Population)

ukraine$prevalence <- round((ukraine$Viremic.Population/ukraine$Total.Population)*100, digits=2)

ukraine <- ukraine[,c(2,3,9,8,23,21,22,20)]
write.csv(ukraine, "data_derived/ukraine.csv")



(quantile(burden, probs = c(0.025,0.5, 0.975)))
