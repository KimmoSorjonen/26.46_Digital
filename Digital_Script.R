

############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

## Loading package

library(MASS)
library(metafor)
library(lavaan)

#################
## Data
## Correlations reported by Che et al. (2026)

rm <- matrix(c(
  
  1.00, 0.67, 0.21, 0.15,
  0.67, 1.00, 0.20, 0.24,
  0.21, 0.20, 1.00, 0.65,
  0.15, 0.24, 0.65, 1.00), nrow=4)

n <-  9163 ## total sample size

df <- data.frame(mvrnorm(n=n, mu=rep(0,4), ## generating data frame with required size and corr.
            Sigma=rm, empirical=T))

names(df) <- c("DMU1","DMU2","CH1","CH2")

## The three models, with dy = Y2-Y1

m1 <- "dy ~ x1 + y1"
m2 <- "dy ~ x1 + y2"
m3 <- "dy ~ x1"

modlist <- list(m1,m2,m3) ## list with all three models

#################
## Figure

panlab <- c("A: DMU → CH", 
            "B: CH → DMU") ## labels on top

slab <- c("4.RMA","","3.no.adj", ## labels for y-axis
          "2.adj.y2","1.adj.y1")

cx <- 0.8 ## sizing factor
r.low <- -0.15 ## range, lower
r.upp <- 0.1 ## range, upper
f.upp <- r.upp+1.8*(r.upp-r.low) ## room for text

if(dev.cur()==2) dev.off() ## removing earlier plots
par(mar=c(1,1,0,0), oma=c(1.5,3,1.4,1.4), mfrow=c(1,2)) ## setting margins and layout

for(i in 1:2){ ## in both directions
  
  if(i==1) names(df) <- c("x1","x2","y1","y2") ## DMU as X and CH as Y
  if(i==2) names(df) <- c("y1","y2","x1","x2") ## CH as X and DMU as Y
  
  dy <- df$y2 - df$y1 ## difference score
  
  numeff <- 5 ## number of rows per panel (one is empty)
  
  plot(c(r.low,f.upp),c(0.5,5.5), type="n",xaxt="n",yaxt="n",xlab="",ylab="") ## empty plot
  axis(1,at=c(-0.1,0,0.1),labels=c(-0.1,0,0.1), cex.axis=cx) ## x-labels
  axis(2,at=1:5,labels=F, las=1, cex.axis=cx) ## y-tics
  if(i==1) axis(2,at=1:5,labels=slab, las=1, cex.axis=cx) ## y-labels
  
  mtext(panlab[i],3, line=0.2, cex=cx) ## label on top
  
  lines(c(0,0),c(-1,15),col="gray") ## vertical gray line at x=0
  
  ball <- vector() ## to be filled with effects below
  seall <- vector() ## to be filled with standard errors below
  
  for(j in 1:3){ ## for the three models
    
    fit <- lm(modlist[[j]], data=df) ## fitting the model
    b <- fit$coefficients[2] ## reg. coefficient
    low <- confint(fit)[2,1] ## CI, low
    upp <- confint(fit)[2,2] ## CI, upp
    se <- summary(fit)$coefficients[2,2] ## standard error
    tx <- paste(round(b,2)," [", round(low,2),"; ", 
                round(upp,2),"]", sep="") ## string with values
    
    ball <- c(ball,b) ## adding effect to object
    seall <- c(seall,se) ## adding se to object
    
    lines(c(-2,r.upp),c(numeff,numeff), col="gray") ## vertical gray line
    arrows(low,numeff,upp,numeff,angle=90,code=3,length=0.05,lwd=2) ## CI
    points(b,numeff,pch=21,col="black",bg="black") ## point for effect
    text(r.upp,numeff,tx,pos=4,cex=cx) ## adding string with values
    
    numeff <- numeff-1 ## next row, please
  }
  
  b.fish <- 0.5*log((1+ball)/(1-ball)) ## Fisher's transformation of effects
  se.fish <- 0.5*log((1+seall)/(1-seall)) ## Fisher's transformation of SE
  
  we <- 1/3 ## weighting the effects (a1, b1, and c1) equally
  
  ###### Following Bartos et al.
  
  random1 <- rma(yi=b.fish, vi=se.fish^2)
  res.ma <- rma(yi=b.fish, vi=(se.fish^2)/we, tau=random1$tau2)
  
  ######
  
  pred <- predict(res.ma, transf=transf.ztor) ## transforms back from Fisher's
  metb <- pred$pred; metlow <- pred$ci.lb; metupp <- pred$ci.ub ## estimate with lower and upper CI limit
  
  lines(c(-2,5),c(2,2),lty=2) ## dashed line
  lines(c(-2,r.upp),c(1,1),col="gray") ## gray vertical line
  
  polygon(c(metlow,metb,metupp,metb), c(1,1.5,1,0.5), ## diamond for RMA 
          col = "black", border = "black", lwd = 1)
  
  tx.ma <- paste(round(metb,2)," [", round(metlow,2),"; ", round(metupp,2),"]", sep="") ## string with values
  text(r.upp,1,tx.ma,pos=4,cex=cx) ## adding string of values
  
}


#############################
## The SPAM

spam <- "

## latent gx, gy

gx =~ b*x1 + b*x2
gy =~ c*y1 + c*y2

## second-order trait-variable

t =~ a*gx+a*gy

## state-correlations

x1 ~~ d*y1
x2 ~~ d*y2

## (error) variances

t ~~ 1*t
gx ~~ 1*gx
gy ~~ 1*gy
x1 ~~ x1
x2 ~~ x2
y1 ~~ y1
y2 ~~ y2

"

names(df) <- c("x1","x2","y1","y2") ## DMU as X and CH as Y

fit <- lavaan(spam, meanstructure=F, data=df) ## fitting the SPAM

summary(fit, fit.measures=T, ## let's have a look 
        standardized=T, ci=T, rsq=T)


########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################

