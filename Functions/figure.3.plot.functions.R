# Plotting functions for figure 3-demographic changes duing transitions

create.turnover.data <- function(mod1, mod2){
  
  mod1.results<-extract.coefs.glmm(mod1)
  mod1.CI <- (confint(mod1))
  mod1.results <- cbind(mod1.results, mod1.CI)
  
  
  mod2.results<-extract.coefs.glmm(mod2)
  mod2.CI <- (confint(mod2))
  mod2.results <- cbind(mod2.results, mod2.CI)
  
  df <- data.frame('transition' = c('Non-Novel:Non-Novel','Non-Novel:Novel',
                                    'Novel:Non-Novel', "Novel:Novel") ,
                   'y'=plogis(c(mod1.results$Estimate[1],
                                mod1.results$Estimate[2]+mod1.results$Estimate[1],
                                mod1.results$Estimate[3]+mod1.results$Estimate[1],
                                mod1.results$Estimate[4]+mod1.results$Estimate[1])),
                   'x'=plogis(c(mod2.results$Estimate[1],
                                mod2.results$Estimate[2]+mod2.results$Estimate[1],
                                mod2.results$Estimate[3]+mod2.results$Estimate[1],
                                mod2.results$Estimate[4]+mod2.results$Estimate[1])),
                   'lower_y'=plogis(c(mod1.results$`2.5 %`[1],
                                      mod1.results$`2.5 %`[2]+mod1.results$`2.5 %`[1],
                                      mod1.results$`2.5 %`[3]+mod1.results$`2.5 %`[1],
                                      mod1.results$`2.5 %`[4]+mod1.results$`2.5 %`[1])),
                   'upper_y'=plogis(c(mod1.results$`97.5 %`[1],
                                      mod1.results$`97.5 %`[2]+mod1.results$`97.5 %`[1],
                                      mod1.results$`97.5 %`[3]+mod1.results$`97.5 %`[1],
                                      mod1.results$`97.5 %`[4]+mod1.results$`97.5 %`[1])),
                   'lower_x'=plogis(c(mod2.results$`2.5 %`[1],
                                      mod2.results$`2.5 %`[2]+mod2.results$`2.5 %`[1],
                                      mod2.results$`2.5 %`[3]+mod2.results$`2.5 %`[1],
                                      mod2.results$`2.5 %`[4]+mod2.results$`2.5 %`[1])),
                   'upper_x'=plogis(c(mod2.results$`97.5 %`[1],
                                      mod2.results$`97.5 %`[2]+mod2.results$`97.5 %`[1],
                                      mod2.results$`97.5 %`[3]+mod2.results$`97.5 %`[1],
                                      mod2.results$`97.5 %`[4]+mod2.results$`97.5 %`[1])),
                   'col.bef'= c('grey', 'grey', 'orange', 'orange'),
                   'col.aft' = c("grey", "orange", "grey","orange"))
  return(df)
}

plot.circles <- function(turnover.data){
  framemat <- rbind(c(0.125,.55,0.7,0.9))
  
  for(n in 1:4){
    shape = "circle"
    
    arrow.shape(x=turnover.data$x[n],
                y=turnover.data$y[n],
                r=0.175, screen=framemat[1,], rads=c(1.5, 0.5),
                col=turnover.data$col.aft[n],
                shape="circle",
                border=NA, plot=TRUE)
    arrow.shape(x=turnover.data$x[n],
                y=turnover.data$y[n],
                r=0.175, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape=shape,
                col=turnover.data$col.bef[n],
                plot=TRUE, border="black", add.arrow=TRUE)
    arrow.shape(x=turnover.data$x[n],
                y=turnover.data$y[n],
                r=0.175, screen=framemat[1,], rads=c(0, 2),
                shape=shape,
                plot=TRUE, border="black", add.arrow=FALSE)
  }
}

figure.3.turnover <- function(emig.mod, ext.mod, immig.mod, orig.mod){
  
  df.ext<-create.turnover.data(emig.mod, ext.mod)
  df.orig<-create.turnover.data(immig.mod, orig.mod)
  par(mfrow=(c(2,1)))
  par(mar=c(5,4,2,4))
  plot(df.ext$y~df.ext$x, xlab='Probability of Extinction', 
       ylab='Probability of Emigration', ylim=c(0,.4), xlim=c(0,.2), type='n')
  # CI's
  for(i in 1:4){
    #vrtical
    arrows(x0=df.ext$x[i],y0=df.ext$y[i],x1=df.ext$x[i],y1=df.ext$upper_y[i], lwd=1, angle=180, length=0.05 )
    arrows(x0=df.ext$x[i],y0=df.ext$y[i],x1=df.ext$x[i],y1=df.ext$lower_y[i], lwd=1, angle=180, length=0.05 )
    #horizontal
    arrows(x0=df.ext$x[i],y0=df.ext$y[i],x1=df.ext$upper_x[i],y1=df.ext$y[i], lwd=1, angle=180, length=0.05 )
    arrows(x0=df.ext$x[i],y0=df.ext$y[i],x1=df.ext$lower_x[i],y1=df.ext$y[i], lwd=1, angle=180, length=0.05 )
    
  }
  plot.circles(df.ext)
  abline(0,1)
  
  par(mar=c(5,4,2,4))
  #plot(2)
  plot(df.orig$y~df.orig$x, xlab='Probability of Origination', 
       ylab='Probability of Immigration', ylim=c(0,.3), xlim=c(0,.2), type='n')
  # CI's
  for(i in 1:4){
    #  vrtical
    arrows(x0=df.orig$x[i],y0=df.orig$y[i],x1=df.orig$x[i],y1=df.orig$upper_y[i], lwd=1, angle=180, length=0.05 )
    arrows(x0=df.orig$x[i],y0=df.orig$y[i],x1=df.orig$x[i],y1=df.orig$lower_y[i], lwd=1, angle=180, length=0.05 )
    #horizontal
    arrows(x0=df.orig$x[i],y0=df.orig$y[i],x1=df.orig$upper_x[i],y1=df.orig$y[i], lwd=1, angle=180, length=0.05 )
    arrows(x0=df.orig$x[i],y0=df.orig$y[i],x1=df.orig$lower_x[i],y1=df.orig$y[i], lwd=1, angle=180, length=0.05 )
    
  }
  plot.circles(df.orig)
  
  abline(0,1)
  
  
}


