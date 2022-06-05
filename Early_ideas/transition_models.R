#### Estimating transition probabilities ####

estimate.observed.expected <- function(prob.model.list,
                                       novel.list,
                                       dist.draws = 1e6){
  
  # extract observed classification probabilities from models
  print("Preparing observed transition data...")
  
  # condense community data into single data-frame.
  comm.dat <- do.call('rbind', lapply(1:4, function(n){
    x <- novel.list[[n]][]
    temp <- do.call("rbind", x)
    temp$BioRealm <- c("PAL", "NEA", "AFRO", "NEO")[n]
    return(temp)
  }))
  
  # WHAT WOULD BE REALLY INTERESTING IS TURNING THIS INTO BASIN INSTEAD OF BIOREALM!
  # set site to be a taxa-specific factor
  comm.dat$site <- as.factor(paste0(comm.dat$BioRealm,":", comm.dat$site))
  
  # observed transition data-frame
  obs.df <- as.data.frame(table(comm.dat$cat.bef,
                                comm.dat$cat,
                                comm.dat$site))
  colnames(obs.df) <- c("cat.bef", "cat.aft", "site", "obs")
  obs.df$taxa <- substr(obs.df$site, 1, regexpr(":", obs.df$site)-1)
  obs.df$taxa.site <- paste0(obs.df$taxa, ":", obs.df$site)
  
  # get all transition combinations
  trans.list <- expand.grid(levels(obs.df$cat.bef), 
                            levels(obs.df$cat.aft),
                            stringsAsFactors = FALSE)
  
  # model observed probabilities of transition occurring
  print("Running observed transition models...")
  trans.preds <- lapply(1:16, function(n){
    
    target.trans <- trans.list[n,]
    
    # subset just frequency of target transition
    trans.df <- do.call("rbind", lapply(split(obs.df, f=obs.df$site), 
                                        function(x){
                                          
                                          x$cat.bef <- as.character(x$cat.bef)
                                          x$cat.aft <- as.character(x$cat.aft)
                                          
                                          data.frame(site=x$site[1],
                                                     taxa=x$taxa[1],
                                                     success = x$obs[x$cat.bef == target.trans[1,1] &
                                                                       x$cat.aft == target.trans[1,2]],
                                                     failure = sum(x$obs[x$cat.bef != target.trans[1,1] |
                                                                           x$cat.aft != target.trans[1,2]]))
                                          
                                        }))
    
    trans.m <- glmer(cbind(success, failure) ~ 1 + (1|taxa),
                     data = trans.df, family = binomial,
                     glmerControl(optimizer ="bobyqa"))
    
    
    return(list(transition = target.trans,
                model = trans.m,
                fixed = summary(trans.m)$coefficients))
    
  })
  
  # calculate expected probabilities from shift probability models
  print("Calculating expected transition probabilities...")
  base.probs <- do.call("rbind", lapply(prob.model.list$random.prob.models,
                                        function(x){x$pred.df}))
  
  rownames(base.probs) = 1:nrow(base.probs)
  base.probs <- as.data.frame(base.probs)
  base.probs$cat = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")
  base.probs <- base.probs[-(1:2),]
  
  print("Drawing from expected and observed distributions...")
  obs.exp.probs.df <- do.call("rbind", lapply(1:16, function(n){
    
    target.trans <- cbind(trans.preds[[n]]$transition,
                          trans.preds[[n]]$fixed)[,1:4]
    colnames(target.trans)[1:2] = c("cat.bef", "cat.aft")
    
    trans1 <- base.probs[base.probs$cat == target.trans[1,1],]
    trans2 <- base.probs[base.probs$cat == target.trans[1,2],]
    
    final.df <- target.trans
    colnames(final.df)[3:4] = paste0("obs.", colnames(final.df)[3:4])
    
    # draw randomly from probability distribution on logit-scale,
    # then calculate combined probabilities, return distribution
    # of expected probabilities
    exp.probs.bef <- rnorm(1e6, trans1$Estimate, trans1$`Std. Error`)
    exp.probs.aft <- rnorm(1e6, trans2$Estimate, trans2$`Std. Error`)
    
    exp.probs.bef.raw = plogis(exp.probs.bef)
    exp.probs.aft.raw = plogis(exp.probs.aft)
    
    exp.probs.comb = exp.probs.bef.raw * exp.probs.aft.raw
    exp.probs.df = data.frame(exp.mean = mean(exp.probs.comb),
                              exp.upper = quantile(exp.probs.comb, 0.975),
                              exp.lower = quantile(exp.probs.comb, 0.025))
    
    obs.probs = rnorm(1e6, final.df$obs.Estimate, final.df$`obs.Std. Error`)
    
    ratio.probs =  plogis(obs.probs) / exp.probs.comb
    
    ratio.probs.df = data.frame(ratio.mean = mean(ratio.probs),
                                ratio.upper = quantile(ratio.probs, 0.975),
                                ratio.lower = quantile(ratio.probs, 0.025))
    
    return(cbind(final.df, exp.probs.df, ratio.probs.df))
    
  }))
  
  obs.exp.probs.df$non.zero <- obs.exp.probs.df$ratio.lower > 1 & obs.exp.probs.df$ratio.upper > 1 |
    obs.exp.probs.df$ratio.lower < 1 & obs.exp.probs.df$ratio.upper < 1
  
  obs.exp.probs.df <- obs.exp.probs.df[order(obs.exp.probs.df$non.zero,
                                             obs.exp.probs.df$ratio.mean),]
  
  return(obs.exp.probs.df)
  
}

figure2.plot <- function(trans.df, plot.name, ylims){
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  light.cols <- rbind(c(205,205,205),
                      c(207,228,237),
                      c(255,179,179),
                      c(255,228,179)) / 255
  light.cols <- rgb(light.cols[,1], light.cols[,2], light.cols[,3])
  
  library(shape)
  library(plotrix)
  
  #pdf(date.wrap(paste0("./plots/transition null model (", 
                       #plot.name, 
                       #")"), ".pdf"), 
      #height=2, width=2.25, useDingbats = FALSE, colormodel = "cmyk")
  
  framemat<-rbind(c(0.19,0.99,0.175,0.99),
                  c(0.8,1,0.25,0.75))
  
  split.screen(framemat)
  
  par(mar=c(0,0,0,0), ps=6, tcl=-0.2, las=1, mgp=c(3,0,0))
  plot(x=NULL, y=NULL, xlim=log(c(2e-4,1)), ylim=ylims,
       axes=FALSE, xlab="", ylab="", yaxs="i")
  
  abline(h=0, lwd=0.75, lty="31")
  
  custom.circle(x=log(1), y=log(1), r=0.5, col=rgb(0.5,0.5,0.5,0.35),
                screen=framemat[1,], border=NA)
  
  draw.ellipse(x=log(0.05), y=log(0.775), a=0.85, b=0.35,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle= 0)
  
  draw.ellipse(x=log(0.0002), y=log(2.8), a=2, b=1,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle = 0)
  
  mtext(side=1, line=0.4, cex=1,
        text="Expected probability of transition")
  mtext(side=1, line=0.8, cex=0.8,
        text="(calculated from occurrence probabilities)")
  
  par(lheight=0.85)
  mtext(side=2, line=1.1, las=0,
        text="Ratio of observed to expected probability")
  
  text(x=log(1), y=log(1.05), pos=3, adj=0.5, offset=0.8,
       labels="Background\ntransitions", col="grey60", font=1)
  
  text(x=log(0.05), y=log(0.75),  adj=0,
       labels="Transitions\nto and\nfrom\nnovelty", col="grey60", font=1)
  
  text(x=log(0.0002), y=log(2.8),  adj=0,
       labels="Transitions\nbetween\nnovelty\ncategories", col="grey60", font=1)
  
  axis(side=1,
       at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
       labels=c(0.0001, 0.001,0.01,0.1,1,10,100),
       mgp=c(3,-0.15,0))
  
  axis(side=1, at=log(c(seq(0.0001,0.001,0.0001),
                        seq(0.001,0.01,0.001),
                        seq(0.01,0.1,0.01),
                        seq(0.1,1,0.1),
                        seq(1,10,1),
                        seq(10,100,10))), tcl=-0.1, labels=NA)
  
  y.locs <- c(seq(1, 15, 0.5),
              1/seq(1, 15, 0.5))
  
  y.locs <- c(seq(1,20,1), seq(10,100,10))
  y.locs <- c(y.locs, 1/y.locs)
  
  axis(side=2, at=log(y.locs), tcl=-0.125, labels=NA)
  
  axis(side = 2, 
       at = log(c(0.5,1,2,3,5,10,20,30,50,100)),
       labels = c(0.5,1,2,3,5,10,20,30,50,100), 
       mgp=c(3,0.4,0))
  
  sapply(1:dim(trans.df)[1], function(x){
    
    if(trans.df$non.zero[x]){
      
      aft.col = c("grey35", cumul.col, "red", "orange")[as.factor(trans.df$cat.aft)[x]]
      bef.col = c("grey35", cumul.col, "red", "orange")[as.factor(trans.df$cat.bef)[x]]
      border.col = "black"
      
    } else {
      
      aft.col = light.cols[as.factor(trans.df$cat.aft)[x]]
      bef.col = light.cols[as.factor(trans.df$cat.bef)[x]]
      border.col = "grey50"
      
    }
    
    segments(y0=log(trans.df$ratio.mean[x]),
             y1=log(trans.df$ratio.mean[x]),
             x0=log(trans.df$exp.upper[x]),
             x1=log(trans.df$exp.lower[x]),
             col=ifelse(trans.df$non.zero[x],
                        "black", "grey50"), lwd=0.75)
    
    segments(y0=log(trans.df$ratio.lower[x]),
             y1=log(trans.df$ratio.upper[x]),
             x0=log(trans.df$exp.mean[x]),
             x1=log(trans.df$exp.mean[x]),
             col=ifelse(trans.df$non.zero[x],
                        "black", "grey50"), lwd=0.75)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,], rads=c(1.5, 0.5),
                col = aft.col,
                shape="circle",
                border=NA)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col = bef.col, 
                border=border.col, add.arrow=TRUE)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,], rads=c(0, 2),
                shape="circle", lwd=0.5,
                plot=TRUE, border=border.col, add.arrow=FALSE)
    
    
  })
  box()
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.375,0.675,len=4))
  
  par(lheight=0.85)
  text(x=0.1, y=1.3, labels="Preceding\nstate", adj=0, cex=0.45,
       lheight=0.5)
  
  text(x=0.85, y=1.135, labels="Succeeding\nstate", adj=1, cex=0.45,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.55, y=1.45,
              r=0.125, screen=framemat[2,], rads=c(1.5,0.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE)
  
  arrow.shape(x=0.45, y=1.45,
              r=0.125, screen=framemat[2,], rads=c(0.5,1.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(0.0725, 0.025, 0.875, 0.925),
           x1=c(0.025, 0.025, 0.925, 0.925),
           y0=c(1.29, 1.29, 1.125, 1.125),
           y1=c(1.29, 1.39, 1.125, 1.35), lwd=0.5)
  
  Arrows(x0=c(0.025, 0.925),
         x1=c(0.275, 0.725),
         y0=c(1.39, 1.35),
         y1=c(1.44, 1.44), lwd=0.5,
         arr.length = 0.05, arr.width = 0.05,
         arr.type = "triangle")
  # 
  rect.pos <- rev(seq(1.175,1.475, len=4))
  
  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[1] + 0.04, ybottom=rect.pos[1] - 0.04, col="grey35")
  
  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[2] + 0.04, ybottom=rect.pos[2] - 0.04, col="red")
  
  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[3] + 0.04, ybottom=rect.pos[3] - 0.04, col=cumul.col)
  
  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[4] + 0.04, ybottom=rect.pos[4] - 0.04, col="orange")
  
  par(lheight=0.7)
  text(x=-2.05, y=rect.pos - 0.01, pos=4, offset=0.75, cex=0.45,
       labels=c("Background", "Instantaneous novelty",
                "Cumulative novelty", "Novel community"))
  par(lheight=1)
  
  par(xpd=FALSE)
  close.screen(2)
  
  
}

custom.circle <- function(x,y,r,rads=c(0,2),nsteps=100,
                          screen=NA, ...){  
  
  plot.window<-dev.size(units="cm")
  
  if(!is.na(screen)[1]){
    plot.window <- dev.size(units="cm") * c(screen[2]-screen[1],
                                            screen[4]-screen[3])
    
  }
  
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  x<-(x-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  y<-(y-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  
  if(rads[2] < rads[1]){
    
    rs <- c(seq(rads[1]*pi, 2*pi, len=nsteps*((2-rads[1]) / (2-rads[1]+rads[2]))),
            seq(0*pi, rads[2]*pi, len=nsteps*(rads[2] / (2-rads[1]+rads[2]))))
    
    width <- 2-max(rads)+min(rads)
  } else {
    rs <- seq(rads[1]*pi,rads[2]*pi,len=nsteps)  
    width <- rads[2]-rads[1]
  }
  
  xc <- x+r*cos(rs) 
  yc <- y+r*sin(rs) 
  
  if(width < 1){
    xc <- c(xc, x)
    yc <- c(yc, y)
  }
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...) 
  
  return(cbind(xc,yc))
}

arrow.shape <- function(x,y,r,nsteps=100, shape="circle",
                        screen=NA, rads=c(0,2), add.arrow=FALSE, ...){
  
  # set aspet ratio to ensure symmetrical shapes
  plot.window<-dev.size(units="cm")
  
  if(!is.na(screen)[1]){
    plot.window <- dev.size(units="cm") * c(screen[2]-screen[1],
                                            screen[4]-screen[3])
    
  }
  
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  x<-(x-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  y<-(y-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  
  if(shape=="circle"){
    
    if(rads[2] < rads[1]){
      
      rs <- c(seq(rads[1]*pi, 2*pi, len=nsteps*((2-rads[1]) / (2-rads[1]+rads[2]))),
              seq(0*pi, rads[2]*pi, len=nsteps*(rads[2] / (2-rads[1]+rads[2]))))
      
      width <- 2-max(rads)+min(rads)
    } else {
      rs <- seq(rads[1]*pi,rads[2]*pi,len=nsteps)  
      width <- rads[2]-rads[1]
    }
    
    xc <- x+r*cos(rs) 
    yc <- y+r*sin(rs) 
    
    if(width < 1){
      xc <- c(xc, x)
      yc <- c(yc, y)
    }
    
  }
  
  if(shape=="triangle"){
    
    xc <- c(x, x-r, x, x+r)
    yc <- c(y+r, y-r, y-r, y-r)
    
  }
  
  if(shape=="square"){
    
    xc <- c(x,   x-r, x-r, x,   x+r, x+r)
    yc <- c(y+r, y+r, y-r, y-r, y-r, y+r)
    
    
  }
  
  if(shape=="diamond"){
    
    xc <- c(x, x-1.2*r, x, x+1.2*r)
    yc <- c(y+1.2*r, y, y-1.2*r, y)
    
  }
  
  if(add.arrow){
    
    if(shape=="circle"){
      shape.top <- rbind(cbind(xc,yc)[yc == max(yc),],
                         cbind(xc,yc)[yc == min(yc),])
      
    } else {
      shape.top <- rbind(cbind(xc,yc)[yc == max(yc) & xc==x,],
                         cbind(xc,yc)[yc == min(yc) & xc==x,])
    }
    
    arrow.dist <- c(max(xc)-min(xc), max(yc)-min(yc))
    if(shape !="circle"){
      arrow.dist <- arrow.dist*c(0.45,1)
    }
    
    arrow.points <- rbind(shape.top[1,] + c(0, - 0.4*arrow.dist[2]),
                          shape.top[1,] + c(0.2*arrow.dist[1],  - 0.4*arrow.dist[2]),
                          shape.top[1,] + c(0.2*arrow.dist[1],  - 0.3*arrow.dist[2]),
                          shape.top[1,] + c(0.6*arrow.dist[1],  - 0.5*arrow.dist[2]),
                          shape.top[2,] + c(0.2*arrow.dist[1],  + 0.3*arrow.dist[2]),
                          shape.top[2,] + c(0.2*arrow.dist[1],  + 0.4*arrow.dist[2]),
                          shape.top[2,] + c(0, + 0.4*arrow.dist[2]))
    
    if(shape=="triangle"){
      arrow.points[,2] <- arrow.points[,2] - 0.2*arrow.dist[2]
    }
    
    xarrow <- c(xc[which(xc == shape.top[1,1] & yc==shape.top[1,2])],
                arrow.points[,1],
                xc[which(xc == shape.top[2,1] & yc==shape.top[2,2]):
                     which(xc == shape.top[1,1] & yc==shape.top[1,2])])
    
    yarrow <- c(yc[which(xc == shape.top[1,1] & yc==shape.top[1,2])],
                arrow.points[,2],
                yc[which(xc == shape.top[2,1] & yc==shape.top[2,2]):
                     which(xc == shape.top[1,1] & yc==shape.top[1,2])])
    
    xc <- xarrow
    yc <- yarrow
    
  }  
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...)
  # return(cbind(xc,yc))
}

figure2.plot.fixed <- function(trans.df, plot.name, xlims, ylims){
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  light.cols <- rbind(c(205,205,205),
                      c(207,228,237),
                      c(255,179,179),
                      c(255,228,179)) / 255
  light.cols <- rgb(light.cols[,1], light.cols[,2], light.cols[,3])
  
  library(shape)
  library(plotrix)
  
  #pdf(date.wrap(paste0("./plots/transition null model (", 
              #         plot.name, 
               #        ")"),
               # ".pdf"), 
    #  height=5.5, width=7.5, useDingbats = FALSE)
  
  framemat<-rbind(c(0.10,0.47,0.525,0.95),
                  c(0.47,0.84,0.525,0.95),
                  c(0.1,0.47,0.125,0.525),
                  c(0.47,0.84,0.125,0.525),
                  c(0.85,1,0.35,0.75))
  
  # framemat<-rbind(c(0.10,0.8,0.135,0.95),
  #                 c(0.8,1,0.25,0.75))
  # 
  split.screen(framemat)
  
  sapply(1:4, function(n){
    
    screen(n)
    par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.2,0))  
    
    sub.trans <- trans.df[trans.df$taxa == levels(trans.df$taxa)[n],]
    
    plot(x=NULL, y=NULL, xlim = xlims, ylim = ylims,
         axes=FALSE, xlab="", ylab="", yaxs="i")
    
    abline(h=0, lwd=2,col="grey60", lty="31")
    
    # custom.circle(x=log(1), y=log(1), r=0.5, col=rgb(0.5,0.5,0.5,0.15),
    #               screen=framemat[n,], border=NA)
    # 
    # draw.ellipse(x=log(0.05), y=log(0.775), a=0.85, b=0.35,
    #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
    #              angle= 0)
    # 
    # draw.ellipse(x=log(0.0002), y=log(2.8), a=2, b=1,
    #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
    #              angle = 0)
    # 
    if(n == 3){
      mtext(side=1, line=2.25, at = par("usr")[2],
            text="Expected probability of transition\n(calculated from occurrence probabilities)")
    }
    
    if(n == 1){
      mtext(side=2, line=1.85, las=0, at = par("usr")[3],
            text="Observed probability / Expected probability")
    }
    
    # text(x=log(1), y=log(1.05), pos=3, adj=0.5, offset=0.8,
    #      labels="Background\nshifts", col="grey60", font=1)
    # 
    # text(x=log(0.05), y=log(0.75),  adj=0,
    #      labels="Shifts\nto and\nfrom\nnovelty", col="grey60", font=1)
    # 
    # text(x=log(0.0002), y=log(2.8),  adj=0,
    #      labels="Shifts\nbetween\nnovelty\ncategories", col="grey60", font=1)
    
    if(n %in% 3:4){
      axis(side=1, 
           at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
           labels=c(0.0001, 0.001,0.01,0.1,1,10,100))
    } else {
      axis(side=1, 
           at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
           labels=NA)
    }
    
    axis(side=1, at=log(c(seq(0.0001,0.001,0.0001),
                          seq(0.001,0.01,0.001),
                          seq(0.01,0.1,0.01),
                          seq(0.1,1,0.1),
                          seq(1,10,1),
                          seq(10,100,10))), tcl=-0.125, labels=NA)
    
    y.locs <- c(seq(1, 10, 1),
                1/seq(1, 10, 0.5))
    
    axis(side=2, at=log(y.locs), tcl=-0.125, labels=NA)
    
    if(n %in% c(1,3)){
      axis(side = 2, 
           at = log(c(0.5,1,2,3,5,10,20)),
           labels = c(0.5,1,2,3,5,10,20), 
           mgp=c(3,0.5,0))
    } else {
      axis(side = 2, 
           at = log(c(0.5,1,2,3,5,10,20)),
           labels = NA, 
           mgp=c(3,0.5,0))
    }
    
    sub.trans <- sub.trans[order(sub.trans$non.zero),]
    
    sapply(1:dim(sub.trans)[1], function(x){
      print(x)
      
      segments(y0=log(sub.trans$ratio.mean[x]),
               y1=log(sub.trans$ratio.mean[x]),
               x0=log(sub.trans$exp.upper[x]),
               x1=log(sub.trans$exp.lower[x]),
               col=ifelse(sub.trans$non.zero[x],
                          "black", "grey70"))
      
      segments(y0=log(sub.trans$ratio.upper[x]),
               y1=log(sub.trans$ratio.lower[x]),
               x0=log(sub.trans$exp.mean[x]),
               x1=log(sub.trans$exp.mean[x]),
               col=ifelse(sub.trans$non.zero[x],
                          "black", "grey70"))
      
      if(sub.trans$non.zero[x]){
        
        aft.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.aft)[x]]
        bef.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.bef)[x]]
        shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
        border.col = "black"
        
      } else {
        
        aft.col = light.cols[as.factor(sub.trans$cat.aft)[x]]
        bef.col = light.cols[as.factor(sub.trans$cat.bef)[x]]
        shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
        border.col = "grey70"
        
      }
      
      arrow.shape(x=log(sub.trans$exp.mean[x]),
                  y=log(sub.trans$ratio.mean[x]),
                  r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
                  col = aft.col,
                  shape=shape,
                  border=NA)
      
      arrow.shape(x=log(sub.trans$exp.mean[x]),
                  y=log(sub.trans$ratio.mean[x]),
                  r=0.2, screen=framemat[1,],
                  rads=c(0.5, 1.5), lwd=0.5,
                  shape=shape,
                  col = bef.col, 
                  border=border.col, add.arrow=TRUE)
      
      arrow.shape(x=log(sub.trans$exp.mean[x]),
                  y=log(sub.trans$ratio.mean[x]),
                  r=0.2, screen=framemat[1,], rads=c(0, 2),
                  shape=shape, border=border.col, add.arrow=FALSE)
      
      
    })
    
    text(x=relative.axis.point(0.02, "x"),
         y = relative.axis.point(0.935, "y"),
         labels = paste0("(",LETTERS[n],") ",
                         c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians")[n]),
         font=2, adj=0)
    box()
    close.screen(n)
  })
  
  screen(5)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.375,0.675,len=4))
  
  text(x=0.575, y=0.75, adj=0.5,
       labels=bquote(bold(underline("Taxa"))), font=2, cex=0.8)
  sapply(taxa.pos, function(y){
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(1.5,0.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(0.5, 1.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE, add.arrow=TRUE, lwd=0.5)
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(0,2),
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)
    
    
  })
  # text(x=0.4, pos=4, y=taxa.pos,
  #      labels=c("Nanno", "Foram", "Radio", "Diatom"), 
  #      cex=0.8, offset=0.75)
  
  text(x=0.13, pos=4, y=taxa.pos,
       labels=sort(c("Nannoplankton", "Foraminifera", "Radiolarians", "Diatoms")), 
       cex=0.8, offset=0.75)
  
  par(lheight=0.85)
  text(x=0.15, y=1.05, labels="Preceding\ncommunity", adj=0, cex=0.8,
       lheight=0.5)
  
  text(x=0.85, y=0.925, labels="Succeeding\ncommunity", adj=1, cex=0.8,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.55, y=1.2,
              r=0.25, screen=framemat[5,], rads=c(1.5,0.5),
              col="white", shape="circle",
              border="black", plot=TRUE)
  
  arrow.shape(x=0.45, y=1.2,
              r=0.25, screen=framemat[5,], rads=c(0.5,1.5),
              col="white", shape="circle",
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(0.125, 0.075, 0.875, 0.925),
           x1=c(0.075, 0.075, 0.925, 0.925),
           y0=c(1.05, 1.05, 0.925, 0.925),
           y1=c(1.05, 1.15, 0.925, 1.15))
  
  Arrows(x0=c(0.075, 0.925),
         x1=c(0.3, 0.7),
         y0=c(1.15, 1.15),
         y1=c(1.2, 1.2),
         arr.length = 0.1, arr.width = 0.1,
         arr.type = "triangle")
  
  rect.pos <- rev(seq(-0.275,0.05, len=4))
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[1] + 0.03, ybottom=rect.pos[1] - 0.03, col="grey35")
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[2] + 0.03, ybottom=rect.pos[2] - 0.03, col="red")
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[3] + 0.03, ybottom=rect.pos[3] - 0.03, col=cumul.col)
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[4] + 0.03, ybottom=rect.pos[4] - 0.03, col="orange")
  
  par(lheight=0.7)
  text(x=0.125, y=rect.pos - c(0, 0.015, 0.015, 0.015), pos=4, offset=0.75, cex=0.8,
       labels=c("Background", "Instantaneous\nnovelty", 
                "Cumululative\nnovelty", "Novel\ncommunity"))
  par(lheight=1)
  
  text(x=0.575, y=0.215, adj=0.5,
       labels=bquote(bold(underline("Community"))), font=2, cex=0.8)
  text(x=0.575, y=0.15, adj=0.5,
       labels=bquote(bold(underline("classification"))), font=2, cex=0.8)
  
  par(xpd=FALSE)
  close.screen(5)
  
  dev.off()
}
