d <- readRDS("correlation-matrix.rds")

library(tidyverse)
library(manipulate)
#library(psychometric)
library(ggsidekick)
predictors <- select(d, ENSO:V5)
responses <- select(d, trend1:IchSW)

names(predictors)=c('ENSO','NPGO','NPI','PDO','Upwell','SST','POLL','PCOD','POP','ARROW')
names(responses)=c('Trend 1','Trend 2', 'Sp Rich','Shannon')

autocorr_pyper<-function(N,tsx,tsy) {
  lag.max<-min(length(tsx),length(tsy))/4
  rxx<-acf(tsx,lag.max=lag.max,plot=F, na.action = na.pass)$acf[-1,1,1]
  ryy<-acf(tsy,lag.max=lag.max,plot=F, na.action = na.pass)$acf[-1,1,1]
  Nstar<-N/(1+2*sum(rxx*ryy))
  ifelse(Nstar<N,Nstar,N)
}

cc <- function(x, y, k, trans = I, pyper = pyper, alpha = 0.9, method = method, nboot = 200L) {
  plyr::ldply(seq(1, nrow(responses) - k + 1), function(i) {
    xdata <-  trans(predictors[i:(k + i),x])
    ydata <-  trans(responses[i:(k + i),y])
    
    if (method == "pearson") {
      m <- cor.test(xdata, ydata, use = "complete.obs", method = method)
      r <- m$estimate[[1]]
      if (!pyper) {
        l <- m$conf.int[[1]]
        u <- m$conf.int[[2]]
      } else {
        n.star <- autocorr_pyper(N = length(xdata), tsx = xdata, tsy = ydata)
        l <- psychometric::CIr(r=r, n = n.star, level = alpha)[1]
        u <- psychometric::CIr(r=r, n = n.star, level = alpha)[2]
      }
    } else { # not pearson
      f <- function(data, i) cor(data[i, 1], data[i, 2], method = method)
      b <- boot::boot(cbind(xdata, ydata), f, nboot)
      ci <- tryCatch({boot::boot.ci(b)}, error = function(e) NA)
      if (!is.na(ci)) {
        l <- ci$bca[[4]]
        u <- ci$bca[[5]]
      } else {
        l <- -1
        u <- 1
      }
      r <- cor.test(xdata, ydata, use = "complete.obs", method = method)$estimate[[1]]
    }
    data.frame(x, y, r, i, l, u)
  })
}


#manipulate({
  k=11
  pyper=TRUE
  method="pearson"
  
  ccnames <- expand.grid(x = names(predictors), y = names(responses), k = k)
  out <- plyr::mdply(ccnames, cc, trans = I, pyper = pyper, method = method)
  
  out <- out %>% mutate(sig = ifelse(u < 0, "low", ifelse(l > 0, "high", "none")))
  
  ggplot(out, aes(i+k/2+1979, r)) + 
    geom_line() +
    #geom_ribbon(aes(x=c(1984,1994),ymin =-1 , ymax =1 ), alpha = 0.2, inherit.aes=F) +
    geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.2) +
    facet_grid(x~y) +
    geom_hline(yintercept = 0, lty = 2, colour = "grey50") +
    #geom_vline(xintercept=c(1994.5,2002.5), lty=2, colour="green")+
    geom_point(data = filter(out, sig != "none"), aes(color = sig)) +
    #theme_light() + 
    guides(colour = FALSE) +
    labs(x="Year", y="Correlation Coefficient") +
    theme_sleek()
  #theme(panel.grid.major = element_blank(),
     # panel.grid.minor = element_blank())
#}, k = slider(4, 33, 11), method = picker("pearson", "spearman", "kendall"), pyper = picker(TRUE, FALSE))


predictors %>% 
  mutate(year = seq_len(nrow(predictors))+1980) %>%
  reshape2::melt(id = "year") %>%
    ggplot(aes(year, value)) +
      geom_line() +
      facet_wrap(~variable, ncol=3, scales = "free_y")
