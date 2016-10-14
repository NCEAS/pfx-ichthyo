d <- readRDS("correlation-matrix.rds")

library(tidyverse)
predictors <- select(d, ENSO:Pollock_F_SSB_thousandt)
responses <- select(d, trend1:IchSW)

autocorr_pyper<-function(N,tsx,tsy) {
  lag.max<-min(length(tsx),length(tsy))/4
  rxx<-acf(tsx,lag.max=lag.max,plot=F, na.action = na.pass)$acf[-1,1,1]
  ryy<-acf(tsy,lag.max=lag.max,plot=F, na.action = na.pass)$acf[-1,1,1]
  Nstar<-N/(1+2*sum(rxx*ryy))
  ifelse(Nstar<N,Nstar,N)
}

cc <- function(x, y, k, trans = I, pyper = FALSE, alpha = 0.9, method = "pearson", nboot = 200L) {
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

library(manipulate)
manipulate({
  ccnames <- expand.grid(x = names(predictors), y = names(responses), k = k)
  out <- plyr::mdply(ccnames, cc, trans = I, pyper = pyper, method = method)
  
  out <- out %>% mutate(sig = ifelse(u < 0, "low", ifelse(l > 0, "high", "none")))
  
  ggplot(out, aes(i+1980, r)) + 
    geom_line() +
    geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.2) +
    facet_grid(x~y) +
    geom_hline(yintercept = 0, lty = 2, colour = "grey50") +
    geom_point(data = filter(out, sig != "none"), aes(color = sig)) +
    theme_light() + guides(colour = FALSE) +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}, k = slider(4, 25, 12), method = picker("pearson", "spearman", "kendall"), pyper = picker(TRUE, FALSE))

d %>% mutate(year = seq_len(nrow(d))) %>%
  reshape2::melt(id = "year") %>%
    ggplot(aes(year, value)) +
      geom_line() +
      facet_wrap(~variable, scales = "free_y")
