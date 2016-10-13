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

cc <- function(x, y, k, trans = I, pyper = FALSE, alpha = 0.9) {
  plyr::ldply(seq(1, nrow(responses) - k + 1), function(i) {
  xdata <-  trans(predictors[i:(k + i),x])
  ydata <-  trans(responses[i:(k + i),y])
  m <- cor.test(xdata, ydata, use = "complete.obs")
  r <- m$estimate[[1]]
  if (!pyper) {
  l <- m$conf.int[[1]]
  u <- m$conf.int[[2]]
  } else {
   n.star <- autocorr_pyper(N = length(xdata), tsx = xdata, tsy = ydata)
   l <- psychometric::CIr(r=r, n = n.star, level = alpha)[1]
   u <- psychometric::CIr(r=r, n = n.star, level = alpha)[2]
  }
  data.frame(x, y, r, i, l, u)
})
}

library(manipulate)
manipulate({
ccnames <- expand.grid(x = names(predictors), y = names(responses), k = k)
out <- plyr::mdply(ccnames, cc, trans = I, pyper = TRUE)

out <- out %>% mutate(sig = ifelse(u < 0, "low", ifelse(l > 0, "high", "none")))

ggplot(out, aes(i, r)) + 
  geom_line() +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.2) +
  facet_grid(x~y) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey50") +
  geom_point(data = filter(out, sig != "none"), aes(color = sig)) +
  theme_light()
}, k = slider(4, 25, 8))


# d %>% mutate(year = seq_len(nrow(d))) %>% 
#   reshape2::melt(id = "year") %>% 
#     ggplot(aes(year, value)) + 
#       geom_line() +
#       facet_wrap(~variable, scales = "free_y")
