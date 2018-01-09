wae <- read.table("WAE Clean.txt", header=T)
attach(wae)
names(wae)

ll <- log10(length)
lw <- log10(weight)
plot(lw~ll)

#linear
l.mod <- lm(lw~ll)
summary(l.mod)
anova(l.mod)

#non-linear
nl.mod <- nls(lw~(alpha*(ll)^2)+(beta*(ll))+charlie, start=list(alpha=-5.73, beta=3.27, charlie=0.01), trace=T)
summary(nl.mod)
#R squared
R2(ll, lw, nl.mod)

#quantile
q75.mod <- nlrq(lw~(alpha*(ll)^2)+(beta*(ll))+charlie, start=list(alpha=-5.73, beta=3.27, charlie=0.01), tau=0.75, trace=T)
q50.mod <- nlrq(lw~(alpha*(ll)^2)+(beta*(ll))+charlie, start=list(alpha=-5.73, beta=3.27, charlie=0.01), trace=T)
summary(q50.mod)
summary(q75.mod)

 plot(lw~ll)
  abline(l.mod, col="red", lty="dashed")
  moda=seq(2.1, 3.0, by=0.001)
  lines(moda, predict(nl.mod, list(ll = moda)), col="red")
  lines(moda, predict(q75.mod, list(ll=moda)), col="purple")