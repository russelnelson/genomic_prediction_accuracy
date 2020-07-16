chrsize <- 10000000
mu <- 1e-8
Ne <- 1000
r <- 2*1e-8

gc <- 2000
chrnum <- 10

theta <- 4*Ne*mu*chrsize

rho <- 4*Ne*r*(chrsize - 1)

options(scipen = 999)

call <- paste0("scrm ", gc, " ", chrnum, " -t ", theta, " -r ", rho, " ", chrsize)
call

