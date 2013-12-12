### Tests for log-concave functions ###
g <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given test function is standard normal
samp01<-ARS(g,n=20000,xlb=-Inf,xub=Inf)
g <- function(x)   exp(-(x-3)^2/2) #unnormalized normal kernel w/ mean 3
samp001<-ARS(g,n=20000,xlb=-1,xub=Inf)
g<-function(x) 0.25*x*exp(-x/2) #gamma(2,2) domain(0,Inf)
samp02<-ARS(g,n=20000,xlb=.01,xub=Inf)
g<-function(x) 1/16*x^2*exp(-x/2) #gamma(3,2)
samp03<-ARS(g,n=20000,xlb=.01,xub=Inf)
g<-function(x) x^2*exp(-x/2) #unnormalized gamma(3,2)
samp003<-ARS(g,n=20000,xlb=.01,xub=Inf)
g<-function(x) x*(1-x)/beta(2,2) #beta(2,2)  domain (0,1)
samp04<-ARS(g,n=20000,xlb=.01,xub=.99)
g<-function(x) x^4*exp(-x) #unnormalized gamma
sample05<-ARS(g,n=20000,xlb=.01,xub=Inf)

### Tests for BAD functions (not log-concave) ###


### PLAYGROUND ###
xval <- seq(-1,1,0.001)

plot(xval, h(xval), pch=".")
lines(xval, u_z)
lines(xval,l_z)
u_z <- vector()
l_z <- vector()
for (i in 1:length(xval)){
  u_z[i]<-compute_u_k(data,xval[i])(xval[i])
  l_z[i]<-compute_l_k(T_k,h_k,xval[i])(xval[i])
}

for (i in 1:length(T_k)){
  
}
points(z_k, u_z)
