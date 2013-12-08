# User input

k # Number of points to sample from the user
g(x) # Density (can be unnormalized)
n # Sample size desired
xlb # lower bound of x
xub # upper bound of x

g<-function(x) exp(-x^2/2)

### OUR PARENT FUNCTION ARS() ###
ARS<-function(k,g,n,xlb,xub){
  #initialization
  # Perform input check
  # ---- n and k need to be a postiive integer
  # ---- g?
  # ---- xlb should be smaller than xub
  
  
  while(length(sample)<n){ # while sample size n not reached, sample+update
    #sampling
    sample_point<-sample_val()
    squeeze<-squeeze_test(sample_point)
    if (squeeze==F){
      reject<-rejection_test(sample_point)
    }
    if(squeeze==T || reject==T){
      sample<-c(sample,sample_point[1])
    }
    
    #updating
    if (squeeze==F){
      #update
    }
  }
  
}
######


# Function
# ---- For now I'm assuming g is not vectorized
h <- function(x) { #Zixiao
  return(log(g(x)))
}

# ---- Calculate Derivative of a function ----
library("numDeriv")
# ----

### Edward's DRAFT functions START ###

k<-20
vec <- c(-2,2) # global Domain
#vec <- c(-1,1) # global D
T_k <- rep(0,k) # global T_k

# func and deriv #
#func <- function(x)  ((2*pi)^-0.5)*exp(-(x)^2/2) # given function as standard normal
func <- function(x)  exp(x)/((1+exp(x))^2) # given function as logistic
body(func) <- deriv(body(func), "x")

# evenly find T_k in D #
compute_T_k <- function(k){ # draw k integers from domain D
  T_k[1] <- vec[1]
  T_k[k] <- vec[2]
  for (i in 2:(k-1)){
    T_k[i] <- T_k[i-1] + (vec[2]-vec[1])/(k-1)
  }
  return(T_k)
}

T_k <- compute_T_k(k) # update global variable T_k
h_x <- sapply(T_k, func) # function evaluated at T_k
h_prime_x <- 0
for (i in 1:length(h_x)){  # derivatives evaluated at T_k
  h_prime_x[i] <- unlist(attributes(func(T_k[i])))
}

# z coordinates #
z <- rep (0, k-1) # tangent line intersections
for (i in 1:k-1){
  z[i] <- (h_x[i+1]-h_x[i]-T_k[i+1]*h_prime_x[i+1]+T_k[i]*h_prime_x[i])/(h_prime_x[i]-h_prime_x[i+1])
}

# upper bound #
# returns a function u_x(x)
# If we want u_x(x) evaluated at x*, call u_x(x*)(x*)
u_x <- function(x) {
  if ((T_k[1] <= x) && (x < z[1])){
    return function(x) (h_x[1] + (x-T_k[1])*h_prime_x[1])
  }
  for (i in 1:(length(T_k)-2)){
    if(z[i] <= x && x < z[i+1]){
      return function(x) (h_x[i+1] + (x-T_k[i+1])*h_prime_x[i+1])
    }
  }
  if ((z[k-1] <= x) && (x <= T_k[k])){
    return function(x) (h_x[k] + (x-T_k[k])*h_prime_x[k])
  }
}

# S_k transformation where x sampled from u_x #
### Need a function for calculating area under u_x(x)
S_k <- function(x) {
  area <- 0 
  area <- area + (z[1]-T_k[1])*(u_x(z[1])(z[1])+h(T_k[1])[1])/2
  area <- area + (T_k[k]-z[k-1])*(u_x(z[k-1])(z[k-1])+h(T_k[k])[1])/2
  for (i in 1:(length(T_k)-2)){
    area <- area + (z[i+1]-z[i])*(u_x(z[i])(z[i])+u_x(z[i])(z[i]))/2
  }
  return(exp(u_x(x)(x))/area)
}

# lower bound #
l_k <- function(x) {
  for (i in 1:(length(T_k)-1)){
    if(T_k[i] <= x && x < T_k[i+1]){
      return(((T_k[i+1]-x)*h_x[i]+(x-T_k[i])*h_x[i+1])/(T_k[i+1]-T_k[i]))
    }
  }
  if (x == T_k[i+1]){
    return(((T_k[k]-x)*h_x[k-1]+(x-T_k[k-1])*h_x[k])/(T_k[k]-T_k[k-1]))
  }
}

### Edward's DRAFT functions END ###


compute_z_k <- function(T_k) { #Zixiao
  z_k <- vector()
  for (i in 2:(length(T_k)) {
    z_k(i) <- (h(T_k[i+1]) - h(T_k[i]) - T_k[i+1] * grad(h, T_k[i+1]) + T_k[i] * grad(h, T_k[i])) / (grad(h, T_k[i]) - grad(h, T_k[i+1]))
  }
  z_k[1] <- T_k[1]
  z_k[(length(T_k))+1] <- tail(T_k)
  return(z_k)
}

### Cindy's DRAFT functions ###

# Assuming u_x(x) looks like:
# u_x <- if(x >= x1 & x <x2) return u1_k(x;x1),
# elseif(x >=x2 & x < x3,  u2_k(x;x2), ...
# So u_x(x) gives me the function
# u_x(x)(x) gives me the value of u_x evaluated at x

# Axillary function
cdf_u<-function(xj,zi,i,temp){
  # depends on h(x); cdf-temp
  ux<-u_x(xj)
  cdf<-ux
  return function(x) (x-zi)*(ux(zi)+ux(x))/2+cumArea[i]-temp
  # I put in a dummy function for mow
  if (xj>=1 && xj<2){
    function(x) x-1-temp
  }else if(xj>=2 && xj<3){
    function(x) x-2-temp
  }else if(xj>=3 && xj<4){
    function(x) x-3-temp
  }else if(xj>=4 && xj<5){
    function(x) x-4-temp
  }
}

#dummy functions for testing
k=5
T_k<-1:5
h <- function (x) x
#area of trapozoid (need to reconstruct with z[i]!!!!!!!!!)
a<-function(i) return((h(T_k[i])+h(T_k[i-1]))*(T_k[i]-T_k[i-1])*0.5)
Area<-unlist(lapply(2:k,a)) 
# add the two trapozoids/triangles on both end with given xub,xlb.
cumArea<-c(0,cumsum(Area))

#dummy cumArea
cumArea<-seq(0.2,1,0.2)
sample_val <- function(T_k,cumArea) {  #Cindy
  # sample x* with p(x) = Uk(x); CDF(x*)=temp_u
  temp<-runif(1)
  # x_star between T_k[1], T_k[k]
  for (i in 1:(k-1)){
    if (temp<cumArea[2]){
      x_star<-uniroot(cdf_u(T_k[1],T_k[1],temp), lower =T_k[1], upper =z[1])[1]
      break
    } else if (temp>cumArea[k]){
       # NEED TO BUILD ON THE AREA FUNCTION      
    } else if(temp>=cumArea[i] && temp<cumArea[i+1]){
      x_star<-uniroot(cdf_u(T_k[i],z[i-1],temp), lower =z[i-1], upper =z[i])[1]
      break
    }
  }
  # sample u* from uniform(0,1)
  u_star<-runif(1)
  return(c(x_star, u_star))
}

squeeze_test <- function(x_star, u_star) { #Cindy
  test<-exp(l_k(x_star)-u_x(x_star)(x_star))
  Boolean<-ifelse(u_star<=test,T,F)
  # T=accept, F=reject
  ### QUESTION: Why do we need to evaluate h(x*) and h'(x*) here???
  return(Boolean)
}

rejection_test <- funcion(x_star, u_star) { #Cindy 
  test<-exp(h(x_star)-u_x(x_star)(x_star))
  Boolean<-ifelse(u_star<=test,T,F)
  return(Boolean)
}

# ---- Function to check if the upper and lower bounds are not actually bounding the density, i.e. the log-concavity is violated
check_concave <- function(u_x, l_k) {
  return(sum((sum(u_x < h(u_x))==0) + (sum(l_k > h(l_k))==0)) == 2)
}

update <- function(x_star) { #Zixiao
  T_k <- sort(append(T_k, x_star))
}

