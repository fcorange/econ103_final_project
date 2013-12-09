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

#################### Edward's 2nd DRAFT functions STARTs #################### 

## Caution: Inputs are assumed reasonable ##
## Also, function must be defined s.t. mode is in D ##

k<-20
xlb <- -1 # lowerbound
xub <- 1 # upperbound
n<-100

######  define g ######  
g <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given test function as standard normal
# h <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given function as standard normal
# body(h) <- deriv(body(h), "x")
# h <- function(x)  exp(x)/((1+exp(x))^2) # given function as standard normal

######  log g = h ######  
composite <- function(f,g) function(...) f(g(...))
log <- function(x) log(x);
h <- composite(log, g)
body(h) <- deriv(body(h), "x")

######  evenly initialize  find T_k in D ######  
compute_T_k <- function(){ # draw k integers from domain D
  T_k <- rep(0,k) # initialize global T_k
  T_k[1] <- xlb
  T_k[k] <- xub
  for (i in 2:(k-1)){
    T_k[i] <- T_k[i-1] + (xub-xlb)/(k-1)
  }
  return(T_k)
}
T_k <- compute_T_k() # update global variable T_k

###### function evaluated at T_k  ######  
h_func <- function(T_k, h){
  return(sapply(T_k, h))
}
h_k <- h_func(T_k, h) # function evaluated at T_k

###### h derivative  ######  
grad <- function(T_k,h){
  h_k_prime <- 0
  for (i in 1:length(T_k)){  # derivatives evaluated at T_k
    h_k_prime[i] <- unlist(attributes(h(T_k[i])))
  }
  return(h_k_prime)
}
h_k_prime <- grad(T_k, h)

###### z_k coordinates ######  
compute_z_k <- function(T_k,h_k,h_k_prime){  # tangent line intersections
  z_k <- rep (0, k-1)
  for (i in 1:k-1){
    z_k[i] <- (h_k[i+1]-h_k[i]-T_k[i+1]*h_k_prime[i+1]+T_k[i]*h_k_prime[i])/(h_k_prime[i]-h_k_prime[i+1])
  }
  return(z_k)
}
z_k <- compute_z_k(T_k,h_k,h_k_prime)

###### upper bound ######  
# returns a function u_x(x)
# If we want u_x(x) evaluated at x*, call u_x(h_k,h_k_prime,z_k,T_k)(x)
# Note: not u_x(x*)(x*)
compute_u_k <- function(h_k,h_k_prime,z_k,T_k,x) {
  if ((T_k[1] <= x) && (x < z_k[1])){
    return (function(x) (h_k[1] + (x-T_k[1])*h_k_prime[1]))
  }
  for (i in 1:(length(T_k)-2)){
    if(z_k[i] <= x && x < z_k[i+1]){
      return (function(x) (h_k[i+1] + (x-T_k[i+1])*h_k_prime[i+1]))
    }
  }
  if ((z_k[k-1] <= x) && (x <= T_k[k])){
    return (function(x) (h_k[k] + (x-T_k[k])*h_k_prime[k]))
  }
}
compute_u_k(h_k,h_k_prime,z_k,T_k,1)(1) # 1st 1 used in if statements, 2nd 1 used as an evaluation point

###### lower bound ######  
compute_l_k <- function(T_k,h_k,x) {
  for (i in 1:(length(T_k)-1)){
    if(T_k[i] <= x && x < T_k[i+1]){
      return(function(x) ((T_k[i+1]-x)*h_k[i]+(x-T_k[i])*h_k[i+1])/(T_k[i+1]-T_k[i]))
    }
  }
  if (x == T_k[i+1]){
    return(function(x) ((T_k[k]-x)*h_k[k-1]+(x-T_k[k-1])*h_k[k])/(T_k[k]-T_k[k-1]))
  }
}
compute_l_k(T_k,h_k,1)(1) # 1st 1 used in if statements, 2nd 1 used as an evaluation point


###### S_k transformation where x sampled from s_x ######  
### Need a function for calculating area under s_x(x)

### Function A computes area of trapezoid indexed between i-1 and i in z_k
A <- function(i, h_k, h_k_prime, z_k, T_k){ 
  composite<-function(f,g) function(...) f(g(...))
  f<-function(x) compute_u_k(h_k,h_k_prime,z_k,T_k,T_k[i])(x) # evaluate u_k at x using the T_k[i] segment 
  g<-function(x) exp(x)
  exp_u_k <- composite(g,f)
  
  if (i==1){
    return (integrate(exp_u_k,T_k[1],z_k[1])$value)
  }
  else if (i==length(T_k)){
    return (integrate(exp_u_k,z_k[length(z_k)],T_k[length(T_k)])$value)
  }
  else{
    return (integrate(exp_u_k,z_k[i-1],z_k[i])$value)
  }
}

###### A_k vector of trapezoid area ###### 
A_k <- 0
for (i in 1:length(T_k)){
  A_k[i] <- A(i,h_k, h_k_prime, z_k, T_k)
}
## cumulative area, starting from zero
cumArea<-c(0,cumsum(A_k))

### S_k returns a function specified in the slides
S_k <- function(x, h_k, h_k_prime, z_k, T_k, A_k) { 
  composite<-function(f,g) function(...) f(g(...))
  f<-function(x) compute_u_k(h_k,h_k_prime,z_k,T_k,x)(x)
  g<-function(x) exp(x)/sum(A_k)
  return(composite(g,f))
}
S_k(.5, h_k, h_k_prime, z_k, T_k, A_k)(.5) # S_k valued at .5

#################### Edward's 2nd DRAFT functions ENDs #################### 


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

# s_k(x) gives me the function
# s_k(x)(x) gives me the value evaluated at x

# Axillary function
cdf_u<-function(xj,hprimej,cumAreaj,temp,h_k, h_k_prime, z_k, T_k, A_k){ #CREATE A DATAFRAME??
  # 1/b*s_k(x)
  cdf<-1/hprimej*S_k(xj,h_k, h_k_prime, z_k, T_k, A_k)
  return function(x) cdf+cumAreaj-temp
}


sample_val <- function(T_k,cumArea,h_k_prime,z_k,...) {  #Cindy
  # sample x* with p(x) = Uk(x); CDF(x*)=temp_u
  temp<-runif(1)
  # x_star between T_k[1], T_k[k]
  for (i in 1:(k-1)){
    if(temp<cumArea[2]){
      x_star<-uniroot(cdf_u(T_k[1],h_k_prime[1],cumArea[1],temp,h_k, h_k_prime, z_k, T_k, A_k), lower =T_k[1], upper =z[1])[1]
      break
    }else if(temp>=cumArea[i] && temp<cumArea[i+1]){
      x_star<-uniroot(cdf_u(T_k[i],h_k_prime[i],cumArea[i],temp,h_k, h_k_prime, z_k, T_k, A_k), lower =z[i-1], upper =z[i])[1]
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
  ### QUESTION: Why do we need to evaluate h(x*) and h'(x*) here?
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

<<<<<<< HEAD
=======
check_input <- function(k,g,n,xlb,xub) {
  if ((k<=0)||(k%%1!=0)) {
     print("Number of points to sample must be a positive integer")
     return(FALSE)
  }
  else {
    if ((n<=0)||(n%%1!=0)) {
      print("Sample size must be a positive integer")
      return(FALSE)
    }
    else {
      if (xlb>=xub) {
        print("Lower bound of domain must be smaller than upper bound")
        return(FALSE)
      }
      else {
        if (typeof(g)!="closure") {
          print("Input density is not a function")
          return(FALSE)
        }
        else {
          return(TRUE)
        }
      }
    }
  }
}




>>>>>>> 2d6e913f354fdd669b10b6c47b424a536734f511
