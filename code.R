# ---- User inputs
k     # Number of points to sample from the user
g     # Density (can be unnormalized)
n     # Sample size desired
xlb   # Lower bound of domain
xub   # Upper bound of domain

# ---- Library for calculating derivatives
library("numDeriv")

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

g<-function(x) x^2*exp(-x/2)

### OUR PARENT FUNCTION ARS() ###
ARS<-function(k,g,n,xlb,xub){
  # Initialization
  h <- function(x) (return(log(g(x))))            # Function h = log(g)
  mode <- optim(0,g, upper=xub, lower =xlb, control=list(fnscale=-1),method="L-BFGS-B")[1]
  if (xub == Inf){
    xub <- mode$par + 10 # arbitary number bounds mode if D specified unbounded
  }
  if (xlb == -Inf){
    xlb <- mode$par - 10
  }
  print(c(mode,xlb,xub))
  sample <- vector()                              # Vector that stores all the sampled points (different from T_k)
  T_k <- compute_T_k(k,xlb,xub)                   # Initialize the evenly spaced x points on domain D
  h_k <- compute_h_k(T_k, h)                      # Obtain the values of h evaluated at T_k 
  h_k_prime <- grad(h, T_k)                       # Obtain the derivative of h evaluated at T_k
  z_k <- compute_z_k(T_k, h_k, h_k_prime)         # Obtain the intersections of tangent lines (k-1 elements)
  z_k <- c(z_k,tail(T_k,n=1))
  A_k <- rep(0,length=k)                          # Initialize the vector of area
  # Create a data frame
  data <- data.frame(T_k,h_k,h_k_prime,z_k,A_k)
  names(data) <- c("T_k","h_k","h_k_prime","z_k","A_k")
  # u_k now takes in the data frame as an input
  #u_k <- vector()
  #for (i in 1:length(T_k)){
  #u_k[i]<-compute_u_k(data,T_k[i])(T_k[i])
  #}
  # Obtain the upper bound on T_k. Note that compute_u_k gives a function
  
  
  #l_k <- vector()
  #for (i in 1:length(T_k)){
  #l_k[i]<-compute_l_k(T_k,h_k,T_k[i])(T_k[i])
  #}
  # Obtain the lower bound on T_k
  
  for (i in 1:length(T_k)) {                    
    A_k[i] <- A(i, data)
  }
  data$A_k <- A_k  # update A_k in data frame
  
  
  
  # Sampling
  while(length(sample) < n) { # While sample size n is ont reached, keep sample & update
    cumArea <- cumsum(data$A_k)    # Cumulative area
    sample_point <- sample_val(data,cumArea)
    x_star<-sample_point[1]
    l_xstar<-compute_l_k(data$T_k,data$h_k,x_star)(x_star)
    u_xstar<-compute_u_k(data,x_star)(x_star)
    squeeze <- squeeze_test(sample_point[1],sample_point[2],l_xstar,u_xstar)
    
    if (squeeze == F){
      h_xstar <- h(x_star) 
      reject <- rejection_test(sample_point[1],sample_point[2],u_xstar,h_xstar)
      if(reject == T){
        sample <- c(sample, x_star)
      }
    }else{
      sample <- c(sample, x_star)
    }
    
    
    ##### Make sure to update the dataframe #########
    
    # Updating
    if (squeeze==F){
      data<-update(sample_point[1], data, u_k, l_k,h)
    }
  }
 
  return(sample)
}
######



#################### Edward's 2nd DRAFT functions STARTs #################### 

## Caution: Inputs are assumed reasonable ##
## Also, function must be defined s.t. mode is in D ##

k<-10
xlb <- -Inf # lowerbound
xub <- Inf # upperbound
n<-100

######  define g ######  
g <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given test function as standard normal
g<-function(x) 0.25*x*exp(-x/2) #gamma(2,2)
g<-function(x) 1/16*x^2*exp(-x/2) #gamma(3,2)
g<-function(x) x*(1-x)/beta(2,2) #beta(2,2)
g<-function(x) (1/768)*x^4*exp(-x/2)
# h <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given function as standard normal
# body(h) <- deriv(body(h), "x")
# h <- function(x)  exp(x)/((1+exp(x))^2) # given function as standard normal

######  log g = h ######  
composite<-function(f,g) function(...) f(g(...))
f <- function(x) log(x);
h <- composite(f, g)

##body(h) <- deriv(body(h), "x", func = TRUE)
##(g) <- deriv(body(g), "x")
##h(1)

######  checking Inf ######  
mode <- optim(0,h, upper=xub, lower =xlb, control=list(fnscale=-1),method="L-BFGS-B")[1]
if (xub == Inf){
  xub <- mode$par + 10 # arbitary number bounds mode if D specified unbounded
}
if (xlb == -Inf){
  xlb <- mode$par - 10
}

######  evenly initialize  find T_k in D ######  
compute_T_k <- function(k,xlb,xub){ # draw k integers from domain D
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
compute_h_k <- function(T_k, h){
  return(h(T_k))
}
h_k <- compute_h_k(T_k, h) # function evaluated at T_k

###### h derivative  ######  (we don't need this function any more)
install.packages("numDeriv")
gradient <- function(T_k,h){
  library(numDeriv)
  h_k_prime <- 0b 
  for (i in 1:length(T_k)){  # derivatives evaluated at T_k
    h_k_prime[i] <- grad(h,T_k[i])
  }
  return(h_k_prime)
}
h_k_prime <- gradient(T_k, h)

###### z_k coordinates ######  
compute_z_k <- function(T_k, h_k, h_k_prime){  # tangent line intersections
  z_k <- rep (0, length(T_k)-1)
  for (i in 1:(length(T_k)-1)){
    z_k[i] <- (h_k[i+1]-h_k[i]-T_k[i+1]*h_k_prime[i+1]+T_k[i]*h_k_prime[i])/(h_k_prime[i]-h_k_prime[i+1])
  }
  return(z_k)
}
z_k <- compute_z_k(T_k,h_k,h_k_prime)

###### upper bound ######  
# returns a function u_x(x)
# If we want u_x(x) evaluated at x*, call u_x(h_k,h_k_prime,z_k,T_k)(x)
# Note: not u_x(x*)(x*)
compute_u_k <- function(data,x) {
  with(data,{
    if ((T_k[1] <= x) && (x < z_k[1])){
      return (function(x) (h_k[1] + (x-T_k[1])*h_k_prime[1]))
    }
    for (i in 1:(length(T_k)-2)){
      if(z_k[i] <= x && x < z_k[i+1]){
        return (function(x) (h_k[i+1] + (x-T_k[i+1])*h_k_prime[i+1]))
      }
    }
    if ((z_k[length(T_k)-1] <= x) && (x <= T_k[length(T_k)])){
      return (function(x) (h_k[length(T_k)] + (x-T_k[length(T_k)])*h_k_prime[length(T_k)]))
    }
  })
}
compute_u_k(data,1)(1) # 1st 1 used in if statements, 2nd 1 used as an evaluation point

###### lower bound ######  
compute_l_k <- function(T_k,h_k,x) {
  for (i in 1:(length(T_k)-1)){
    if((T_k[i] <= x) && (x < T_k[i+1])){
      return(function(x) ((T_k[i+1]-x)*h_k[i]+(x-T_k[i])*h_k[i+1])/(T_k[i+1]-T_k[i]))
    }
  }
  if (x == T_k[length(T_k)]){
    return(function(x) ((T_k[length(T_k)]-x)*h_k[length(T_k)-1]+(x-T_k[length(T_k)-1])*h_k[length(T_k)])/(T_k[length(T_k)]-T_k[length(T_k)-1]))  }
}
compute_l_k(T_k,h_k,1)(1) # 1st 1 used in if statements, 2nd 1 used as an evaluation point


###### S_k transformation where x sampled from s_x ######  
### Need a function for calculating area under s_x(x)

### Function A computes area of trapezoid indexed between i-1 and i in z_k
A <- function(i, data){ 
  with(data,{
    a<-h_k_prime[i]
    b<-h_k[i]-T_k[i]*a
    if(a==0){
      return(exp(b)*(z_k[i]-z_k[i-1]))
    }else if (i==1){
      return (-(exp(a*T_k[1]+b)-exp(a*z_k[1]+b))/a)
    }else{
      return (-(exp(a*z_k[i-1]+b)-exp(a*z_k[i]+b))/a)
    }
  })
}   


#################### Edward's 2nd DRAFT functions ENDs #################### 


compute_z_k2 <- function(T_k) { #Zixiao
  z_k <- vector()
  for (i in 2:(length(T_k))) {
    z_k[i] <- (h(T_k[i+1]) - h(T_k[i]) - T_k[i+1] * grad(h, T_k[i+1]) + T_k[i] * grad(h, T_k[i])) / (grad(h, T_k[i]) - grad(h, T_k[i+1]))
  }
  z_k[1] <- T_k[1]
  z_k[(length(T_k))+1] <- tail(T_k,n=1)
  return(z_k)
}

### Cindy's DRAFT functions ###

sample_val <- function(data,cumArea) {  #Cindy
  with(data, {
    # sample x* with p(x) = Uk(x); CDF(x*)=temp_u
    temp<-runif(1)*sum(A_k)
    k<-length(T_k)
    for (i in 1:k){
      a<-h_k_prime[i]
      b<-h_k[i]-T_k[i]*a
      if(i==1 && temp<cumArea[1]){
        x_star<-(log(exp(a*T_k[1]+b)+temp*sum(A_k)*a)-b)/a
        break
      }else if(temp>=cumArea[i-1] && temp<cumArea[i]){
        if (a==0){
          x_star<-(temp-cumArea[i-1])*sum(A_k)/exp(b)+z_k[i-1]
          break
        }
        
        x_star<-(log(exp(a*z_k[i-1]+b)+(temp-cumArea[i-1])*sum(A_k)*a)-b)/a
        break
      }
    }
    # sample u* from uniform(0,1)
    u_star<-runif(1)
    return(c(x_star, u_star))
  })
}

squeeze_test <- function(x_star, u_star,l_xstar,u_xstar) { #Cindy
  test<-exp(l_xstar-u_xstar)
  Boolean<-ifelse(u_star<=test,T,F)
  # T=accept, F=reject
  return(Boolean)
}

rejection_test <-function(x_star, u_star,u_xstar,h_xstar) { #Cindy 
  test<-exp(h_xstar-u_xstar)
  Boolean<-ifelse(u_star<=test,T,F)
  return(Boolean)
}


# ---- Function to check if the upper and lower bounds are not actually bounding the density, i.e. the log-concavity is violated
check_concave <- function(u_k, l_k) {
  return(sum((sum(u_x < h(u_x))==0) + (sum(l_k > h(l_k))==0)) == 2)
}

update <- function(x_star, data, u_k, l_k,h) { #Zixiao
  with(data,{
    T_k <- sort(append(data$T_k, x_star))
    position <- (which(T_k == x_star) - 1)
    h_k <- append(data$h_k, compute_h_k(x_star, h), after = position)
    h_k_prime <- append(data$h_k_prime, grad(h, x_star), after = position)
    z_k <- compute_z_k(T_k, h_k, h_k_prime)
    z_k <- c(z_k,tail(T_k,n=1))
    #u_k <- append(u_k, compute_u_k(data,x_star)(x_star), after = position)
    #l_k <- append(l_k, compute_l_k(data$T_k,data$h_k,x_star)(x_star), after = position)
    A_k <- rep(0,length=length(T_k))
    data <- data.frame(T_k,h_k,h_k_prime,z_k,A_k)
    for (i in 1:length(T_k)) {                    
      A_k[i] <- A(i, data)
    }
    data$A_k <- A_k
    
    return(data)
  })
} 

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
