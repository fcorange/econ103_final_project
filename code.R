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


###### function evaluated at T_k  ######  
compute_h_k <- function(T_k, h){
  return(h(T_k))
}



###### z_k coordinates ######  
compute_z_k <- function(T_k, h_k, h_k_prime){  # tangent line intersections
  z_k <- rep (0, length(T_k)-1)
  for (i in 1:(length(T_k)-1)){
    z_k[i] <- (h_k[i+1]-h_k[i]-T_k[i+1]*h_k_prime[i+1]+T_k[i]*h_k_prime[i])/(h_k_prime[i]-h_k_prime[i+1])
  }
  return(z_k)
}


###### upper bound ######  
# returns a function u_x(x)
# If we want u_x(x) evaluated at x*, call u_x(data,T_k)(x)
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


### Function A computes area between i-1 and i in z_k
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



sample_val <- function(data,cumArea) {
  with(data, {
    # sample x*: CDF(x*)=temp-cumArea
    temp<-runif(1)*sum(A_k)
    k<-length(T_k)
    for (i in 1:k){
      a<-h_k_prime[i]
      b<-h_k[i]-T_k[i]*a
      if(i==1 && temp<cumArea[1]){
        x_star<-(log(exp(a*T_k[1]+b)+temp*a)-b)/a
        break
      }else if(temp>=cumArea[i-1] && temp<cumArea[i]){
        if (a==0){
          x_star<-(temp-cumArea[i-1])/exp(b)+z_k[i-1]
          break
        }
        x_star<-(log(exp(a*z_k[i-1]+b)+(temp-cumArea[i-1])*a)-b)/a
        break
      }
    }
    # sample u* from uniform(0,1)
    u_star<-runif(1)
    return(c(x_star, u_star))
  })
}

squeeze_test <- function(x_star, u_star,l_xstar,u_xstar) {
  test<-exp(l_xstar-u_xstar)
  Boolean<-ifelse(u_star<=test,T,F)
  # T=accept, F=reject
  return(Boolean)
}

rejection_test <-function(x_star, u_star,u_xstar,h_xstar) {
  test<-exp(h_xstar-u_xstar)
  Boolean<-ifelse(u_star<=test,T,F)
  return(Boolean)
}


# ---- Function to check if the upper and lower bounds are not actually bounding the density, i.e. the log-concavity is violated
check_concave <- function(data,h) {
  sample_points <- vector()
  for (i in 1:(length(data$T_k)-1)) {
    sample_points[i] <- mean(c(data$T_k[i],data$T_k[i+1]))
  }
  u_k <- vector()
  for (i in 1:length(sample_points)){
    u_k[i]<-compute_u_k(data,sample_points[i])(sample_points[i])
  } 
  l_k <- vector()
  for (i in 1:length(sample_points)){
    l_k[i]<-compute_l_k(data$T_k,data$h_k,sample_points[i])(sample_points[i])
  }
  return(sum((sum(u_k < h(sample_points))==0) + (sum(l_k > h(sample_points))==0)) == 2)
}

update <- function(x_star, data, u_k, l_k,h) { 
  with(data,{
    T_k <- sort(append(data$T_k, x_star))
    position <- (which(T_k == x_star) - 1)
    h_k <- append(data$h_k, compute_h_k(x_star, h), after = position)
    h_k_prime <- append(data$h_k_prime, grad(h, x_star), after = position)
    z_k <- compute_z_k(T_k, h_k, h_k_prime)
    z_k <- c(z_k,tail(T_k,n=1))
    A_k <- rep(0,length=length(T_k))
    data <- data.frame(T_k,h_k,h_k_prime,z_k,A_k)
    for (i in 1:length(T_k)) {                    
      A_k[i] <- A(i, data)
    }
    data$A_k <- A_k
    
    return(data)
  })
} 

check_input <- function(g,n,xlb,xub) {
  if ((n<5)||(n%%1!=0)) {
      warning("Sample size must be a bigger integer")
      return()
    }
  if(xlb>=xub) {
      warning ("Lower bound must be smaller than the upper bound")
      return()
  }
  if (typeof(g)!="closure") {
      warning("Input density is not a function")
      return()
  }
}



### OUR PARENT FUNCTION ARS() ###
ars<-function(g,n,xlb,xub){
  library("numDeriv")
  check_input(g,n,xlb,xub)
  # Initialization
  k<-5
  h <- function(x) (return(log(g(x))))            # Function h = log(g)
  mode <- try(optim(0,h, upper=xub, lower =xlb, control=list(fnscale=-1),method="L-BFGS-B")[1],silent=T)
  if (class(mode)=="try-error"){
    warning("Mode not found in range. Please check your kernel/density or change xlb/xub.")
    return()
  }
  if (xub == Inf){
    xub <- mode$par + 10 # arbitary number bounds mode if D specified unbounded
  }
  if (xlb == -Inf){
    xlb <- mode$par - 10
  }
  
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
  if(any(is.na(z_k))){
    warning("Your kernel is uniform, please sample with runif() instead.")
    return()
  }
  for (i in 1:length(T_k)) {                    
    A_k[i] <- A(i, data)
  }
  data$A_k <- A_k # update A_k in data frame
  
  
  ## Sampling
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
    
    # Updating
    if (squeeze==F){
      data<-update(sample_point[1], data, u_k, l_k,h)
      # Check concavity
      if (check_concave(data,h)==FALSE) {
        warning("Input function not concave!")
        return()
      }
    }
  }
  hist(sample,main="")
  return(sample)
}
######

test<-function(){
  ######  define g ######  
  print("Testing on a standard normal density.")
  g <- function(x)   (((2*pi)^-0.5)*exp(-(x)^2/2)) # given test function is standard normal
  samp01<-ars(g,n=20000,xlb=-Inf,xub=Inf)
  title(main="Sampling from a Standard Normal")
  summary(samp01)
  print("Testing on a standard normal kernel with mean=3.")
  g <- function(x)   exp(-(x-3)^2/2) #unnormalized normal w/ mean 3
  samp001<-ars(g,n=20000,xlb=-1,xub=Inf)
  title(main="Sampling from a Normal mean=3")
  summary(samp001)
  print("Testing on a gamma(2,2) density.")
  g<-function(x) 0.25*x*exp(-x/2) #gamma(2,2)
  samp02<-ars(g,n=20000,xlb=.01,xub=Inf)
  title(main="Sampling from a gamma density.")
  summary(samp02)
  print("Testing on an unmormalized gamma kernel.")
  g<-function(x) x^2*exp(-x/2) #unnormalized gamma(3,2)
  samp003<-ars(g,n=20000,xlb=.01,xub=Inf)
  title(main="Sampling from an unmormalized gamma kernel.")
  summary(samp003)
  print("Testing on a beta(2,2) kernel.")
  g<-function(x) x*(1-x)/beta(2,2) #beta(2,2)  domain (0,1)
  samp04<-ars(g,n=20000,xlb=.01,xub=.99)
  title(main="Sampling from a beta(2,2) kernel.")
  summary(samp04)
  
  ## Non log-concave example ##
  print("Testing on a kernel/density that is not log concave.")
  g <- function(x) 2*(((2*pi)^-0.5)*exp(-(x-4)^2/2)) + (((2*pi)^-0.5)*exp(-(x)^2/2)) 
  sample05<-ars(g,n=20000,xlb=.01,xub=Inf)
  
  ## Bad lower bound and upper bound input test ##
  print("Testing on undefined domain of a kernel/density: beta with xlb and xub out of range.")
  g<-function(x) x*(1-x)/beta(2,2) #beta(2,2)  domain (0,1)
  samp04d<-ars(g,n=20000,xlb=-1,xub=1)
  
  ## Given density is uniform ##
  print("Testing on a uniform kernel/density.")
  g<-function(x) 1
  sampu<-ars(g,n=20000,xlb=-1,xub=1)
  
  ## Bad input ##
  print("Testing on bad input: n not postive integer.")
  g<-function(x) x*(1-x)/beta(2,2) #beta(2,2)  domain (0,1)
  samp04d<-ars(g,n=2.5,xlb=-1,xub=1)
}

