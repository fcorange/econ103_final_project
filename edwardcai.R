## Edward Cai's directory: /Users/edwsurewin/econ103_final_project ##
## setwd("/Users/edwsurewin/econ103_final_project")

## Caution: Inputs are assumed reasonable ##
## Also, function must be defined s.t. mode is in D ##

k<-20
vec <- c(-2,2) # global Domain
#vec <- c(-1,1) # global D
T_k <- rep(0,k) # global T_k

# func and deriv #
#func <- function(x)  ((2*pi)^-0.5)*exp(-(x)^2/2) # given function as standard normal
func <- function(x)  exp(x)/((1+exp(x))^2) # given function as standard normal
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
u_k <- function(x) {
  if ((T_k[1] <= x) && (x < z[1])){
    return(h_x[1] + (x-T_k[1])*h_prime_x[1])
  }
  for (i in 1:(length(T_k)-2)){
    if(z[i] <= x && x < z[i+1]){
      return(h_x[i+1] + (x-T_k[i+1])*h_prime_x[i+1])
    }
  }
  if ((z[k-1] <= x) && (x <= T_k[k])){
    return(h_x[k] + (x-T_k[k])*h_prime_x[k])
  }
}

# S_k transformation where x sampled from u_k #
S_k <- function(x) {
    area <- 0 
    area <- area + (z[1]-T_k[1])*(u_k(z[1])+func(T_k[1])[1])/2
    area <- area + (T_k[k]-z[k-1])*(u_k(z[k-1])+func(T_k[k])[1])/2
    for (i in 1:(length(T_k)-2)){
      area <- area + (z[i+1]-z[i])*(u_k(z[i])+u_k(z[i]))/2
    }
    return(exp(u_k(x))/area)
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



############## testing below, not part of the submission

# solve for a root of function G at G=-5
G <- function(k)  (2*(1 - k)*(1 + 2*k)^.5)/(1+3*k) 
GZero <-function(k,G0) G(k)-G0 
plot(1:100,G(1:100)) 
uniroot(GZero,c(0,100),G=-5) 
uniroot(GZero,c(0,Inf),G=-5) 
uniroot(GZero,c(0,10^20),G=-5) 

##################

# optimize() to find mode

func <- function(x) { # given function as standard normal
  ((2*pi)^-0.5)*exp(-(x)^2/2)
}

vec <- c(-10^(17),300) # vector to indicate D
vec <- c(-10^(17),100) # vector to indicate D

compute_T_k <- function(k){ # user k int
  T_k <- rep(0,k)
  x_max <- optimize(func, interval=vec, maximum=TRUE)[1]
  x_max
}
compute_T_k(3)


######################

compute_T_k <- function(k){ # user k int
  T_k <- rep(0,k)
  x_max <- optimize(func, interval=vec, maximum=TRUE)[1]
  x_max
}


compute_T_k <- function(func_expr,vec,k){ # function expression, domain vector, user k int
  T_k <- rep(0,k)
  
  func <- function(x) {
    func_expr
  }
  x_max <- optimize(func, interval=c(-1,1), maximum=TRUE)[1]
  ## int <- integrate(eval(func),2,3)
  ## names(x_max) <- NULL
  
  if (vec[1]==(-Inf)){
    if (vec[2]==Inf){
      
    }
    else{
      
      T_k[k] <- vec[2]
    }
  }
  else{
    T_k[1] <- vec[1]
    T_k[k] <- vec[2]
  }
  
  return (x_max)
}

#D(log(x^(3-1)*(1-x)^(6-1)), 'x') # e.g. log beta pdf
#eval(func)
## func_expr <- expression(((2*pi)^-0.5)*exp(-(x)^2/2))
