# User input

k # Number of points to sample from the user
g(x) # Density
n # Sample size desired
xlb # lower bound of x
xub # upper bound of x

### OUR PARENT FUNCTION ARS() ###
ARS<-function(k,g,n,xlb,xub){
  #initialization
  
  while(length(sample)<n){ # while sample size n not reached, sample+update
    #sampling
    sample_point<-sample_val()
    squeeze<-squeeze_test(sample_point)
    if (squeeze==F){
      reject<-rejection_test(sample_point)
    }
    if(squeeze==T || reject==T){
      sample<-c(sample,sample_point)
    }
    
    #updating
    if (squeeze==F){
      #update
    }
  }
  
}
######


# Function
h(x)


T_k <- vector()


u_k <- function(x) { #Edward
  
}

s_k <- function(x) { #Edward
  
}

l_k <- function() { #Edward
  
}

get_z <- function(T_k) { #Zixiao
  return(z)
}

### Cindy's DRAFT functions ###

# Assuming u_k(x) looks like:
# u_k <- if(x >= x1 & x <x2) return u1_k(x;x1),
# elseif(x >=x2 & x < x3,  u2_k(x;x2), ...
# So u_k(x) gives me the function
# u_k(x)(x) gives me the value of u_k evaluated at x

# Axillary function
cdf_u<-function(xj,temp){
  # depends on h(x); cdf-temp
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
#area of trapozoid
a<-function(i) return((h(T_k[i])+h(T_k[i-1]))*(T_k[i]-T_k[i-1])*0.5)
Area<-unlist(lapply(2:k,a))
cumArea<-c(0,cumsum(Area))

#dummy cumArea
cumArea<-seq(0.2,1,0.2)
sample_val <- function(T_k,cumArea) {  #Cindy
  # sample x* with p(x) = Uk(x); CDF(x*)=temp_u
  temp<-runif(1)
  # x_star between T_k[1], t_k[k]
  for (i in 1:(k-1)){
    if(temp>=cumArea[i] && temp<cumArea[i+1]){
      x_star<-uniroot(cdf_u(T_k[i],temp), lower =T_k[i], upper =T_k[i+1])[1]
      break
    }
  }
  # sample u* from uniform(0,1)
  u_star<-runif(1)
  return(c(x_star, u_star))
}

squeeze_test <- function(x_star, u_star) { #Cindy
  test<-exp(l_k(x_star)(x_star)-u_k(x_star)(x_star))
  Boolean<-ifelse(u_star<=test,T,F)
  # T=accept, F=reject
  ### QUESTION: Why do we need to evaluate h(x*) and h'(x*) here???
  return(Boolean)
}

rejection_test <- funcion(x_star, u_star) { #Cindy 
  test<-exp(h(x_star)(x_star)-u_k(x_star)(x_star))
  Boolean<-ifelse(u_star<=test,T,F)
  return(Boolean)
}

update <- function(x_star) { #Zixiao
  
}





