
#Reading data
dat <- read.table("http://home.iitk.ac.in/~dootika/assets/course/MixG_data/181075.txt",
                  header = F, col.names = c("X", "Y"))
dat<-as.matrix(dat)

#Function to calculate multivariate normal density value
multivariate_normal<-function(x,mu,sigma_matrix){
  
  p<-ncol(x)
  n<-nrow(x)
  
  value<-rep(1,n)
  
  for(i in 1:n){
    log_value<-(-p/2)*log(2*pi)-(0.5*log(det(sigma_matrix)))-(0.5*t((x[i,]-mu))%*%solve(sigma_matrix)%*%(x[i,]-mu))
    value[i]<-exp(log_value)
  }
  value
}

#Function for estimating mean vectors, covariance matrices and mixture probabilities
#for given number of classes c
em_estimates<-function(x,c){
  
  n<-nrow(x) #Getting number of data points.
  p<-ncol(x) #Getting dimension of each data points. Here it will be 2.
  
  sample_mean<-colSums(x)/n
  sample_variance<-var(x)
  
  pie_vector<-rep(1/c,c)
  mean_vectors<-sample_mean+matrix(rnorm(p*c),nrow=p) #Creating c different mean vectors
  #Here mean_vectors is a matrix, where ith column is mean vector for ith class
  sigma_square<-list(1,2) #Creating an  empty list.
  for(i in 1:c){
    sigma_square[[i]]<-sample_variance
  }
  #ith element of sigma_square list is covariance matrix of ith class
  theta<-c(pie_vector,mean_vectors,unlist(sigma_square))
  
  d<-1
  ep<-matrix(1,nrow=n,ncol=c)
  i<-0
  while(d>10^-5){
    
    i<-i+1
    old<-theta
    for(i in 1:c){
      ep[,i]<-pie_vector[i]*multivariate_normal(x,mu=mean_vectors[,i],sigma=sigma_square[[i]])
    }
    ep<-ep/rowSums(ep) #Posterior probabilities.
    
    pie_vector<-colSums(ep)/n #Updated prior probabilities
    for(i in 1:c){
      mean_vectors[,i]<-colSums(ep[,i]*x)/sum(ep[,i]) #Updated mean vectors
    }
    for(i in 1:c){
      sigma_square[[i]]<-cov.wt(x,ep[,i])$cov #Updated covariance matrices
    }
    theta<-c(pie_vector,mean_vectors,unlist(sigma_square))
    d<-max(abs(theta-old))
  }
  
  list("mu.est"=t(mean_vectors),"sig.est"=sigma_square,"mix.est"=pie_vector,"posterior_probabilities"=ep)
  
}

#Creating function to calculate log likelihood values.
log_likelihood_value<-function(x,estimates){
  
  #Extracting values
  pie_vector<-estimates$mix.est
  mean_vectors<-t(estimates$mu.est)
  sigma_square<-estimates$sig.est
  
  n<-nrow(x)
  p<-ncol(x)
  
  classes<-length(pie_vector)
  
  s<-rep(0,n)
  for(i in 1:classes){
    s<-s+(pie_vector[i]*multivariate_normal(x,mu=mean_vectors[,i],sigma=sigma_square[[i]]))
  }
  value<-(-sum(log(s)))
  
  value
}

#Creating function to calculate AIC values.
aic_value<-function(x,estimates){
  
  #Extracting values
  pie_vector<-estimates$mix.est
  mean_vectors<-t(estimates$mu.est)
  sigma_square<-estimates$sig.est
  
  n<-nrow(x)
  p<-ncol(x)
  
  classes<-length(pie_vector) #Classes contains number of classes.
  number_of_parameters<-(classes-1)+(classes*p)+(classes*p*(p+1)/2)
  
  s<-rep(0,n)
  for(i in 1:classes){
    s<-s+(pie_vector[i]*multivariate_normal(x,mu=mean_vectors[,i],sigma=sigma_square[[i]]))
  }
  value<-( -2*sum(log(s)))+(2*number_of_parameters)
  
  value
}

#Conducting 5-Fold Cross-Validation
set.seed(5)
n<-nrow(dat)
permutation <- sample(1:n, replace = FALSE)
K <- 5

test.index <- split(permutation, rep(1:K, length = n, each = n/K))
pot_c<-2:6
cv_log_likelihood_values<-1:length(pot_c)
for(c in 1:length(pot_c))
{
  foo<-0
  for(k in 1:K)
  {
    x_train <- dat[-test.index[[k]], ]
    x_test <- dat[test.index[[k]], ]
    
    estimates<-em_estimates(x_train,pot_c[c])
    
    foo<-foo+log_likelihood_value(x_test,estimates)
  }
  cv_log_likelihood_values[c]<-foo/n
}

cv_log_likelihood_values

#Determining AIC values for different models with dfferent c values.
aic_values<-1:length(pot_c)
for(c in 1:length(pot_c)){
  estimates<-em_estimates(dat,pot_c[c])
  aic_values[c]<-aic_value(dat,estimates)
}

aic_values

#Minimum AIC is obtained for the model with C=5. But we will choose C=4.
estimates<-em_estimates(dat,4)
mle.est<-estimates[-4]
mle.est

#Assigning data to the class for which posterior probability is maximum
posterior_probabilities<-estimates$posterior_probabilities
assigned_class<-apply(posterior_probabilities,1,which.max)

#Ploting data according to their assigned class
palette(rainbow(4))
plot(dat,col=c(assigned_class),pch=19)