#Reading Data
dat <- read.table("http://home.iitk.ac.in/~dootika/assets/course/Log_data/181075.txt",
                  header = F)
dat<-as.matrix(dat)

#Extracting data matrix and response values
x<-dat[,-1]
y<-dat[,1]

#Function to calculate log of posterior distribution value
log_posterior<-function(x,y,beta_vector){
  
  sum(beta_vector^2)/200 + sum(y*(x%*%beta_vector)) - sum(log(1+exp(x%*%beta_vector)))
  
}

#Function to get samples
#Arguments of our function:
#n: number of MCMC samples we want to generate, x: Data matrix
#y: response values, starting_point: A vector of starting values
#Here we are using N(mu,sigma) as proposal distribution. Here mu is current sample
#and sigma is a diagonal matrix, which has diagonal elements as h's.
#Length of h should be same as number of columns of data matrix.

mcmc_logistic_regression_sample<-function(n,x,y,starting_point,h){
  p<-ncol(x)
  
  sampled_beta<-matrix(1,nrow=n,ncol=p) #Creating a blank matrix to store samples
  sampled_beta[1,]<-starting_point
  
  acceptance<-0 #acceptance counts the number of times the proposed sample is accepted
  
  for(i in 2:n){
    
    sampled_beta[i,]<-sampled_beta[i-1,]#If acceptance criteria is not met we will not move from current sample
    #Thus we can remove else condition
    
    proposed_sample<-sampled_beta[i-1,]+rnorm(p,sd=h)
    posterior_ratio<-exp(log_posterior(x,y,proposed_sample)-log_posterior(x,y,sampled_beta[i-1,]))
    alpha<-min(1,posterior_ratio)
    
    u<-runif(1)
    if(u<alpha){
      sampled_beta[i,]<-proposed_sample
      acceptance<-acceptance+1
    }
    
  }
  acceptance_probability<-acceptance/n
  posterior_mean<-colSums(sampled_beta)/n
  list("acc.prob"=acceptance_probability,"chain"=sampled_beta,"beta.est"=posterior_mean)
}

set.seed(5)
#Using logistic regression estimates to determine coefficients and step sizes.
dat <- read.table("http://home.iitk.ac.in/~dootika/assets/course/Log_data/181075.txt",
                  header = F)
required_data<-dat
colnames(required_data)<-c("y","x1","x2","x3","x4","x5")
logistic_regression<-glm(y~.,data=required_data[,-2],family=binomial)
estimates<-logistic_regression$coefficients
estimates #initial choice of beta.

summary(logistic_regression)
standard_error<-c(0.74,0.743,0.38,0.2,0.55) #Tuned h values
required_sample_and_outcome<-mcmc_logistic_regression_sample(100000,x,y,estimates,standard_error)

chain<-required_sample_and_outcome$chain
first_component<-chain[,1]
second_component<-chain[,2]
third_component<-chain[,3]
fourth_component<-chain[,4]
fifth_component<-chain[,5]

par(mfrow=c(2,3))
acf(first_component)
acf(second_component)
acf(third_component)
acf(fourth_component)
acf(fifth_component)

par(mfrow=c(1,2))
plot(density(first_component),main="Density Plot for First Component") 
plot(density(fifth_component),main="Density Plot for Fifth Component") 

acc.prob<-required_sample_and_outcome$acc.prob
beta.est<-required_sample_and_outcome$beta.est

print(acc.prob)
print(beta.est)