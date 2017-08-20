

# MCMC fit for a logistic regression model
ls()
rm(list=ls())
gc()
cols=4;# one intercept + 3 attributes
rows=10000;
sigma2=0.05;

datax <- read.csv(file="R-simulated-logist-data-x.csv",head=FALSE,sep=",")
x1=datax$V1
x2=datax$V2
x3=datax$V3
beta_curr=seq(0.5,cols)/1.0
beta_new=seq(0.5,cols)/1.0
datay <- read.csv(file="R-simulated-logist-data-y.csv",head=FALSE,sep=",")
y=datay$V1
d1 <- data.frame(1.5, t(beta_curr))                
write.table(d1, "R-mcmc_logistic_version_2.txt", row.names = FALSE)
m=5000 # the length of the MCMC samples
for (i in 1:m)
{
  if (i%%500 ==0){
    print(i)
  }
  for (k in 1:cols)
  {
    temp1=beta_curr[k] +sigma2*rnorm(1, mean = 0, sd = 1)
    beta_new[k]=temp1;
    likelihood=0.0;
    for (h in 1:rows)
    {
      temp_new=beta_new[1] +beta_new[2]*x1[h] +beta_new[3]*x2[h]+beta_new[4]*x3[h];
      temp_new=exp(temp_new);
      temp_curr=beta_curr[1] +beta_curr[2]*x1[h] +beta_curr[3]*x2[h]+beta_curr[4]*x3[h];
      temp_curr=exp(temp_curr);
      likelihood=likelihood+ log( temp_new /(1.0 +temp_new))*y[h];
      likelihood=likelihood+ log(1.0- temp_new /(1.0 +temp_new))*(1.0-y[h]);
      likelihood=likelihood- log(temp_curr /(1.0 +temp_curr))*y[h];
      likelihood=likelihood- log(1.0- temp_curr /(1.0 +temp_curr))*(1.0-y[h]);
    }
    u=runif(1, min = 0, max = 1);
    if ( u<min(1.0, exp(likelihood)))
    {
      beta_curr[k]=beta_new[k];
    }
  }
  d1 <- data.frame(likelihood, t(beta_curr))
  write.table(d1, "R-mcmc_logistic_version_2.txt", row.names = FALSE, col.names = FALSE, append = TRUE)
}






 
 