## single level GLM.
# Tom Wallis wrote it. tsawallis at gmail dot com

data {
  int<lower=0> N; // number of rows across all subjects
  int<lower=0> S; // number of subjects
  int<lower=1,upper=S> ss[N]; // vector associating each row with a subject.
  int<lower=0,upper=1> y[N]; 
  int<lower=0> D; // number of coefficients to estimate.
  vector[D] X[N]; // design matrix.
}
  
transformed data {

}
  
parameters {
  # subject level params:
  vector<lower=-5, upper=5>[D] beta[S]; # beta for each subject.

  # collect samples of the prior:
  vector<lower=-5, upper=5>[D] prior_beta[S];
}
  
transformed parameters {
}
  
model {
  vector[N] prob;
  real gamma;

  for (s in 1:S) {
    beta[s] ~ normal(0,2);
    prior_beta[s] ~ normal(0,2);
  }

  gamma <- 0.25;

  for (n in 1:N) {
    prob[n] <- gamma + (1.0 - gamma) * inv_logit( dot_product(X[n],beta[ss[n]]) );
  }
  y ~ bernoulli(prob);  

  # Print statements for debugging:
  #print("beta = ", beta);
}
  
