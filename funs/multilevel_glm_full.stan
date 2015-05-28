## multi level GLM.
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
  vector[D] e_beta[S]; # beta offset for each subject.

  # population level params:
  vector<lower=-5, upper=5>[D] mu_beta; 
  vector<lower=0,upper=5>[D] sigma_beta;

  # collect samples of the prior:
  vector[D] prior_e_beta[S];
  vector<lower=-5, upper=5>[D] prior_mu_beta; 
  vector<lower=0,upper=5>[D] prior_sigma_beta;
}
  
transformed parameters {
  vector[D] beta[S]; # beta is a beta estimate for each subject. 
  vector[D] prior_beta[S];

  # convert e_beta into betas:
  for (s in 1:S) {
     beta[s] <- mu_beta + e_beta[s] .* sigma_beta; 
     prior_beta[s] <- prior_mu_beta + prior_e_beta[s] .* prior_sigma_beta; 
  }
}
  
model {
  vector[N] prob;
  real gamma;

  mu_beta ~ normal(0,1);
  sigma_beta ~ cauchy(0,1);

  prior_mu_beta ~ normal(0,1);
  prior_sigma_beta ~ cauchy(0,1);

  for (s in 1:S) {
    e_beta[s] ~ normal(0,1);
    prior_e_beta[s] ~ normal(0,1);
  }

  gamma <- 0.25;

  for (n in 1:N) {
    prob[n] <- gamma + (1.0 - gamma) * inv_logit( dot_product(X[n],beta[ss[n]]) );
  }
  y ~ bernoulli(prob);  

  # Print statements for debugging:
  #print("beta = ", beta);
}
  
