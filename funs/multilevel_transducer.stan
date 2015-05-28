## multilevel transducer model.
# Tom Wallis wrote it. tsawallis at gmail dot com

data {
  int<lower=0> N; // number of rows across all subjects
  int<lower=0> S; // number of subjects
  int<lower=1,upper=S> ss[N]; // vector associating each row with a subject.
  int<lower=0,upper=1> y[N]; 
  vector[N] ped; // pedestal contrast.
  vector[N] ped_plus_inc; // increment contrast.
}
  
transformed data {

}
  
parameters {
  # subject level params:
  vector<lower=0>[S] p; 
  vector<lower=0>[S] q; 
  vector<lower=0, upper=1>[S] z;
  vector<lower=0, upper=100>[S] rmax;

  # population level params:
  real<lower=0> mu_p;
  real<lower=0> mu_q;
  real<lower=0, upper=1> mu_z;
  real<lower=0, upper=100> mu_rmax;

  real<lower=0> sigma_p;
  real<lower=0> sigma_q;
  real<lower=0> sigma_z;
  real<lower=0> sigma_rmax;

  # store prior samples (don't update with data)
  vector<lower=0>[S] prior_p; 
  vector<lower=0>[S] prior_q; 
  vector<lower=0, upper=1>[S] prior_z;
  vector<lower=0, upper=100>[S] prior_rmax;

  real<lower=0> prior_mu_p;
  real<lower=0> prior_mu_q;
  real<lower=0, upper=1> prior_mu_z;
  real<lower=0, upper=100> prior_mu_rmax;

  real<lower=0> prior_sigma_p;
  real<lower=0> prior_sigma_q;
  real<lower=0> prior_sigma_z;
  real<lower=0> prior_sigma_rmax;

}
  
transformed parameters {
}
  
model {
  vector[N] prob;
  real r_ped;
  real r_pedInc;
  real delta_r;
  real denom_ped;
  real denom_ped_plus_inc;
  real gamma;
  real raw_pc;
  vector[2] max_check;
  int subj;

  # population-level priors:
  mu_p ~ normal(2.0, 1);
  mu_q ~ normal(0.4, 1);
  mu_z ~ uniform(0, 1);
  mu_rmax ~ lognormal(log(40), 0.5);

  sigma_p ~ cauchy(0, 1);
  sigma_q ~ cauchy(0, 1);
  sigma_z ~ cauchy(0, 1); 
  sigma_rmax ~ cauchy(0, 1); // log units if rmax is lognormal. 

  p ~ normal(mu_p, sigma_p);
  q ~ normal(mu_q, sigma_q);
  z ~ normal(mu_z, sigma_z);
  rmax ~ lognormal(log(mu_rmax), sigma_rmax);

  # collect samples from the priors:
  prior_mu_p ~ normal(2.0, 1);
  prior_mu_q ~ normal(0.4, 1);
  prior_mu_z ~ uniform(0, 1);
  prior_mu_rmax ~ lognormal(log(40), 0.5);

  prior_sigma_p ~ cauchy(0, 1);
  prior_sigma_q ~ cauchy(0, 1);
  prior_sigma_z ~ cauchy(0, 1); 
  prior_sigma_rmax ~ cauchy(0, 1); // log units if rmax is lognormal. 

  prior_p ~ normal(mu_p, sigma_p);
  prior_q ~ normal(mu_q, sigma_q);
  prior_z ~ normal(mu_z, sigma_z);
  prior_rmax ~ lognormal(log(mu_rmax), sigma_rmax);

  gamma <- 0.25; // chance performance rate.
  max_check[2] <- 0.0;

  # Trial loop begins ----------------
  for (n in 1:N) {
    subj <- ss[n];

    # calculate denominator :
    denom_ped <- pow(z[subj],p[subj]) + pow(ped[n],p[subj]);
    denom_ped_plus_inc <- pow(z[subj],p[subj]) + pow(ped_plus_inc[n],p[subj]);
  
    r_ped <- (pow(ped[n],(p[subj] + q[subj]))) / denom_ped;
  
    r_pedInc <- (pow(ped_plus_inc[n],(p[subj] + q[subj]))) / denom_ped_plus_inc;
    
    max_check[1] <- r_pedInc - r_ped;
    delta_r <- max(max_check);

    # manual cumulative weibull with parameters fit to dprime function:
    raw_pc <- 1.0 - exp( - pow( (delta_r / 1.545903), 1.481270) );
    prob[n] <- gamma + (1.0 - gamma) * raw_pc;
  }
  y ~ bernoulli(prob);  

  # Print statements for debugging:
  #print("p = ", p, " ","q = ", q, " ");
  #print("mu_p =", mu_p, ", sigma_p = ", sigma_p, 
  #      ", mu_q = ", mu_q, ", sigma_q = ", sigma_q);
  #print("mu_z =", mu_z, ", sigma_z = ", sigma_z);
  #print("mean prob = ", mean(prob));
}
  
