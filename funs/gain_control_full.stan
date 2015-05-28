## full gain control model. 
# Tom Wallis wrote it. tsawallis at gmail dot com

data {
  int<lower=0> N; // number of rows across all subjects
  int<lower=0> S; // number of subjects
  int<lower=0> SF; // number of spatial frequencies.
  int<lower=1,upper=S> ss[N]; // vector associating each row with a subject.
  int<lower=1,upper=SF> sf[N]; // vector associating each row with a spatial frequency.
  int<lower=0,upper=1> y[N]; 
  vector[N] ped; // pedestal contrast.
  vector[N] ped_plus_inc; // increment contrast.
  matrix[N,SF] norm_pool; // normalisation pool. 
}
  
transformed data {

}
  
parameters {
  # subject level params:
  vector<lower=0>[S] p; 
  vector<lower=0>[S] q; 
  matrix<lower=0>[SF,S] z;   
  matrix<lower=0, upper=1>[SF,S] w;  

  # population level params:
  real<lower=0> mu_p;  // for some reason, upper bound on p breaks stan...
  real<lower=0> mu_q;
  vector<lower=0>[SF] mu_z;
  vector<lower=0, upper=1>[SF] mu_w;

  real<lower=0> sigma_p;
  real<lower=0> sigma_q;
  vector<lower=0>[SF] sigma_z;
  vector<lower=0, upper=1>[SF] sigma_w;
}
  
transformed parameters {
}
  
model {
  vector[N] prob;
  real r_ped;
  real r_pedInc;
  real delta_r;
  real w_norm_pool;
  real denom_ped;
  real denom_ped_plus_inc;
  real gamma;
  real raw_pc;
  vector[SF] summands;
  vector[2] max_check;
  int subj;
  int targ_sf;  

  mu_p ~ normal(2.0, 1); 
  mu_q ~ normal(.4, 1); 
  mu_z ~ lognormal(log(.01),2);
  mu_w ~ normal(0.5,2);

  sigma_p ~ cauchy(0, 0.5);
  sigma_q ~ cauchy(0, 0.25);
  sigma_z ~ cauchy(0, 0.2); // log units if z is lognormal.
  sigma_w ~ cauchy(0, 0.25);

  # each subject mu is derived from the population mu:
  p ~ normal(mu_p, sigma_p);
  q ~ normal(mu_q, sigma_q);

  for (i in 1:SF){
    z[i] ~ lognormal(log(mu_z[i]), sigma_z[i]);
    w[i] ~ normal(mu_w[i], sigma_w[i]);
  }
  
  gamma <- 0.25; // chance performance rate.
  max_check[2] <- 0.0;

  # Trial loop begins ----------------
  for (n in 1:N) {
    subj <- ss[n];
    targ_sf <- sf[n];

    # calculate norm pool for each trial.
    for (i in 1:SF){
      summands[i] <- pow(norm_pool[n,i],p[subj]) * w[i,subj];
    }

    w_norm_pool <- sum(summands);

    # calculate denominator from norm pool and semi-saturation:
    denom_ped <- w_norm_pool + pow(z[targ_sf,subj],p[subj]);
    denom_ped_plus_inc <- w_norm_pool + pow(z[targ_sf,subj],p[subj]);
  
    r_ped <- (pow(ped[n],(p[subj] + q[subj]))) / denom_ped;
  
    r_pedInc <- (pow(ped_plus_inc[n],(p[subj] + q[subj]))) / denom_ped_plus_inc;
    
    max_check[1] <- r_pedInc - r_ped;
    delta_r <- max(max_check);

    # manual cumulative weibull with parameters fit to real dprime function:
    raw_pc <- 1.0 - exp( - pow( (delta_r / 1.545903), 1.481270) );
    prob[n] <- gamma + (1.0 - gamma) * raw_pc;
  }
  y ~ bernoulli(prob);  

  # Print statements for debugging:
  #print("p = ", p, " ","q = ", q, " ");
  #print("mu_p =", mu_p, ", sigma_p = ", sigma_p, 
  #      ", mu_q = ", mu_q, ", sigma_q = ", sigma_q);
  #print("mu_z =", mu_z, ", sigma_z = ", sigma_z, 
  #      ", mu_w = ", mu_w, ", sigma_w = ", sigma_w);
  #print("mean prob = ", mean(prob));
}
  
