## single level transducer model.
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

  vector<lower=0>[S] prior_p; 
  vector<lower=0>[S] prior_q; 
  vector<lower=0, upper=1>[S] prior_z;
  vector<lower=0, upper=100>[S] prior_rmax;
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

  p ~ normal(2.0, 1); // sd is mean /2
  q ~ normal(0.4, 0.2); // sd is mean /2
  z ~ uniform(0, 1);
  #rmax ~ lognormal(log(40), 1);

  prior_p ~ normal(2.0, 1);
  prior_q ~ normal(0.4, 0.2);
  prior_z ~ uniform(0, 1);
  #prior_rmax ~ lognormal(log(40), 1);

  gamma <- 0.25; // chance performance rate.
  max_check[2] <- 0.0;

  # Trial loop begins ----------------
  for (n in 1:N) {
    subj <- ss[n];

    # calculate denominator :
    denom_ped <- pow(z[subj],p[subj]) + pow(ped[n],p[subj]);
    denom_ped_plus_inc <- pow(z[subj],p[subj]) + pow(ped_plus_inc[n],p[subj]);
  
    r_ped <- rmax[subj] * (pow(ped[n],(p[subj] + q[subj]))) / denom_ped;
  
    r_pedInc <- rmax[subj] * (pow(ped_plus_inc[n],(p[subj] + q[subj]))) / denom_ped_plus_inc;
    
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
  
