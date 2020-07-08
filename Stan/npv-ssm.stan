//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic.
// We assume that the sensitivity of the diagnostic decays with time.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// We have RT-PCR data from two studies where patients eventually developed
// antibodies (IgG or IgM) to covid-19.
// 'N' is the number of rows in our input data, which has a row for each number of
// days since symptom onset. 'T_max' is the number of days that we would like to
// predict into the future. 'test_n' is the number of people who were tested
// and 'test_pos' is the number of people who tested positive. 't_symp_test' is
// the number of days since symptom onset. 'exposed_n' and 'exposed_pos' are
// numbers to estimate the attack rate of covid-19 from other papers.
data {
    int<lower=1> N;
    int<lower=1> J; // number of studies
    int<lower=1> T_max;
    int<lower=1> test_n[N];
    int<lower=0> test_pos[N];
    int<lower=1> study_idx[N];
    int<lower=1,upper=T_max> t[N];
    int zeroDay;
    real zeroMean;
    real zeroSd;
}

// the beta terms are the coefficients for the cubic polynomial for log-time.
// 'attack_rate' is the probability of infection given exposure.
parameters{
    real<lower=0> sigma;
    real<lower=0> daySigma;
    real zeroDayRate;
    real rawDayDiff[T_max-1];
    vector[J] beta_j;
}

transformed parameters{
    vector[T_max] dayRate;
    vector[N] mu;
    vector[T_max-1] dayDiff;
    dayRate[1]=zeroDayRate;
    dayDiff[1]=rawDayDiff[1]*daySigma;
    for(ii in 2:(T_max-1)) dayDiff[ii]=dayDiff[ii-1]+rawDayDiff[ii]*daySigma;
    for(ii in 2:T_max) dayRate[ii]=dayRate[ii-1]+dayDiff[ii-1];
    mu = beta_j[study_idx]+dayRate[t];
}

model {
    beta_j ~ normal(0,sigma);
    rawDayDiff~normal(0,1);
    zeroDayRate~normal(zeroMean,zeroSd);
    sigma~gamma(1,1);
    daySigma~gamma(1,1);
    test_pos ~ binomial_logit(test_n, mu);
}

