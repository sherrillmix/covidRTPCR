#  Analysis of SARS-CoV-2 RT-qPCR Detection Rates Prior to and After Symptom Onset

Code and data for the preprint "Variation in the Probability of SARS-CoV-2 Detection Using RT-qPCR Testing Prior to and After Symptom Onset". 

The figures in [rtpcr.pdf](rtpcr.pdf) can be regenerated by running [manuscript/rt-pcr-analysis_markov.Rmd](manuscript/rt-pcr-analysis_markov.Rmd). Stan code for the Bayesian model is in [Stan/npv-ssm.stan](Stan/npv-ssm.stan).

For processing the previously published data, the package `tidyverse` is required. MCMC fitting is accomplished using the package `rstan`. The packages can be obtained by e.g.:

```{R}
install.packages(c('tidyverse','rstan'))
```
