library(rstan)
library(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Stan settings
n_iter <- 1e4
n_warmup <- 1e4-5e3
p_adapt_delta <- 0.99
n_max_treedepth <- 20
## the max number of days after exposure to estimate
T_max <- 21
exposed_n <- 686
exposed_pos <- 77

source("R/utils.R")

## read in raw data
raw_data <- read_csv("data/antibody-test-data.csv") %>% 
    filter(grepl("RT_PCR", test),
           study != "Danis_no_4")

pcr_dat <- raw_data %>% 
    ## add non-quantified positives to other positives for Danis et al.
    mutate(n_adj=n+nqp,
           test_pos_adj=test_pos+nqp) %>% 
    ## remove estimates without observations
    filter(n_adj > 0,
           ## days needs to be above -5
           day > -5,
           ## only use the nasal swabs from Kujawski, not throat swabs
           !(study == "Kujawski" & test == "RT_PCR_oro")) %>% 
    mutate(study_idx=paste(study, test, sep="_") %>% as.factor() %>% as.numeric(),
           pct_pos=test_pos_adj/n_adj)

## create orthogonal polynomials for days since exposure
day_poly <- poly(log(pcr_dat$day+5), degree=3)
poly_predict <- predict(day_poly, log(1:T_max))


npv_onset_model <- stan_model("Stan/npv-fixed-onset.stan")
main_analysis <- make_analysis_data(stan_model=npv_onset_model,
                                    dat=pcr_dat,
                                    T_max=T_max,
                                    poly_est=as.matrix(day_poly),
                                    poly_pred=poly_predict,
                                    exposed_n=exposed_n,
                                    exposed_pos=exposed_pos,
                                    spec=1,
                                    iter=n_iter,
                                    warmup=n_warmup,
                                    chains=20,
                                    control=list(adapt_delta=p_adapt_delta,
                                                 max_treedepth=n_max_treedepth),
                                    save_warmup=F,
                                    save_stan=T)


newDat<-read.csv('../newData/newData.csv',stringsAsFactors=FALSE)
newDat$n_adj<-newDat$FALSE.+newDat$TRUE.
newDat$test_pos_adj<-newDat$TRUE.
newDat$study_idx<-max(pcr_dat$study_idx)+7
newDat<-newDat[newDat$day<=max(pcr_dat$day),]

dat<-rbind(pcr_dat[,c('study_idx','day','n_adj','test_pos_adj')],newDat[,c('study_idx','day','n_adj','test_pos_adj')])
dat$dayAdjust<-dat$day-min(dat$day)+1

npv_onset_model3 <- stan_model("Stan/model.stan")
stan_sample <- sampling(npv_onset_model3,
                        data=list(N=nrow(dat),
                                  J=max(dat$study_idx),
                                  T_max=max(dat$dayAdjust),
                                  t=dat$dayAdjust,
                                  test_n=dat$n_adj,
                                  test_pos=dat$test_pos_adj,
                                  study_idx=dat$study_idx,
                                  zeroDay=unique(dat[dat$day==0,'dayAdjust',drop=TRUE])
                                  ),
                            iter=n_iter,
                            warmup=n_warmup,
                            control=list(adapt_delta=.98,
                                          max_treedepth=15),
                            chains=40
                        )

zz<-as.matrix(main_analysis$stan_sample)
invLogit<-function(xx)1/(1+exp(-xx))
days<-min(pcr_dat$day):max(pcr_dat$day)#sort(unique(pcr_dat$day))
dayP<-predict(day_poly,log(days+5)) #do.call(rbind,lapply(days,function(xx)day_poly[which(pcr_dat$day==xx)[1],]))
mat<-invLogit(zz[,'beta_0']+zz[,'beta_j[5]']+zz[,c('beta_1','beta_2','beta_3')] %*% t(dayP[,1:3]))
predsM<-apply(mat,2,mean)
predsU<-apply(mat,2,quantile,.975)
predsL<-apply(mat,2,quantile,.025)
preds<-invLogit(mean(zz[,'beta_0'])+dayP[,1]*mean(zz[,'beta_1'])+dayP[,2]*mean(zz[,'beta_2'])+dayP[,3]*mean(zz[,'beta_3']))
tests<-data.frame('pos'=tapply(pcr_dat$test_pos,pcr_dat$day,sum),'n'=tapply(pcr_dat$n,pcr_dat$day,sum))
tests$day<-as.numeric(rownames(tests))
tests2<-data.frame('pos'=tapply(dat$test_pos_adj,dat$day,sum),'n'=tapply(dat$n_adj,dat$day,sum))
tests2$day<-as.numeric(rownames(tests2))
posStudy<-tapply(pcr_dat$test_pos,list(pcr_dat$day,pcr_dat$study),sum)
nStudy<-tapply(pcr_dat$n,list(pcr_dat$day,pcr_dat$study),sum)
nStudy[is.na(nStudy)]<-0
studyDays<-as.numeric(rownames(nStudy))
nStudy<-cbind('all'=tests$n,nStudy)
posStudy<-cbind('all'=tests$pos,posStudy)
#zz2<-as.matrix(main_analysis2$stan_sample)
#preds2<-invLogit(mean(zz2[,'beta_0'])+dayP[,1]*mean(zz2[,'beta_1'])+dayP[,2]*mean(zz2[,'beta_2'])+dayP[,3]*mean(zz2[,'beta_3']))
zz3<-as.matrix(stan_sample)
predDays<-min(dat$day):max(dat$day)
preds3<-invLogit(apply(zz3[,grep('dayRate',colnames(zz3))]+zz3[,'beta_j[5]'],2,mean))
preds3Upper<-invLogit(apply(zz3[,grep('dayRate',colnames(zz3))]+zz3[,'beta_j[5]'],2,quantile,.975))
preds3Lower<-invLogit(apply(zz3[,grep('dayRate',colnames(zz3))]+zz3[,'beta_j[5]'],2,quantile,.025))
tests2$origPos<-tests[as.character(tests2$day),'pos']
tests2$origN<-tests[as.character(tests2$day),'n']
tests2[is.na(tests2$origPos),'origPos']<-0
tests2[is.na(tests2$origN),'origN']<-0
tests2$origNeg<-tests2$origN-tests2$origPos
#
origCols<-c('#ffe9e9','#b89999') #,'#b30000'
newCols<-c('#e9e9ff','#9999b8')#,'#045a8d'
#prevent pdf bugginess at exact rect border meetings
pdfFill<-.0005
#origCols<-c('#fef0d9','#fc8d59') #,'#b30000'
#newCols<-c('#f1eef6','#74a9cf')#,'#045a8d'
pdf('modelsCompare.pdf',width=4.5,height=4.5)
par(mar=c(3.2,3.2,.4,3.5))
plot(predDays,preds3,xlab='Day after onset of symptoms',ylab='Probability of detection',ylim=c(0,1),type='n',yaxs='i',las=1,mgp=c(2.25,.8,0),xlim=c(-9,30))
#lines(predDays,preds3Upper,xpd=NA)
#lines(predDays,preds3Lower,xpd=NA)
lines(days,predsM,col='red')
polygon(c(days,rev(days)),c(predsU,rev(predsL)),col='#FF000033',border=NA)
lines(predDays,preds3,col='blue')
polygon(c(predDays,rev(predDays)),c(preds3Upper,rev(preds3Lower)),col='#0000FF33',border=NA)
abline(v=-.5,lty=2)
rectHeight<-diff(par('usr')[3:4])*.35
rectMin<-par('usr')[3]
perTest<-rectHeight/max(tests2$n)
rect(tests2$day-.4,rectMin+tests2$origPos*perTest+pdfFill,tests2$day+.4,rectMin,col=origCols[2],border=NA)
rect(tests2$day-.4,rectMin+(tests2$origNeg+tests2$pos)*perTest+pdfFill,tests2$day+.4,rectMin+tests2$pos*perTest,col=origCols[1],border=NA)
rect(tests2$day-.4,rectMin+tests2$pos*perTest+pdfFill,tests2$day+.4,rectMin+tests2$origPos*perTest,col=newCols[2],border=NA)
rect(tests2$day-.4,rectMin+tests2$n*perTest+pdfFill,tests2$day+.4,rectMin+(tests2$origNeg+tests2$pos)*perTest,col=newCols[1],border=NA)
rect(tests2$day-.4,rectMin,tests2$day+.4,rectMin+tests2$n*perTest,col=NA,border='#00000088')
segments(tests2$day-.4,rectMin+tests2$pos*perTest,tests2$day+.4,rectMin+tests2$pos*perTest,col='#00000088',lwd=.5)
axPos<-pretty(1:max(tests$n))
axis(4,rectMin+axPos*perTest,as.character(axPos),las=1,mgp=c(2.5,.8,0))
text(dnar::convertLineToUser(3,4),50*perTest+rectMin,'Number of tests',xpd=NA,srt=270,mgp=c(2.25,.8,0))
#legend('topright',c('This study','Previous estimate'),lwd=1,col=c('red','blue'),bty='n',y.intersp=.75)
#legend('topright',c('This study','Previous'),lwd=1,fill=c('#FF000033','#0000FF33'),bty='n',y.intersp=.75)
dev.off()





