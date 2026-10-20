#######################################
#######################################
## Updated Longitudinal Code for:    ##
## Thompson et al. (2018)            ##           
## Wes Thompson (wes.stat@gmail.com) ##
## Apr 19, 2023                      ##
#######################################
#######################################

# Note this is to run the Stan code only - subsequent analyses (and many of the data checks) are left out! #

####################
####################
## Load libraries ##
####################
####################

library(mvtnorm) # version 1.1.3
library(parallel) # version 4.1.0
library(rstan) # version 2.21.8
library(loo) # version 2.6.0
library(gamm4) # version 0.2.6
library(nFactors) # version 2.4.1.1
library(reshape2) # version 1.4.4
library(psych) # version 2.3.3
library(data.table) # version 1.14.2
library(abind) # version 1.4.5
library(cvTools) # version 0.3.2

####################################
####################################
## Load and manipulate abcd data  ##
####################################
####################################

## Read Rds file from DEAP (ABCD NDA version 4.0 release)

data_dir <- '~/'
setwd(data_dir)
# Load in the file containing 5.0 data, already cleaned and subsetted 
abcd <- readRDS('Release5-0_SubsetData.rds')

#########################################
#########################################
## Select & process data for analyses  ##
#########################################
#########################################

## Select demographics
ind_demog = c(which(names(abcd)=="age"),which(names(abcd)=="female"),which(names(abcd)=="race.6level"),
              which(names(abcd)=="high.educ.bl"),which(names(abcd)=="married.bl"),which(names(abcd)=="household.income.bl"))
#names(abcd)[ind_demog]
#summary(abcd[,ind_demog])

## Select nesting variables
ind_nest = c(which(names(abcd)=="abcd_site"),which(names(abcd)=="rel_family_id"));summary(abcd[,ind_nest])

abcd$abcd_site = as.character(abcd$abcd_site)
abcd$abcd_site[abcd$abcd_site=="site22"] = "site07"
abcd$abcd_site = factor(abcd$abcd_site)

## Select neuropsychological measures
#par(mfrow=c(1,2))
#hist(abcd$pea_ravlt_ld)
#hist(abcd$lmt_scr_perc_correct)
abcd$lmt_scr_perc_correct[abcd$lmt_scr_perc_correct>1 & !is.na(abcd$lmt_scr_perc_correct)] = abcd$lmt_scr_perc_correct[abcd$lmt_scr_perc_correct>1 & !is.na(abcd$lmt_scr_perc_correct)]/100
#hist(abcd$lmt_scr_perc_correct)

# Specify variables to select
ind_np = c(which(names(abcd)=="nihtbx_picvocab_uncorrected"),which(names(abcd)=="nihtbx_flanker_uncorrected"),
           which(names(abcd)=="nihtbx_list_uncorrected"),which(names(abcd)=="nihtbx_cardsort_uncorrected"),
           which(names(abcd)=="nihtbx_pattern_uncorrected"),which(names(abcd)=="nihtbx_picture_uncorrected"),
           which(names(abcd)=="nihtbx_reading_uncorrected"),which(names(abcd)=="pea_ravlt_ld"),
           which(names(abcd)=="lmt_scr_perc_correct"))

#summary(abcd[abcd$eventname=="baseline_year_1_arm_1",ind_np])
#summary(abcd[abcd$eventname=="2_year_follow_up_y_arm_1",ind_np])

# NOW FURTHER REMOVE LMT (9) in addition to RAVLT (8) in addition to "list" and "cardsort" (3 and 4)
ind_np = ind_np[-c(3,4,8,9)]
#summary(abcd[abcd$eventname=="baseline_year_1_arm_1",ind_np])
#summary(abcd[abcd$eventname=="2_year_follow_up_y_arm_1",ind_np])

###################################
###################################
## Subset variables for analyses ##
###################################
###################################

# Select visits to use (baseline, year 2, year 4)
visits.to.use = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1","4_year_follow_up_y_arm_1")

# Subset data
data = abcd[abcd$eventname %in% visits.to.use ,
            c(ind_nest,ind_np, 
              which(names(abcd)=="src_subject_id"), 
              which(names(abcd)=="eventname"))]

#names(data);dim(data)

# Take complete cases
data = data[complete.cases(data),]
#table(data$event_name)

# Define `subjectid` as src_subject_id var
data$subjectid <- data$src_subject_id
data$n_visits = 0
for(i in unique(data$subjectid)){
  data$n_visits[data$subjectid==i] = length(data$n_visits[data$subjectid==i])
}
# Now let's just use those that have all 3 visits 
data = data[data$n_visits==3,]
#data = data[data$n_visits==2,]

data$site_num = as.numeric(substr(data$abcd_site,5,6))
data$fam_num = 0
ind=0
for(i in unique(data$rel_family_id)){
  ind = ind+1
  data$fam_num[data$rel_family_id==i & !is.na(data$rel_family_id)] = ind
}
data = data[order(data$site_num,data$fam_num,data$rel_family_id),]
data$id_num = 0
ind=0
for(i in sort(unique(data$subjectid))){
  ind=ind+1
  data$id_num[data$subjectid==i] = ind
}

# Check number of families and subjects
length(unique(data$fam_num)) 
length(unique(data$subjectid)) 

####################
####################
## Fit usual PCA  ## 
####################
####################

## PCA on unimputed data 
ind_Y = c(3:7); names(data)[ind_Y] 
Y = as.matrix(scale(data[complete.cases(data[,c(ind_Y)]),ind_Y]))

#dim(Y[data$event_name==visits.to.use[1],])
#dim(Y[data$event_name==visits.to.use[2],])

ev = eigen(cor(Y))
ap = parallel(subject=nrow(Y),var=ncol(Y),rep=100,cent=.05)
nS = nScree(x=ev$values,aparallel=ap$eigen$qevpea)
#plotnScree(nS)
ncomp = 3
#y.pca = psych::principal(Y, rotate="promax", nfactors=ncomp, scores=TRUE)
#y.pca$loadings
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
#y.pca$loadings

#########################################
#########################################
## Fit 5-variable PCA Model using Stan ##
#########################################
#########################################

set.seed(314)
#keep = unique(sample(data$fam_num,120,replace=FALSE))
#data1 = data[data$fam_num %in% keep,]
#dim(data1)
data1=data

data1$site_num = 0
ind=0
for(i in unique(data1$abcd_site)){
  ind = ind+1
  data1$site_num[data1$abcd_site==i] = ind
}
data1$fam_num = 0
ind=0
for(i in unique(data1$rel_family_id)){
  ind = ind+1
  data1$fam_num[data1$rel_family_id==i & !is.na(data1$rel_family_id)] = ind
}
data1 = data1[order(data1$site_num,data1$fam_num,data1$rel_family_id),]
data1$id_num = 0
ind=0
for(i in sort(unique(data1$subjectid))){
  ind=ind+1
  data1$id_num[data1$subjectid==i] = ind
}

Part = data1$id_num
Site = data1$site_num
Fam = data1$fam_num
np = length(unique(Part))
nf = length(unique(Fam))
ns = length(unique(Site))
Y = (as.matrix(scale(data1[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################
## Run stan model ##
####################

Nsamples = 1000
Nchains = 3
D_max = 5

model_file = "bppcaLong2.stan"
smod = stan_model(model_file)

# Set cores to run in parallel
options(mc.cores = parallel::detectCores())
# Print out number of cores just to confirm
print(paste("number of cores:", parallel::detectCores()))

D_max = 5
sa.list = list()
log_lik.list = list()
looic.list = list()
for(d in 1:D_max){
  print(d)
  pca_data <- list(Y = Y, N = N, P = P, D = d, Fam = Fam, Site = Site, Part = Part, ns = ns, nf = nf, np = np)
  #set.seed(6)
  sa.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains, init="random",
                          cores = Nchains) # init="random"
  #log_lik.list[[d]] <- extract_log_lik(sa.list[[d]])
  log_lik.list[[d]] <- extract(sa.list[[d]],"log_lik_marg")[[1]]
  looic.list[[d]] = loo(log_lik.list[[d]])
  save(sa.list,log_lik.list,looic.list,file="bppca2_results_Year4.RData")
  print("###############################")
}

#################################
## Model selection using LOOIC ##
#################################

load("bppca2_results_Year4.RData")
# Check LOOIC values
show(looic.list)


d=3
# Check convergence
# print(sa.list[[d]], pars=c("Q"), probs=c(.025,.5,.975))
# print(sa.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma2_f","sigma_c","sigma_d","sigma_g"), probs=c(.025,.5,.975))

#######################################
## Extract parameters from the model ##
#######################################

sa = sa.list[[d]]
Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma2_f = extract(sa,"sigma2_f",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_g = extract(sa,"sigma_g",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

par(mfrow=c(2,3))
hist(sigma2_a,main="Var(a)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_a,.5),3)))
hist(sigma2_b,main="Var(b)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_b,.5),3)))
hist(sigma2_f,main="Var(f)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_f,.5),3)))
hist(sigma_c^2,main="Var(c)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_c^2,.5),3)))
hist(sigma_d^2,main="Var(d)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_d^2,.5),3)))
hist(sigma_g^2,main="Var(g)",xlim=c(0,1),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_g^2,.5),3)))

### Save variance components and get estimated random effect variance
vc_tab = array(0, dim=c(6,3))
vc_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))
vc_tab[5,] = quantile(sigma2_f,c(0.25,.5,.975))
vc_tab[6,] = quantile(sigma_g^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
R_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
S_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
W_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
Theta_old = array(0, dim=c(d,N,Nchains*Nsamples/2))
Theta_new = array(0, dim=c(d,N,Nchains*Nsamples/2))
Lambda_old = array(0, dim=c(P,d,Nchains*Nsamples/2))
Lambda_new = array(0, dim=c(P,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind+ 1
    Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
    Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
    Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
    W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
    R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
    S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
    Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
  }
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
  for(j in 1:(Nsamples*Nchains/2)){
    if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
      ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
    ll_mx[k] = max(abs(Lambda_new[,k,j]))
  }
}
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    for(k in 1:d){
      Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
    }
    Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
      Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
  }
}

## Varimax & Promax Rotations
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
Lambda_pmx = Lambda_new
Theta_pmx = Theta_new
Rot_pmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    tmp = varimax(Lambda_new[,,ind])
    Rot_vmx[,,ind] = tmp$rot
    Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
    Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
    tmp = psych::Promax(Lambda_new[,,ind])
    Rot_pmx[,,ind] = tmp$rot
    Lambda_pmx[,,ind] = Lambda_new[,,ind]%*%Rot_pmx[,,ind]
    Theta_pmx[,,ind] = t(Rot_pmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
  }
}
## Permute Varimax factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
  tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    ll=rep(0,d)
    for(d1 in 1:d){
      tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
      tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
      if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
      if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
    }
    Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
    Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
    Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
    for(d1 in 1:d){
      Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
      Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
      Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
    }
  }
}

## Permute Promax factors
tr=principal(Y,nfactors=d,rotate="promax")$loadings[,1:d]
for(d1 in 1:d){
  tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    ll=rep(0,d)
    for(d1 in 1:d){
      tmp1=apply((Lambda_pmx[,d1,ind] - tr[,])^2,2,sum)
      tmp2=apply((Lambda_pmx[,d1,ind] + tr[,])^2,2,sum)
      if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
      if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
    }
    Lambda_pmx[,,ind] = Lambda_pmx[,ll,ind]
    Theta_pmx[,,ind] = Theta_pmx[ll,,ind]
    Rot_pmx[,,ind] = Rot_pmx[ll,ll,ind]
    for(d1 in 1:d){
      Theta_pmx[d1,,ind] = Theta_pmx[d1,,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
      Rot_pmx[,d1,ind] = Rot_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
      Lambda_pmx[,d1,ind] = Lambda_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
    }
  }
}
## Variance explained by retained factors
var_new = array(0, dim=c(d,Nchains*Nsamples/2))
var_vmx = array(0, dim=c(d,Nchains*Nsamples/2))
var_pmx = array(0, dim=c(d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
    var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
    var_pmx[,ind] = diag(t(Lambda_pmx[,,ind])%*%(Lambda_pmx[,,ind]))
  }
}


## Communalities & Uniquenesses
comm_uniq_new = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_vmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_pmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
  for(j in 1:dim(Q)[2]){
    ind = ind + 1
    comm_uniq_new[,,ind] = cbind(diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])),1-diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])))
    comm_uniq_vmx[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
    comm_uniq_pmx[,,ind] = cbind(diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])),1-diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])))
  }
}
## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
  for(p2 in 1:P){
    Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
  }
}
comm_uniq_new_m = array(0,dim=c(P,2,3))
comm_uniq_vmx_m = array(0,dim=c(P,2,3))
comm_uniq_pmx_m = array(0,dim=c(P,2,3))
Lambda_new_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
Lambda_pmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
  comm_uniq_new_m[p,1,] = quantile(comm_uniq_new[p,1,],probs=c(.025,.5,.975))
  comm_uniq_vmx_m[p,1,] = quantile(comm_uniq_vmx[p,1,],probs=c(.025,.5,.975))
  comm_uniq_pmx_m[p,1,] = quantile(comm_uniq_pmx[p,1,],probs=c(.025,.5,.975))
  comm_uniq_new_m[p,2,] = quantile(comm_uniq_new[p,2,],probs=c(.025,.5,.975))
  comm_uniq_vmx_m[p,2,] = quantile(comm_uniq_vmx[p,2,],probs=c(.025,.5,.975))
  comm_uniq_pmx_m[p,2,] = quantile(comm_uniq_pmx[p,2,],probs=c(.025,.5,.975))
  for(k in 1:d){
    Lambda_new_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
    Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
    Lambda_pmx_m[p,k,] = quantile(Lambda_pmx[p,k,],probs=c(.025,.5,.975))
  }
}
Theta_new_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
Theta_pmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
  for(i in 1:N){
    Theta_new_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
    Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
    Theta_pmx_m[k,i,] = quantile(Theta_pmx[k,i,],probs=c(.025,.5,.975))
  }
}		
var_new_m = array(0,dim=c(d,3))
var_vmx_m = array(0,dim=c(d,3))
var_pmx_m = array(0,dim=c(d,3))
for(k in 1:d){
  var_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
  var_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
  var_pmx_m[k,] = quantile(var_pmx[k,],c(.025,.5,.975))
}
var_expl_new = apply(var_new_m,2,sum)/P
var_expl_new 
var_expl_vmx = apply(var_vmx_m,2,sum)/P
var_expl_vmx 
var_expl_pmx = apply(var_pmx_m,2,sum)/P
var_expl_pmx 

## Add PC scores into dataset
data1$pc1_new = Theta_vmx_m[1,,2]
data1$pc2_new = Theta_vmx_m[2,,2]
data1$pc3_new = Theta_vmx_m[3,,2]

loadings_new = Lambda_vmx_m[,,2]
# This should be helpful to figure out if the factor structure is the same
loadings_new 


# Now save all the objects listed above that are needed for subsequent analyses
save(Y, ind_Y, p, d, Q_m, comm_uniq_new_m, comm_uniq_vmx_m, comm_uniq_pmx_m, Lambda_new_m, Lambda_vmx_m, Lambda_pmx_m, 
     Theta_new_m, Theta_vmx_m, Theta_pmx_m, var_new_m, var_vmx_m, var_pmx_m, 
     vc_tab,
     file = 'Saved_Objects_StanPCA_Year4.Rdata')
