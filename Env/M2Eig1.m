clear all
cd '/Volumes/TOSHIBA/Env'


addpath(genpath('./release2.17'))
ntrain = 500;
ntest = 500;
nsim = 100;
model =2; dimx = 4:1:10; dimy = 1:1:6;evals = sort((0.1:0.15:9),'descend')/2; 


[FLR ,ENV, Oracle]=simfun(model,dimx,dimy,evals,ntrain,ntest,nsim);



FLRmean = mean(FLR(2:nsim+1,:), 1);
csvwrite('M2Eig1FLRmean.csv',reshape(FLRmean,[6,7]))
FLRstd = std(FLR(2:nsim+1,:), 1)./sqrt(nsim);
csvwrite('M2Eig1FLRstd.csv',reshape(FLRstd,[6,7]))



ENVmean = mean(ENV(2:nsim+1,:), 1);
csvwrite('M2Eig1ENVmean.csv',reshape(ENVmean,[6,7]))
ENVstd = std(ENV(2:nsim+1,:), 1)./sqrt(nsim);
csvwrite('M2Eig1ENVstd.csv',reshape(ENVstd,[6,7]))


Oracleout= vertcat(mean(Oracle(2:nsim+1)), std(Oracle(2:nsim+1))./sqrt(nsim));
csvwrite('M2Eig1Oracle.csv',Oracleout)


