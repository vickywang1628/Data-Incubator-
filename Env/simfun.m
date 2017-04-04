function[FLR ,ENV, Oracle]=simfun(model,dimx,dimy,evals,ntrain,ntest,nsim)


n = ntrain+ntest;
vary = 1;
varx = 0.1;
lint_x = 10;lint_y = 10;
Sgrid = 50;  % y
Tgrid = 60;  % x
regular = 2; % regular =  0 means creating sparse data
criteria= 1; %kernel smoothing parameter selection criteria;1 is AIC; 2 is BIC


FIT = 0;
alpha = 0.05;
param_X = setOptions('selection_k','AIC_R','verbose','on', 'ngrid',Tgrid);
param_Y = setOptions('selection_k','FVE','FVE_threshold',0.9,'verbose','on','ngrid',Sgrid);


tic
FLR = zeros(1,(length(dimx)*length(dimy)));
ENV = zeros(1,(length(dimx)*length(dimy)));
Oracle = 0;
for ss = 1:nsim
rng(ss)
[t_x, s_y, x, y, t, s, beta_st, Ey] = getFunData(n, lint_x, lint_y, Sgrid, Tgrid, regular , varx, vary, evals, model);

ytr = y(:,1:ntrain);
yte = y(:,(ntrain+1):n);
Eyte = Ey(:,(ntrain+1):n);
[WUx, WUt_x, WUy, WUs_y, ~, newx, new_tx, ~, new_ty, ~] = pre(x, t_x, ytr, s_y, ntest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testerFLR=0;
testerEnv=0;
for i= 1:length(dimx)
    for j = 1:length(dimy)
        K_x = dimx(i); K_y = dimy(j); uux=dimx(i); uuy=dimy(j);
        [res,FPCxx,FPCyy]= FPCreg(x, t_x, ytr, s_y,param_X, param_Y,FIT,K_x,K_y,ntest,alpha); 
        newy = getVal(res,'newy');  
        temp1=mean((cell2mat(newy)-cell2mat(yte)).^2);
        testerFLR=[testerFLR temp1];


        [Esimulout] = Env(x, t_x, ytr, s_y, Sgrid, Tgrid, ntest,uux,uuy, lint_x, lint_y);
        newyEnv = Esimulout.newy;
        temp2=mean((cell2mat(newyEnv)-cell2mat(yte)).^2);
        testerEnv=[testerEnv,temp2];
    end
end
FLR = [FLR;testerFLR(2:(length(dimx)*length(dimy))+1)];
ENV = [ENV;testerEnv(2:(length(dimx)*length(dimy))+1)];



FPCmu_xx=getVal(FPCxx, 'mu');
FPCmu_yy=getVal(FPCyy, 'mu');
[newyOracle] = PredicFUN(ntest, newx, FPCmu_xx, FPCmu_yy, beta_st*sqrt(10), new_tx, new_ty, lint_x, lint_y, Tgrid, Sgrid);
testerOracle=mean((cell2mat(newyOracle)-cell2mat(yte)).^2);
Oracle = [Oracle,testerOracle];
end
toc