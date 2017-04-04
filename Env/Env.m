function [X] = Env(x, t_x, ytr, s_y, Sgrid, Tgrid, ntest,uux,uuy, lint_x, lint_y)


[x, t_x, y, s_y, ~, newx, new_tx, ~, new_ty, ~] = pre(x, t_x, ytr, s_y, ntest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_x = mean(cell2mat(x'),1);
mu_y = mean(cell2mat(y'),1);
[covxy, xcov, ycov] = getRawcovariance(x, y, mu_x, mu_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%========================================================================
%MX
out21 = linspace(0,lint_x,Tgrid);
out1 = out21;
K = 2;
d = 1;
R = covxy;
T = covxy;
for k=2:K
        for j = 1:d
            for r=1:length(out21); %out21 : 1 * 51 double
                A(r,j)=trapz(out21, xcov(r,:).*T(:,j)');
            end
        end
        T = A;
        R= [R,A];
    end
MX = R;
W = MX * MX';
W = (W+W')/2;
[Wlam, ~, Wphi21, ~] = getEigens(W, out1, out21, 20);
[~, I] = sort(Wlam, 'descend');
PsiphiX = Wphi21(:,I);  %Lphi21 : 51 * 20 double

clear out21 Lphi21 V R T A K dd L out1 W Wlam Wphi21 
%-----------------------------------------------------------------------
%MY
out21 = linspace(0,lint_y,Sgrid);
out1 = out21;
K = 2;
d = 1;
R = covxy';
T = covxy';
    for k=2:K
        for j = 1:d
            for r=1:length(out21); %out21 : 1 * 51 double
                A(r,j)=trapz(out21, ycov(r,:).*T(:,j)');
            end
        end
        T = A;
        R= [R,A];
    end
MY = R;

W = MY * MY';
W = (W+W')/2;
[Wlam, ~, Wphi21, ~] = getEigens(W, out1, out21, 20);
[~, I] = sort(Wlam, 'descend');
PsiphiY  = Wphi21(:,I);  %Lphi21 : 51 * 20 double


clear out21 Lphi21 V R T A K dd L out1 W Wlam Wphi21 


% Project X onto e.d.r. space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(x,2);
PSIprojX = zeros(size(PsiphiX,2), n);
for i = 1:size(PsiphiX,2)
    for j = 1:n
        PSIprojX(i,j) = trapz(t_x{j}, PsiphiX(:,i)'.*x{j});
    end
end
clear tmp t_beta


% Project Y onto e.d.r. space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSIprojY = zeros(size(PsiphiY,2), n);
for i = 1:size(PsiphiY,2)
    for j = 1:n
        PSIprojY(i,j) = trapz(s_y{j}, PsiphiY(:,i)'.*y{j});
    end
end
clear tmp t_beta


beta = mvregress(PSIprojX(1:uux,:)',PSIprojY(1:uuy,:)');
betafun = PsiphiX(:,1:uux)*beta*PsiphiY(:,1:uuy)';


[newy] = PredicFUN(ntest, newx, mu_x, mu_y, betafun, new_tx, new_ty, lint_x, lint_y, Tgrid, Sgrid);

%----------------------------------------------------------------------

% Output
X.PsiphiX = PsiphiX;
X.PsiphiY = PsiphiY;
X.projX = PSIprojX;
X.projY = PSIprojY; 
X.beta =  betafun;
X.b = beta;
X.newy = newy;
X.new_ty=new_ty;
X.new_tx=new_tx;
X.uux=uux;
X.uuy=uuy;
end