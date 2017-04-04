function[newy] = PredicFUN(ntest, newx, mu_x, mu_y, betafun, new_tx, new_ty, lint_x, lint_y, Tgrid, Sgrid)

%==================================================================
%compute score
%mu_newx = mean(cell2mat(newx'),1);
n = length(newx);
temp.xx = reshape(cell2mat(newx),length(newx{1}), n)';
temp.mux = repmat(mu_x,n,1);
%temp.mux = repmat(mu_newx,n,1);
temp = (temp.xx-temp.mux);



cellnewx = cell(1,n); 
for i = 1:n
    cellnewx{i}=temp(i,:);
end

newxscore = zeros(size(betafun,2), n);
t_beta = linspace(0,lint_x,Tgrid);
for i = 1:size(betafun,2)
    for j = 1:n
        tmp = interp1(t_beta, betafun(:,i), new_tx{j}, 'splines');
        newxscore(i,j) = trapz(new_tx{j}, tmp.*cellnewx{j});
    end
end
newxscore = newxscore';
clear tmp t_beta
%==================================================================
yout1 = linspace(0,lint_y,Sgrid);
ymu = mu_y;
new_ty = new_ty{1};
newymu = interp1(yout1,ymu,new_ty,'spline');
newymat = repmat(newymu,ntest,1) + newxscore;
newy = mat2cell(newymat, ones(1,ntest))';
end