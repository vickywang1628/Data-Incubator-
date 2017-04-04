
function [ccov, xcov, ycov] = getRawcovariance(x, y, mu_x, mu_y)
tt_x = []; tt_y = [];
xx = []; yy = [];
count = [];
n = length(x);
xx = reshape(cell2mat(x),length(x{1}), n)';
mux = repmat(mu_x,n,1);
yy = reshape(cell2mat(y),length(y{1}), n)';
muy = repmat(mu_y,n,1);
ccov = (xx-mux)'*(yy-muy)/n;
xcov = (xx-mux)'*(xx-mux)/n;
ycov = (yy-muy)'*(yy-muy)/n;


xcov = (xcov + xcov')/2; % xcov: 51 *51 double
ycov = (ycov + ycov')/2; % xcov: 51 *51 double
end