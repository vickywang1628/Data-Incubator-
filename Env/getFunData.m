function [t_x, s_y, x, y, t, s, beta_st, Ey, phi_x, phi_y] = getFunData(n, lint_x, lint_y, Sgrid, Tgrid, regular , varx, vary, evals, model)


% Create regression coefficient functions
if model == 1
    b = [3 0 0 0; 0 3 0 0; 0 0 3 0; 0 0 0 3];
elseif model == 2
    btemp = ones(10,10)/5;
    btemp(7:10,7:10) = [4 0.2 0.2 0.2; 0.2 4 0.2 0.2; 0.2 0.2 4 0.2; 0.2 0.2 0.2 4];
    b=btemp;
elseif model == 3
    btemp = zeros(10,10);
    btemp(7:10,7:10) =  [2 1 1 4 ;  .2 5 2 2;  2 1 1 .2; 1 2 1 2]*.5+diag([2 1.5 1 2.5])*1.5;
    b = btemp; 
end


numeig_x = size(b,1); numeig_y = size(b,2);


s = linspace(0, lint_y, Sgrid); %y
t = linspace(0, lint_x, Tgrid); %x
t_x = cell(1,n); 
s_y = cell(1,n); 
Ex = cell(1,n); 
Ey = cell(1,n);
x = cell(1,n); 
y = cell(1,n);
phi_x = xeig(t,lint_x,numeig_x)/sqrt(10);      % 30 eigenfunctions for X
phi_y = xeig(s,lint_y,numeig_y)/sqrt(10);      % 30 eigenfunctions for Y 
mu_x = mu_true(t,Tgrid);     
mu_y = mu_true(s,Sgrid); 

pc_x = randn(n, numeig_x)*diag(sqrt(evals(1:numeig_x)));   % FPC score for X

% calculate the true beta function
beta_st = zeros(length(t),length(s));


gamma = zeros(numeig_y,1);
for i=1:numeig_x
    for j=1:numeig_y
        beta_st = beta_st+b(i,j)*phi_x(i,:)'*phi_y(j,:);
        gamma(j) = sum(b(:,j).^2*evals(i));
    end
end

%plot(gamma)


% Create data
for i = 1:n
    % Sparse case:
    if regular == 0
        % Get number of time points for subject
        nwu = floor(random('unif',5,11,1,1));
        t_x{i} = sort(random('unif',0,lint,1,nwu));
    % Dense Case:
    elseif regular == 2
        % Create dense time points
        t_x{i} = t;
        s_y{i} = s;
    end
        
    % Create eigenfunctions
    phi = xeig(t, lint_x, size(pc_x, 2)); 
    
    % True x
    Ex{i} = mu_x+ pc_x(i,:)*phi;
    tmp = interp1(t, Ex{i}, t_x{i}, 'splines');
    
    % Observed x
    x{i} = tmp + sqrt(varx)*randn(1,Tgrid);
    clear tmp
    
    
    Ey{i} = mu_y+(pc_x(i,:)*b)*phi_y;
    y{i} = interp1(s,Ey{i},s_y{i},'spline')+sqrt(vary)*randn(1,Sgrid);
    
end


end