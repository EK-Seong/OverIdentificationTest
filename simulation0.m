
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters -                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=[0;0];
rho = 0.1; % Represents the degree of endogeneity
sigma=[[1,rho];[rho,1]];

samplesize = 1000;

pi1 = 1;
pi2 = 0.8;
pi3 = 0.5;
pi4 = 0.1;
pi = [pi1;pi2;pi3;pi4];

alpha = 0.1; % Represents Exogeneity of the IV

beta = 1;
beta_constant = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model: Y=X'beta+Z'alpha+e
%            X=Zpi_1+Z^2pi_2+...+u
%            (e,u)~N(mu,sigma)
%            Z~N(0,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_vector = mvnrnd(mu,sigma,samplesize); % let the first column be e and the second be u
e = error_vector(:,1);
u = error_vector(:,2);
instrument = normrnd(0,1,samplesize,1);
Z = [instrument,instrument.^2,instrument.^3,instrument.^4];

X = zeros(samplesize,1);
X = Z*pi+u;

Y = zeros(samplesize,1);
Y = beta_constant+X.*beta+instrument.*alpha+e;

% %% OLS
% x_ols = [ones(samplesize,1) X];
% beta_ols = (x_ols'*x_ols)\x_ols'*Y;
% x_axis = -2:0.1:8;
% Y_ols = beta_ols(1,1)+x_axis.*beta_ols(2,1);
% min(Y)
% max(Y)
% 
% %% IV

%% Overidentification Test
test_size = 0.05; % test size alpha
reps = 5000;
r = 1;
Jstat_matrix = zeros(reps,3);
test_matrix = zeros(reps,3); % the column size is 3 because the model is up to 4th polynomial of Z
while r <= reps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The model: Y=X'beta+Z'alpha+e
    %            X=Zpi_1+Z^2pi_2+...+u
    %            (e,u)~N(mu,sigma)
    %            Z~N(0,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    error_vector = mvnrnd(mu,sigma,samplesize); % let the first column be e and the second be u
    e = error_vector(:,1);
    u = error_vector(:,2);
    instrument = normrnd(0,1,samplesize,1);
    Z = [instrument,instrument.^2,instrument.^3,instrument.^4];
    
    X = zeros(samplesize,1);
    X = Z*pi+u;
    
    Y = zeros(samplesize,1);
    Y = beta_constant+X.*beta+instrument.*alpha+e;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OverIdentification tests              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 2:4
        df = i-1;
        Zn = [ones(samplesize,1) Z(:,1:i)];
        % 2SLS
        % First stage:
        pi_hat = (Zn'*Zn)\Zn'*X;
        X_fit = Zn*pi_hat;
        % Second stage:
        Xn = [ones(samplesize,1) X_fit];
        beta_2sls = (Xn'*Xn)\Xn'*Y;
    
        % 2-step efficient GMM
        Omega_hat = diag((Y-Xn*beta_2sls).^2);
        beta_egmm = ((Xn'*Zn)/(Zn'*Omega_hat*Zn)*(Zn'*Xn))\((Xn'*Zn)/(Zn'*Omega_hat*Zn)*(Zn'*Y));
    
        % Test statistic
        nJn = (Zn'*(Y-Xn*beta_egmm))'/(Zn'*Omega_hat*Zn)*(Zn'*(Y-Xn*beta_egmm));

        % Save the count of rejection and J-stat to each matrices
        test_matrix(r,i-1) = int8(test_size > chi2cdf(nJn,df,"upper"));
        Jstat_matrix(r,i-1) = nJn;
    end

    r = r + 1;
end
test = mean(test_matrix,1)
J1 = Jstat_matrix(:,1);
J2 = Jstat_matrix(:,2);
J3 = Jstat_matrix(:,3);

chi2 = chi2rnd(3,reps,1);
% J1 = sort(J1,1,"ascend","ComparisonMethod","real");
% J2 = sort(J2,1,"ascend","ComparisonMethod","real");
% J3 = sort(J3,1,"ascend","ComparisonMethod","real");

% figure(1)
% histogram(J1,50)
% figure(2)
% histogram(J2,50)
% figure(3)
% histogram(J3,50)
% figure(4)
% histogram(chi2,50)
% chi2inv(1.2426e-299,3)

% F-Test
R = [[0,0,0];[0,0,0];eye(3,3,"double")];
c = [0;0;0];
W0 = ((R'*pi_hat-c)'/(R'/(Zn'*Zn)*R)*(R'*pi_hat-c))/(((X-X_fit)'*(X-X_fit))/(samplesize-5));
chi2cdf(W0,1,"upper")