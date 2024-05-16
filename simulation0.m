clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters -                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=[0;0];
rho = 0.5; % Represents the degree of endogeneity
sigma=[[1,rho];[rho,1]];

samplesize = 100;

pi0 = 1;
pi1 = 1;
pi2 = 0.5;
pi3 = 0.2;
pi4 = 0.1;
pi5 = 0.05;
pi6 = 0.01;
pi = [pi1;pi2;pi3;pi4;pi5;pi6];

alpha = [0;0;0;0;0;0]; % Represents Exogeneity of the IV

beta = 1;
beta_constant = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model: Y=X'beta+Z'alpha+e
%            X=Zpi_1+Z^2pi_2+...+u
%            (e,u)~N(mu,sigma)
%            Z~N(0,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error_vector = mvnrnd(mu,sigma,samplesize); % let the first column be e and the second be u
% e = error_vector(:,1);
% u = error_vector(:,2);
% instrument = normrnd(0,1,samplesize,1);
% Z = [instrument,instrument.^2,instrument.^3,instrument.^4, instrument.^5, instrument.^6];
% 
% X = zeros(samplesize,1);
% X = Z*pi+u;
% 
% Y = zeros(samplesize,1);
% Y = beta_constant+X.*beta+instrument.*alpha+e;

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
reps = 1000;
r = 1;
K = 8; % We expand up to K-th order of polynomial
Jstat_matrix = zeros(reps,K-1);
test_matrix = zeros(reps,K-1); % the column size is 3 because the model is up to 4th polynomial of Z
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
    instrument = mvnrnd(0,1,samplesize);
    Z = [instrument,instrument.^2,instrument.^3,instrument.^4, instrument.^5, instrument.^6];
    
    X = zeros(samplesize,1);
    X = pi0+Z*pi+u;
    
    Y = zeros(samplesize,1);
    Y = beta_constant+X.*beta+Z*alpha+e;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OverIdentification tests              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Zmat = [ones(samplesize,1) instrument.^1 instrument.^2 instrument.^3 instrument.^4 instrument.^5 instrument.^6 instrument.^7 instrument.^8]; %Preallocation for speed
    for i = 2:K
        df = i-1; % degree of freedom of the limiting distribution
        % This part is For Abstraction - now it is referenced for speed
%         Zn = [ones(samplesize,1)];
%         for j = 1:K
%             Zn = [Zn instrument^j];
%         end
        
        Zn = Zmat(:,1:i+1);
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

        if(i==2)
            R = [zeros(2,i-1);eye(i-1,i-1,"double")];
            c = zeros(i-1,1);
            W0_1 = ((R'*pi_hat-c)'/(R'/(Zn'*Zn)*R)*(R'*pi_hat-c))/(((X-X_fit)'*(X-X_fit))/(samplesize-5));
        end
    end

    r = r + 1;
end
test = mean(test_matrix,1)
J1 = Jstat_matrix(:,1);
J2 = Jstat_matrix(:,2);
J3 = Jstat_matrix(:,3);
J4 = Jstat_matrix(:,4);
J5 = Jstat_matrix(:,5);
J6 = Jstat_matrix(:,6);
J7 = Jstat_matrix(:,7);

chi2 = chi2rnd(5,reps,1);
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
R = [zeros(2,K-1);eye(K-1,K-1,"double")];
c = zeros(K-1,1);
W0 = ((R'*pi_hat-c)'/(R'/(Zn'*Zn)*R)*(R'*pi_hat-c))/(((X-X_fit)'*(X-X_fit))/(samplesize-5));
chi2cdf(W0,1,"upper")
chi2cdf(W0_1,1,"upper");

histogram(chi2);hold on;histogram(Jstat_matrix(:,5));hold off;