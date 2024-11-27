clear

rho = 0.5;
pi = [1;0.1;0.01];
alpha = [1;0];
beta = [1;1];
n = 100;
rep = 1000;
testsize=0.05;

b2sls = zeros(rep,size(beta,1));
b2step = zeros(rep,size(beta,1));
Jstat = zeros(rep,1);
rejection = zeros(rep,1);

for i=1:rep
    for K=100:100:rep
        if i==K
            i
        end
    end
    
    u=mvnrnd([0;0],[[1,rho];[rho,1]],n);
    z=mvnrnd([0;0],[[1,rho];[rho,1]],n);
    x=[ones(n,1),z]*pi+u(:,2);
    y=[ones(n,1),x]*beta+z*alpha+u(:,1);
    
    Z=[ones(n,1),z];
    X=[ones(n,1),x];
    Pz=Z/(Z'*Z)*Z';
    b=(X'*Pz*X)\X'*Pz*y;
    b2sls(i,:)=b';

    omega = diag((y-X*b).^2);
    S=(Z'*omega*Z/n);
    Pw = Z/S*Z';
    % Sc=(Z'*omega*Z/n)-(Z'*(y-X*b)/n)*(Z'*(y-X*b)/n)';
    % Pw = Z/Sc*Z';
    b_egmm = (X'*Pw*X)\X'*Pw*y;
    b2step(i,:)=b_egmm';
    nJn=n*(Z'*(y-X*b_egmm)/n)'/S*(Z'*(y-X*b_egmm)/n);
    % nJn=n*(Z'*(y-X*b_egmm)/n)'/Sc*(Z'*(y-X*b_egmm)/n);
    Jstat(i,1)=nJn;
    rejection(i,1)=chi2cdf(nJn,1,"upper")<testsize;

end
% mean(b2sls,1)
% std(b2sls)
% mean(b2step,1)
% std(b2step)
% mean(Jstat,1)
% std(Jstat)

chi = chi2rnd(1,rep,1);
mean(rejection,1)

histogram(chi);hold on;histogram(Jstat);hold off;