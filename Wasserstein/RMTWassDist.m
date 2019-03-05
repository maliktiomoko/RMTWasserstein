function [est,esthat] = RMTWassDist(X,Y)
%Function that compute the Wasserstein distance between Gaussian centered
%distribution based on the article Random  Matrix-Improved  Estimation  of  the  Wasserstein  Distance
%between  two  Centered  Gaussian  Distribution (Malik TIOMOKO & Romain Couillet)
%Input Need the samples from the first class X of dimension p*n and the
%samples from the second class Y of size p*n
%Return the estimate est proposed in the article and the classical esthat
%Define the dimensions
p=size(X,1);
n1=size(X,2);
n2=size(Y,2);
c1=p/n1;c2=p/n2;
%Sample covariance estimate
hatC1=X*X'/n1;hatC2=Y*Y'/n2;
lambda=sort(eig(hatC1*hatC2));
m=@(z) mean(1./(lambda-z));
phi=@(z) z./(1-c1-c1.*z.*m(z));
psi=@(z) 1-c2-c2*z.*m(z);
f=@(z) sqrt(z);
eta=sort(real(eig(diag(lambda)-(1/n1)*sqrt(lambda)*sqrt(lambda)')));
zeta=sort(real(eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
phi_test=@(z) z;
psi_test=@(z) 1;
phipsi=@(z) sqrt(z)/(c2);
for i=1:length(lambda)
    phi_test=@(z) phi_test(z).*((z-lambda(i))./(z-eta(i)));
    psi_test=@(z) psi_test(z).*(z-zeta(i))./(z-lambda(i));
    phipsi=@(z) phipsi(z).*sqrt((z-zeta(i))./(z-eta(i)));
end
% Distinguish the case where n1<n2 to the case where n1>n2
if eta(1)<zeta(1)
    my_eta=zeta;
    my_zeta=eta;
else
    my_zeta=zeta;
    my_eta=eta;
end
other=@(z) 2*sum(1./(z-zeta))-2*sum(1./(z-lambda));
integrand_real=@(z) (1/(2*pi))*2*f(-(phi(z)./psi(z))).*other(z).*(psi(z)/c2);
%Computing the second term (real_integral)
real_integral=0;
for i=1:length(my_zeta)
    real_integral=real_integral+integral(integrand_real,my_zeta(i),my_eta(i));
end
%Computing the first term (pole in lambda)
pole=2*(sqrt(c2/c1))*sum(sqrt(lambda))/c2;
esty=pole+real_integral;
est=(1/p)*trace(hatC1+hatC2)-2*esty;
%Distinguish the case n1=n2
if n1==n2
    est=(1/p)*trace(hatC1+hatC2)-2*(sum(sqrt(lambda))-sum(sqrt(zeta)))*(2*n1/p);
end
%Classical estimate
esthat=(1/p)*(trace(hatC1)+trace(hatC2)-2*trace((hatC1^(1/2)*hatC2*hatC1^(1/2))^(1/2)));
end

