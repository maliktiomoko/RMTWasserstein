clear all
clc
close all
p=64;n1=1024;n2=2048;
C1=toeplitz(0.2.^(0:p-1));
C2=toeplitz(0.4.^(0:p-1));
n_simulation=10;
for i=1:n_simulation
    X=zeros(p,n1);
    Y=zeros(p,n2);
for k=1:n1
    X(:,k) = mvnrnd(zeros(1,p),C1);
end
for k=1:n2
    Y(:,k) = mvnrnd(zeros(1,p),C2);
end
%Proposed Wasserstein distance
[est(i),esthat(i)] = RMTWassDist(X,Y);
%Real Wasserstein distance
est_vrai(i)=(1/p)*(trace(C1)+trace(C2)-2*trace((C1^(1/2)*C2*C1^(1/2))^(1/2)));
end
est_mean=mean(est)
esthat_mean=mean(esthat)
est_vrai_mean=mean(est_vrai)