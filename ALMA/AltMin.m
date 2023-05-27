function [Q,W]=AltMin(A,params)
% alternating minimization

K=params.K;
L=params.L;
M=params.M;
n=params.n;

Kmax=max(K);

%% initialization ??

Uhat=squeeze(sum(A,1));

DEG=squeeze(sum(A,2));
DEG=sum(DEG,1);
delta=2*sqrt(Kmax)*max(DEG)/sqrt(sum(DEG.^2));


% I=eye(n);
% delta=0;
% for i=1:n
%     delta=max(delta,norm(I(i,:)*Uhat));
% end

Uhat1=Uhat;

for i=1:n
    d=norm(Uhat(1,:));
    Uhat1(i,:)=Uhat(1,:)*min(delta,d)/d;
end

[U,~,~]=svd(Uhat1); % sum(A,1) returns sum of horizontal slices
% [U,~,~]=svd(Uhat);
Qini=zeros(L,Kmax,Kmax);

for i=1:L
    [U1,~,~]=svd(U'*squeeze(A(i,:,:))*U);
    Qini(i,:,:)=U1(:,1:Kmax)'*squeeze(A(i,:,:))*U1(:,1:Kmax);
end


D=zeros(L,Kmax*Kmax);

for i=1:L
    temp=Qini(i,:,:);
    temp=temp(:);
    D(i,:)= temp;
end

label=kmeans(D,M);

Z = zeros(L,M); % clustering matrix for networks

for i=1:L
    Z(i,label(i))=1;
end

% W=orth(W);

W=Z*(Z'*Z)^(-0.5);
 

%%


orthoprojector = @(X) X*(X'*X)^(-0.5);

iterMax=200;
etol=1e-3;
err = Inf;
idx = 0;

Amatrix=reshape(A,[L,n^2]); % mode-1 unfolding of tensor A

while err>etol && idx<iterMax
    idx = idx+1;
    Wold=W;
    
    Qmatrix=W'*Amatrix;
    Q=reshape(Qmatrix,[M,n,n]);
    Q=Kprojector(Q,K);
    Qmatrix=reshape(Q,[M,n^2]);
    W=orthoprojector(Amatrix*Qmatrix');
    err=norm(Wold-W,'fro');
end


fprintf('Objective error: %f \n',err)
fprintf('# of iterations: %d \n',idx)


end





%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function Y = Kprojector(X,K)
% projectionto the nearest rank-Km matrix for each horizontal slice

M=size(X,1);
Y=zeros(size(X));

for i=1:M
    k=K(i);
    [U,S,V] = svd(squeeze(X(i,:,:)));
    Y(i,:,:) = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
end


end