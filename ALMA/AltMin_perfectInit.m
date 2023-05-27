function [Q,W]=AltMin_perfectInit(A,Z,params)
% alternating minimization

K=params.K;
L=params.L;
M=params.M;
n=params.n;

Kmax=max(K);

%% initialization ??
% 
% Uhat=squeeze(sum(A,1));
% 
% I=eye(n);
% delta=0;
% for i=1:n
%     delta=max(delta,norm(I(i,:)*Uhat));
% end
% 
% Uhat1=Uhat;
% 
% for i=1:n
%     d=norm(Uhat(1,:));
%     Uhat1(i,:)=Uhat(1,:)*min(delta,d)/d;
% end
% 
% [U,~,~]=svd(Uhat1); % sum(A,1) returns sum of horizontal slices
% % [U,~,~]=svd(Uhat);
% Qini=zeros(L,Kmax,Kmax);
% 
% for i=1:L
%     [U1,~,~]=svd(U'*squeeze(A(i,:,:))*U);
%     Qini(i,:,:)=U1(:,1:Kmax)'*squeeze(A(i,:,:))*U1(:,1:Kmax);
%     
% end
% 
% D=zeros(L,Kmax*Kmax);
% 
% for i=1:L
%     temp=Qini(i,:,:);
%     temp=temp(:);
%     D(i,:)= temp;
% end
% 
% label=kmeans(D,M);
% 
% W = zeros(L,M); % clustering matrix for networks
% 
% for i=1:L
%     W(i,label(i))=1;
% end


% W=orth(Z);

W=Z*(Z'*Z)^(-0.5);

%%


orthoprojector = @(X) X*(X'*X)^(-0.5);

iterMax=1;
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
    err=norm(Wold-W);
    
    W;
    
end


fprintf('Objective error: %f \n',err)
fprintf('# of iterations: %d \n',idx)


end