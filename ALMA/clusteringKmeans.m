function [Zhat,Thetahat]=clusteringKmeans(Q1,W1,params)


K=params.K;
L=params.L;
M=params.M;
n=params.n;

labelhat=kmeans(W1,M,'Replicates',100);
Zhat = zeros(L,M); % clustering matrix for networks

for i=1:L
    Zhat(i,labelhat(i))=1;
end


% err=clusteringErr2(Zhat,Z); % network error


label_comhat=cell(1,M);
Thetahat=cell(1,M);


err_com=zeros(1,M);
deg=sum(Zhat,1);

for i=1:M

    Qslice=squeeze(Q1(i,:,:));
    k=K(i);
    [U,S,V]=svd(Qslice);
    label_comhat{i}=kmeans(U(:,1:k),k,'Replicates',100);
%     label_comhat{i}=kmeans(Qslice,k,'Replicates',100);
    
    Thetahat{i}=zeros(n,k);
    for j=1:n
        Thetahat{i}(j,label_comhat{i}(j))=1;
    end
    
%     err_com(i)=clusteringErr2(Thetahat{i},Theta{i});
    
end



end