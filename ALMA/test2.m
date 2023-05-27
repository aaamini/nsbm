clc
clear

rng(1)

% outerP=0.4*innerP;

for q=1:30
    
    clearvars -except ERR q
    
    q
    
    L = 60; % L layers(networks)
    M = 2; % number of network classes
    innerP=0.6;
    outerP=0.5;
    n = q*5; % number of vertices
    maxrank = 2; % maximum of communities of each layer class
    K = randi(maxrank,[1,M])+1; % Km, number of communities for each layer class
    
    params.K=K;
    params.L=L;
    params.M=M;
    params.n=n;
    
    %% true data construction
    label = randi(M,[1,L])';
    
    % make sure each class of network has one layer
    while length(unique(label))~=M
        label = randi(M,[1,L])';
    end
    
    Z = zeros(L,M); % clustering matrix for networks
    
    for i=1:L
        Z(i,label(i))=1;
    end
    
    Theta = cell(1,M); % clustering matrix for communities of each layer class
    label_com = cell(1,M); % community labels for each layer class
    B = cell(1,M); % connectivity matrix for each layer class
    Q = zeros(M,n,n); % stack of probability matrices for layer classes
    
    for m=1:M
        Km=K(m);
        label_com{m} = randi(Km,[1,n])';
        
        while length(unique(label_com{m}))~=Km
            label_com{m} = randi(Km,[1,n]);
        end
        
        
        Theta{m}=zeros(n,Km);
        for i=1:n
            Theta{m}(i,label_com{m}(i))=1;
        end
        
        B{m}=outerP*ones(Km)+(innerP-outerP)*eye(Km); % inner probability: 0.9, outer: 0.1
        Q(m,:,:)=Theta{m}*B{m}*Theta{m}';
        
    end
    
    P = zeros(L,n,n); % stack of probability matrices for all layers
    for i=1:L
        l=label(i);
        P(i,:,:)=Q(l,:,:);
    end
    
    RAND = rand([L,n,n]); % generate L*n*n random matrix with element in (0,1)
    A = zeros(L,n,n); % observation tensor
    Index = find(RAND<P);
    A(Index) = 1;
    
    %%
    
        [Q1,W1]=AltMin(A,params); % recover tensor Q and matrix W
%     [Q1,W1]=AltMin_perfectInit(A,Z,params);
    
    
    %% clustering based on Q1 and W1
    
    labelhat=kmeans(W1,M,'Replicates',100);
    Zhat = zeros(L,M); % clustering matrix for networks
    
    for i=1:L
        Zhat(i,labelhat(i))=1;
    end
    
    
    err=clusteringErr2(Zhat,Z); % network error
    
    
    label_comhat=cell(1,M);
    Theta_hat=cell(1,M);
    
    
    err_com=zeros(1,M);
    deg=sum(Zhat,1);
    
    for i=1:M
        
        Qslice=squeeze(Q1(i,:,:));
        k=K(i);
        [U,S,V]=svd(Qslice);
        label_comhat{i}=kmeans(U(:,1:k),k,'Replicates',100);
        %     label_comhat{i}=kmeans(Qslice,k,'Replicates',100);
        
        Theta_hat{i}=zeros(n,k);
        for j=1:n
            Theta_hat{i}(j,label_comhat{i}(j))=1;
        end
        
        err_com(i)=clusteringErr2(Theta_hat{i},Theta{i});
        
    end
    
    ERR(q)=sum(err_com)/M;
    
    
end

plot((1:30)*10,ERR,'LineWidth',2)
xlabel('n')
ylabel('error rate')
title('L=60,M=2')

