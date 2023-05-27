clc
clear

rng(1)
L = 40; % L layers(networks)
M = 3; % number of network classes
innerP=0.6;
% outerP=0.45;
outerP=0.5;

n = 100; % number of vertices
maxrank = 3; % maximum of communities of each layer class
% K = randi(maxrank,[1,M])+1; % Km, number of communities for each layer class

K=maxrank*ones(1,M);

params.K=K;
params.L=L;
params.M=M;
params.n=n;
params.innerP=innerP;
params.outerP=outerP;

%%

for j=1:10
    
    j
    params.L=5*j;
    %     params.outerP=j/10;
    %     params.n=j*20;
    
    [Z,Theta,A]=GenMatrices(params);
    [Q1,W1]=AltMin(A,params); % Xia's Initialization
    [Q2,W2]=AltMin(A,params); %
    
    
    %     W=Z*(Z'*Z)^(-0.5);
    [Zhat,Thetahat]=clusteringKmeans(Q1,W1,params);
    
    [networkError,order]=clusteringErr(Zhat,Z);
    
    communityError=zeros(1,M);
    
    for i=1:M
        temp=zeros(1,M);
        for k=1:M
            temp(k) = clusteringErr2(Thetahat{i},Theta{k});
        end
        communityError(i)=min(temp);
    end
    
    %     fprintf('Network clustering error: %f \n',networkError);
    %     fprintf('Community clustering error: %f \n',sum(communityError)/M);
    %     fprintf('Community Errors: ');
    %     fprintf('%g ', communityError);
    
    
    Err1(j)=sum(communityError)/M;
    
    [Zhat,Thetahat]=clusteringKmeans(Q2,W2,params);
    
    [networkError,order]=clusteringErr(Zhat,Z);
    
    communityError=zeros(1,M);
    
    for i=1:M
        temp=zeros(1,M);
        for k=1:M
            temp(k) = clusteringErr2(Thetahat{i},Theta{k});
        end
        communityError(i)=min(temp);
    end
    
    Err2(j)=sum(communityError)/M;
    
end

%%

figure;
plot((1:10)*5,Err1,'-o','LineWidth',2)
hold on
plot((1:10)*5,Err2,'-d','LineWidth',2)

xlabel('L')
ylabel('Community clustering error')
% title('L=40,M=3,K=3')
title('n=100,M=3,K=3')

legend('Xia Init.','Previous Init.')

