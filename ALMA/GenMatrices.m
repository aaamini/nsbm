function [Z,Theta,A]=GenMatrices(params)
% generate matrices, Z: clustering matrix for networks, Theta: clustering
% matrix for communities of each layer class, A: observation tensor

K=params.K;
L=params.L;
M=params.M;
n=params.n;
innerP=params.innerP;
outerP=params.outerP;

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


end