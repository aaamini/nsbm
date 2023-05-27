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