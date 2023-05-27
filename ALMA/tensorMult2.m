function Y=tensorMult2(X,U)
% 1-mode product of tensor
[n1,~,n3]=size(X);
[m2,~]=size(U);
Y=zeros(n1,m2,n3);

for i=1:n3
    Y(:,:,i)=U*X(:,:,i);
end

end