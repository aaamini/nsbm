function Y=tensorMult1(X,U)
% 1-mode product of tensor
[~,n2,n3]=size(X);
[m1,~]=size(U);
Y=zeros(m1,n2,n3);

for i=1:n3
    Y(:,:,i)=U*X(:,:,i);
end

end