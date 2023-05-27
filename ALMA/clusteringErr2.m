function err=clusteringErr2(Zhat,Z)

[L,M]=size(Z);

Order=perms(1:M); % all possible orders of labels
n=length(Order);

diff=zeros(n,1);
for i=1:n
%     diff(i)=norm(Zhat(:,Order(i,:))-Z,1);
    diff(i)=sum(sum(abs(Zhat(:,Order(i,:))-Z)));
end

err=min(diff)/2/L;




end