function [err,order]=clusteringErr(Zhat,Z)

[L,M]=size(Z);

Order=perms(1:M); % all possible orders of labels
n=length(Order);

diff=Inf;

for i=1:n
    
    temp=sum(sum(abs(Zhat(:,Order(i,:))-Z)));
    
    if temp<=diff
        diff=temp;
        order=Order(i,:);
    end
    
end

err=diff/2/L;



end