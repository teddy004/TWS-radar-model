function [costMat] = gattting1(range,d,m)

n=numel(range);
k=numel(m);
costMat=zeros(n,k);
    
for i=1:n
    for j=1:k
        if and((range(i)-d)<= m(j) , m(j) <=(range(i)+d))
            costMat(i,j)=m(j);
            
        else
            costMat(i,j)= inf;
            
        end
        
        costMat(i,j)=costMat(i,j);
    end
    
end

