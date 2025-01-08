function [L,U]=croutsAlgorithm(A)

[n,n]=size(A);
L=zeros(n,n);
U=zeros(n,n);

for i=1:n
    L(i,1)=A(i,1);
end
for j=2:n
    U(1,j)=A(1,j)/L(1,1);
end

for j=2:n
    for i=j:n
        sum=0;
        for k=1:j-1
            sum=sum+L(i,k)*U(k,j);
        end
        L(i,j)=A(i,j)-sum;
    end
        
    for k=j+1:n
        sum=0;
        for i=1:j-1
            sum=sum+L(j,i)*U(i,k);
        end
        U(j,k)=(A(j,k)-sum)/L(j,j);
    end
end

for i=1:n
   U(i,i)=1.0; 
end
    
end     