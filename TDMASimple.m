function [X] = TDMASimple(A,d)
%Solves TriDigonal matrix using Thomas TDMA Alghorithm 
% A: Coe. Matrix, d=RHS of the solutiion Returns X=Solution
%   ai*xi-1+bi*xi+ci*x=di
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%First check the matrix
s=size(A); 
l=size(d);
if s(1)~=s(2)
    fprintf('A must be a Square matrix\n')
    return
end
 if s(1)~=l(1)
         fprintf('Size A and d Matrices must agree\n')
    return
 end
n=s(1);
a=zeros(1,n);
b=zeros(1,n);
c=zeros(1,n);
%First fill coe.
for k=1:n
    b(k)=A(k,k);
    if k~=1
        a(k)=A(k,k-1);
    else
        a(1)=0;
    end
    if k~=n
        c(k)=A(k,k+1);
    else
        c(k)=0;
    end 
end

%Forward elimination phase
for k=2:n
    m=a(k)/b(k-1);
    b(k)=b(k)-m*c(k-1);
    d(k)=d(k)-m*d(k-1);   
end

%Backward substitution phase
X(n)=d(n)/b(n);
for k=n-1:-1:1
    X(k)=(d(k)-c(k)*X(k+1))/b(k);
end
% for k=2:n
%     m=A(k,k-1)/A(k-1,k-1);
%     A(k,k)=A(k,k)-m*A(k-1,k);
%     d(k)=d(k)-m*d(k-1);   
% end
% 
% %Backward substitution phase
% X(n)=d(n)/A(n,n);
% for k=n-1:1
%     X(k)=(d(k)-A(k,k+1)*X(k+1))/A(k,k);
% end

end

