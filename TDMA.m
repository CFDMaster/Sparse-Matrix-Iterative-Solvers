function [X] = TDMA(A,d)
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
    if s(1)~=l(2)
        fprintf('Size A and d Matrices must agree\n')
        return
    end
    
end
if l(1)~=1 && l(2)~=1
    fprintf('Second Matix must be vector or row\n')
    return
end
% display('Please wait..... Calculating the solution by TDMA Solver')
n=s(1);
X=zeros(n,1);
for k=2:n
    m=A(k,k-1)/A(k-1,k-1);
    A(k,k)=A(k,k)-m*A(k-1,k);
    d(k)=d(k)-m*d(k-1);
end

%Backward substitution phase
X(n)=d(n)/A(n,n);
for k=n-1:-1:1
    X(k)=(d(k)-A(k,k+1)*X(k+1))/A(k,k);
end

end

