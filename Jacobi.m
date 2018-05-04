function [x] = Jacobi(A,X,b,MaxIT,tol )
%Jacobi Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%Revison Modifications 1:Making Faster by removing Slow SUM Function
%Over 4 times Faster over Scaler Version 4 Simple FD Test
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged

Xold=zeros(size(X));
res=norm(A*X-b); %Compute Residual
IT=1;
n=size(A,1);
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    Xold=X;
    for i=1:n
        %jj=[1:i-1 i+1:n]; %----------Removed in R1
        %X(i)=(b(i)-sum(A(i,jj)*Xold(jj)))/A(i,i);%----------Replaced in R1
        X(i)=(b(i)-( A(i,:)*Xold-A(i,i)*Xold(i)) )/A(i,i);
    end
    res=norm(A*X-b); %Compute Residual
    IT=IT+1;
    %fprintf('Jacobi Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Jacobi Solver Diverged\n');
    fprintf('Jacobi Iteration=%i\tResidual=%2.6f\n',IT,res);
    x=X;
else
    fprintf('\nJacobi Solver Converged');
    fprintf('\nJacobi Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

