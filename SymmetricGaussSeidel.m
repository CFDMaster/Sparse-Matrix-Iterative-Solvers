function [x] = SymmetricGaussSeidel(A,X,b,MaxIT,tol)
%Symmetric Gauss-Seidel Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
% X must contain initial guess for X
% Symmetric verson applies forward and backward sweep to accelerate iteration
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged

res=norm(A*X-b); %Compute Residual
IT=1;
n=size(A,1);
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    for i=1:n %Forward Sweep
        jj=[1:i-1 i+1:n];
        X(i)=(b(i)-sum(A(i,jj)*X(jj)))/A(i,i);
    end
    for i=n:-1:1  %Backword Sweep
        jj=[n:-1:i+1 i-1:-1:1];
        X(i)=(b(i)-sum(A(i,jj)*X(jj)))/A(i,i);
    end
res=norm(A*X-b); %Compute Residual
IT=IT+1;
%     fprintf('Symmetric Gauss-Seidel Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Gauss-Seidel Solver Diverged\n');
    fprintf('Symmetric Gauss-Seidel Iteration=%i\tResidual=%2.6f\n',IT,res);
    x=X;
else
    fprintf('\nSymmetric Gauss-Seidel Solver Converged');
    fprintf('\nSymmetric Gauss-Seidel Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

