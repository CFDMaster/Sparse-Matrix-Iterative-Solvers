function [x] = GaussSeidelMatrix( A,X,b,MaxIT,tol )
%Gauss-Seidel Matrix Iteration For Solotion of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%Revison Modifications 1:Making Faster by removing Slow Inverse Rearranging
%WOW!!Over 40 times Faster over Privous Matrix & Scaler Version 4 Simple D FD Test
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
% D=diag(diag(A));---------Removed in R1
% L=-tril(A,-1);---------Removed in R1
% U=-triu(A,1);---------Removed in R1
TrilA=tril(A); % Added in R1
%res=norm(A*X-b); %Compute Residual---------Replaced in R1
r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    %X=inv(D-L)*U*X+inv(D-L)*b; ---------Replaced in R1
    X=X+TrilA\(r);
    %res=norm(A*X-b);---------Replaced in R1
    r=b-A*X;
    res=norm(r); %Compute Residual
    IT=IT+1;
    %fprintf('Gauss-Seidel MatrixIteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Gauss-Seidel Solver Diverged\n');
    fprintf('Gauss-Seidel Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
    x=X;
else
    fprintf('\nGauss-Seidel Matrix Solver Converged');
    fprintf('\nGauss-Seidel Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

