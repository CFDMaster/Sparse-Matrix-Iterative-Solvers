function [x] = SymmetricGaussSeidelMatrix(A,X,b,MaxIT,tol)
%Symmetric Gauss-Seidel Matrix Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%By Mohammad Aghakhani
%Usesd Forward and Backword Sweep to Improve Convergence
%Revison Modifications 1:Making Faster by removing Slow Inverse Rearranging
%WOW!!Over 60 times Faster over Privous Matrix & Scaler Version 4 Simple D FD Test
% D=diag(diag(A));%---------Removed in R1
L=-tril(A,-1);
U=-triu(A,1);
TriLa=tril(A)\eye(size(A));
TriUa=triu(A)\eye(size(A));
TriLaU=TriLa*U;
TriUaL=TriUa*L;
% Bggs=(triu(A)\eye(size(A)))*diag(diag(A))*(triu(A)\eye(size(A))); %
% =inv(D+U)*D*(D+L)----Added in R1--Not implemented Due to lower Accuracy
res=norm(A*X-b); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    X=TriLaU*X+TriLa*b; %Forward Sweep---------Replaced in R1
    X=TriUaL*X+TriUa*b; %Backrward Sweep-------Replaced in R1   
    %X=X+Bggs\(r);----Not implemented Due to lower Accuracy
    res=norm(A*X-b);
    IT=IT+1;
    %     fprintf('Symmetric Gauss-Seidel Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Gauss-Seidel Solver Diverged\n');
    fprintf('Symmetric Gauss-Seidel Iteration=%i\tResidual=%2.6e\n',IT,res);
    x=X;
else
    fprintf('\nSymmetric Gauss-Seidel Solver Converged');
    fprintf('\nSymmetric Gauss-Seidel Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

