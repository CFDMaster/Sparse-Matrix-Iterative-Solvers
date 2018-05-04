function [x] = SymmetricSOR(A,X,b,Omega,MaxIT,tol)
%Symmetric SOR Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X, Omega is over relaxation Coe.
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%   Omega: Over Relaxation Coe.

res=norm(A*X-b); %Compute Residual
IT=1;
n=size(A,1);
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    Xold=X;
    for i=1:n
        jj=[1:i-1 i+1:n];
        X(i)=(1-Omega)*Xold(i)+Omega*(b(i)-sum(A(i,jj)*X(jj)))/A(i,i);
    end
    for i=n:-1:1
        jj=[1:i-1 i+1:n];
        X(i)=(1-Omega)*Xold(i)+Omega*(b(i)-sum(A(i,jj)*X(jj)))/A(i,i);
    end
    
    res=norm(A*X-b); %Compute Residual
    IT=IT+1;
%         fprintf('Symmetric SOR Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Symmetric SOR Solver Diverged\n');
    fprintf('Symmetric SOR Iteration=%i\tResidual=%2.6f\n',IT,res);
    x=X;
else
    fprintf('\nSymmetric SOR Solver Converged');
    fprintf('\nSymmetric SOR Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

