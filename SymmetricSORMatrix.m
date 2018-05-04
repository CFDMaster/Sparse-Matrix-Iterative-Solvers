function [x] = SymmetricSORMatrix(A,X,b,Omega,MaxIT,tol)
%Symmetric SOR Matrix  Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%Revison Modifications 1:Making Faster by removing Slow Inverse Rearranging
%WOW!!Over 60 times Faster over Privous Matrix & Scaler Version 4 Simple FD Test
%Best Iterative Solver(Expet Than CG) with Optimal Relaxation Coe.
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%   Omega: Over Relaxation Coe.

D=diag(diag(A));
L=tril(A,-1);%---------Sign Changed in R1
U=triu(A,1);%---------Sign Changed in R1
Bssor=Omega*(2-Omega)*((D+Omega*U)\eye(size(A)))*D*((D+Omega*L)\eye(size(A))); %-------Added in R1
% res=norm(A*X-b); %Compute Residual---------Replaced in R1
r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    %Xhalf=inv(D-Omega*L)*(Omega*U+(1-Omega)*D)*X+inv(D-Omega*L)*Omega*b;
    %%Symmetric SOR Iteration Forward;---------Removed in R1
    %X=inv(D-Omega*U)*(Omega*L+(1-Omega)*D)*Xhalf+inv(D-Omega*U)*Omega*b;
    %%Symmetric SOR Iteration Backward;---------Removed in R1
    X=X+Bssor*r;
    r=b-A*X;
    res=norm(r); %Compute Residual
    IT=IT+1;
    %     fprintf('Symmetric SOR Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Symmetric SOR Solver Diverged\n');
    fprintf('Symmetric SOR(Matrix) Iteration=%i\tResidual=%2.6e\n',IT,res);
    x=X;
else
    fprintf('\nSymmetric SOR(Matrix) Solver Converged');
    fprintf('\nSymmetric SOR(Matrix) Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

