function [x] = SORMatrix( A,X,b,Omega,MaxIT,tol)
%SOR Matrix Iteration For Solotion of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%Revison Modifications 1:Making Faster by removing Slow Inverse Rearranging
%WOW!!Over 10 times Faster over Privous Matrix & Scaler Version 4 Simple FD Test
%Among the Best Iterative Solvers with Optimal Relaxation Coe.
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%   omega: Over Relaxation Coe
D=diag(diag(A));
L=tril(A,-1); % ---------Sign Changed From - to + in R1
Bsor=Omega*(D+Omega*L)\eye(size(A));
% U=-triu(A,1);---------Not required for Forward SOR----Removed in R1
% res=norm(A*X-b); %Compute Residual---------Replaced in R1
r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    X=X+Bsor*r;
    %X=inv(D-Omega*L)*(Omega*U+(1-Omega)*D)*X+inv(D-Omega*L)*Omega*b; %SOR Iteration ---------Replaced in R1
    %res=norm(A*X-b);---------Replaced in R1
    r=b-A*X;
    res=norm(r); %Compute Residual
    IT=IT+1;
    %fprintf('SOR Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. SOR Solver Diverged\n');
    fprintf('SOR Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
    x=X;
else
    fprintf('\nSOR Matrix Solver Converged');
    fprintf('\nSOR Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

