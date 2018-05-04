function [x] = PCG(A,X,b,B,MaxIT,tol )
%Perconditioned Conjugate Gradient Method For Solotion of a linear system Ax=b
%B is the Perconditioning Matrix
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
r=b-A*X;
IT=1;
rho=1;
display('Please wait..... Calculating the solution')
while(sqrt(rho)>=tol) && (IT<=MaxIT)
    Br = B*r;
    rho=dot(r,Br);
    if IT==1
        p=Br;
    else
        beta = rho/rho_old;
        p = Br + beta*p;
    end
    Ap=A*p;
    alpha=rho/dot(Ap,p);
    X=X+alpha*p;
    r=r-alpha*Ap;
    rho_old=rho;
    IT=IT+1;
    %fprintf('PCG Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT || sqrt(rho)>tol
    fprintf('\nMaximum Iteratons Reached. PCG Solver Diverged\n');
    fprintf('PCG Iteration=%i\tResidual=%2.6f\n',IT,sqrt(rho));
    x=X;
else
    fprintf('\nPCG Solver Converged');
    fprintf('\nPCG Iteration=%i\tResidual=%2.6e\n',IT,sqrt(rho));
end
x=X;
end
