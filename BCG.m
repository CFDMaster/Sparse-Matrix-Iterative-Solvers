function [x] = BCG(A,X,b,MaxIT,tol )
%BiConjugate Gradient Method For Solution of a linear system Ax=b-Matrix Form
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
tol = tol*norm(b);
r=b-A*X;
r_star=10.*r.*rand(size(r));
dot(r,r_star)
p=r;
p_star=r_star;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    Ap=A*p;
    alpha=dot(r,r_star)/dot(Ap,p_star);
    X=X+alpha*p;
    r_old=r;
    rstar_old=r_star;
    r=r-alpha*Ap;
    r_star=r_star-alpha*A'*p_star;  
    betha=dot(r,r_star)/dot(r_old,rstar_old);
    p=r+betha*p;
    p_star=r_star+betha*p_star;
    res=norm(r); %Compute Residual
    %fprintf('BiConjugate Gradien Iteration=%i\tResidual=%2.6e\n',IT,res);
    IT=IT+1;
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Gauss-Seidel Solver Diverged\n');
    fprintf('BiConjugate Gradient Iteration=%i\tResidual=%2.6f\n',IT-1,res);
    x=X;
else
    fprintf('\nBiConjugate Gradient Solver Converged');
    fprintf('\nBiConjugate Gradient Iteration=%i\tResidual=%2.6e\n',IT-1,res);
end
x=X;
end
