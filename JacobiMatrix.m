function [x,resi] = JacobiMatrix( A,X,b,MaxIT,tol )
%Jacobi Matrix Iteration For Solution of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X, resi is residial vector for plotting
%Revison Modifications 1:Making Faster by removing Slow Inverse Rearranging
%WOW!!Over 40 times Faster over Privous Matrix & Scaler Version 4 Simple FD Test
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
D=diag(diag(A));
% L=-tril(A,-1);---------Removed in R1
% U=-triu(A,1);---------Removed in R1
% res=norm(A*X-b); %Compute Residual---------Replaced in R1
r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
% if nargout==2, res1=zeros(1,20000); end % If output >1,set residial vector
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    X=X+D\(r);
    %X=inv(D)*(L+U)*X+inv(D)*b; %Jacobi Iteration
    %  res=norm(A*X-b);---------Replaced in R1
    r=b-A*X;
    res=norm(r); %Compute Residual
    if nargout==2, res1(IT)=res; end % If output >1,set residial vector
    IT=IT+1;
    %     fprintf('Jacobi Matrix Iteration=%i\tResidual=%2.6f\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. Jacobi Solver Diverged\n');
    fprintf('Jacobi Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
    x=X;
    if nargout==2
        resi=res1;
    end
else
    fprintf('\nJacobi Matrix Solver Converged');
    fprintf('\nJacobi Matrix Iteration=%i\tResidual=%2.6e\n',IT,res);
    if nargout==2
        resi=res1;
    end
end
x=X;
if nargout==2
    resi=res1;
end
end

