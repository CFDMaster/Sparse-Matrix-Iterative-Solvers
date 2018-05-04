function [x] = SOR4Tri(A,X,b,Omega,MaxIT,tol)
%SOR Iteration For Solution of a linear system Ax=b
%Optimized for Tri-Diagonal System, Exploiting Sparsity of the Matrix
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X, Omega is over relaxation Coe.
%Revison Modifications 1:Making Faster by removing Slow SUM Function
%Over 7 times Faster over Scaler Version 4 Simple FD Test
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%   A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
%   b: Right Hand Side Matix
%   X: Unkown Matrix-Also contains the initial Guess
%   tol: Tolerence for Residual Convergense
%   MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
%   omega: Over Relaxation Coe
res=norm(A*X-b); %Compute Residual
IT=1;
n=size(A,1);
display('Please wait..... Calculating the solution')
while(res>=tol) && (IT<=MaxIT)
    %Calcualte First & Last Element
    X(1)=(1-Omega)*X(1)+Omega*(b(1)-( A(1,2)*X(2) ))/A(1,1);
    X(n)=(1-Omega)*X(n)+Omega*(b(n)-( A(n,n-1)*X(n-1) ))/A(n,n);
    for i=2:n-1
        %switch  i%-------Replaced in R1
        %      case 1
        %        jj=[i+1];
        %      case n
        %        jj=[i-1];
        %      otherwise
        %        jj=[i-1 i+1];
        %end
        % X(i)=(1-Omega)*Xold(i)+Omega*(b(i)-sum(A(i,jj)*X(jj)))/A(i,i);--- Replaced in R1             
        X(i)=(1-Omega)*X(i)+Omega*(b(i)-( A(i,i-1)*X(i-1)+A(i,i+1)*X(i+1)) )/A(i,i);
    end
    res=norm(A*X-b); %Compute Residual
    IT=IT+1;
    %fprintf('SOR Iteration=%i\tResidual=%2.6e\n',IT,res);
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. SOR Solver Diverged\n');
    fprintf('SOR4Tri Iteration=%i\tResidual=%2.6f\n',IT,res);
    x=X;
    return
else
    fprintf('\nSOR4Tri Solver Converged');
    fprintf('\nSOR4Tri Iteration=%i\tResidual=%2.6e\n',IT,res);
end
x=X;
end

