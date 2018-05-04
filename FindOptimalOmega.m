function [OptimalOmega] = FindOptimalOmega(A)
%Finds the Optimal relaxation coe. for SOR method in the range 0< Omga <2
% For Symmetric Positive Matrix with Ro(A)<1
% A is the solution Matrix Ax=b
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
n=size(A,1);
D=diag(diag(A));
Bj=eye(n)-D\A; 
if ~isreal(eig(Bj))
    fprintf('Bad Shaped Matrix, Imaginary Eigenvalues');
    OptimalOmega=-1;
    return 
end

RoBj=max(abs(eig(Bj))); %Jacobian Specteral 
OptimalOmega=2/(1+sqrt(1-RoBj^2));

end

