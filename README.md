# Sparse-Matrix-Iterative-Solvers
This Repo Contains Series of Functions Coded in Matlab for System of  Spare Linear  Equations Ax=B solution, Specially those arise in Computational Fluid Mechanics(CFD)
Description:
Let Ax=B be a system of algeraic equation in which A is Left Hand Side or Coefficient Matrix, B is right hand side or knows matrix and x is the unkown matrix. For most methods to work properly A must be a Sparse Matrix. Although some methods like Bi-Conjugate Gradient work well with dense Marices too.
Parameters:
  A : Left Hand Side(Coefficient) Matrix(Usually Sparse)
  b: Right Hand Side Matix
  X: Unkown Matrix-Also contains the initial Guess
  tol: Tolerence for Residual Convergense
  MaxIT: Maximum allowed Iteration-If soluiton does not converges after MaxIT iterations, it considered diverged
Function List:
  BCG(A,X,b,MaxIT,tol ): BiConjugate Gradient Method
  CG(A,X,b,MaxIT,tol ):  Conjugate Gradient Method-Matrix Form
  FindOptimalOmega(A): Function to Finds the Optimal relaxation coe. for SOR method in the range 0< Omga <2
  GaussSeidel(A,X,b,MaxIT,tol ):  Gauss-Seidel Iteration-Scaler Form
  GaussSeidelMatrix(A,X,b,MaxIT,tol ):  Gauss-Seidel Iteration-Matrix Form 
  GaussSeidelRedBlack(A,X,b,MaxIT,tol ):  Gauss-Seidel Iteration-Red-Black Decomposition-1D Case
  Jacobi(A,X,b,MaxIT,tol ): Jacobi Iteration-Scaler Form
  JacobiMatrix(A,X,b,MaxIT,tol ): Jacobi Iteration-Matrix Form
  PCG(A,X,b,B,MaxIT,tol ): Perconditioned Conjugate Gradient Method-B is the Perconditioning Matrix
  SOR(A,X,b,Omega,MaxIT,tol) :  SOR(Successive Over Relaxation) Iteration, Scaler Form-Omega: Over Relaxation Coe.
  SOR4Tri(A,X,b,Omega,MaxIT,tol ) : SOR(Successive Over Relaxation) Iteration Optimized for Tri-Diagonal System, Exploiting Sparsity of       the Matrix--Omega: Over Relaxation Coe.
  SORMatrix( A,X,b,Omega,MaxIT,tol) : SOR(Successive Over Relaxation) Iteration, Matrix Form-Omega: Over Relaxation Coe.
  SymmetricGaussSeidel(A,X,b,MaxIT,tol) : Symmetric Gauss-Seidel Iteration, Scaler Form- Gauss-Seidel Optimized for symmetric matrices
  SymmetricGaussSeidelMatrix(A,X,b,MaxIT,tol) : Symmetric Gauss-Seidel Iteration, Matrix Form- Gauss-Seidel Optimized for symmetric           matrices
  SymmetricSOR(A,X,b,Omega,MaxIT,tol):  Symmetric SOR Iteration, Explointing advantage of symmetric Coe. Marix A-Scaler Form-Omega: Over     Relaxation Coe
  SymmetricSORMatrix(A,X,b,Omega,MaxIT,tol):  Symmetric SOR Iteration, Explointing advantage of symmetric Coe. Marix A-Scaler Form-         Omega: Over Relaxation Coe
  TDMA(A,d) : Solves TriDigonal matrix using Thomas TDMA Alghorithm, A: Coe. Matrix, d=RHS of the solutiion Returns X=Solution(Works         only for Tridiagonal matrices likes 2d convection-diffusion discritizations by FVM)-- A: Coe. Matrix, d=RHS of the solutiion Returns     X=Solution
  TDMASimple(A,d): Another version of TDMA Solver- Same as above
  
  Usage :
  Copy Corresponding m File to the working directory or add the directory to the Matlab path and use in the way as Matlab Build-in Functions.
  Example:
  A=...  % Assemble Coe. matrix-Sparse nxm matix
  B=...% Assembel RHS matrix 1xm matrix
  X=...% Set intial guess 1xn
  tol=1e-5; % Set tolerance for system of algerbaric equations solution
  MaxIT=1e3; % Set Maximum allowed iterations number
  x = GaussSeidelMatrix(A,X,b,MaxIT,tol); % Solve equatations and store result in x
  % Do somthing with result(x) lik plot or etc
  
  

  
  
  
  
  
  
  
  

