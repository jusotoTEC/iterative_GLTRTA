% Numerical Experiment 1

% Reference:
%   Paper   = "A modified proximal point algorithm for solving a 
%              generalized low-tubal-rank tensor approximation problem 
%              based on the t-product"
%   Author = Soto-Quiros, Pablo (jusoto@tec.ac.cr)

clc; clear; close all

%Define tensors A, B, and C
p=5; q=4; m=3; n=4;  s=2; r=2; 
A=zeros(5,4,2);
A(:,:,1)=[-1  -2   1  -11; -38 -27  17  34; 43  24 -13 -35; -105 -47  18  29; 76  33 -13 -38];
A(:,:,2)=[-16  7  -11  -8; -43 -24  13  35; 38  27 -17 -34; -125 -35   2  33; 71  36 -17 -37];
B=zeros(5,3,2);
B(:,:,1)=[2   1   0; 3   2  -2; 0  -1   2; 2   3  -1; 0  -1   3];
B(:,:,2)=[-2  -2  -2; -2   1   3; 0  -2  -2; 3  -1  -2; 1  -2   1];
C=zeros(4,4,2);
C(:,:,1)=[1   1  -2   0; 3   1  -2  -1; 3  -1   3   1; 0   2   1  -1];
C(:,:,2)=[1   2   3  -1; 0   2   2  -2; 1   2   0   2; 3  -1  -2   0];

%Define initial parameters in the PGS method in Algorithm 3
iterMax=1000; tol=1e-10;
X0=zeros(3,2,2);
X0(:,:,1)=[-2   1;  0  -1; -1  -1]; 
X0(:,:,2)=[ -1  -2; -1  -1;  0   2];
Y0=zeros(2,4,2);
Y0(:,:,1)=[ 1   1   0  -2; -1   0   1   1]; 
Y0(:,:,2)=[-2  -2   2   0;  0   2  -2   2];

%Compute tensors Xk and Yk using the PGS method in Algorithm 3
[Xk,Yk,k,erk_v]=pgsMethod(A,B,C,X0,Y0,tol,iterMax);
semilogy(1:k,erk_v,'-o','LineWidth',1.25)
xlabel('Iterations (k)')
ylabel('Error')
grid  on

%Results:
%Error
disp(['Number of iterations = ', num2str(k)])
disp(['Error of the method = ', num2str(erk_v(end))])
disp(['Tensor Xk'])
Xk
disp(['Tensor Yk'])
Yk
disp(['Tensor Zk=Xk*Yk'])
Zk=tprod(Xk,Yk)
D=tprod(tprod(B,Zk),C);
obj_fun=tFrobNorm(A-D);
disp(['Value in objective function ||A-B*Zk*C||_{fr}^2 = ', num2str(obj_fun)])



