% Numerical Experiment 3

% Reference:
%   Paper   = "A modified proximal point algorithm for solving a 
%              generalized low-tubal-rank tensor approximation problem 
%              based on the t-product"
%   Author = Soto-Quiros, Pablo (jusoto@tec.ac.cr)

clc; clear; close all

%Source Image
A=im2double(imread('Nico.jpg'));
figure
imshow(A)
%Noisy Image
C=im2double(imread('Nico_Noisy.jpg'));
figure
imshow(C)

[m,q,s]=size(A); p=size(C,1); n=m;
B=teye(m,s);

r=64; 
cr=(m+n)*r/(m*n);
disp(['Using r = ', num2str(r), ' and a compression ratio = ', num2str(cr)])


%Proposed Method
iterMax=10000; tol=1e-1;
[X0,Y0] = initial_values(A,B,C,r);
[Xk,Yk,~,~]=pgsMethod(A,B,C,X0,Y0,tol,iterMax);
At_pgs=tprod(tprod(Xk,Yk),C);
figure
imshow(At_pgs)


ssim_R_pgs=ssim(A(:,:,1),At_pgs(:,:,1));
ssim_G_pgs=ssim(A(:,:,2),At_pgs(:,:,2));
ssim_B_pgs=ssim(A(:,:,3),At_pgs(:,:,3));
ssimColor_pgs=(ssim_R_pgs+ssim_G_pgs+ssim_B_pgs)/3;
disp(['SSIM using PGS method = ', num2str(ssimColor_pgs)])


%Method IMP
At_imp=zeros(size(A));
%Red Channel
A1=A(:,:,1); C1=C(:,:,1);
Z1=impMethod(A1,C1,r);
At_imp(:,:,1)=Z1*C1;
ssim_R_imp=ssim(A(:,:,1),At_imp(:,:,1));
%Red Green
A2=A(:,:,2); C2=C(:,:,2);
Z2=impMethod(A2,C2,r);
At_imp(:,:,2)=Z2*C2;
ssim_G_imp=ssim(A(:,:,2),At_imp(:,:,2));
%Red Blue
A3=A(:,:,3); C3=C(:,:,3);
Z3=impMethod(A3,C3,r);
At_imp(:,:,3)=Z3*C3;
ssim_B_imp=ssim(A(:,:,3),At_imp(:,:,3));
figure
imshow(At_imp)

ssimColor_imp=(ssim_R_imp+ssim_G_imp+ssim_B_imp)/3;
disp(['SSIM using IMP method = ', num2str(ssimColor_imp)])