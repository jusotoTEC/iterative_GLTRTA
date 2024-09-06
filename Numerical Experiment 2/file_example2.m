% Numerical Experiment 2

% Reference:
%   Paper   = A modified proximal point algorithm for solving a generalized 
%             low-tubal-rank tensor approximation problem based on the t-product
%   Author = Soto-Quiros, Pablo

clc; clear; close all

%Dimension of tensors and initial parameters
iterMax=1000; tol=1e-5;
p=8; q=7; m=6; n=4; s=2; r=2; 
numSimulations=5000;
vect_k=zeros(1,numSimulations);
vect_er=zeros(1,numSimulations);

for i=1:numSimulations
    disp(['No. Simulation = ' num2str(i), ' [Process completed = ', num2str(i/numSimulations*100), '%]'])
    %Define tensors A, B, and C
    A=rand(p,q,s); B=rand(p,m,s); C=rand(n,q,s);   
    %Using initial values given by Section 3.3.1
    [X0,Y0] = initial_values(A,B,C,r);
    %PGS method
    [Xk,Yk,k,erk_v]=pgsMethod(A,B,C,X0,Y0,tol,iterMax);
    vect_k(i)=k; vect_er(i)=erk_v(end);
    clc
end

mean_k=mean(vect_k);
disp(['Mean of iterations = ', num2str(mean_k)])
std_k=std(vect_k);
disp(['STD of iterations = ', num2str(std_k)])
mean_er=mean(vect_er);
disp(['Mean of errors = ', num2str(mean_er)])
std_er=std(vect_er);
disp(['STD of errors = ', num2str(std_er)])

%Historgram of frequency of number of iterations
histogram(vect_k,max(vect_k)-min(vect_k)+1);
xlabel('Iterations (k)')
ylabel('Frecuency')
set(gcf, 'Position', [680 558 800 200])



