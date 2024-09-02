function Xp=tpinv(X)
% Xp=tpinv(X) computes the tensor pseudoinverse of third-order tensor X
%
% Input:
%       X       -   m*n*p tensor
% Ouput:
%       Y       -   n*m*p tensor Y=tensor pseudoinverse of X
%
% References:
% Jin, H., Bai, M., Ben?tez, J., & Liu, X. (2017). 
% The generalized inverses of tensors and an application to linear models. 
% Computers & Mathematics with Applications, 74(3), 385-397.
%
% Written by Pablo Soto-Quiros (jusoto@tec.ac.cr)

[m,n,p]=size(X);
A=fft(X,[],3);
B=zeros(n,m,p);

for j=1:p
    B(:,:,j)=fastPinv(A(:,:,j));        
end

Xp=real(ifft(B,[],3));

end