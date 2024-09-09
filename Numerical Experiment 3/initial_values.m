function [X0,Y0] = initial_values(A,B,C,r)   
    Bp=tpinv(B); Cp=tpinv(C);
    Z=tprod(tprod(Bp,A),Cp);
    [U,S,V]=tsvd(Z);
    X0=tprod(U(:,1:r,:),S(1:r,1:r,:));
    Y0=tTranspose(V(:,1:r,:));        
end

function [U,S,V]=tsvd(X)
% Written by Pablo Soto-Quiros (jusoto@tec.ac.cr)
[m,n,p]=size(X);
A=fft(X,[],3);
U1=zeros(m,m,p); S1=zeros(m,n,p); V1=zeros(n,n,p);
for j=1:p
    [U1(:,:,j),S1(:,:,j), V1(:,:,j)] =svd(A(:,:,j));        
end
U=real(ifft(U1,[],3)); S=real(ifft(S1,[],3)); V=real(ifft(V1,[],3));
end

