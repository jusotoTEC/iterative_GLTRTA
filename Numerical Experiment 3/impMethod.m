function Z=glrma(A,C,r)
    Cp=pinv(C);
    T=A*Cp*C;
    [U,S,V]=svd(T);
    Z=U(:,1:r)*S(1:r,1:r)*(V(:,1:r)).'*Cp;
end