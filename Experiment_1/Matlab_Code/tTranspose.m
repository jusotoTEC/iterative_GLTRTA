function Y=tTranspose(X)

    [m,n,s]=size(X);
    Y=zeros(n,m,s);
    Y(:,:,1)=(X(:,:,1))';    
    for k=2:s
        Y(:,:,k)=(X(:,:,s-k+2))';
    end
    
end

