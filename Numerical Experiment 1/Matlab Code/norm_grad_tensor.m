function z = norm_grad_tensor(A,B,C,X,Y)
    
    aux=A-tprod(tprod(B,tprod(X,Y)),C);
    
    aux1=tprod(Y,C);
    dX=-2*tprod(tprod(tTranspose(B),aux),tTranspose(aux1));
        
    aux2=tprod(B,X);
    dY=-2*tprod(tprod(tTranspose(aux2),aux),tTranspose(C));
    
    z=sqrt(tFrobNorm(dX)^2+tFrobNorm(dY)^2);
    
    
end

