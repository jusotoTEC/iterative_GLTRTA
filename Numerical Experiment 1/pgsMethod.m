function [Xk,Yk,k,erk_v]=pgsMethod(A,B,C,X0,Y0,tol,iterMax)

    r=size(X0,2); s=size(X0,3);
    Xk=X0; Yk=Y0; q=0.5;           
    Bp=tpinv(B); Cp=tpinv(C); Ir=teye(r,s);     
    erk_v=zeros(1,iterMax);
    
    for k=1:iterMax
        %Compute X^(k+1)
        alphak=q^k;
        At=[A tprod(alphak^2*B,Xk)];
        Ctp=tpinv([tprod(Yk,C) alphak^2*Ir]);
        Xk=tprod(Bp,tprod(At,Ctp));  

        %Compute Y^(k+1)
        betak=q^k;
        At=[A; tprod(betak^2*Yk,C)];
        Btp=tpinv([tprod(B,Xk); betak^2*Ir]);
        Yk=tprod(Btp,tprod(At,Cp));

        %Compute error
        erk = norm_grad_tensor(A,B,C,Xk,Yk);
        erk_v(k)=erk;        
        if erk<tol
            break
        end
    end
    erk_v=erk_v(1:k);
end