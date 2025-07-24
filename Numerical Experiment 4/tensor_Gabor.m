function [A,C]=tensor_Gabor(a,b,h1,h2)
    x = a:h1:b;
    y = a:h1:b;
    z = 0:h2:pi;
    s=length(z);
    A=zeros(length(x),length(y),s);
    for k=1:s
        A(:,:,k) = generar_matriz(x, y, z(k));
    end
    C=A+0.05*randn(size(A));      
end



function H = generar_matriz(x,y,p)
    % Parámetros
    f = 0.1; 
    [X, Y] = meshgrid(x, y);

    % Evaluación de la función h(x,y)
    H = exp(-(X.^2 + Y.^2)/2) .* cos(2*pi*f*X + p);

    %H = sin(2*pi*f1*X + p1)+ cos(2*pi*f1*X + p1);
end