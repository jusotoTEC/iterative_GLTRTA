function 



function H = generar_matriz_h()
    % Parámetros
    f = 0.1; %Entre 0.05 y 0.25
    p = 1; %Entre 0 y pi

    % Dominio
    x = -5:0.1:5;
    y = -5:0.1:5;
    [X, Y] = meshgrid(x, y);

    % Evaluación de la función h(x,y)
    H = exp(-(X.^2 + Y.^2)/2) .* cos(2*pi*f*X + p)+randn(length(x),length(y));
    surf(x, y, H)
end