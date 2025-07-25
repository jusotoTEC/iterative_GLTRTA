function file_example3()

    % Numerical Experiment 4
    
    % Reference:
    %   Paper   = "A proximal Gauss-Seidel algorithm for solving a 
    %              generalized low-tubal-rank tensor approximation problem 
    %              based on the t-product"
    %   Author = Soto-Quiros, Pablo (jusoto@tec.ac.cr)

    clc; clear; close all


    % Numerical Solution
    a=-3; b=3;
    h1=0.1;
    h2=0.025;
    z = 0:h2:pi;

    [A,C]=tensor_Gabor(a,b,h1,h2);

    dimTensor=length(a:h1:b);
    s=size(A,3);
    B=teye(dimTensor,s);
    r=dimTensor;

    X0=rand(dimTensor,r,s);
    Y0=rand(r,dimTensor,s);

    tol=1e-5; iterMax=1000;
    [Xk,Yk,k,~]=pgsMethod(A,B,C,X0,Y0,tol,iterMax);
    display(['Number of iterations k=', num2str(k)])

    A1=tprod(tprod(Xk,Yk),C);

    display(['Error given by = ', num2str(tFrobNorm(A-A1))])

    % Video
    video = VideoWriter('enhanced_signal.mp4', 'MPEG-4');
    video.FrameRate = 10;
    open(video);
    fig = figure;

    for k=1:size(A,3)
        subplot(1,3,1)
        surf(a:h1:b,a:h1:b,A(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$f(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Source Signal $f$', 'Interpreter', 'latex')
        subplot(1,3,2)
        surf(a:h1:b,a:h1:b,C(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$h(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Noisy Signal $h$', 'Interpreter', 'latex')
        subplot(1,3,3)
        surf(a:h1:b,a:h1:b,A1(:,:,k))
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex')
        zlabel(['$\widehat{f}(x,y,', num2str(z(k)), ')$'], 'Interpreter', 'latex');
        title('Reconstructed Signal $\widehat{f}$', 'Interpreter', 'latex')
        pause(0.025)

        drawnow;  % Actualiza la figura
        frame = getframe(fig);  % Captura el frame
        writeVideo(video, frame);  % Escribe el frame en el video
    end
    close(video);

end
