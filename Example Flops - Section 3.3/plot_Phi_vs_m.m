function plot_Phi_vs_m()   


    % Numerical example in Section 3.3
    
    % Reference:
    %   Paper   = "A proximal Gauss-Seidel algorithm for solving a 
    %              generalized low-tubal-rank tensor approximation problem 
    %              based on the t-product"
    %   Author = Soto-Quiros, Pablo (jusoto@tec.ac.cr)

    % Dimensions
    n = 50;
    p = 30;
    q = 40;
    r = 5;
    s = 10;
    k = 15;
    
    m_values = 10:10:1000;
        
    Phi_values = zeros(size(m_values));
       
    for i = 1:length(m_values)
        m = m_values(i);
        Phi_values(i) = calculate_Phi(m, n, p, q, r, s, k);
    end
        
    figure;
    plot(m_values, Phi_values, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Dimension (m)', 'FontSize', 12);
    ylabel('Flops (\Phi)', 'FontSize', 12);

end

function Phi = calculate_Phi(m, n, p, q, r, s, k)
    % Compute Φ1
    term1 = m.^2 * p;
    term2 = n^2 * q;
    term3 = 3 * r * (m*p + n*q);
    term4 = (p + 1) * r^3;
    term5 = (q + m + r^2) * r^2;
    term6 = (m + r) .* (p + r) * q;
    
    Phi1 = s * (term1 + term2 + term3 + term4 + term5 + term6);
    
    % Compute Φ2
    termA = m * (2*p + r + 3*s);
    termB = n * (3*q + r);
    termC = p * (2*r + s);
    termD = q * (2*r + s);
    termE = r * (2*r + 3*s);
    
    Phi2 = s * (termA + termB + termC + termD + termE);
    
    % Compute Φ 
    Phi = k * (Phi1 + Phi2 * log(s));
end