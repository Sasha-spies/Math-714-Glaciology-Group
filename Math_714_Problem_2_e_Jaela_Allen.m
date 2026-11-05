%% Second-Order Accuracy Test: Solve the BVP with q(x) = sin(\pi x)


% Set grid resolution
m_values= [20, 40 , 80 , 160];
n = length(m_values);
E_values = zeros(n,1);
h_values= zeros(n,1);
 

for idx = 1:n
    m = m_values(idx);
    h = 1/m;
    x = linspace(0, 1, m+1)';

    % Generate the matrix A and function F
    F=zeros(m,1);
    A = zeros(m,m);

    %Row 1: u'(0) = 0 ; 4U_1 - U_2 = 0
    A(1,1) = 4;
    A(1,2) = -1;
    F(1) = 0;

    %Row 2: u"(0) = 0; -4U_1+6U_2 - 4U_3 +U_4 = h^4 q(x_2)
    A(2,1) = -4;
    A(2,2) = 6;
    A(2,3) = -4;
    A(2,4) = 1;
    F(2)  = h^4 * sin(pi*x(3)); %x_3 because of the DANG indexing so, x_2 = x_3


    % Interior rows (j = 3..m-2)
    for j = 3:m-2
        A(j,j-2) = 1;
        A(j,j-1) = -4;
        A(j,j) = 6;
        A(j,j+1) = -4;
        A(j,j+2) = 1;
        F(j) = h^4 * sin(pi*x(j+1)); %There MATLAB goes again where x_j+1 = x_j


    end

    %Row m-1: 2U_m - 5U_m-1 + 4U_m-2 -U_m_3 = 0
    A(m-1,m-3) = -1;
    A(m-1, m-2) = 4;
    A(m-1, m-1) = -5;
    A(m-1, m) = 2;
    F(m-1) = 0;

    %Row m: u'''(1) = 0 ; 5U_m - 18U_m-1 + 24U_m-2 - 14U_m-3 + 3U_m-4 = 0
    A(m, m-4) = 3;
    A(m, m-3) = -14;
    A(m, m-2)  = 24;
    A(m, m-1) = -18;
    A(m, m)  = 5;

    F(m) = 0;

    %Solve the linear system
    U_interior_rows = A\F;
    U = [0; U_interior_rows];


    %% Print the exact solutions, and compute the infinity norm
    % error
    uex= 1/(pi^4) * sin(pi * x) - x.^3/(6*pi) + x.^2/(2*pi) - x/(pi^3) ;

    %weight
    w  = ones(m+1,1);
    w(1) = 0.5;
    w(end) = 0.5;

    Em = sqrt(h * sum(w .* (U - uex).^2));

    h_values(idx) = h;
    E_values(idx) = Em;
    fprintf('m = %4d, h = %8.6, E_m = %.6e\n', m, h, Em);

end


%plotting
figure(1);
loglog(h_values, E_values, 'b-o', 'LineWidth', 1.5);
xlabel('h');
ylabel ('E_m');
title('Convergence of Euler Bernoulli');
grid on;