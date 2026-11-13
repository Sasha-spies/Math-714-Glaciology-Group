%% Mass Balance Equation, Final Math 714 Project

%elliptic equation, dH/dt + d(uH)/dx = mass balance (adot)

% k = 0.5; %so half a year
% t = 200; %in years
% nsteps = t/k;


%% Parameters (constants)
rho_i = 917;                % density of ice (kg/m^3)
n = 3;                      % Glen's law stress exponent
g = 9.81;                   % gravitational acceleration (m/s^2)
A = 3e-24;                  % flow rate parameter Pa^-3 s^-1
beta = 0.007;               % mass balance gradient yr^-1
ELA  = 1300;                % ELA in meters
beta_s = beta / 3.1536e7;   % beta in seconds
adot = beta_s * (s - ELA);

%% Set up for grid and computation 
% x-vector
dx = 100;              % grid spacing in meters 
x = 0:dx:20000;        % grid points from 0 - 20 km 
nx = length(x);        % number of grid points 

% Bed elevation computation
z_0  = 1600;           % assumed in meters asl 
alpha = 0.05;          % slope angle 
z_b = z_0 - alpha*x;   % given equation

% Initialiation of surface and thickness
H = zeros(1,nx);      %no ice initially 
s = z_b + H;          %surface here is equal to bed 

u =                 ; %placeholder for Alex's function, need an actual function 

%In time
% dsdt = zeros(1,nx);
% dsdt(2:end-1) = bdot + (q(2:end-1) - q(1:end-1) /dx);

dHdt = dHdt_spatial(H,u,adot,dx);

%% ODE Function
function dHdt = dHdt_spatial(H,u,adot,dx)
    F = u .* H; 

    % Initialize derivative vector for below 
    dFdx = zeros(size(H));
    %centered difference scheme
    dFdx(2:end-1) = (F(3:end) - F(1:end-2)) / (2*dx); 
    dFdx_stag=  0.5* dFdx(2:end-1); 
    %boundary conditions (dirichlet)
    dFdx_stag(1) = 0;
    dFdx_stag(end) = 0;

    dHdt = adot - dFdx_stag; 

end 



%plotting?

% separate questions
% ask marianne about the definition of flux (why not density x u -- would this be for 2/3 dimensions?) 
% ask marianne about the staggered grid and what scheme I need to receive
% why are we incorporating longitudinal stress -- it's generally negligible
% for steady-state and thus can be omitted 
