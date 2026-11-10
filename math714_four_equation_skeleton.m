%% Math 714 final project skeleton -- Alex notes 
% -------------------------------------------------------------------------
clear; clc; close all;
%% We need to decide a domain and grid (group meeting November 5)
% For now, we should have every equation on identacle 1-D x grid.
% we need to discuss if any of our methods utilize a staggered grid or
% spectral discretization, because then it might change how the data is
% seen on shared grid
x_l  = 0;           % this is the left end of domain [m]
x_r = 10000;       % same  but right end [m]
num_space_grid     = 101;         % number of spatial grid poiN_ts
x  = linspace(x_l, x_r, num_space_grid)';   % Column vector, size num_space_grid x 1
dx   = x(2) - x(1);        % Spatial grid spacing [m]

%% Time stepping (We should also have this decided during our november 5th meeting)
t_i = 0;  % Initial time [s]
t_f   = 24*3600;      % end time [s]
dt     = 3600;   % time step [s]
N_t     = round((t_f-t_i)/dt) + 1;
t      = t_i + (0:N_t-1)*dt;

%% Physical Paramters
p.g       = 9.81;     % Gravitational acceleration [m/s^2]
p.rho_i   = ???;      % Ice density [kg/m^3]
p.rho_w   = ???;      % Water density [kg/m^3]
p.A       = ???;      % Glen flow-rate factor [units depend on n]
p.n       = ???;      % Glen exponeN_t [-]
%  basal water/till parameters
p.s0      = ???;      % CoefficieN_t s0 from H*s0*(de/dt); verify units
p.e_ref   = ???;      % Reference till void ratio [-]
p.C       = ???;      % Void-ratio/log-pressure coefficieN_t [-]
p.N_ref   = ???;      % Reference effective pressure [Pa]
p.Lh      = ???;      % LateN_t heat [J/kg]
p.q_geo   = ???;      % Geothermal heat flux [W/m^2]
p.k       = ???;      % Thermal conductivity [W/(m K)]
p.deltaT  = ???;      % Temperature difference [K]
%GIA carson stuff 
p.eta     = ???;      % Viscosity parameter [Pa s] (verify model definition)
p.rho_r   = ???;      % Rock/maN_tle density [kg/m^3]
p.D       = ???;      % Flexural rigidity [N m] (verify dimensions/model)

% FORCING / CONSTITUTIVE RELATIONS -- Need to decided these tg 

%% BOUNDARY CONDITIONS!!!!!!!!!! (this needs to be decided before ppl work on their solvers)
bc.H.left    = struct('A','B','C',???); % thickness x=x_l
bc.H.right   = struct('A','B','C',???); % thickness x=x_r
bc.u.left    = struct('A','B','C',???); % velocity x=x_l
bc.u.right   = struct('A','B','C',???); % velocity x=x_r
bc.water.left  = struct('A','B','C',???); % e, N, or q BC?
bc.water.right = struct('A','B','C',???); % decide variable
bc.bed.left  = struct('A','B','C',???); % bed displacement
bc.bed.right = struct('A','B','C',???); % bed displacement
%% INITIAL CONDITIONS!!!!!!!!!!! This is another thing we have to decide as a group
H  = ???(num_space_grid,N_t); % Ice thickness H(x,t) [m]
u  = ???(num_space_grid,N_t);  % Ice velocity u(x,t) [m/s]
N  = ???(num_space_grid,N_t);% Effective pressure N(x,t)[Pa]
e  = ???(num_space_grid,N_t);  % Till void ratio e(x,t)[-]
db = ???(num_space_grid,N_t); % Change in bed elevation Delta b(x,t)[m]
q  = ???(num_space_grid,N_t);% water flux? does kit want this? 

%% Main time lop should go here -- this is where the coupling happens 
%------------------------------------------------------------------------
% Known at time j: H, u, N, e, db
%  (i)  surface el.   s = b0 + db + H
%  (ii) Solve vel from H, N, and surface slope
%  (iii) Update thickness from mass conservation and u
%  (iv) Update void ratio/effective pressure from water balance,
%        melting and the void-ratio/pressure relationship
%  (v)  Update bed displacement from carsons GIA equation
%  (vi) store vals 
p.maxCouplingIters = 1; % need to choose iterations amd convergence tol

for j = 1:N_t-1
    % Extract fields at the old time t(j).
    H_old  = H(:,j);
    u_old  = u(:,j);
    N_old  = N(:,j);
    e_old  = e(:,j);
    db_old = db(:,j);
    % Initial guesses for the new time t(j+1).
    H_guess  = H_old;
    u_guess  = u_old;
    N_guess  = N_old;
    e_guess  = e_old;
    db_guess = db_old;
    % inner loop for iterating coupled solvers maybe? ~~~~~~~~
    for iter = 1:p.maxCouplingIters
     
    end
    % Save the new stuff after the coupling iteration
    H(:,j+1)  = H_new;
    u(:,j+1)  = u_new;
    N(:,j+1)  = N_new;
    e(:,j+1)  = e_new;
    db(:,j+1) = db_new;
    q(:,j+1)  = q_new;

end

%% Plotting Stuff Will go here _________________________________
%%------------------------------------------------------------------------
% I hate plotting were going to have to ask Yue and Marianne some questions


%%LOCAL FUNCTIONSSSSSSSSSSSSSS
%% ------------------------------------------------------------------------

function u_new = vel_solve(H,N,s,u_guess,dx,p,bc_u)
% Alex's Eq: Velocity and Stress Balance 
%   d/dx[ 2 A^(-1/n) H |du/dx|^(1/n-1) (du/dx) ]
%      - tau_b(N,u) = rho_i*g*H*(ds/dx)
end

function H_new = mass_bal(H_old,u_new,x,time,dt,dx,p,bc_H)
% Jaela's Eq: Mass Balance!! 
%   dH/dt + d(u*H)/dx = a_dot(x,t)
end

function [e_new,N_new,q_new] = water_bal(e_old,N_old,H,u,x,time,dt,dx,p,bc_water)
%Kit's Eq: Watwer Balance and Eff Press
%   H*s0*de/dt + d(q)/dx = m_dot
%   m_dot = [tau_b(N,u)*u + q_geo - k*DeltaT/H]/(rho_w*Lh)
%   e = e_ref - C*log(N/N_ref)
end

function db_new = gia(db_old,H,u,x,time,dt,dx,p,bc_bed)
% Carson's eq: Glacial isostatic adjustmeN_t 
%   2*eta*|nabla|*d(db)/dt + rho_r*g*db + D*nabla^4(u) = sigma_zz
end
