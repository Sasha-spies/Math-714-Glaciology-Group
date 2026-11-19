clear all

indices=1; 
print_all=1;
plot_all=1;
istart=0;

addpath ./
addpath('/Users/kross29/Downloads/model_petsc')
set(0,'DefaultFigureWindowStyle','docked')

%% Define variables

denw  = 1000;      % kg/m3, density of water
Lh    = 330000;    % J/kg, latent heat
h0    = 1;         % m, drained hydraulically active sediment thickness
qgeo  = 0.065;     % W/m2, geothermal heat flux
Tmelt = 273;       % K
Tsurf = -20 + Tmelt;
delT  = -20;       % K, difference between surface and melt
mu    = 0.5;       % sliding param
k     = 2;         % W/mK, thermal conductivity
Ce    = 0.0345;    % coeff of compressibility

er = 0.78;         % void reference parameter
Nr = 1000;         % Pa, reference value for effective pressure

dx = 1;
x  = 0:dx:1000;    % inclusive of 1000, like np.arange(0, 1000+dx, dx)
disp(numel(x))

Hi = max(2000 - 0.006*x, 100);   % m, from Haseloff 2025
ui = 1*x;
Ni = 20.85 - 1e-6*x;             % MPa

hsi = h0;
u0  = 1;                         % m/yr, for calculating taub
taub = mu*Ni.*ui ./ sqrt(ui.^2 + u0^2);   % basal shear resistance

Nfrome = @(ei) Nr .* exp((er - ei) ./ Ce);

efromN = @(Ni) er - Ce .* log(Ni ./ Nr);

dedt = @(Hi, ui, ei) (1 ./ (hsi*denw*Lh)) .* ...
    ( (Nfrome(ei) .* mu .* ui ./ sqrt(ui.^2 + u0^2)) .* ui ...
    + qgeo - k*delT ./ Hi );

%% RK4

dt      = 0.01;
tottime = 100;
time    = 0:(tottime/dt - 1);       

ei   = efromN(Ni);                 
x    = 0:numel(ei)-1;              
eold = ei;

earr_rk4        = zeros(numel(time), numel(ei));
earr_rk4(1,:)   = ei;

figure; hold on;
colors = jet(numel(time));

for i = 1:numel(time)-1
    t  = time(i);
    k1 = dedt(Hi, ui, eold);
    k2 = dedt(Hi, ui, eold + dt/2*k1);
    k3 = dedt(Hi, ui, eold + dt/2*k2);
    k4 = dedt(Hi, ui, eold + dt*k3);
    enew = eold + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    eold = enew;
    earr_rk4(i+1,:) = eold;

    cidx = int32(t) + 1;            

    if i == 1
        plot(x, eold, 'Color', colors(cidx,:), 'DisplayName', 'initial');
    elseif i == numel(time)-1
        plot(x, eold, 'Color', colors(cidx,:), 'DisplayName', 'final');
    elseif mod(i-1, 1000) == 0
        plot(x, eold, 'Color', colors(cidx,:), 'HandleVisibility', 'off');
        disp(i-1)
    end
end

legend show;
xlabel('Distance (x)');
ylabel('Till Void Ratio (e)');
title("Till Void Ratio Using RK4")


%% 3 step multistep using RK4

dt      = 0.01;
tottime = 100;
time    = 0:(tottime/dt - 1);       % 0, 1, ..., 9999

ei   = efromN(Ni);                  % must come BEFORE x is defined from numel(ei)
x    = 0:numel(ei)-1;
eold = ei;

earr_ab3      = zeros(numel(time), numel(ei));
earr_ab3(1,:) = ei;

figure; hold on;
colors = jet(numel(time));

for i = 1:numel(time)-1
    t = time(i);

    if i <= 3   % RK4 startup for the first three steps (Python i = 0, 1, 2)
        k1 = dedt(Hi, ui, eold);
        k2 = dedt(Hi, ui, eold + dt/2*k1);
        k3 = dedt(Hi, ui, eold + dt/2*k2);
        k4 = dedt(Hi, ui, eold + dt*k3);
        enew = eold + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    else        % Adams-Bashforth 3rd order
        enew = eold ...
            + 23/12*dt*dedt(Hi, ui, eold) ...
            - 4/3*dt*dedt(Hi, ui, earr_ab3(i-1,:)) ...
            + 5/12*dt*dedt(Hi, ui, earr_ab3(i-2,:));
    end

    eold = enew;
    earr_ab3(i+1,:) = eold;

    cidx = int32(t) + 1;

    if i == 1
        plot(x, eold, 'Color', colors(cidx,:), 'DisplayName', 'initial');
    elseif i == numel(time)-1
        plot(x, eold, 'Color', colors(cidx,:), 'DisplayName', 'final');
    elseif mod(i-1, 1000) == 0
        plot(x, eold, 'Color', colors(cidx,:), 'HandleVisibility', 'off');
        disp(i-1)
    end
end

legend show;
xlabel('Distance (x)');
ylabel('Till Void Ratio (e)');
title("Till Void Ratio Using Three-Step Adams-Bashforth")

%% Define functions


