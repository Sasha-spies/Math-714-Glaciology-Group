function [u_face,H_face] = alex_staggered_u(x_face,H,N,b,A,mu)
% Constants
n = 3;
rho = 910;
g = 9.81;
u_0 = 1/(365.25*24*3600);  
eps_0 = 1e-12;% Strain-rate regularization [1/s]
%staggered grid!!!!!!!!!!!!!!!!!!!!
x_face = x_face(:);
H = H(:); N = N(:); b = b(:);
m = length(H);
dx = x_face(2)-x_face(1);
% Ice thickness @ velocity faces 
H_face = [H(1); (H(1:end-1)+H(2:end))/2; 0];
% Surface elevation @ cell centers
s = H+b;
% Initial face velocities
u_face = zeros(m+1,1);
for k = 1:400
    % longitudinal coefficient at cell centers; remark -- this uses old u
    u_x_center = diff(u_face)/dx;
    K = 2*A^(-1/n).*H.*(u_x_center.^2+eps_0^2).^((1/n-1)/2);
    % basal drag coefficient at ~~interior~~ velocity faces
    N_face = (N(1:end-1)+N(2:end))/2;
    C = mu.*N_face./sqrt(u_face(2:end-1).^2+u_0^2);
    % new face velocity linear system
    M = sparse(m+1,m+1);
    rhs = zeros(m+1,1);
    M(1,1) = 1;           % u(0)=0
    M(end,end) = 1;       % u(xL)=0
    % Face j is between centers j-1 and j     also j=2,...,m
    for j=2:m
        left=K(j-1)/dx^2;
        right=K(j)/dx^2;
        M(j,j-1) = -left;
        M(j,j)=left+right+C(j-1);
        M(j,j+1) = -right;
        rhs(j)=-rho*g*H_face(j)*(s(j)-s(j-1))/dx;
    end
    u_trial = M\rhs;
    u_new = 0.5*u_face + 0.5*u_trial; %relaxation
    if norm(u_new-u_face,inf) < 1e-6*max(norm(u_new,inf),u_0)
        u_face = u_new;
        return
    end
    u_face = u_new;
end
warning('picard not converging in these number of iterations.');
end
