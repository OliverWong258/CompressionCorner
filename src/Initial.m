function U = Initial(X,Y)
computational_X = 0.5*(X(1:end-1,1) + X(2:end,1));
computational_Y = 0.5*(Y(1:end-1,1) + Y(2:end,1));
p = 99719;
T = 293.15;
uu = 686.47;
R = 287.14;
gamma = 1.4;
Nx = size(computational_X,1);
Ny = size(computational_Y,1);
rho = zeros(Nx,Ny);
u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
P = zeros(Nx,Ny);
rho(1:Nx,1:Ny) = p/(R*T);
u(1:Nx,1:Ny) = uu;
v(1:Nx,1:Ny) = 0;
P(1:Nx,1:Ny) = p;
% 初始条件全场统一
U = zeros(Nx,Ny,4);
U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) = P./(gamma-1) + 0.5.*rho.*(u.^2+v.^2);

end

