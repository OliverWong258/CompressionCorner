function [F_,G_] = Flux_Splitting(U,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type)

Nx = size(U,1);
Ny = size(U,2);
U_p=U./J(4:Nx+3,4:Ny+3);
[Uxp,Uyp] = BoundaryCondition(U_p);

% x方向
Uhx = 0.5*(Uxp(3:Nx+3,:,:) + Uxp(4:Nx+4,:,:));
% 重构通量函数
Wx=Uxp.*J(:,4:3+Ny);
Fp = Flux_X(Uxp,gamma); 
Gp = Flux_Y(Uxp,gamma);
F=J(:,4:3+Ny).*(px1(:,4:3+Ny).*Fp+py1(:,4:3+Ny).*Gp);

F_pos = zeros(size(F));
F_neg = zeros(size(F));
switch splitting_type
    case 'LF'
        lamda_x = MaxEigX(Uhx,gamma);
        F_pos = 0.5*( F + lamda_x*Wx );    
        F_neg = F - F_pos;
    case 'SW'
        F_pos = FluxSW_X(Uxp,gamma);
        F_pos = 0.95.*J(:,4:3+Ny).*F_pos;
        F_neg = F - F_pos;
    case 'VL'
        F_pos = FluxSW_X(Uxp,gamma);
        F_pos = J(:,4:3+Ny).*F_pos;
        F_neg = F - F_pos;
end
% 重构五阶WENO所需的通量
Fp_2 = F_pos(1:Nx+1,:,:); Fp_1 = F_pos(2:Nx+2,:,:); Fp0 = F_pos(3:Nx+3,:,:); Fp1 = F_pos(4:Nx+4,:,:); Fp2 = F_pos(5:Nx+5,:,:);
F_pos_ = Weno5_p(Fp_2,Fp_1,Fp0,Fp1,Fp2,rec,epsilon);

Fn_2 = F_neg(2:Nx+2,:,:); Fn_1 = F_neg(3:Nx+3,:,:); Fn0 = F_neg(4:Nx+4,:,:); Fn1 = F_neg(5:Nx+5,:,:); Fn2 = F_neg(6:Nx+6,:,:);
F_neg_ = Weno5_n(Fn_2,Fn_1,Fn0,Fn1,Fn2,rec,epsilon);

F_ = F_pos_ + F_neg_;


% y方向
Uhy = 0.5*(Uyp(:,3:Ny+3,:) + Uyp(:,4:Ny+4,:));
% 重构通量函数
Fp = Flux_X(Uyp,gamma); 
Gp = Flux_Y(Uyp,gamma);
G=J(4:3+Nx,:).*(px2(4:3+Nx,:).*Fp+py2(4:3+Nx,:).*Gp);

G_pos = zeros(size(G));
G_neg = zeros(size(G));
switch splitting_type
    case 'LF'
        %采用LF重构
        Wy=Uyp.*J(4:3+Nx,:);
        lamda_y = MaxEigY(Uhy,gamma);
        G_pos = 0.5*( G + lamda_y*Wy );    
        G_neg = 0.5*( G - lamda_y*Wy ); 
    case 'SW'
        %采用SW重构
        [G_pos,G_neg] = FluxSW_Y(Uyp,gamma);
        G_pos = 0.95.*G_pos.*J(4:3+Nx,:);
        %G_neg = G_neg.*J(4:3+Nx,:);
        G_neg = G - G_pos;
    case 'VL'
        [G_pos, G_neg] = FluxSW_Y(Uyp,gamma);
        G_pos = J(4:3+Nx,:).*G_pos;
        %G_neg = G_neg.*J(4:3+Nx,:);
        G_neg = G - G_pos;
end

% 重构五阶WENO所需的通量
Gp_2 = G_pos(:,1:Ny+1,:); Gp_1 = G_pos(:,2:Ny+2,:);  Gp0 = G_pos(:,3:Ny+3,:); Gp1 = G_pos(:,4:Ny+4,:); Gp2 = G_pos(:,5:Ny+5,:);
G_pos_ = Weno5_p(Gp_2,Gp_1,Gp0,Gp1,Gp2,rec,epsilon);

Gn_2 = G_neg(:,2:Ny+2,:); Gn_1 = G_neg(:,3:Ny+3,:); Gn0 = G_neg(:,4:Ny+4,:); Gn1 = G_neg(:,5:Ny+5,:);  Gn2 = G_neg(:,6:Ny+6,:);
G_neg_ = Weno5_n(Gn_2,Gn_1,Gn0,Gn1,Gn2,rec,epsilon);

G_ = G_pos_ + G_neg_;
end

% 计算X方向通量
function F = Flux_X(U,gamma)
rho = U(:,:,1);
E   = U(:,:,4);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;

F(:,:,1) = U(:,:,2);
F(:,:,2) = 0.5*rho.*( (3-gamma)*u.^2 + (1-gamma)*v.^2 ) + (gamma-1)*E;
F(:,:,3) = U(:,:,2).*v;
F(:,:,4) = gamma*E.*u - 0.5*(gamma-1)*rho.*u.*(u.^2+v.^2);

end

% 计算Y方向通量
function G = Flux_Y(U,gamma)
rho = U(:,:,1);
E   = U(:,:,4);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;

G(:,:,1) = U(:,:,3);
G(:,:,2) = U(:,:,2).*v;
G(:,:,3) = 0.5*rho.*( (1-gamma)*u.^2 + (3-gamma)*v.^2 ) + (gamma-1)*E;
G(:,:,4) = gamma*E.*v - 0.5*(gamma-1)*rho.*v.*(u.^2+v.^2);

end
