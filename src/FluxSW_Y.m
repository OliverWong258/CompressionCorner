function [g_pos,g_neg] = FluxSW_Y(U,gamma)
rho = U(:,:,1);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;
E = U(:,:,4);
e = E./rho - 0.5.*(u.^2+v.^2);
P = (gamma-1).*rho.*e;
c = sqrt(gamma*P./rho);
h = gamma.*e;

lamda1 = v;
%lamda2 = v;
lamda3 = v-c;
lamda4 = v+c;

lamda1_pos = 0.5.*(lamda1+abs(lamda1));
%lamda2_pos = 0.5.*(lamda2+abs(lamda2));
lamda3_pos = 0.5.*(lamda3+abs(lamda3));
lamda4_pos = 0.5.*(lamda4+abs(lamda4));

%lamda1_neg = 0.5.*(lamda1-abs(lamda1));
%lamda2_neg = 0.5.*(lamda2-abs(lamda2));
%lamda3_neg = 0.5.*(lamda3-abs(lamda3));
%lamda4_neg = 0.5.*(lamda4-abs(lamda4));

g_pos = zeros(size(U));
g_neg = zeros(size(U));

g_pos(:,:,1) = (rho./(2*gamma)).*(2.*(gamma-1).*lamda1_pos + lamda3_pos + lamda4_pos);
g_pos(:,:,3) = (rho./(2*gamma)).*(2.*(gamma-1).*v.*lamda1_pos + (v-c).*lamda3_pos + (v+c).*lamda4_pos);
g_pos(:,:,2) = (rho./(2*gamma)).*(2.*(gamma-1).*u.*lamda1_pos + u.*lamda3_pos + u.*lamda4_pos);
g_pos(:,:,4) = (rho./(2*gamma)).*((gamma-1).*(u.^2+v.^2).*(lamda1_pos) + (h-c.*u).*lamda3_pos + (h+c.*u).*lamda4_pos);

%g_neg(:,:,1) = (rho./(2*gamma)).*(2.*(gamma-1).*lamda1_neg + lamda3_neg + lamda4_neg);
%g_neg(:,:,3) = (rho./(2*gamma)).*(2.*(gamma-1).*v.*lamda1_neg + (v-c).*lamda3_neg + (v+c).*lamda4_neg);
%g_neg(:,:,2) = (rho./(2*gamma)).*(2.*(gamma-1).*u.*lamda1_neg + u.*lamda3_neg + u.*lamda4_neg);
%g_neg(:,:,4) = (rho./(2*gamma)).*((gamma-1).*(u.^2+v.^2).*(lamda1_neg) + (h-c.*u).*lamda3_neg + (h+c.*u).*lamda4_neg);
end



