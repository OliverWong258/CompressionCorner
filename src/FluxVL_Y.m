function [f_pos, f_neg] = FluxVL_Y(U,gamma)
rho = U(:,:,1);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;
E = U(:,:,4);
V = sqrt(u.^2+v.^2);
e = E./rho - 0.5.*(u.^2+v.^2);
P = (gamma-1).*rho.*e;
c = sqrt(gamma*P./rho);
M = sqrt(u.^2+v.^2)./c;

f_mass_pos = rho.*c.*(M+1).^2./4;
%f_mass_neg = -rho.*c.*(M-1).^2./4;

mask = M < -1;  % 创建一个逻辑掩码
f_mass_pos(mask) = 0;  % 使用掩码设置值
%mask = M > 1;
%f_mass_neg(mask) = 0;

f_pos = zeros(size(U));
f_neg = zeros(size(U));

f_pos(:,:,1) = f_mass_pos;
f_pos(:,:,3) = f_mass_pos.*((-V+2.*c)./gamma+v);
f_pos(:,:,2) = f_mass_pos.*u;
f_pos(:,:,4) = f_mass_pos.*(((gamma-1).*V+2.*c).^2/(2.*(gamma.^2-1)));
%f_neg(:,:,1) = f_mass_neg;
%f_neg(:,:,3) = f_mass_neg.*((-V-2.*c)./gamma+v);
%f_neg(:,:,2) = f_mass_neg.*u;
%f_neg(:,:,4) = f_mass_neg.*(((gamma-1).*V-2.*c).^2/(2.*(gamma.^2-1)));

end

