function  fmx = Weno5_n(f_2,f_1,f,f1,f2,rec,epsilon)
% 重构通量的负向部分.

f0 = (2 * f_2 - 7 * f_1 + 11 * f) / 6;
f1 = (-f_1 + 5 * f + 2 * f1) / 6;
f2 = (2 * f + 5 * f1 - f2) / 6;

beta0 = 13 / 12 * ((f_2 - 2 * f_1 + f).^2) + 1 / 4 .* ((f_2 - 4 * f_1 + 3 * f).^2);
beta1 = 13 / 12 * ((f_1 - 2 * f + f1).^2) + 1 / 4 .* ((f_1 - f1).^2);
beta2 = 13 / 12 * ((f - 2 * f1 + f2).^2) + 1 / 4 .* ((3 * f - 4 * f1 + f2).^2);

alpha0 = rec(1)./((epsilon + beta0).^2);
alpha1 = rec(2)./((epsilon + beta1).^2);
alpha2 = rec(3)./((epsilon + beta2).^2);
alphasum = alpha0 + alpha1 + alpha2;

w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

fmx = w0.*f0 + w1.*f1 + w2.*f2;
end
