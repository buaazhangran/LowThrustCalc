%% system calc for differential algebra equations
syms p f g h k L
syms miu F Lm
% two most used variable
w = 1 + f * cos(L) + g * sin(L);
s2 = 1+h*h+k*k;
% f and g don't have an analytic expression
dp_max = 2*p/(1-sqrt(f^2+g^2))*sqrt(p/miu)*F;
df_max = sqrt(p/miu)*F/w*sqrt(w^2*sin(Lm)*sin(Lm)+((w+1)*cos(Lm)+f)^2+(h*sin(Lm)-k*cos(Lm))^2*g^2);
dg_max = sqrt(p/miu)*F/w*sqrt(w^2*cos(Lm)*cos(Lm)+((w+1)*sin(Lm)+g)^2+(h*sin(Lm)-k*cos(Lm))^2*f^2);
dh_max = -sqrt(p/miu)*s2*sqrt(1-g^2)/2/(1-f*sqrt(1-g^2)-g^2);
dk_max = -sqrt(p/miu)*s2*sqrt(1-f^2)/2/(1-g*sqrt(1-f^2)-f^2);

syms pT fT gT hT kT 
syms wp wf wg wh wk
% Lyapunov function
V = wp*((p-pT)/dp_max)^2+wf*((f-fT)/df_max)^2+wg*((g-gT)/dg_max)^2+wh*((h-hT)/dh_max)^2+wk*((k-kT)/dk_max)^2;
tic
% the result of dV/dp, dV/df, dV/dg
dVdp = diff(V,f);
% dVdp_simple = simplify(dVdp)
% pretty(dVdp_simple)
pretty(dVdp)
toc