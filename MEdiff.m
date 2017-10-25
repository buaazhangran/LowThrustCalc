%% this is the differential equations of the modified equinotical elements
% the Lyapunov feedback control are used in this script.
function dX = MEdiff(X,Const);
% X are modified equinoctial elements.
p = X(1);
f = X(2);
g = X(3);
h = X(4);
k = X(5);
L = X(6);
m = X(7);
% Const = [Tmax,pT,fT,gT,hT,kT,wp,wf,wg,wh,wk];
Tmax = Const(1);
pT = Const(2);
fT = Const(3);
gT = Const(4);
hT = Const(5);
kT = Const(6);
wp = Const(7);
wf = Const(8);
wg = Const(9);
wh = Const(10);
wk = Const(11);

miu = 3.986004415e14; 
F = Tmax/m;
w = 1 + f * cos(L) + g * sin(L);
s2 = 1+h*h+k*k;
% we consider de_max as constant variables
% maximum rates for different elements
dp_max = 2*p/(1-sqrt(f^2+g^2))*sqrt(p/miu)*F;
% we should find df_max that makes df_max the largest first numerically.
df_max = sqrt(p/miu)*F*finddfmax(f,g,h,k);
dg_max = sqrt(p/miu)*F*finddgmax(f,g,h,k);
dh_max = -sqrt(p/miu)*s2*sqrt(1-g^2)/2/(1-f*sqrt(1-g^2)-g^2);
dk_max = -sqrt(p/miu)*s2*sqrt(1-f^2)/2/(1-g*sqrt(1-f^2)-f^2);
% partial differential functions
Kp = 2*wp*(p-pT)/dp_max^2;
Kf = 2*wf*(f-fT)/df_max^2;
Kg = 2*wg*(g-gT)/dg_max^2;
Kh = 2*wh*(h-hT)/dh_max^2;
Kk = 2*wk*(k-kT)/dk_max^2;
% calc the flight angles
tanalpha = w*(Kf*sin(L)-Kg*cos(L))/(2*p*Kp+Kf*((1+w)*cos(L)+f)+Kg*((1+w)*sin(L)+g));
alpha = atan(tanalpha);




