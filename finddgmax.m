%% find dgmax for df which makes the rate dg the largest
function dgmax = finddgmax(f,g,h,k);
L = linspace(0,2*pi,100);
w = 1+ f*cos(L) + g*sin(L);
dg = 1./w.*sqrt(w.*w.*cos(L).*cos(L)+((1+w).*sin(L)+g).^2+f^2*(h*sin(L)-k*cos(L)).^2);
dgmax = max(dg);