%% find Lm for df which makes the rate df the largest
function Lm = findLm4df(f,g,h,k);
L = linspace(0,2*pi,100);
w = 1+ f*cos(L) + g*sin(L);
df = 1./w.*sqrt(w.*w.*sin(L).*sin(L)+((1+w).*cos(L)+f).^2+g^2*(h*sin(L)-k*cos(L)).^2);

