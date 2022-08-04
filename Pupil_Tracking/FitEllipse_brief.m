function params = FitEllipse_brief(ix,iy)
params = EllipseDirectFit([iy ix]);

if params(1)<0
    params = -params;
end
a = params(1);
b = params(2)/2;
c = params(3);
d = params(4)/2;
e = params(5)/2;
f = params(6);
clear params;
% checking if the ellipse is a proper real ellipse
emat = [a b d; b c e; d e f];
lmat = [a b; b c];
good_ellipse = 0;
if det(emat)~=0 && det(lmat)>0 && det(emat)/(a+c)<0
    good_ellipse = 1;
end

% find minor/major axis and center
xc = (c*d - b*e)/(b^2 - a*c);
yc = (a*e - b*d)/(b^2 - a*c);
u0 = 2*(a*e^2 + c*d^2 + f*b^2 - 2*b*d*e - a*c*f);
ra = sqrt(u0/((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
rb = sqrt(u0/((b^2-a*c)*(-1*sqrt((a-c)^2+4*b^2)-(a+c))));
% angle
% keyboard
% alpha=asin(b*2/(1/(rb^2)-1/(ra^2)))/2;

params.ra = ra;
params.rb = rb;
params.xc  = xc;
params.yc  = yc;
% params.alpha=alpha;
params.isgood = good_ellipse;
% params.coefs=[a b c d e f];
