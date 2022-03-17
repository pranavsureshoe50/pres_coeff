
%% INPUT

export_body = 0;
path = '.';

%% BODY
angle = 0; 
alfa  = angle * pi/180;
AR = 4;
n  = 700;
b  = 0.5;
shift = 0;
bo = 2; %outer domain

%% CALCULATION
% BODY

a=AR*b;
ao=AR*bo;

x = linspace(-a,a,n);
xo = linspace(-ao,ao,n); %outer domain

%non uniform body definition
stb = @(x) sin(x.*pi/(2*a));
%plot(x,a*sin(x*pi/(2*a)))
stbo = @(x) sin(x.*pi/(2*ao));

x=a*stb(x);
xo=ao*stbo(xo);

zp = sqrt(b.^2-AR^(-2)*x.^2);
zn = -sqrt(b.^2-AR^(-2)*x.^2);
%plot(x,zp,x,zn);

zpo = sqrt(bo.^2-AR^(-2)*xo.^2); %outer domain
zno = -sqrt(bo.^2-AR^(-2)*xo.^2);

xpl=cos(alfa)*x - sin(alfa)*zp; %local reference frame
zpl=sin(alfa)*x + cos(alfa)*zp +shift;
xnl=cos(alfa)*x - sin(alfa)*zn;
znl=sin(alfa)*x + cos(alfa)*zn +shift;

xpol=cos(alfa)*xo - sin(alfa)*zpo; %local reference frame, outer domain
zpol=sin(alfa)*xo + cos(alfa)*zpo +shift;
xnol=cos(alfa)*xo - sin(alfa)*zno;
znol=sin(alfa)*xo + cos(alfa)*zno +shift;
