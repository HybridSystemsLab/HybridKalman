function xplus = g(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Function  Author: Ricardo Sanfelice (Revised by Giampiero Campa)
%
% Project: Simulation of a hybrid system (Bouncing ball)
%
% Name: g.m
%
% Description: Jump map
%
% Version: 1.0
% Required files: - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K NOISE1 tnoise V sigmam
xsys = x(1:6);
ysys = x(7:12);
xhat = x(13:18);
PO1 = x(19:24);
PO2 = x(25:30);
PO3 = x(31:36);
PO4 = x(37:42);
PO5 = x(43:48);
PO6 = x(49:54);
tau = x(55);
h = x(56);
tn = x(57);

tmp = abs(tnoise-tn);
[row col] = min(tmp);
Ynoise = NOISE1(col);


Pkma = [PO1';PO2';PO3';PO4';PO5';PO6'];

%----to avoid an possible imaginary matrix elements------------%
Pkm = 0.5*(Pkma+transpose(Pkma));

%---- Angle only feedback for 2 states LVLH no outof plane z----%
% LOS output estimation ----- to avoid singularities with atan() function

% Pkma = [PO1';PO2';PO3';PO4'];
% Pkm = 0.5*(Pkma+transpose(Pkma));
% rho = [xhat(1);xhat(2)];
% nrho = norm(rho);
% H1 = 1/nrho *eye(2)-rho*rho'/nrho^3;
% H2 = zeros(2,2);
% Hkm = [H1 H2];
% ys = [ysys(1)/sqrt(ysys(1)^2+ysys(2)^2); ysys(2)/sqrt(ysys(1)^2+ysys(2)^2)]+[(Ynoise); (Ynoise)];
% yhat = [xhat(1)/sqrt(xhat(1)^2+xhat(2)^2); xhat(2)/sqrt(xhat(1)^2+xhat(2)^2)];   
% 
% Rk  = [((0.001))^2 0;0 ((0.001))^2];    
%%Kk  = Pkm*Hkm'/(Hkm*Pkm*Hkm'+Rk); % use transpose() instead of '
%Kk  = Pkm*transpose(Hkm)/(Hkm*Pkm*transpose(Hkm)+Rk);
%    xhatp = xhat+Kk*(ys-yhat);
%----regular update for the covariance matrix--------%    
%    Ptp   = (eye(4)-Kk*Hkm)*Pkm;
%--------Joseph form to avoid numerical issues with the Covariance matrix
%collapse or go to zero------------------%
%Ptp   = (eye(4)-Kk*Hkm)*Pkm*transpose(eye(4)-Kk*Hkm) +Kk*Rk*transpose(Kk) ;    
%    Ptk = reshape(Ptp,[16,1]);
%------------------------------------------------------------------------------------%


H1 = eye(3);
H2 = zeros(3,3);

Hkm = [H1 H2];


ys = [ysys(1); 
      ysys(2);
      ysys(3)]+[(Ynoise); (Ynoise);(Ynoise)];
yhat = [xhat(1);
        xhat(2);
        xhat(3)];   


Rk  = sigmam^2*[1 0 0;
            0 1 0;
            0 0 1];




% Information Matrix Propogation
Ptp = Pkm + (transpose(Hkm)/(Rk) *Hkm);
Kk = Ptp\transpose(Hkm)/(Rk);
xsysp = xsys;
ysysp = ysys;
xhatp = xhat+Kk*(ys-yhat);
%----regular update for the covariance matrix--------%
%Ptp   = (eye(6)-Kk*Hkm)*Pkm;
%--------Joseph form to avoid numerical issues with the Covariance matrix
%collapse or go to zero------------------%
%Ptp   = (eye(6)-Kk*Hkm)*Pkm*transpose(eye(6)-Kk*Hkm) +Kk*Rk*transpose(Kk) ;
Ptk = reshape(Ptp,[36,1]);
hp = 1;
taup = 0;

Tmin = 2;
Tmax = 6;
V = Tmin+(Tmax-Tmin)*rand(1);

xplus =[xsysp;ysysp;xhatp;Ptk;taup;hp;tn];

end

