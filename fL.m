function xdot = fL(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Function  Author: Ricardo Sanfelice (Revised by Giampiero Campa)
%
% Project: Simulation of a hybrid system (Bouncing Ball)
%
% Name: f.m
%
% Description: Flow map
%
% Version: 1.0
% Required files: - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K NOISEP tnoise
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
mu = 3.98600444*10^14; 
ro = 7100*1000; rd = sqrt((ro+ysys(1))^2+ysys(2)^2+ysys(3)^2);
n = sqrt(mu/ro^3);
A = [0     0  0  1 0 0;
     0     0  0  0 1 0;
     0     0  0  0 0 1;
     3*n^2 0  0  0 2*n 0;
     0     0  0 -2*n 0 0;
     0     0  -n^2 0 0 0];
m = 1*500;  % needs to be updated 
B = [0  0 0;0 0  0;0 0 0; 1/m  0 0;0  1/m 0;0  0 1/m];

%   


% FNL =  [0;0;0;(n^2-(mu/rd^3))*(ro+ysys(1));n^2*ysys(2)-(mu/rd^3)*ysys(2);-(mu/rd^3)*ysys(3)];
FNL =  [0;0;0;-2*n^2*ysys(1)+(mu/ro^2)-((mu/rd^3)*(ro+ysys(1)));n^2*ysys(2)-(mu/rd^3)*ysys(2);n^2*ysys(3)-(mu/rd^3)*ysys(3)];

% %-------------------%
 FNLhat =  [0;0;0;-2*n^2*xhat(1)+(mu/ro^2)-((mu/rd^3)*(ro+xhat(1)));n^2*xhat(2)-(mu/rd^3)*xhat(2);n^2*xhat(3)-(mu/rd^3)*xhat(3)];
 FNLs =  [0;0;0;-2*n^2*x(1)+(mu/ro^2)-((mu/rd^3)*(ro+x(1)));n^2*x(2)-(mu/rd^3)*x(2);n^2*x(3)-(mu/rd^3)*x(3)];

%SOM =  (mu/ro^4)*[0;0;0;-3*x(1)^2+(3/2*x(2)^2)+(3/2*x(3)^2);3*x(1)*x(2);3*x(1)*x(3)];
%SOMhat =  (mu/ro^4)*[0;0;0;-3*xhat(1)^2+(3/2*xhat(2)^2)+(3/2*xhat(3)^2);3*xhat(1)*xhat(2);3*xhat(1)*xhat(3)];

%------------Noise----------%
    tmp = abs(tnoise-tn);
    [row col] = min(tmp);
    %col
    pnoise = NOISEP(col);  % process noise
%---------Input-----------%

%up = -1*B*K*xhat-B*m^2*transpose(B)*FNLhat;
% up = -1*B*K*xhat-B*m^2*B'*FNLhat;
%up = -1*B*K*xhat-B*m^2*B'*SOMhat;
up = -1*B*K*xhat;
inp = up;
nrinf = norm(up,inf);
if nrinf > 0.02
    inp = 0.02*(up/nrinf);
end
% % %----------------------FULL NONLINEAR MODEL---------------%
% 
 
 xdot = A*xsys+[0 0 0 pnoise pnoise pnoise]'+inp+FNLs;
 
 ydot = A*ysys+inp+FNL;
 
 xhatdot = A*xhat+inp+FNLhat;
 
 
 
 
 AdFNL = [zeros(3,6);
         -2*n^2-(mu/rd^3) 0 0 0 0 0;
         0 n^2-(mu/rd^3) 0 0 0 0;
         0 0 n^2-(mu/rd^3) 0 0 0];
 Ft = A+AdFNL;

% %----------------------SECOND ORDER MODEL---------------%
% 
% 
% xdot = A*xsys+SOM+[0 0 0 pnoise pnoise pnoise]'+inp;
% 
% ydot = A*ysys+SOM+inp;
% 
% xhatdot = A*xhat+SOMhat+inp;
% 
% Ad = (mu/ro^4)*[0 0 0 0 0 0;0  0 0 0 0;0 0 0 0 0 0;-6*xhat(1) 3*xhat(2) 3*xhat(3) 0 0 0;...
%       3*xhat(2) 3*xhat(1) 0 0 0 0;3*xhat(3) 0 3*xhat(1) 0 0 0];
% 
% Ft = A+Ad;

%-------------------------------------------------%
%Pt = [transpose(PO1);transpose(PO2);transpose(PO3);transpose(PO4);transpose(PO5);transpose(PO6)];

Pta = [PO1';PO2';PO3';PO4';PO5';PO6'];

Pt = 0.5*(Pta+Pta');
Q = 10^-8*[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];

%Ptdot = Ft*Pt+Pt*Ft'+10^-14*[.00001 0 0 0 0 0;0 .00001 0 0 0 0;0 0 .00001 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
Ptdot = -Ft'*Pt - Pt*Ft - Pt*Q*Pt; 

POdot = reshape(Ptdot,[36,1]);


taudot = 1;
tndot = 1;


xdot = [xdot;ydot;xhatdot;POdot;taudot;0;tndot];

