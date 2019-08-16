%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file       Project: HyEQ Toolbox  @ Hybrid Dynamics and Control
% Lab, http://www.u.arizona.edu/~sricardo/index.php?n=Main.Software
%
% Filename: run_switch.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NOISE1 tnoise K NOISEP V sigmap sigmam

V = 4;
format long
mu = 3.98600444*10^14; ro = 7100*1000;
Tc = 6000;%4*pi*sqrt(ro^3/mu);
n = sqrt(mu/ro^3); 
A = [0     0  0  1 0 0;
     0     0  0  0 1 0;
     0     0  0  0 0 1;
     3*n^2 0  0  0 2*n 0;
     0     0  0 -2*n 0 0;
     0     0  -n^2 0 0 0];
m = 1*500;  % needs to be updated 
B = [0  0 0;0 0  0;0 0 0; 1/m  0 0;0  1/m 0;0  0 1/m];


% xint = [-10000 0 0 .5 .5 .5];
% yint = [-10000 0 0 .5 .5 .5];
% xhatint = [-10000 0 0 .5  .5 .5]+[-1000 1000 1000 0 0 0];
% xint = [0 -10000 0 0.5 .5 .5];
% yint = [0 -10000 0 .5 .5 .5];
% xhatint = [0 -10000 0 .5 .5 .5]+[-1000 1000 1000 0 0.0 0.0];
% xint = [0 0 -10000 0.5 .5 .5];
% yint = [0 0 -10000 .5 .5 .5];
% xhatint = [0 0 -10000 .5 .5 .5]+[-1000 1000 1000 0 0.0 0.0];
xint = [7071 -7071 7071  .5 .5 .5];
yint = [7071 -7071 7071 .5 .5 .5];
xhatint = [7071 -7071 7071 .5 .5 .5]+[-1000 1000 1000 0.0 0.0 0.0];
% xint = [0 100 -100 0 0 0];
% yint = [0 100 -100 0 0 0];
% xhatint = [0 100 -100 0 0 0]+[-1 1 1 0 0.0 0.0];
% po = 4*10^10*([1*10^5 0 0 0 0 0;
%       0 1*10^4 0 0 0 0;
%       0 0 1*10^5 0 0 0;
%       0 0 0 1*10^-1 0 0;
%       0 0 0 0 1*10^-2 0;
%       0 0 0 0 0 1*10^-2]);
%---------------TEST-------------%

%po = 0.192*eye(6);

po = 4*([1*10^8 0 0 0 0 0;
      0 1*10^8 0 0 0 0;
      0 0 1*10^8 0 0 0;
      0 0 0 0.8*10^3 0 0;
      0 0 0 0 0.8*10^3 0;
      0 0 0 0 0 .8*10^4]);
  
  
poi = reshape(po,[1,36]);
tau = 0;
hin = 1;
tnN = 0;   

x0 = [xint yint xhatint poi tau hin tnN];
 
%-----------------NOISE-----------%
N=150000;
Fs = 10;
tnoise = (0:N-1)/Fs;
sigma = 1;
sigmap=sqrt(10^-8);
sigmam = 20;
NOISE1 = (20)*sigma*randn(size(tnoise));
NOISEP = sigmap*sigma*randn(size(tnoise));

% Q1 = 1.5e10*eye(6);

Q1 = 1*[1.5e-1 0 0 0 0 0;
      0 2.5e-1 0 0 0 00
      0 0 1.5e-1 0 0 0;
      0 0 0 1.5e-2 0 0;
      0 0 0 0 1.5e-2 0;
      0 0 0 0 0 1.5e-2];

 %R1 = 20e4*eye(3); 
 R1 = 1*[20e4 0 0;
       0 15e4 0;
       0 0 99e3];
 [K,s,e] = lqr(A,B,Q1,R1);


% simulation horizon
T = [0 Tc];                                                                
J = [0 10e80];
% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);
%options = odeset('RelTol',1e-6,'MaxStep',1);

% simulate
[t, j, x] = HyEQsolver( @fL,@g,@C,@D,x0',T,J,rule,options);

% figure(200)
% plotflows(t,j,x(:,55))
% hold on
% plot(t,2*ones(length(t)))

%%
% %plot hybrid arc   
figure(1)
plot(t,x(:,1)/1000); 
hold on
plot(t,x(:,13)/1000,'r'); 
hold off
grid on
xlabel('t')                    
ylabel('x')                    

%%
% %plot hybrid arc   
figure(2)
plot(t,x(:,2)/1000); 
hold on
plot(t,x(:,14)/1000,'r'); 
hold off
grid on
xlabel('t')                    
ylabel('y')                    
%zlabel('x1')   

%%
%plot hybrid arc   
figure(3)
plot(t,x(:,3)/1000); 
hold on
plot(t,x(:,15)/1000,'r'); 
hold off
grid on
xlabel('time (sec)')                    
ylabel('Z (KM)')                    
%zlabel('x1')    
%%
clear xnom
for i= 1:1:length(t)
  xnomervel(i,:)= norm([(x(i,4)-x(i,16)) (x(i,5)-x(i,17)) (x(i,6)-x(i,18))],2);
  xnomerfnl(i,:)= norm([(x(i,1)-x(i,13)) (x(i,2)-x(i,14)) (x(i,3)-x(i,15))],2);
  xfnl(i,:)= ((x(i,1)-x(i,13)));
  yfnl(i,:)= ((x(i,2)-x(i,14)));
  zfnl(i,:)= ((x(i,3)-x(i,15)));
  Lfnl(i,:)= norm([(x(i,4)-x(i,16)) (x(i,5)-x(i,17)) (x(i,6)-x(i,18))],inf)/norm([(x(i,1)-x(i,13)) (x(i,2)-x(i,14)) (x(i,3)-x(i,15))],inf);
end
%%
figure(12)  
plot(t,xnomerfnl,'b')
% hold on
% plot(t,xnomerfnl,'r')
xlabel('time(sec)')
ylabel('error $e = \eta-\eta_c$ (m)')
grid on
hold on
%%
figure(500)
set(gca,'Xdir','reverse')
sp = 1:1000:length(t);
angle = linspace(0,2*pi,360);
xc = 1000*cos(angle);
yc = 1000*sin(angle);
plot(xc/100,yc/100,'k') 
plot(x(sp,14)/1000,x(sp,13)/1000,'r')
hold on
plot(x(sp,2)/1000,x(sp,1)/1000,'k')
hold on
plot(x(1,14)/1000,x(1,13)/1000,'r*','linewidth',5)
hold on
plot(x(1,2)/1000,x(1,1)/1000,'b*','linewidth',1)
angle = linspace(0,2*pi,360);
xc = 1000*cos(angle);
yc = 1000*sin(angle);
plot(xc/1000,yc/1000,'g')
xlabel('True and Estimated y-position in Km')
ylabel('True and Estimated x-position in Km')
grid on
hold on
plot(xc/100,yc/100,'k')
%%
figure(800)
set(gca,'Xdir','reverse')
sp = 1:1000:length(t);
angle = linspace(0,2*pi,360);
xc = 1000*cos(angle);
yc = 1000*sin(angle);
plot(xc/100,yc/100,'k') 
plot(x(sp,14)/1000,x(sp,15)/1000,'g')
hold on
plot(x(sp,2)/1000,x(sp,3)/1000,'b')
hold on
plot(x(1,14)/1000,x(1,15)/1000,'r*','linewidth',5)
hold on
plot(x(1,2)/1000,x(1,3)/1000,'b*','linewidth',1)
angle = linspace(0,2*pi,360);
xc = 1000*cos(angle);
yc = 1000*sin(angle);
plot(xc/1000,yc/1000,'g')
xlabel('True and Estimated y-position in Km')
ylabel('True and Estimated z-position in Km')
grid on
hold on
plot(xc/100,yc/100,'k')
%%
figure(600)
view(73,10)
[x1 y1 z1] = sphere;
 plot3(10*x1, 10*y1, 10*z1,'color',[0 0 0],'LineStyle',':'); 
 hold on
 contour3(10*x1, 10*y1, 10*z1,'k:');
grid on
 hold on
plot3(x(sp,14)/1000,x(sp,15)/1000,x(sp,13)/1000,'b','linewidth',2)
hold on
plot3(x(sp,2)/1000,x(sp,3)/1000,x(sp,1)/1000,'r','linewidth',2)
xlabel('Local horizontal - y axis')
h=get(gca,'xlabel');
set(h,'rotation',-28)
zlabel('Local verical - x axis')
ylabel('Out of plane - z axis')
h=get(gca,'ylabel');
set(h,'rotation',0)
%%
for i= 1:1:length(t)
   % NL = B*m^2*B'*(mu/ro^4)*[0;0;-3*x(i,9)^2+(3/2*x(i,10)^2);3*x(i,9)*x(i,10)];
    unom(i,1)= norm((B*K*[x(i,13) x(i,14) x(i,15) x(i,16) x(i,17) x(i,18)]'),inf);
    if unom(i,1) >0.02
        unom(i,1) = 0.02;
    end
end

figure(5)
%plot(t,x(:,6))
grid on
hold on
plot(t,0.02*ones(size(t)),'r')
hold on
plot(t,unom)
xlabel('time(sec)')
ylabel('$\frac{1}{m}\|u\|_{\infty}$ (m/sec$^2$)')
legend('Max input','calculated input')
hold on
%CF = trapz(t,unom)

%%
figure(6)
plot(t,atan2(x(:,8),x(:,7)),'r')
grid on
hold on
plot(t,atan2(x(:,14),x(:,13)))
hold on
xlabel('time(sec)')
ylabel('angle (rad)')
%legend('Max input','calculated input')
hold on
%CF = trapz(t,unom)

%%


% Plot Lyapunov function
% Pm = vec2mat(x(:,(19:54)),6);
for i= 1:1:length(t)
    Pm = vec2mat(x(i,(19:54)),6); 
    PMC = inv(Pm);
    sig1(i,1) = PMC(1,1);
    sig2(i,1) = PMC(2,2);
    sig3(i,1) = PMC(3,3);
    clear PMC
    e = x(i,(1:6))' - x(i,(13:18))'; 
    emn(i) = max(eig(Pm));
    emx(i) = min(eig(Pm));
    tm(i,1) = t(i);
    Lf(i) = exp(-1*x(i,55))* norm(e'* Pm* e,2);
  P(i) = norm(Pm,2);
end
%%
figure(13)
SIG = 3;
subplot(3,1,1)
plot(t,xfnl,'b')
hold on
plot(t,SIG*sig1,'r')
hold on
plot(t,-SIG*sig1,'r')
grid on
ylabel('$e_x$ (m)')
hold off
subplot(3,1,2)
plot(t,yfnl,'b')
hold on
plot(t,SIG*sig2,'r')
hold on
plot(t,-SIG*sig2,'r')
grid on
ylabel('$e_y$ (m)')
hold off
subplot(3,1,3)
plot(t,zfnl,'b')
hold on
plot(t,SIG*sig3,'r')
hold on
plot(t,-SIG*sig3,'r')
grid on
ylabel('$e_z$ (m)')
xlabel('time(sec)')
hold off
%% PLOT THE LYAPUNOV FUNCTION
figure(6)
plot(t,Lf)

 Q = 10^-8*[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
a = 1+2*n
L = 2*n;%10^-7; %2*n; %
Tmin = 2;
%2*n; 
%L = 1/ro^3 %;%2.79*10e-21% 1/|ro^3|
uq = max(eig(Q));
lq = min(eig(Q));
r = 20^2;
k1 = min(eig(po))
k2 = (4*a)/lq% + (r^2*1^2*lq/(4*L)) - (4*a*r*1)/lq
k1t = min(emx)
c = 1/(k2*(0.9*k1t +r*Tmin))
k1tc_ub = sqrt(4*L*(k2)/lq)
k1tc_lb = 1/c/k2-r*Tmin
theta = 10*4*L*k2/(k1tc_lb^2*lq)
Ktime = lq*k1tc_lb-4*L*k2/(theta*k1tc_lb)
sigma = Ktime/2
T_max = Tmin*(1+(Ktime/sigma))
%%
figure(28)
plot(t(1:i),emx)
xlabel('time(sec)')
ylabel('Lowerbound on \phi_P(t,j), \tilde k_1')
hold on
plot(t(1:i),k1*ones(size(t(1:i))))
 hold on
 plot(t(1:i),k1tc_lb*ones(size(t(1:i))))
grid on
hold off
% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes

%indexOfInterest = (t < 100) & (t > 80); % range of t near perturbation

plot(t(1:i),emx(1:i)) % plot on new axes
hold on
%plot(t(1:i),k1tc*ones(size(t(1:i))),'b') 
%hold on
plot(t(1:i),k1tc_lb*ones(size(t(1:i))))
grid on
axis tight
%%
figure(29)
plot(t(1:i),emn)
hold on
plot(t(1:i),k2*ones(size(t(1:i))))
hold on
plot(t(1:i),k1*ones(size(t(1:i))))
grid on
xlabel('time(sec)')
ylabel('Upperbound on $phi_P(t,j)$, k_2')
hold off
% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes

%indexOfInterest = (t < 100) & (t > 80); % range of t near perturbation

plot(t(1:i),emn(1:i)) % plot on new axes
hold on
plot(t(1:i),k1*ones(size(t(1:i))))
% grid on
%axis tight

%%
figure(100)
plotflows(t(1:5500),j(1:5500),x((1:5500),55))
hold on
plot(t(1:5500),2*ones(length(t(1:5500))),'r','Linewidth',2)
hold on
plot(t(1:5500),6*ones(length(t(1:5500))),'r','Linewidth',2)
xlabel('time(sec)')
ylabel('Timer, $\tau$')
grid on
hold off


%%
clear xnom
for i= 1:1:length(t)
    for k = 1:1:length(j)
  xnomervel(i,k,:)= norm([(x(i,k,4)-x(i,k,16)) (x(i,k,5)-x(i,k,17)) (x(i,k,6)-x(i,k,18))],2);
%   xnomerfnl(i,:)= norm([(x(i,1)-x(i,13)) (x(i,2)-x(i,14)) (x(i,3)-x(i,15))],2);
%   xfnl(i,:)= norm((x(i,1)-x(i,13)),2);
%   yfnl(i,:)= norm((x(i,2)-x(i,14)),2);
%   zfnl(i,:)= norm((x(i,3)-x(i,15)),2);
%   Lfnl(i,:)= norm([(x(i,4)-x(i,16)) (x(i,5)-x(i,17)) (x(i,6)-x(i,18))],inf)/norm([(x(i,1)-x(i,13)) (x(i,2)-x(i,14)) (x(i,3)-x(i,15))],inf);
     end
end