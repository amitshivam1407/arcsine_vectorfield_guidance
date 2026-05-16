%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%--------------- declaration of constant parameters ----------------------%
global  Vgd kv k_kai k_kaidot k2 figure_type_good
Vgd = 10;
kv = 5;
k_kai = 50;  
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
figure_type_good = 1;

%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%

tspan = [0 25];

x0 = 100;
y0 = -100;
k2 = 0.002;
if x0 < 0
   kaid0 =  asin(1./(1+ k2 *(x0).^2)) ;
%    kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
   kaid0 = pi - asin(1./(1 + k2 *(x0).^2)) ;  
%   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
kai0 = kaid0;
kaidot0 = -(2*k2*Vgd.*x0)./((1+k2*(x0).^2).^2);

% k2 = abs(14.92597/x0^2);

x_initial = [x0;y0;kai0;kaidot0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -110:10:110;
[X,Y] = meshgrid(range);
[kaid,kaid_dot] = vf_proposed(X);
xdot  = Vgd*cos(kaid);
ydot  = Vgd*sin(kaid);

kaid_mod = vf_proposed_mod(X);
xdot_mod  = Vgd*cos(kaid_mod);
ydot_mod  = Vgd*sin(kaid_mod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@(t,x)fun_stline_upd(t,x) ,tspan, x_initial,options);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vg_act = sqrt(x(:,5).^2 + wx^2 + wy^2 + 2*wx*x(:,5).*cos(psi) + 2*wy*x(:,5).*sin(psi));
x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;
for i = 1:length(t)
if x(i,1) < 0
   kai_des(i,1) =  asin(1./(1+ k2*(x(i,1)).^2)) ;
   kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
   kai_des(i,1) = pi - asin(1./(1 + k2*(x(i,1)).^2)) ;  
   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
end
kaidot_actual = x(:,4);
kappa_des = kaidot_des./Vgd;
kappa_actual = x(:,4)./Vgd ;
% [kai_des,kaidot_des] = vf_proposed(x(:,1));

% Va_des = sqrt((Vgd*cos(kai_des) - wx).^2 + (Vgd*sin(kai_des) - wy).^2) ;
% psi_des = kai_des + asin((wx*sin(kai_des) - wy*cos(kai_des))./ Va_des );
% % psiddot_des = kaidot_des.*((Va_des.*cos(psi_des - kai_des) + (wx*cos(kai_des) + wy*sin(kai_des)))./(Va_des.*cos(psi_des - kai_des)));
% c1 = (1./x(:,5)).*kaidot_des.*(wx*cos(kai_des)+wy*sin(kai_des));
% c2 = (wx*sin(kai_des)-wy*cos(kai_des)).*((1./x(:,5)).^2).*(Vgd./Va_des).*kaidot_des.*(-wx*sin(kai_des)+wy*cos(kai_des));
% psiddot_des = kaidot_des + (1./sqrt(1 - ((1./x(:,5)).*(wx*sin(kai_des)-wy*cos(kai_des))).^2)).*(c1 + c2);
% psi_actual = x(:,3);
% psidot_actual = x(:,4);
% Vg_actual = sqrt((x(:,5).*cos(x(:,3)) + wx).^2 + (x(:,5).*sin(x(:,3)) + wy).^2) ;
% % Vg_actual = (wx.*cos(x(:,3)) + wy.*sin(x(:,3))) + sqrt(x(:,5).^2 - (wx.*sin(x(:,3)) - wy.*cos(x(:,3))).^2);
% kai_actual = x(:,3) + asin((-wx*sin(x(:,3)) + wy*cos(x(:,3)))./ Vg_actual );
% kaidot_actual = x(:,4).*((Vg_actual.*cos(kai_actual - x(:,3)) - (wx*cos(x(:,3)) + wy*sin(x(:,3))))./(Vg_actual.*cos(kai_actual - x(:,3))));
% kappa_actual = kaidot_actual./Vg_actual;
%---------------------- plotting figures ---------------------------------%

figure_plot = fun_stlineplot(t,X,Y,xdot,ydot,xdot_mod,ydot_mod,x(:,1),x(:,2),x_ini,y_ini,x_end,y_end,...
    kaidot_des,kaidot_actual,kappa_des,kappa_actual,kai_des,x(:,3));



function out = fun_stline_upd(t,x)
global Vgd k2 wx wy  kv k_kai k_kaidot


if x(1) < 0
   kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
   kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
else
   kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
   kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
end
% 
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
% out(5,1) = kv*(Vad - x(5));
end

% function out = fun_stline_upd(t,x)
% global Vgd k2 wx wy  kv k_kai k_kaidot
% 
% 
% if x(1) < 0
%    kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
%    kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
% else
%    kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
%    kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
% end
% 
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% % psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% c1 = (1/x(5))*kaid_dot*(wx*cos(kaid)+wy*sin(kaid));
% c2 = (wx*sin(kaid)-wy*cos(kaid))*((1/x(5))^2)*(Vgd/Vad)*kaid_dot*(-wx*sin(kaid)+wy*cos(kaid));
% psiddot = kaid_dot + (1/sqrt(1 - ((1/x(5))*(wx*sin(kaid)-wy*cos(kaid)))^2))*(c1 + c2);
% % psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% % Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));
% 
% 
% out(1,1) = x(5)*cos(x(3))  + wx ;
% out(2,1) = x(5)*sin(x(3))  + wy ;
% out(3,1) = x(4) ;
% out(4,1) = k_kai*(psid - x(3)) + k_kaidot*(psiddot - x(4))  ;
% out(5,1) = kv*(Vad - x(5));
% end

function [kaid,kaid_dot] = vf_proposed(x)
global Vgd  k2 
for i = 1:length(x)
    for j = 1:length(x)
if x(i,j) < 0
   kaid(i,j) = asin(1./(1+k2*(x(i,j)).^2)) ;
   kaid_dot(i,j) =  2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
else
   kaid(i,j) = pi - asin(1./(1+k2*(x(i,j)).^2)) ;
   kaid_dot(i,j) = - 2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
end
    end
end
end



function [kaid_mod,kaid_dot] = vf_proposed_mod(x)
global Vgd  k2 
q = 1.5;
for i = 1:length(x)
    for j = 1:length(x)
if x(i,j) < 0
   % kaid(i,j) = asin(1./(1+k2*(x(i,j)).^2)) ;
   kaid_mod(i,j) = pi/2 - asin(500*k2.*abs(x(i,j)).^q./(1 + 500*k2.*abs(x(i,j)).^q));
   % kaid_dot(i,j) =  2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
else
   % kaid(i,j) = pi - asin(1./(1+k2*(x(i,j)).^2)) ;
   % kaid_dot(i,j) = - 2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
kaid_mod(i,j) = pi/2 + asin(500*k2.*abs(x(i,j)).^q./(1 + 500*k2.*abs(x(i,j)).^q));
end
    end
end
end