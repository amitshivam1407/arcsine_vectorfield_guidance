%-------------------------------------------------------------------------%
%---------------------     15th December 2023      -----------------------%
%---- Sinusoidal path following using Adaptive Sliding Mode Control ------%
%-------------------------------------------------------------------------%

close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global A_amp omega tau k_s Va m
A_amp = 100;
Va = 20;
omega = 0.003;
tau = 100;
k_s = 0.09;
m = 0;
xd = -200:5:200;
% t = linspace(0,50,length(xd));
omega1 = omega*exp(-xd./tau);
% A_amp = -50*exp(-t./tau);
% yd =  m*xd + A_amp.*sin(omega1.*xd);
yd =  m*xd + A_amp.*sin(omega1.*xd);

range      = -200:5:200;
% range_y = linspace(-50,50,length(xd));
% range_x = linspace(-50,50,length(t));
[X,Y]      = meshgrid(range);
% tt = linspace(0,50,length(X));
omega11 = omega*exp(-X./tau);
% Yp =  m*X + A_amp.*sin(omega11.*X);
% dy_ddx_des =  m*1 -(A_amp/tau).*omega11.*cos(omega11.*X);
% Yp =  m*X + A_amp.*sin(omega11.*X);
% Yp = A_amp.*sin(omega11.*X);
Yp = A_amp.*sin(omega11.*X);
% dy_ddx_des =  m*1 -(A_amp/tau).*omega11.*cos(omega11.*X);
% dy_ddx_des =   -(A_amp/tau).*omega11.*cos(omega11.*X);
dy_ddx_des =  (-A_amp/tau)*omega11*cos(omega11.*X);
kai_p_des  = atan(dy_ddx_des);
e    =  Yp - Y;
% c = 0.3; a = 0.2; b = 0.8 ; 
% ks = a*e.^2 + b*e + c ;
% tspan = 
% x0 = 2;
% options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% [t,x] = ode45(@(t,x) get_fun(t,x) ,tspan, x0, options);
pp = (k_s*e)./sqrt(1 + (k_s*e).^2); 
qq =  1./sqrt(1 + (k_s*e).^2);
k1 = 10*k_s ;
rho = 1;
for i = 1:length(X)
    for j = 1:length(X)
        % kai_o(i,j) = pi/2 - asin(1./(1 + (k_s*e(i,j)).^2));
        kai_o(i,j) = pi/2 - asin(1./(1 + k_s*(e(i,j)).^2));

%           if e(i,j)< 0 
%           kaid(i,j) = (kai_p_des(1,j)) - kai_o(i,j)  ;
           
               kaid(i,j) = (1-rho)*sech(k_s*e(i,j))*(kai_p_des(1,j)) + rho*tanh(k_s*e(i,j))*kai_o(i,j)  ;

             % kaid(i,j) = (1-rho)*qq(i,j)*(kai_p_des(1,j)) + rho*pp(i,j)*kai_o(i,j)  ;
%          else
%              kaid(i,j) = (1-rho)*(kai_p_des(1,j)) + rho*(2/pi)*atan(k_s*e(i,j)).*kai_o(i,j) ;
%              kaid(i,j) = (1-rho)*(kai_p_des(1,j)) + rho*(2/pi)*atan(k_s*e(i,j)).*kai_o(i,j) ;
%             kaid(i,j) = (kai_p_des(1,j)) + kai_o(i,j)  ;
           
%           end
    end
end

Xdot       = Va*cos(kaid);
Ydot       = Va*sin(kaid);
Xdot_norm = Xdot./sqrt(Xdot.^2 + Ydot.^2) ;
Ydot_norm = Ydot./sqrt(Xdot.^2 + Ydot.^2) ;
% end
% [kaid]     = vf_sinusoidal(X,Y);

% Plot format control variables
    lw = 3;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 17;        % Axis font size
    lbl_fnt = ax_fnt+2; % Label font size
    leg_fnt = ax_fnt-1; % Legend font size
    ax_wdth = 3;        % Axis line width


figure(1)
quiver(X,Y,Xdot_norm,Ydot_norm,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
plot(xd,yd,'r','linewidth',3);
% plot(x(:,1),x(:,2),'m','linewidth',3);hold on;
% plot(x_ini,y_ini,'-o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','green'); hold on;
% plot(x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_wdth) ;
xlabel(ax1,' $ x, $ m','Fontsize',23);
ylabel(ax1,'$ y, $ m','Fontsize',23);
% legend(ax1,'Vector field','Desired path','UAV trajectory','','','Fontsize',leg_fnt);
axis(ax1, 'equal')

function [kaid] = vf_sinusoidal(X,Y)
global Va  k_s A_amp omega 
dydx = A_amp*omega*cos(omega.*X);
kai_p_des = atan(dydx);
e =  A_amp*sin(omega.*X) - Y;
kai_o = pi/2 - asin(1./(1 + k_s*(e).^2));

for i = 1:length(X)
    for j = 1:length(X)
        if e(i,j)< 0
            
            kaid(i,j) = (kai_p_des(1,j)) - kai_o(i,j)  ;
        else
            kaid(i,j) = (kai_p_des(1,j)) + kai_o(i,j) ;
        end
    end
end
end

function out = get_fun(t,x)
global c1 epsilon
out(1,1) = -c1*x(1) + epsilon
end