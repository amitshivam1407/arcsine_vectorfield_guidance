%-------------------------------------------------------------------------%
%---------------------     13th December 2023      -----------------------%
%------------------- Compare Vector field method -------------------------%
%-------------------------------------------------------------------------%

close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global A_amp omega  Va k_s
omega = 0.06; 
A_amp = 5; 
Va = 15;
k_s = 0.002;
%-------------------------------------------------------------------------%
xd = -105:.005:110;
yd = A_amp*sin(omega.*xd);
%----------------  vector field construction  ----------------------------%
range      = -100:8:100;
[X,Y]      = meshgrid(range);
dy_ddx_des = A_amp*omega*cos(omega.*X);
kai_p_des  = atan(dy_ddx_des);
e          =  A_amp*sin(omega.*X)-Y;
[kaid]     = vf_sinusoidal(X,Y);
Xdot       = Va*cos(kaid);
Ydot       = Va*sin(kaid);

[kaid_mod] = vf_sinusoidal_mod(X,Y);
Xdot_mod       = Va*cos(kaid_mod);
Ydot_mod       = Va*sin(kaid_mod);

% plot parameter
ax_fnt = 16;
lbl_fnt = 19;
ax_wdth = 3;
lgd_fnt = 14;

%% Colors
    blue = [0 0.4470 0.7410];
    red = [0.8500 0.3250 0.0980];
    orange = [0.9290 0.6940 0.1250];
    violet = [0.4940 0.1840 0.5560];
    green = [0.4660 0.6740 0.1880];
    cyan = [0.3010 0.7450 0.9330];
    maroon = [0.6350 0.0780 0.1840];
    black = [0 0 0];
    color = [red;blue;maroon];



figure(1)
% h13 = plot(xd,yd,'k','linewidth',3,'DisplayName','Straight line');hold on;
h13 = plot(xd,yd,'k','linewidth',3,'DisplayName','Sinusoidal path');hold on;
h11 = quiver(X,Y,Xdot,Ydot,'color',[0.75  0.75   0.75],'linewidth',1,'DisplayName','Vector field');hold on;
% h12 = quiver(X,Y,Xdot_mod,Ydot_mod,'color',[1  0   0],'linewidth',1,'DisplayName','Vector field Modified');hold on;

ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_wdth) ;
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
legend(ax1,'Fontsize',lgd_fnt);
axis(ax1, 'equal')

figure(2)
% h11 = quiver(X,Y,Xdot,Ydot,'color',[0.5  0.5   0.5],'linewidth',1,'DisplayName','Vector field');hold on;
h22 = quiver(X,Y,Xdot_mod,Ydot_mod,'color',[1  0   0],'linewidth',1,'DisplayName','Vector field Modified');hold on;
h23 = plot(xd,yd,'k','linewidth',3,'DisplayName','Desired curved path');
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',ax_wdth) ;
xlabel(ax2,' $ x, $ m','Fontsize',23);
ylabel(ax2,'$ y, $ m','Fontsize',23);
legend(ax2,'Fontsize',lgd_fnt);
axis(ax2, 'equal')


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

function [kaid_mod] = vf_sinusoidal_mod(X,Y)
global Va  k_s A_amp omega 
fx = A_amp*omega*cos(omega.*X);
fy = 1;
% kai_p_des = atan(fx);
kai_p_des = atan2(fx,fy);
e =  A_amp*sin(omega.*X) - Y;
kai_o = pi/2 - asin(1./(1 + k_s*(e).^2));

% del = (2/pi).*asin(1./(1 + k_s*(e).^2));

for i = 1:length(X)
    for j = 1:length(X)
        if e(i,j)< 0
            
%             kaid_mod(i,j) = del(1,j).*(kai_p_des(1,j)) - (1-del(1,j)).*kai_o(i,j)  ;
             kaid_mod(i,j) = (kai_p_des(1,j)) - kai_o(i,j)  ;
%              kaid_mod(i,j) = -sech(e(i,j)).*(kai_p_des(1,j)) + tanh(e(i,j)).*kai_o(i,j)  ;
        else
%             kaid_mod(i,j) = del(1,j).*(kai_p_des(1,j)) + (1-del(1,j)).*kai_o(i,j) ;
             kaid_mod(i,j) = (kai_p_des(1,j)) + kai_o(i,j)  ;
%             kaid_mod(i,j) = sech(e(i,j)).*(kai_p_des(1,j)) + tanh(e(i,j)).*kai_o(i,j)  ;
        end
    end
end
end