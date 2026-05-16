%-------------------------------------------------------------------------%
%------------            11th September 2022            ------------------%
%------------         Autopilot design with no wind     ------------------%
%------------         Straight line path following      ------------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%--------------- declaration of constant parameters ----------------------%
global  Vgd kv k_kai k_kaidot k2 figure_type_good
Vgd = 10;
kv = 5;
k_kai = 30;
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

% plot parameter
ax_fnt = 23;
lbl_fnt = 25;
ax_wdth = 3;
lgd_fnt = 17;

 hFig1 = figure('Position', [0 0 1920 1080]);
%  hAxes1 = gca;
 ax1 = gca;
 h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
xl = xline(ax1,0,'k','linewidth',3);hold(ax1,'on');
h1 = plot(ax1,x(1,1),x(1,2),'ob','MarkerSize',10,'MarkerFaceColor','blue');hold(ax1,'on');
 hPlot1 = plot(ax1,NaN,NaN,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4);
 ax1.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
% l= legend(ax1,[h xl h11],'Vector field','Desired path','UAV trajectory');hold on;
hold(ax1,'on');grid(ax1,'on');
axis(ax1,'equal');
ax1.XLim = [-100 100];
ax1.YLim = [-110 100];

 hFig2 = figure('Position', [0 0 1920 1080]);
%  hAxes2 = gca;
 ax2 = gca;
 hPlot2 = plot(ax2,NaN,NaN,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4);
% ax2.Position(2) = [1010 550 900 450];
% h21 = plot(ax2,t,x(:,1),'-r','LineWidth',4);hold on; grid on;
ax2.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax2,' $ t, $ s','Fontsize',lbl_fnt);
ylabel(ax2,'$ x, $ m','Fontsize',lbl_fnt);
hold(ax2,'on');grid(ax2,'on');
ax2.XLim = [0 t(end,1)];
ax2.YLim = [x(end,1) x(1,1)];
axis(ax2,'tight');


fig3 = figure('Position', [0 0 1920 1080]);
% subplot(2,2,3)
% f3.Position = [1010 50 900 450];
ax3 = gca;
% ax3.Position(3) = [1010 50 900 450];
% h31 = plot(ax3,t,x(:,3),'-r','LineWidth',4);hold on; grid on;
hPlot3= plot(ax3,NaN,NaN,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4);
ax3.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax3,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax3,'$$ \chi, $$ rad.','Fontsize',lbl_fnt);
hold(ax3,'on');grid(ax3,'on');
ax3.XLim = [0 t(end,1)];
ax3.YLim = [x(end,3) x(1,3)];
axis(ax3,'tight');

 for k = 1:length(x)
     xdata1 = get(hPlot1,'XData');
     ydata1 = get(hPlot1,'YData');
     set(hPlot1,'XData',[xdata1 x(k,1)],'YData',[ydata1 x(k,2)]);
     xdata = get(hPlot2,'XData');
     ydata = get(hPlot2,'YData');
     set(hPlot2,'XData',[xdata t(k,1)],'YData',[ydata x(k,1)]);
     ax2.XLim = [0 t(end,1)];
     ax2.YLim = [x(end,1) x(1,1)];
     xdata = get(hPlot3,'XData');
     ydata = get(hPlot3,'YData');
     set(hPlot3,'XData',[xdata t(k,1)],'YData',[ydata x(k,3)]);
     ax3.XLim = [0 t(end,1)];
     ax3.YLim = [x(end,3) x(1,3)];
     pause(0.005);
 end

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