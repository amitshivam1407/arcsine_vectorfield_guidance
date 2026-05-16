%-------------------------------------------------------------------------%
%---------------------     31st March 2024      --------------------------%
%----------------  Sinusoidal path following with wind -------------------%
%-------------------------------------------------------------------------%

close all;clear all;clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%--------------- declaration of constant parameters ----------------------%
global  Vgd  kv wx wy k_kai k_kaidot k2 omega A
omega = 0.06; 
A = 5; 
Vgd = 10;
wx = -1;
wy = 2;
kv = 20;
k_kai = 50;  
k_kaidot = 10;  % 20 for better curvature profile  %%
k2 = .005;

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xd = -300:10:300;
yd = A*sin(omega.*xd);
%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range      = -300:10:300;
[X,Y]      = meshgrid(range);
dy_ddx_des = A*omega*cos(omega.*X);
kai_p_des  = atan(dy_ddx_des);
epsilon    =  A*sin(omega.*X)-Y;
[kaid]     = vf_sinusoidal(X,Y);
Xdot       = Vgd*cos(kaid);
Ydot       = Vgd*sin(kaid);
%-------------------  initial conditions   -------------------------------%
%% Initialization   ------------------------------%
tspan = [0 30];  % run time of ode
x0   = [100 -100 -250 -200];
y0   = [-150 110 200 -110];

% Va0 = 10;

for i = 1:length(x0)
% kai0(i) = fun_getkai(x0(i),y0(i));
[kai0(i),kaidot0(i)] = get_initialvalue(x0(i),y0(i));
kai0(i) = wrapToPi(kai0(i));

Va0 = 8;
x_init(:,i) = [x0(i);y0(i);kai0(i);kaidot0(i);Va0];

% tspan = [0 30];
% 
% xa_0 = -100;
% ya_0 = -90;
% 
% % [kai0, kaidot0] = fun_propkaid(xa_0,ya_0);
% kai0_deg = rad2deg(wrapToPi(kai0));
% 
% 
% x_initial = [xa_0;ya_0;(kai0);kaidot0; Va0];




%--------------- ode solver without stopping condition -------------------%
% 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode45(@(t,x)fun_sinusoidal(t,x) ,tspan, x_init(:,i),options);

% [t1,y] = ode45(@(t,y)fun_sinusoidal(t,y) ,tspan, x_initial,options);
%----------------  ode solver with stopping condition  -------------------%

% option_x = odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(t,x) stopping_sinusoidal(t,x));
% 
% [t,x, te, ze, ie] = ode45(@(t,x) fun_sinusoidal(t,x), tspan, x_initial, option_x);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------- Proposed parameter -------------------------------%
x_ini(:,i)     = x(1,1);
y_ini(:,i)     = x(1,2);
x_end(:,i)     = x(end,1);
y_end(:,i)     = x(end,2);

path_error =   A*sin(omega.*x(:,1)) - x(:,2);

% [kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x(:,1),x(:,2),x(:,3),x(:,5));
% kappa_des = kaidot_des./Vgd ;
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vad = sqrt(Vgd^2  + wx^2 + wy^2 - 2.*(Vgd*wx.*cos(kai_des) + Vgd*wy.*sin(kai_des)));
% Vg = sqrt((x(:,5).*cos(psi) + wx).^2 + (x(:,5).*sin(psi) + wy).^2) ;
% kaidot_actual = x(:,4);
% kappa_actual = kaidot_actual./Vg;
%---------------------- plotting figures ---------------------------------%

    %% Colors
    blue = [0 0.4470 0.7410];
    red = [0.8500 0.3250 0.0980];
    orange = [0.9290 0.6940 0.1250];
    violet = [0.4940 0.1840 0.5560];
    green = [0.4660 0.6740 0.1880];
    cyan = [0.3010 0.7450 0.9330];
    maroon = [0.6350 0.0780 0.1840];
    black = [0 0 0];

     color = [blue; red; violet; green];
    %%

% Plot format control variables
    lw = 3;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 17;        % Axis font size
    lbl_fnt = 20; % Label font size
    leg_fnt = 16; % Legend font size
    ax_lw = 3;        % Axis line width
%---------    
figure(1)
if i == 1
h11 = quiver(X,Y,Xdot,Ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12 = plot(xd,yd,'k','linewidth',lw);
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
h13 = plot(x(:,1),x(:,2),'Color',color(i,:),'linewidth',lw);hold on;
% plot(x(:,1),x(:,2),'--','Color',red,'linewidth',lw);hold on;
h14 =plot(x_ini(:,i),y_ini(:,i),'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold on;
h14.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
else
h13 = plot(x(:,1),x(:,2),'Color',color(i,:),'linewidth',3);hold on ;
% h1.Color = red;
% h11 =plot(x_ini(:,i),y_ini(:,i),'-o','Color',color(i,:),'MarkerSize',10); hold on;
h14 =plot(x_ini(:,i),y_ini(:,i),'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold on;
h14.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h22 = plot(x_end(:,i),y_end(:,i),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
% h22.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_lw) ;
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
% legend(ax1,'Vector field','Desired path','UAV trajectory','','','Fontsize',leg_fnt);
axis(ax1, 'equal')

figure(2)
plot(t,path_error,'Color',color(i,:),'LineWidth',lw);hold on;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',ax_lw) ;
xlabel(ax2,' $t$, s','Fontsize',lbl_fnt);
ylabel(ax2,'$$(f(x) - y)$$, m','Fontsize',lbl_fnt);
ax2.XLim = [0 30];
ax2.YLim = [-250 250];

% figure(3)
% plot(t,wrapToPi(kai_des)*(180/pi),'Color',blue,'linewidth',lw);hold on; grid on;
% plot(t,wrapToPi(x(:,3))*(180/pi),'--','Color',red,'linewidth',lw);hold on; grid on;
% plot(t,wrapToPi(kai_p)*(180/pi),'Color',maroon,'linewidth',lw);hold on; grid on;
% plot(t,wrapToPi(kai_o)*(180/pi),'-k','linewidth',lw);hold on; grid on;
% ax3 = gca;
% ax3.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',ax_lw) ;
% xlabel(ax3,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax3,' $ \chi$, deg.','Fontsize',lbl_fnt);
% legend(ax3,'Commanded','Achieved','Slope along the path','Shaping function','Fontsize',leg_fnt);
% 
% figure(4)
% plot(path_error,wrapToPi(kai_des)*(180/pi),'m','linewidth',lw);hold on; grid on;
% plot(path_error,wrapToPi(x(:,3))*(180/pi),'--m','linewidth',lw);
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',ax_lw) ;
% xlabel(ax4,'Tracking error,  m','Fontsize',lbl_fnt);
% ylabel(ax4,'$ \chi $, deg.','Fontsize',lbl_fnt);
% legend(ax4,'Commanded','Achieved','Fontsize',leg_fnt);
% grid on;
% 
% figure(5)
% % plot(t,kaidot_des,'b','LineWidth',lw);hold on;
% plot(t,kaidot_des,'Color',blue,'LineWidth',lw);hold on;
% % plot(t,x(:,4),'--m','LineWidth',lw);hold on; grid on;
% plot(t,x(:,4),'--','Color',red,'LineWidth',lw);hold on; grid on;
% ax5 = gca;
% ax5.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax5.XColor = 'black';         % Box horizontal lines' color
% ax5.YColor = 'black';         % Box vertical lines' color
% set(ax5,'linewidth',ax_lw) ;
% xlabel(ax5,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax5,' $\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt);
% legend(ax5,'Commanded','Achieved','Fontsize',leg_fnt);
% 
% 
% figure(6)
% % plot(t,kappa_des,'m','LineWidth',2);hold on; grid on;
% plot(t,kappa_actual,'Color',red,'LineWidth',lw);hold on; grid on;
% ax6 = gca;
% ax6.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax6.XColor = 'black';         % Box horizontal lines' color
% ax6.YColor = 'black';         % Box vertical lines' color
% set(ax6,'linewidth',ax_lw) ;
% xlabel(ax6,' $ t, $  s','Fontsize',lbl_fnt);
% ylabel(ax6,' $$\kappa$$, m $ ^{-1} $ ','Fontsize',lbl_fnt);
% % legend(ax6,'Commanded','Achieved','Fontsize',leg_fnt);
% grid on;
% 
% figure(7)
% % plot(path_error,kaidot_des,'m','LineWidth',2);hold on;
% plot(path_error,x(:,4),'-m','LineWidth',lw);hold on; grid on;
% ax7 = gca;
% ax7.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax7.XColor = 'black';         % Box horizontal lines' color
% ax7.YColor = 'black';         % Box vertical lines' color
% set(ax7,'linewidth',ax_lw) ;
% xlabel(ax7,'Tracking  error, m','Fontsize',lbl_fnt);
% ylabel(ax7,'$\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt);
% % legend(ax7,'Commanded','Achieved','Fontsize',leg_fnt);
% 
% 
% figure(8)
% % plot(path_error,kappa_des,'m','LineWidth',2);hold on; grid on;
% plot(path_error,kappa_actual,'-m','LineWidth',lw);hold on; grid on;
% ax8 = gca;
% ax8.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax8.XColor = 'black';         % Box horizontal lines' color
% ax8.YColor = 'black';         % Box vertical lines' color
% set(ax8,'linewidth',ax_lw) ;
% xlabel(ax8,'Tracking error,  m','Fontsize',lbl_fnt);
% ylabel(ax8,'Curvature $$\kappa$$,  m $ ^{-1} $ ','Fontsize',lbl_fnt);
% % legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt);
% grid on;
% 
% figure(9)
% plot(t,Vgd*ones(size(t)),'Color',blue,'linewidth',lw);hold on;
% plot(t,Vg,'--','Color',red,'linewidth',lw);
% ax9 = gca;
% ax9.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax9.XColor = 'black';         % Box horizontal lines' color
% ax9.YColor = 'black';         % Box vertical lines' color
% set(ax9,'linewidth',ax_lw) ;
% xlabel(ax9,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax9,' $ V_{\mathrm{g}}, $ m/s ','Fontsize',lbl_fnt);
% axis(ax9,[0 t(end,1) 0 15]);
% legend(ax9,'Commanded','Achieved','Fontsize',leg_fnt);
% grid on;
% % 
% figure(10)
% plot(t,Vad,'Color',blue,'linewidth',lw);hold on;grid on;
% plot(t,x(:,5),'--','Color',red,'linewidth',lw);grid on;
% ax10 = gca;
% ax10.FontSize = ax_fnt;
% box on                      % Switch on the box around the axis
% ax10.XColor = 'black';         % Box horizontal lines' color
% ax10.YColor = 'black';         % Box vertical lines' color
% set(ax10,'linewidth',ax_lw) ;
% xlabel(ax10,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax10,' $ V_{\mathrm{a}}, $ m/s ','Fontsize',lbl_fnt);
% axis(ax10,[0 t(end,1) 0 15]);
% legend(ax10,'Commanded','Achieved','Fontsize',leg_fnt);

end
Legend1 = cell(length(x0),1);
Legend2 = cell(length(x0),1);

for j = 1:length(x0)   
    Legend1{j} = sprintf('$$(x_{0},y_{0}) = (%d $m$,%d $m$)$$',x0(1,j),y0(1,j)) ;   
end
legend(ax1, Legend1,'Fontsize',leg_fnt,'NumColumns',1)

for j = 1:length(x0)
    Legend2{j} = sprintf('$$(x_{0},y_{0}) = (%d $m$,%d $m$)$$',x0(1,j),y0(1,j)) ;     
end
legend(ax2,Legend2,'Fontsize',leg_fnt,'NumColumns',1)


function out = fun_sinusoidal(t,x)
global Vgd k2 kv wx wy k_kai k_kaidot A omega

dydx = A*omega*cos(omega.*x(1));
kai_p =  atan(dydx);
epsilon =  A*sin(omega.*x(1)) - x(2);
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2)) ;



%  psid = kai + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
 
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));

Vg = sqrt((x(5).*cos(psi) + wx).^2 + (x(5).*sin(psi) + wy).^2) ;

out(1,1) = Vg*cos(x(3)) ;
out(2,1) = Vg*sin(x(3)) ;
if epsilon < 0
    kaid = (kai_p) - kai_o;
    Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = -(A*omega^2 *sin(omega*x(1))*xdot)./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = A*omega*cos(omega*x(1))*xdot - ydot;
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = (kai_p) + kai_o ;
    Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = -(A*omega^2 *sin(omega*x(1))*xdot)./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = A*omega*cos(omega*x(1))*xdot - ydot;
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end
% [kaid, kaid_dot] = fun_propkaid(x(1),x(2));


out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
out(5,1) = kv*(Vad - x(5));
end



function [kaid] = vf_sinusoidal(X,Y)
global Vgd  k2 A omega 
dydx = A*omega*cos(omega.*X);
kai_p_des = atan(dydx);
epsilon =  A*sin(omega.*X) - Y;
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2));

for i = 1:length(X)
    for j = 1:length(X)
        if epsilon(i,j)< 0
            
            kaid(i,j) = (kai_p_des(1,j)) - kai_o(i,j)  ;
        else
            kaid(i,j) = (kai_p_des(1,j)) + kai_o(i,j) ;
        end
    end
end
end

function [kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x1,x2,x3,x5)
global  Vgd k2 A wx wy omega 

for i = 1:length(x1)
dydx(i,:) = A.*omega*cos(omega*x1(i,:));
kai_p(i,:) = atan(dydx(i,:));
epsilon(i,:) = A.*sin(omega*x1(i,:)) - x2(i,:);
kai_o(i,:) = pi/2 - asin(1./(1 + k2*(epsilon(i,:)).^2)) ;

psi = x3(i,:) + asin((wx*sin(x3(i,:))-wy*cos(x3(i,:)))/x5(i,:));
 
Vg = sqrt(x5(i,:)^2  + wx^2 +wy^2 + 2*x5(i,:)*(wx*cos(psi) + wy*sin(psi)));


if epsilon(i,:) < 0
    kai_des(i,:)     = kai_p(i,:) - kai_o(i,:);
    xdot(i,:)        = Vg*cos(x3(i,:));
    ydot(i,:)        = Vg*sin(x3(i,:));
    kaip_dot(i,:)    = -(A*omega^2*sin(omega*x1(i,:)).*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
    factor_1(i,:)    = (2*k2)./((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
    epsilon_dot(i,:) = A*omega*cos(omega*x1(i,:))*xdot(i,:) - ydot(i,:);
    kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
    kaidot_des(i,:)  = kaip_dot(i,:) - kaio_dot(i,:);
else
    kai_des(i,:)     = kai_p(i,:) + kai_o(i,:);
    xdot(i,:)        = Vg*cos(x3(i,:));
    ydot(i,:)        = Vg*sin(x3(i,:));
    kaip_dot(i,:)    = -(A*omega^2*sin(omega*x1(i,:)).*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
    factor_1(i,:)    = (2*k2)/((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
    epsilon_dot(i,:) =  A*omega*cos(omega*x1(i,:))*xdot(i,:) - ydot(i,:);
    kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
    kaidot_des(i,:)  = kaip_dot(i,:) + kaio_dot(i,:);
end
end

end

function [value,isterminal,direction] = stopping_sinusoidal(t,x)
global  omega A

value(1) = ( (A*sin(omega*x(1))) - x(2)) + 0.02; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1) = 0; % negative direction(as R decreases from positive to zero d=-1;If R increases from negative to zero d=+1;d=0 implies no need of direction )

end

% function out = fun_griffiths_sinusoidal(t,y)
% global Vgd k2  k_kai k_kaidot A omega
% k = sqrt(k2*y(1).^2 + 2*k2);
% dydx = A*omega*cos(omega.*y(1));
% kai_p = atan(dydx);
% epsilon =  A*sin(omega.*y(1)) - y(2);
% epsilon_x = A*omega*cos(omega.*y(1));
% epsilon_y = -1;
% 
% kai_o =  atan(k.*(epsilon)) ;
% if epsilon < 0
%     kaid = (kai_p) - kai_o;
%     xdot = Vgd*cos(kaid);
%     ydot = Vgd*sin(kaid);
%     kaip_dot = (k.*(epsilon_x*xdot + epsilon_y.*ydot))./(1 + (tan(kai_p)).^2);
%     factor_1 = (k./(1 + (tan(kai_o)).^2));
%     epsilon_dot = A*omega*cos(omega*y(1))*xdot -ydot;
%     kaio_dot = factor_1*epsilon_dot;
%     kaid_dot = kaip_dot - kaio_dot;
% else
%     kaid = (kai_p) + kai_o ;
%     xdot = Vgd*cos(kaid);
%     ydot = Vgd*sin(kaid);
%      kaip_dot = (k.*(epsilon_x*xdot + epsilon_y.*ydot))./(1 + (tan(kai_p)).^2);
%     factor_1 = (k./(1 + (tan(kai_o)).^2));
%     epsilon_dot = A*omega*cos(omega*y(1))*xdot -ydot;
%     kaio_dot = factor_1*epsilon_dot;
%     kaid_dot = kaip_dot + kaio_dot;
% end
% 
% out(1,1) = Vgd*cos(y(3)) ;
% out(2,1) = Vgd*sin(y(3)) ;
% out(3,1) = y(4) ;
% out(4,1) = k_kai*(kaid - y(3)) + k_kaidot*(kaid_dot - y(4) )  ;
% 
% end
function [kai0,kaidot0] = get_initialvalue(xa_0,ya_0)
global Vgd k2 omega A

epsilon0 =  A*sin(omega*xa_0) - ya_0;
dydx0 = A*omega*cos(omega.*xa_0);
kai_p0 = atan(dydx0);

kai_o0 =  pi/2 - asin(1./(1 + k2*(epsilon0).^2));

if epsilon0 < 0
    kai0          = (kai_p0) - kai_o0;
    xdot0         = Vgd*cos(kai0);
    ydot0         = Vgd*sin(kai0);
    kaip_dot0     = -(A*omega^2*sin(omega*xa_0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k2)/((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0  =  A*omega*cos(omega*xa_0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaidot0       = kaip_dot0 - kaio_dot0;
else
    kai0          = kai_p0 + kai_o0 ;
    xdot0         = Vgd*cos(kai0);
    ydot0         = Vgd*sin(kai0);
    kaip_dot0     = -(A*omega^2*sin(omega*xa_0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k2)/((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0  = A*omega*cos(omega*xa_0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaidot0       = kaip_dot0 + kaio_dot0;
end
end