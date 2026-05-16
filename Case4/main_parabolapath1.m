%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Parabola path following    ----------------------%
%-------------------------------------------------------------------------%

close all;clear all; clc
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%--------------- declaration of constant parameters ----------------------%
global  Vgd k_kai k_kaidot k2 a h k
Vgd = 10;
k_kai = 100;
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
k2 = 0.002;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = .015;
h = 0;
k = -50;
xd = -100:.005:100;
yd = a*(xd - h).^2 + k;
%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 25];

xa_0 = -100;
ya_0 = -80;
epsilon0 = -((ya_0 - k) - a*(xa_0 - h)^2);
dydx0 = 2*a*(xa_0 - h);
kai_p0 = atan(dydx0);
kai_o0 = pi/2 - asin(1./(1 + k2*(epsilon0).^2));

if epsilon0 < 0
    kai0 = kai_p0 -  kai_o0 ;
    xdot0 = Vgd*cos(kai0);
    ydot0 = Vgd*sin(kai0);
    kaip_dot0 = (2*a*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10 = (2*k2)./((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0 = -(ydot0 - (2*a*(xa_0-h)*xdot0));
    kaio_dot0 = factor_10*epsilon_dot0;
    kaidot0 = kaip_dot0 - kaio_dot0;
else
    kai0 = kai_p0 + kai_o0 ;
    xdot0 = Vgd*cos(kai0);
    ydot0 = Vgd*sin(kai0);
    kaip_dot0 = (2*a*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10 = (2*k2)./((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0 = -(ydot0 - (2*a*(xa_0-h)*xdot0));
    kaio_dot0 = factor_10*epsilon_dot0;
    kaidot0 = kaip_dot0 + kaio_dot0;
end
kai0_deg = rad2deg(wrapToPi(kai0));
% Va0 = sqrt(Vgd^2  + wx^2 + wy^2 - 2.*(Vgd*wx.*cos(kai0) + Vgd*wy.*sin(kai0)));
% Va0 = 8;
x_initial = [xa_0;ya_0;kai0;kaidot0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -100:5:100;
[X,Y] = meshgrid(range);
[kaid] = vf_parabola(X,Y);
Xdot  = Vgd*cos(kaid);
Ydot  = Vgd*sin(kaid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%----------------- ode solver   ------------------------------------------%
[t,x] = ode45(@(t,x)fun_parabola(t,x) ,tspan, x_initial,options);
%
% %----------------------  parameters to plot   ----------------------------%
% % % psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% % % Vg_act = sqrt(x(:,5).^2 + wx^2 + wy^2 + 2*wx*x(:,5).*cos(psi) + 2*wy*x(:,5).*sin(psi));
x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;
path_error =  -((x(:,2) - k) - a.*(x(:,1) - h).^2);
[kai_p, kai_o, kai_des, kaidot_des] = fun_propkaid(x(:,1),x(:,2));
kappa_des = kaidot_des./Vgd ;
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vad = sqrt(Vgd^2  + wx^2 + wy^2 - 2.*(Vgd*wx.*cos(kai_des) + Vgd*wy.*sin(kai_des)));
% Vg = sqrt((x(:,5).*cos(psi) + wx).^2 + (x(:,5).*sin(psi) + wy).^2) ;
kappa_actual = x(:,4)./Vgd ;
%
% %---------------------- plotting figures ---------------------------------%
% Plot format control variables
lw = 3;            % Line width
ms = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = ax_fnt+2; % Label font size
leg_fnt = ax_fnt-1; % Legend font size
ax_lw = 3;        % Axis line width

figure(1)
quiver(X,Y,Xdot,Ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
plot(xd,yd,'k','linewidth',2);
plot(x(:,1),x(:,2),'m','linewidth',3);hold on;
plot(x_ini,y_ini,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green'); hold on;
plot(x_end,y_end,'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold on;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3) ;
xlabel(ax1,' $ x, $ m','Fontsize',23);
ylabel(ax1,'$ y, $ m','Fontsize',23);
legend(ax1,'Vector field','Desired path','UAV trajectory','','','Fontsize',leg_fnt);
axis(ax1, 'equal')

figure(2)
plot(t,path_error,'m','LineWidth',2);hold on;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3) ;
xlabel(ax2,' $t$, s','Fontsize',lbl_fnt);
ylabel(ax2,'Tracking error $$(f(x(t) - y(t))$$, m','Fontsize',lbl_fnt);

figure(3)
plot(t,wrapToPi(kai_des)*(180/pi),'-g','linewidth',3);hold on; grid on;
plot(t,wrapToPi(x(:,3))*(180/pi),'--b','linewidth',3);hold on; grid on;
plot(t,wrapToPi(kai_p)*(180/pi),'-r','linewidth',3);hold on; grid on;
plot(t,wrapToPi(kai_o)*(180/pi),'-k','linewidth',3);hold on; grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3) ;
xlabel(ax3,' $ t, $ s','Fontsize',23);
ylabel(ax3,' $ \chi$, deg.','Fontsize',23);
legend(ax3,'Commanded','Achieved','Tangential angle','Shaping function angle','Fontsize',leg_fnt);

figure(4)
plot(path_error,wrapToPi(kai_des)*(180/pi),'m','linewidth',2);hold on; grid on;
plot(path_error,wrapToPi(x(:,3))*(180/pi),'--m','linewidth',2);
ax4 = gca;
ax4.FontSize = 16;
box on                      % Switch on the box around the axis
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',1.5) ;
xlabel(ax4,'Tracking error,  m','Fontsize',lbl_fnt);
ylabel(ax4,'$ \chi $, deg.','Fontsize',lbl_fnt);
legend(ax4,'Commanded','Achieved','Fontsize',leg_fnt);
grid on;

figure(5)
plot(t,kaidot_des,'m','LineWidth',2);hold on;
plot(t,x(:,4),'--m','LineWidth',2);hold on; grid on;
ax5 = gca;
ax5.FontSize = 16;
box on                      % Switch on the box around the axis
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',1.5) ;
xlabel(ax5,' $ t, $ s','Fontsize',lbl_fnt);
ylabel(ax5,' $\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt);
% legend(ax5,'Commanded','Achieved','Fontsize',leg_fnt);
grid on;

figure(6)
% plot(t,kappa_des,'-m','LineWidth',2);hold on; grid on;
plot(t,kappa_actual,'-m','LineWidth',3);hold on; grid on;
ax6 = gca;
ax6.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',3) ;
xlabel(ax6,' $ t, $  s','Fontsize',23);
ylabel(ax6,'Curvature $$\kappa$$, m $ ^{-1} $ ','Fontsize',23);
% legend(ax6,'Commanded','Achieved','Fontsize',leg_fnt);
grid on;

figure(7)
plot(path_error,kaidot_des,'-m','LineWidth',2);hold on;
plot(path_error,x(:,4),'--m','LineWidth',2);hold on; grid on;
ax7 = gca;
ax7.FontSize = 16;
box on                      % Switch on the box around the axis
ax7.XColor = 'black';         % Box horizontal lines' color
ax7.YColor = 'black';         % Box vertical lines' color
set(ax7,'linewidth',1.5) ;
xlabel(ax7,'Tracking  error, m','Fontsize',lbl_fnt);
ylabel(ax7,'$\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt);
% legend(ax7,'Commanded','Achieved','Fontsize',leg_fnt);
grid on;
%
figure(8)
plot(path_error,kappa_des,'m','LineWidth',2);hold on; grid on;
plot(path_error,kappa_actual,'--m','LineWidth',2);hold on; grid on;
ax8 = gca;
ax8.FontSize = 16;
box on                      % Switch on the box around the axis
ax8.XColor = 'black';         % Box horizontal lines' color
ax8.YColor = 'black';         % Box vertical lines' color
set(ax8,'linewidth',1.5) ;
xlabel(ax8,'Tracking error,  m','Fontsize',lbl_fnt);
ylabel(ax8,'Curvature $$\kappa$$,  m $ ^{-1} $ ','Fontsize',lbl_fnt);
% legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt);
grid on;
%
% figure(9)
% plot(t,Vgd*ones(size(t)),'r','linewidth',2);hold on;
% plot(t,Vg,'b','linewidth',2);
% ax9 = gca;
% ax9.FontSize = 16;
% box on                      % Switch on the box around the axis
% ax9.XColor = 'black';         % Box horizontal lines' color
% ax9.YColor = 'black';         % Box vertical lines' color
% set(ax9,'linewidth',1.5) ;
% xlabel(ax9,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax9,' $ V_{\mathrm{g}}, $ m/s ','Fontsize',lbl_fnt);
% % axis(ax9,[0 t(end,1) 0 30]);
% legend(ax9,'Desired','Actual','Fontsize',leg_fnt);
% grid on;
%
% figure(10)
% plot(t,Vad,'r','linewidth',2);hold on;grid on;
% plot(t,x(:,5),'b','linewidth',2);grid on;
% ax10 = gca;
% ax10.FontSize = 16;
% box on                      % Switch on the box around the axis
% ax10.XColor = 'black';         % Box horizontal lines' color
% ax10.YColor = 'black';         % Box vertical lines' color
% set(ax10,'linewidth',1.5) ;
% xlabel(ax10,' $ t, $ s','Fontsize',lbl_fnt);
% ylabel(ax10,' $ V_{\mathrm{a}}, $ m/s ','Fontsize',lbl_fnt);
% % axis(ax10,[0 t(end,1) 0 30]);
% legend(ax10,'Desired','Actual','Fontsize',leg_fnt);


function out = fun_parabola(t,x)
global Vgd k2 k_kai k_kaidot a h  k

dydx = 2*a.*(x(1) - h);
kai_p = atan(dydx);
epsilon = -((x(2) - k) - a.*(x(1) - h)^2);
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2)) ;
if epsilon< 0
    kaid = kai_p - kai_o;
    xdot = Vgd*cos(kaid);
    ydot = Vgd*sin(kaid);
    kaip_dot = (2*a*xdot)./(1 + (tan(kai_p))^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = -(ydot - 2*a*(x(1)-h)*xdot);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = kai_p + kai_o;
    xdot = Vgd*cos(kaid);
    ydot = Vgd*sin(kaid);
    kaip_dot = (2*a*xdot)./(1 + (tan(kai_p))^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = -(ydot - 2*a*(x(1)-h)*xdot);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end

% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%
%  psid = kai + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
%
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3)) ;
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
% out(5,1) = kv*(Vad - x(5));
end



function [kaid] = vf_parabola(X,Y)
global Vgd  k2 a h k
% dydx = X./(2*a);
dydx = 2*a*(X-h);
kai_p_des = atan(dydx);
% epsilon =  Y - (X.^2)./(4*a);
epsilon =  -((Y-k) - a*(X-h).^2) ;
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2));
for i = 1:length(X)
    for j = 1:length(X)
        if epsilon(i,j)< 0
            kaid(i,j) = (kai_p_des(i,j)) - kai_o(i,j)  ;
        else
            kaid(i,j) = (kai_p_des(i,j)) + kai_o(i,j)  ;
        end
    end
end
end

function [kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x1,x2)
global  Vgd k2 a h k
dydx = 2*a.*(x1 - h);
kai_p = atan(dydx);
epsilon = -((x2 - k) - a.*(x1 - h).^2);
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2)) ;
for i = 1:length(x1)
    if epsilon(i,:)< 0
        kai_des(i,:)     = kai_p(i,:) - kai_o(i,:);
        xdot(i,:)        = Vgd*cos(kai_des(i,:));
        ydot(i,:)        = Vgd*sin(kai_des(i,:));
        kaip_dot(i,:)    = (2*a*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
        factor_1(i,:)    = (2*k2)/((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
        epsilon_dot(i,:) = -(ydot(i,:) - 2*a*(x1(i,:) - h)*xdot(i,:));
        kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
        kaidot_des(i,:)  = kaip_dot(i,:) - kaio_dot(i,:);
    else
        kai_des(i,:)     = kai_p(i,:) + kai_o(i,:);
        xdot(i,:)        = Vgd*cos(kai_des(i,:));
        ydot(i,:)        = Vgd*sin(kai_des(i,:));
        kaip_dot(i,:)    = (2*a*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
        factor_1(i,:)    = (2*k2)/((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
        epsilon_dot(i,:) = -(ydot(i,:) - 2*a*(x1(i,:) - h)*xdot(i,:));
        kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
        kaidot_des(i,:)  = kaip_dot(i,:) + kaio_dot(i,:);
    end
end

end
