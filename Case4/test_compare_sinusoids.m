%-------------------------------------------------------------------------%
%---------------------     15th December 2023      -----------------------%
%---- Sinusoidal path following using Adaptive Sliding Mode Control ------%
%-------------------------------------------------------------------------%

close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%--------------- declaration of constant parameters ----------------------%
global A_amp omega  Va k_s W psi_w A_w A_phiw epsilon nabla zeta_3 zeta_1 zeta_2 c0 c1 c2 c3 ...
    wxs wys tex sx ux tey sy uy tez sz uz kai_ddot mu variance
omega = 0.06; 
A_amp = 5; 
W = 4;
A_w = 3;
A_phiw = pi;
Va = 10;
psi_w = 230*pi/180 ;
wxs = 4;
wys = 0;
mu  = 0;
variance = 0.5;
k_s = 0.005;
epsilon = 0.00002;
nabla = 10;
c1 = 5;

zeta_1 = 0.01;
zeta_2 = 0.01;
zeta_3 = 0.01;

c0 = 1;
c2 = 2;
c3 = 3;

%-------------------------------------------------------------------------%
xd = -100:.005:100;
yd = A_amp*sin(omega.*xd);
%----------------  vector field construction  ----------------------------%
range      = -100:5:100;
[X,Y]      = meshgrid(range);
dy_ddx_des = A_amp*omega*cos(omega.*X);
kai_p_des  = atan(dy_ddx_des);
e    =  A_amp*sin(omega.*X)-Y;
[kaid]     = vf_sinusoidal(X,Y);
Xdot       = Va*cos(kaid);
Ydot       = Va*sin(kaid);

%-------------------  initial conditions   -------------------------------%
a_arr  = [] ;
b_arr  = [] ;
c_arr  = [] ;
a_tran  = [] ;
b_tran  = [] ;
c_tran  = [] ;

tspan = [0 25];
t0=0;
tf=tspan(end);
dt=0.01;

% tf = tspan(end);
[u,v,Time] = DrydenDatagen(tf);


x_0 = -90;
y_0 = -60;
e10 = pi/2;
e20 = 0;
e10_hat = e10;
e20_hat = e20;
e30_hat = 0.05;
err_10 = e10 - e10_hat;
err_20 = e20 - e20_hat;

x60 = 0.1;
x70 = 0.2;
x80 = 0.3;
% [kai0,kaidot0] = get_initial(xa_0,ya_0);
% x_initial = [xa_0;ya_0;(kai0);kaidot0];

x0 = [x_0;y_0;e10;e20;e10_hat;e20_hat;err_10;err_20;x60;x70];
y0 = [x_0;y_0;e10;e20;e10_hat;e20_hat;err_10;err_20;x60;x70;x80];
z0 = [x_0;y_0;e10;e20;e10_hat;e20_hat;err_10;err_20;x60;x70];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode45(@(t,x) fun_sinusoidal_case1(t,x,u,v,Time) ,tspan, x0,options);
[t1,y] = ode45(@(t,y) fun_sinusoidal_case2(t,y,u,v,Time) ,tspan, y0, options);
[t2,z] = ode45(@(t,z) fun_sinusoidal_case3(t,z,u,v,Time) ,tspan, y0, options);

x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;

[kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x(:,1),x(:,2));
kai = x(:,3) + kai_des ;
% Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
V_g = get_Vg(t,kai,u,v,Time);
kappa_des = kaidot_des./V_g ;
e_path_x = A_amp*sin(omega.*x(:,1)) - x(:,2);

for i = 1:length(t)
    if abs(e_path_x(i,1))  <=0.1
%         a(i) = abs(e_path_x(i,1)) ;
        a_arr = [a_arr abs(e_path_x(i,1))];
    else
        a_tran =[a_tran abs(e_path_x(i,1))];
%         c_arr_col(j,:) = c_arr ;
    end
    
end

% rmse_z = rms(c);
tran_x = mean(a_tran);
ss_x = mean(a_arr);

e_path_y = A_amp*sin(omega.*y(:,1)) - y(:,2);
for i = 1:length(t1)
    if abs(e_path_y(i,1))  <=0.1
%         b(i) = abs(e_path_y(i,1)) ;
        b_arr = [b_arr abs(e_path_y(i,1))];
        else
        b_tran =[b_tran abs(e_path_y(i,1))]; 
    end
end
tran_y = mean(b_tran);
% mean_tran_y(gg,j) = tran_y ;
ss_y = mean(b_arr);
% rmsey(gg,j) = ss_y ;

%--- D = k0 + k1s -----

e_path_z = A_amp*sin(omega.*z(:,1)) - z(:,2);
for i = 1:length(t2)
    if abs(e_path_z(i,1))  <=0.1
%         c(i) = abs(e_path_z(i,1)) ;
        c_arr = [c_arr abs(e_path_z(i,1))];
        else
        c_tran =[c_tran abs(e_path_z(i,1))]; 
    end
end

tran_z = mean(c_tran);
% mean_tran_z(gg,j) = tran_z ;
ss_z = mean(c_arr);
% rmsez(gg,j) = ss_z ;

% Plot format control variables
    lw = 3;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 17;        % Axis font size
    lbl_fnt = ax_fnt+2; % Label font size
    leg_fnt = ax_fnt-1; % Legend font size
    ax_wdth = 3;        % Axis line width

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
h11 = quiver(X,Y,Xdot,Ydot,'color',[0.75  0.75   0.75],'linewidth',1,'DisplayName','Vector field');hold on;
h12 = plot(xd,yd,'k','linewidth',3,'DisplayName','Desired curved path');
h13 = plot(x(:,1),x(:,2),'-','linewidth',3,'Color',red,'DisplayName','$$d = k_{0} + k_{1}\hat{e}_{1} + k_{2}\hat{e}_{2}$$ ');hold on;
h14 = plot(y(:,1),y(:,2),'-','linewidth',3,'Color',blue,'DisplayName','$$d_{1} = k_{0} + k_{1}\hat{s} + k_{2}\hat{s}^{2}$$ ');hold on;
h15 = plot(z(:,1),z(:,2),'-','linewidth',3,'Color',maroon,'DisplayName','$$d_{2} = k_{0} + k_{1}\hat{s} $$' );hold on;
h16 = plot(x_ini,y_ini,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green'); hold on;
h16.Annotation.LegendInformation.IconDisplayStyle = 'off' ;
h17 = plot(x_end,y_end,'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold on;
h17.Annotation.LegendInformation.IconDisplayStyle = 'off' ;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_wdth) ;
xlabel(ax1,' $ x, $ m','Fontsize',23);
ylabel(ax1,'$ y, $ m','Fontsize',23);
legend(ax1,'Fontsize',leg_fnt);
axis(ax1, 'equal')

figure(2)
plot(t,e_path_x,'-','LineWidth',3,'Color',red,'DisplayName','$$d = k_{0} + k_{1}\hat{e}_{1} + k_{2}\hat{e}_{2}$$ ');hold on;grid on;
plot(t1,e_path_y,'-','LineWidth',3,'Color',blue,'DisplayName','$$d_{1} = k_{0} + k_{1}\hat{s} + k_{2}\hat{s}^{2}$$ ');hold on;grid on;
plot(t2,e_path_z,'-','LineWidth',3,'Color',maroon,'DisplayName','$$d_{2} = k_{0} + k_{1}\hat{s} $$' );hold on;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
% Outer box setup
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',ax_wdth)
xlabel(ax2,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax2,'Tracking error $$(f(x(t) - y(t))$$, m','Fontsize',lbl_fnt);
legend(ax2,'Fontsize',leg_fnt);


% figure(3)
% plot(t,wrapToPi(kai_des)*(180/pi),'-g','linewidth',3);hold on; grid on;
% plot(t,wrapToPi(kai)*(180/pi),'--b','linewidth',3);hold on; grid on;
% plot(t,wrapToPi(kai_p)*(180/pi),'-r','linewidth',3);hold on; grid on;
% plot(t,wrapToPi(kai_o)*(180/pi),'-k','linewidth',3);hold on; grid on;
% ax3 = gca;
% ax3.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',3) ;
% xlabel(ax3,' $ t, $ s','Fontsize',23);
% ylabel(ax3,' $ \chi$, deg.','Fontsize',23);
% legend(ax3,'Commanded','Achieved','Path-tangential angle','Offset angle','Fontsize',leg_fnt);
% 
% figure(4)
% % plot(t2,kai,'linewidth',3,'DisplayName',['Case ', num2str(j)]);grid on;hold on;
% % plot(t2,z(:,3),'linewidth',3,'DisplayName',['$$e_{\chi}$$ ', num2str(j)]);grid on;hold on;
% % plot(t2,z(:,5),'linewidth',3,'DisplayName',['$$\hat{e}_{\chi}$$  ', num2str(j)]);grid on;hold on;
% plot(t,x(:,3),'linewidth',3,'DisplayName','$$e_{\chi}$$ ');grid on;hold on;
% plot(t,x(:,5),'linewidth',3,'DisplayName','$$\hat{e}_{\chi}$$ ');grid on;hold on;
% % plot(t2,kaid,'r','linewidth',3,'DisplayName',' Commanded');grid on;hold on;
% % plot(t2,kai,'b','linewidth',3,'DisplayName',' Achieved');grid on;hold on;
% % plot(t2,z(:,7),'g','linewidth',3,'DisplayName','Observed');grid on;hold on;
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% % Outer box setup
% box on                        % Switch on the box around the axis
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',ax_wdth)
% xlabel(ax4,' $$ t, $$ s','Fontsize',lbl_fnt);
% ylabel(ax4,' $$ \hat{e}_{1}, $$ rad.','Fontsize',lbl_fnt);
% legend(ax4,'Commanded','Achieved','Observed','Fontsize',15)
% 
% figure(5)
% % plot(t2,z(:,4),'linewidth',3,'DisplayName',['$$e_{\dot{\chi}}$$ ', num2str(j)]);grid on;hold on;
% % plot(t2,z(:,6),'linewidth',3,'DisplayName',['$$\hat{e}_{\dot{\chi}}$$ ', num2str(j)]);grid on;hold on;
% plot(t,x(:,4),'linewidth',3,'DisplayName','$$e_{\dot{\chi}}$$ ');grid on;hold on;
% plot(t,x(:,6),'linewidth',3,'DisplayName','$$\hat{e}_{\dot{\chi}}$$ ');grid on;hold on;
% ax5 = gca;
% ax5.FontSize = ax_fnt;
% % Outer box setup
% box on                        % Switch on the box around the axis
% ax5.XColor = 'black';         % Box horizontal lines' color
% ax5.YColor = 'black';         % Box vertical lines' color
% set(ax5,'linewidth',ax_wdth)
% xlabel(ax5,'  $$ t, $$ s','Fontsize',lbl_fnt);
% ylabel(ax5,' $$ \hat{e}_{2}, $$ rad/s','Fontsize',lbl_fnt);
% legend(ax5,'Commanded','Achieved','Fontsize',15)
% % 

figure(6)
plot(tex,sx,'linewidth',3,'DisplayName','$$d = k_{0} + k_{1}\hat{e}_{1} + k_{2}\hat{e}_{2}$$');grid on;hold on;
plot(tey,sy,'linewidth',3,'DisplayName','$$d_{1} = k_{0} + k_{1}\hat{s} + k_{2}\hat{s}^{2}$$ ');grid on;hold on;
plot(tez,sz,'linewidth',3,'DisplayName','$$d_{2} = k_{0} + k_{1}\hat{s} $$');grid on;hold on;
ax6 = gca;
ax6.FontSize = ax_fnt;
% Outer box setup
box on                        % Switch on the box around the axis
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',ax_wdth)
xlabel(ax6,'  $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax6,' $$ s = c_{1}\hat{e}_{1} + \hat{e}_{2} $$ ','Fontsize',lbl_fnt);
legend(ax6,'Fontsize',leg_fnt)

figure(10)
plot(t,x(:,10),'linewidth',3,'DisplayName','$$k_{2x}$$ ');grid on;hold on;
plot(t1,y(:,11),'linewidth',3,'DisplayName','$$k_{2y}$$ ');grid on;hold on;
plot(t2,z(:,11),'linewidth',3,'DisplayName','$$k_{2z}$$ ');grid on;hold on;
ax10 = gca;
ax10.FontSize = ax_fnt;
% Outer box setup
box on                        % Switch on the box around the axis
ax10.XColor = 'black';         % Box horizontal lines' color
ax10.YColor = 'black';         % Box vertical lines' color
set(ax10,'linewidth',ax_wdth)
xlabel(ax10,'  $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax10,' $$ k_{2} $$ ','Fontsize',lbl_fnt);
legend(ax10,'$$d_{1}$$','$$d_{2}$$','$$d$$','Fontsize',15)
% 
figure(11)
plot(t,x(:,9),'linewidth',3,'DisplayName','$$k_{0x}$$ ');grid on;hold on;
plot(t1,y(:,9),'linewidth',3,'DisplayName','$$k_{0y}$$ ');grid on;hold on;
plot(t2,z(:,9),'linewidth',3,'DisplayName','$$k_{0z}$$ ');grid on;hold on;
ax11 = gca;
ax11.FontSize = ax_fnt;
% Outer box setup
box on                        % Switch on the box around the axis
ax11.XColor = 'black';         % Box horizontal lines' color
ax11.YColor = 'black';         % Box vertical lines' color
set(ax11,'linewidth',ax_wdth)
xlabel(ax11,'  $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax11,' $$ k_{0} $$ ','Fontsize',lbl_fnt);
legend(ax11,'$$d_{1}$$','$$d_{2}$$','$$d$$','Fontsize',15)
 
figure(12)
plot(t,x(:,10),'linewidth',3,'DisplayName','$$k_{1x}$$ ');grid on;hold on;
plot(t1,y(:,10),'linewidth',3,'DisplayName','$$k_{1y}$$ ');grid on;hold on;
plot(t2,z(:,10),'linewidth',3,'DisplayName','$$k_{1z}$$ ');grid on;hold on;
ax12 = gca;
ax12.FontSize = ax_fnt;
% Outer box setup
box on                        % Switch on the box around the axis
ax12.XColor = 'black';         % Box horizontal lines' color
ax12.YColor = 'black';         % Box vertical lines' color
set(ax12,'linewidth',ax_wdth)
xlabel(ax12,'  $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax12,' $$ k_{1} $$ ','Fontsize',lbl_fnt);
legend(ax12,'$$d_{1}$$','$$d_{2}$$','$$d$$','Fontsize',15)

function out = fun_sinusoidal_case1(t,x,u,v,Time)
global Va A_amp omega W psi_w A_w A_phiw c0 c2 c1 epsilon tex ux zeta_1 zeta_2 sx k_s wxs wys nabla kai_ddot

if t==0
    sx = [];
    ux = [];
    tex = [];
    kai_ddot = [];
end

A = A_w*cos(0.01*t);
psi_a = A_phiw*sin(0.01*t);

Wtx = A*cos(psi_a) ;
Wty = A*sin(psi_a) ;
Wcx = W*cos(psi_w) ;
Wcy = W*sin(psi_w) ;

wxgb=interp1(Time,u,t);
wygb=interp1(Time,v,t);
% 
% % % wg=[cos(x(3)),-sin(x(3));sin(x(3)),cos(x(3))]*[wxgb;wygb];
% 
wg = [wxgb;wygb];
% wg(1) = 0;
% wg(2) = 0;

Wgx = wxs+wg(1);
Wgy = wys+wg(2);

Wx = Wcx + Wtx + Wgx ;
Wy = Wcy + Wty + Wgy ;

dydx = A_amp*omega*cos(omega.*x(1));
kai_p =  atan(dydx);
e =  A_amp*sin(omega.*x(1)) - x(2);
kai_o = pi/2 - asin(1./(1 + k_s*(e).^2)) ;

if e < 0
    kaid = (kai_p) - kai_o;
    kai = x(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*x(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*x(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = (kai_p) + kai_o ;
    kai = x(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*x(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*x(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end

%---  system error dynamics -----
z_e1 = x(7)./(abs(x(7)) + epsilon);
k1 = 5; k2 = 10;
z1 = k1*(abs(x(7))^(1/2))*z_e1 ;
z2 = k2*z_e1 ; 



s = c1*x(5) + x(6);

z_sat = s./(abs(s) + epsilon);

del_bar = c0 + c2*abs(s)   ;

rho_est = x(9) + x(10)*abs(s) ;

u_x = -c1*x(6) - z2 -nabla*s - rho_est*z_sat ;

% kaiddot = ((Vg.^2).*(2*k_s*(3*k_s*e.^2 - 1))./(1 + k_s*e.^2).^3).*(cos(kai));
kaiddot = 0;

out(3,1) = x(4);
out(4,1) = u_x + del_bar - kaiddot ;

%--- observer error dynamics ----


out(5,1) = x(6) + z1;
out(6,1) = u_x + del_bar + z2 ;

%----- error between observer and system dynamics 

out(7,1) = x(8) - z1 ;
out(8,1) = -z2 - kaiddot;

out(9,1) = abs(s)- zeta_1*x(9);
out(10,1) = abs(s).^2 -  zeta_2*x(10);

kaid_ddot = u_x + del_bar - out(4,1);

ux = [ux u_x];
tex = [tex t];
sx = [sx s];
kai_ddot = [kai_ddot kaid_ddot];

end

function out = fun_sinusoidal_case2(t,y,u,v,Time)
% function out = fun_proposed_case2(t,y)
global Va A_amp omega k_s W psi_w epsilon nabla zeta_1 zeta_2 zeta_3 c0 c1 c2 c3 wxs wys sy tey uy kai_ddoty A_w A_phiw

if t==0
    sy = [];
    uy = [];
    tey = [];
    kai_ddoty = [];
end

A = A_w*cos(0.01*t);
psi_a = A_phiw*sin(0.01*t);

Wtx = A*cos(psi_a) ;
Wty = A*sin(psi_a) ;
Wcx = W*cos(psi_w) ;
Wcy = W*sin(psi_w) ;

wxgb=interp1(Time,u,t);
wygb=interp1(Time,v,t);

% % wg=[cos(x(3)),-sin(x(3));sin(x(3)),cos(x(3))]*[wxgb;wygb];

wg = [wxgb;wygb];

Wgx = wxs + wg(1);
Wgy = wys + wg(2);

Wx = Wcx + Wtx + Wgx ;
Wy = Wcy + Wty + Wgy ;


dydx = A_amp*omega*cos(omega.*y(1));
kai_p =  atan(dydx);
e =  A_amp*sin(omega.*y(1)) - y(2);
kai_o = pi/2 - asin(1./(1 + k_s*(e).^2)) ;

if e < 0
    kaid = (kai_p) - kai_o;
    kai = y(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*y(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*y(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = (kai_p) + kai_o ;
    kai = y(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*y(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*y(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end


kai = y(5) + kaid ;
Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai)).^2) ;

out(1,1) = Vg*cos(kai);
out(2,1) = Vg*sin(kai);


%---  system error dynamics -----
y_e1 = y(7)./(abs(y(7)) + epsilon);
k1 = 5; k2 = 10;
z1 = k1*(abs(y(7))^(1/2))*y_e1 ;
z2 = k2*y_e1 ; 

s = c1*y(5) + y(6);

% s = y(5);

z_sat = s./(abs(s) + epsilon);

del_bar = c0  + c2*abs(y(5)) + c3*abs(y(5)).^2;

rho_est = y(9) + y(10)*abs(y(5)) + y(11)*abs(y(5)).^2 ;

u_y = -c1*y(8) - z2 - nabla*s - rho_est*z_sat ;

kaiddot = 0;

% kaiddot = ((Vg.^2).*(2*k_s*(3*k_s*e.^2 - 1))./(1 + k_s*e.^2).^3).*(cos(kai));

out(3,1) = y(4);
out(4,1) = u_y + del_bar - kaiddot ;

%--- observer error dynamics ----


out(5,1) = y(6) + z1;
out(6,1) = u_y + del_bar + z2 ;

%----- error between observer and system dynamics 

out(7,1) = y(8) - z1 ;
out(8,1) = -z2 - kaiddot ;

out(9,1) = abs(y(5)) - zeta_1*y(9);
out(10,1) = abs(y(5)).^2 -  zeta_2*y(10);
out(11,1) = abs(y(5)).^3 - zeta_3*y(11) ;
 
kaid_ddot = u_y + del_bar - out(4,1);


tey = [tey t];
sy = [sy s];
uy = [uy u_y];
kai_ddoty = [kai_ddoty kaid_ddot];

end

function out = fun_sinusoidal_case3(t,z,u,v,Time)
global Va A_amp omega k_s W psi_w epsilon nabla zeta_3 zeta_1 zeta_2 c0 c1 c2 c3 wxs wys sz tez uz kai_ddotx  A_w A_phiw
if t==0
    sz = [];
    uz = [];
    tez = [];
    kai_ddotx = [];
end

A = A_w*cos(0.01*t);
psi_a = A_phiw*sin(0.01*t);

Wtx = A*cos(psi_a) ;
Wty = A*sin(psi_a) ;
Wcx = W*cos(psi_w) ;
Wcy = W*sin(psi_w) ;

wxgb=interp1(Time,u,t);
wygb=interp1(Time,v,t);



wg = [wxgb;wygb];
Wgx = wxs+wg(1);
Wgy = wys+wg(2);

Wx = Wcx + Wtx + Wgx ;
Wy = Wcy + Wty + Wgy ;

dydx = A_amp*omega*cos(omega.*z(1));
kai_p =  atan(dydx);
e =  A_amp*sin(omega.*z(1)) - z(2);
kai_o = pi/2 - asin(1./(1 + k_s*(e).^2)) ;

if e < 0
    kaid = (kai_p) - kai_o;
    kai = z(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*z(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*z(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = (kai_p) + kai_o ;
    kai = z(3) + kaid ;
    Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;
    out(1,1) = Vg*cos(kai);
    out(2,1) = Vg*sin(kai);
    kaip_dot = -(A_amp*omega^2 *sin(omega*z(1))*out(1,1))./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k_s)/((1 + k_s*(e).^2).*(sqrt(2*k_s + (k_s*e).^2)));
    epsilon_dot = A_amp*omega*cos(omega*z(1))*out(1,1) -out(2,1);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end

kai = z(5) + kaid ;
Vg = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai)).^2) ;

out(1,1) = Vg*cos(kai);
out(2,1) = Vg*sin(kai);

%---  system error dynamics -----
x_e1 = z(7)./(abs(z(7)) + epsilon);
k1 = 5; k2 = 10;
z1 = k1*(abs(z(7))^(1/2))*x_e1 ;
z2 = k2*x_e1 ; 

s = c1*z(5) + z(6);

z_sat = s./(abs(s) + epsilon);



del_bar = c0  + c2*abs(z(5) )  + c3*abs(z(6) );

rho_est = z(9) + z(10).*abs(z(5) ) + z(11).*abs(z(6) ) ;

u_z = -c1*z(6) - z2 -nabla*s - rho_est*z_sat ; % proportional-rate reaching law

% kaiddot = ((Vg.^2).*(2*k_s*(3*k_s*e.^2 - 1))./(1 + k_s*e.^2).^3).*(cos(kai));
kaiddot = 0;

out(3,1) = z(4);
out(4,1) = u_z + del_bar - kaiddot;

%--- observer error dynamics ----


out(5,1) = z(6) + z1;
out(6,1) = u_z + del_bar + z2 ;


%----- error between observer and system dynamics 

out(7,1) = z(8) - z1 ;
out(8,1) = -z2 - kaiddot ;

out(9,1) = abs(s)- zeta_1*z(9);
out(10,1) = abs(z(5))*abs(s) -  zeta_2*z(10);
out(11,1) = abs(z(6))*abs(s) -  zeta_3*z(11);

kaid_ddot = u_z + del_bar - out(4,1);

tez = [tez t];
sz = [sz s];
uz = [uz u_z];
kai_ddotx = [kai_ddotx kaid_ddot];

end

function [kai0,kaidot0] = get_initial(xa_0,ya_0)
global Va  k_s A_amp omega 
epsilon0 =  A_amp*sin(omega*xa_0) - ya_0;
dydx0 = A_amp*omega*cos(omega.*xa_0);
kai_p0 = atan(dydx0);
kai_o0 =  pi/2 - asin(1./(1 + k_s*(epsilon0).^2));

if epsilon0 < 0
    kai0          = (kai_p0) - kai_o0;
    xdot0         = Va*cos(kai0);
    ydot0         = Va*sin(kai0);
    kaip_dot0     = -(A_amp*omega^2*sin(omega*xa_0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k_s)/((1 + k_s*(epsilon0).^2).*(sqrt(2*k_s + (k_s*epsilon0).^2)));
    epsilon_dot0  =  A_amp*omega*cos(omega*xa_0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaidot0       = kaip_dot0 - kaio_dot0;
else
    kai0          = kai_p0 + kai_o0 ;
    xdot0         = Va*cos(kai0);
    ydot0         = Va*sin(kai0);
    kaip_dot0     = -(A_amp*omega^2*sin(omega*xa_0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k_s)/((1 + k_s*(epsilon0).^2).*(sqrt(2*k_s + (k_s*epsilon0).^2)));
    epsilon_dot0  = A_amp*omega*cos(omega*xa_0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaidot0       = kaip_dot0 + kaio_dot0;
end
end

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

function [u,v,Time]= DrydenDatagen(tend)
global mu variance
% This function generates the u and v component of velocity
% Input:
%      tend:- Time for which the wind turbulence data is required
%Output:-
%       u = x component of velocity
%       v = y component of velocity
%       Time = time vector

Gu=tf(20,[20 1]);
Gv=tf([1 (sqrt(3)/60)],[1 (1/400) (1/10)]);
del=0.01;
% ntot=101;
% tend=350;                            %  total time of simulation
t=0:del:tend;
ntot=(tend/del)+1;
% w=normrnd(0,1/sqrt(del),1,ntot);
w=normrnd(mu,variance,1,ntot);
% t=0:del:tend;
u=lsim(Gu,w,t);
v=lsim(Gv,w,t);
Time=t;

end


function [z,t] = rk45(f,z0,t0,tf,dt) %  function f, intial values x0 a vector 
                                     %  at t0, final time tf, step size dt
                              
 t = t0:dt:tf;    %vector of times
 nt = numel(t);   %no of elements in the time vector
 
 nx = numel(z0);  %no of states, (no of elemets in initial vector) 
 z = nan(nx,nt);  %matrix with nx rows and nt columns (not a number)


 z(:,1) = z0;     %first column of the matrix is initial x0
 for k = 1:nt-1   % At each step in the loop below, changed x(i) to x(:,i) 

     k1 = dt*f(t(k),z(:,k));
     k2 = dt*f(t(k) + dt/2, z(:,k) + k1/2);
     k3 = dt*f(t(k) + dt/2, z(:,k) + k2/2);
     k4 = dt*f(t(k) + dt, z(:,k) + k3);
     
     dx=(k1+2*k2+2*k3+k4)/6;
     z(:,k+1)=z(:,k)+dx;
     
 end
end

function [kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x1,x2)
global  Va k_s A_amp omega 

for i = 1:length(x1)
dydx(i,:) = A_amp.*omega*cos(omega*x1(i,:));
kai_p(i,:) = atan(dydx(i,:));
epsilon(i,:) = A_amp.*sin(omega*x1(i,:)) - x2(i,:);
kai_o(i,:) = pi/2 - asin(1./(1 + k_s*(epsilon(i,:)).^2)) ;
if epsilon(i,:) < 0
    kai_des(i,:)     = kai_p(i,:) - kai_o(i,:);
    xdot(i,:)        = Va*cos(kai_des(i,:));
    ydot(i,:)        = Va*sin(kai_des(i,:));
    kaip_dot(i,:)    = -(A_amp*omega^2*sin(omega*x1(i,:)).*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
    factor_1(i,:)    = (2*k_s)./((1 + k_s*(epsilon(i,:)).^2).*(sqrt(2*k_s + (k_s*epsilon(i,:)).^2)));
    epsilon_dot(i,:) = A_amp*omega*cos(omega*x1(i,:))*xdot(i,:) - ydot(i,:);
    kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
    kaidot_des(i,:)  = kaip_dot(i,:) - kaio_dot(i,:);
else
    kai_des(i,:)     = kai_p(i,:) + kai_o(i,:);
    xdot(i,:)        = Va*cos(kai_des(i,:));
    ydot(i,:)        = Va*sin(kai_des(i,:));
    kaip_dot(i,:)    = -(A_amp*omega^2*sin(omega*x1(i,:)).*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
    factor_1(i,:)    = (2*k_s)/((1 + k_s*(epsilon(i,:)).^2).*(sqrt(2*k_s + (k_s*epsilon(i,:)).^2)));
    epsilon_dot(i,:) =  A_amp*omega*cos(omega*x1(i,:))*xdot(i,:) - ydot(i,:);
    kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
    kaidot_des(i,:)  = kaip_dot(i,:) + kaio_dot(i,:);
end
end

end

function V_groundspeed = get_Vg(t,kai,u,v,Time)
global  Va k_s W psi_w A_w A_phiw wxs wys


A = A_w*cos(0.01*t);
psi_a = A_phiw*sin(0.01*t);

Wtx = A.*cos(psi_a) ;
Wty = A.*sin(psi_a) ;
Wcx = W*cos(psi_w) ;
Wcy = W*sin(psi_w) ;

wxgb=interp1(Time,u,t);
wygb=interp1(Time,v,t);
% 
% % % wg=[cos(x(3)),-sin(x(3));sin(x(3)),cos(x(3))]*[wxgb;wygb];
% 
wg = [wxgb;wygb];
% wg(1) = 0;
% wg(2) = 0;

Wgx = wxs+wg(1);
Wgy = wys+wg(2);

Wx = Wcx + Wtx + Wgx ;
Wy = Wcy + Wty + Wgy ;


V_groundspeed = Wx.*cos(kai) + Wy.*sin(kai) + sqrt(Va^2 - (Wx.*sin(kai) - Wy.*cos(kai))) ;

end