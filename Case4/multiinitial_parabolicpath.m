%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Parabola path following    ----------------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%--------------- declaration of constant parameters ----------------------%
global  Vgd kv wx wy k_kai k_kaidot k2 a h k
Vgd = 12;
wx = -1;
wy = 2;
kv = 20;
k_kai = 100;
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
k2 = 0.001;

%-------------------------------------------------------------------------%
a = .01;
h = 0;
k = -100;
xd = -170:.005:170;
yd = a*(xd - h).^2 + k;
%------------------ initial conditions   ---------------------------------%

tspan = [0 30];

x0   = [70 -50 10 -180];
y0   = [-170 130 5 -50];

Va0 = 10;

for i = 1:length(x0)
% kai0(i) = fun_getkai(x0(i),y0(i));
[kai0(i),kaidot0(i)] = get_initialvalue(x0(i),y0(i));
kai0(i) = wrapToPi(kai0(i));


x_init(:,i) = [x0(i);y0(i);kai0(i);kaidot0(i);Va0];

%-------------------------------------------------------------------------%


%----------------  vector field construction  ----------------------------%
range_x  = -200:5:200;
range_y  = linspace(-10,150,length(range_x));
[X,Y] = meshgrid(range_x,range_x);
[kaid] = vf_parabola(X,Y);
Xdot  = Vgd*cos(kaid);
Ydot  = Vgd*sin(kaid);


%-------------------------------------------------------------------------%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%----------------- ode solver   ------------------------------------------%
[t,x] = ode45(@(t,x)fun_parabola_withwind(t,x) ,tspan, x_init(:,i),options);
%
% %----------------------  parameters to plot   ----------------------------%
% % % psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% % % Vg_act = sqrt(x(:,5).^2 + wx^2 + wy^2 + 2*wx*x(:,5).*cos(psi) + 2*wy*x(:,5).*sin(psi));

x_ini(:,i)     = x(1,1);
y_ini(:,i)     = x(1,2);
x_end(:,i)     = x(end,1);
y_end(:,i)     = x(end,2);


path_error =  -((x(:,2) - k) - a.*(x(:,1) - h).^2);

% [kai_p, kai_o, kai_des, kaidot_des] = fun_propkaid(x(:,1),x(:,2),x(:,3),x(:,5));
% kappa_des = kaidot_des./Vgd ;
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vad = sqrt(Vgd^2  + wx^2 + wy^2 - 2.*(Vgd*wx.*cos(kai_des) + Vgd*wy.*sin(kai_des)));
% Vg = sqrt((x(:,5).*cos(psi) + wx).^2 + (x(:,5).*sin(psi) + wy).^2) ;
% kappa_actual = x(:,4)./Vg ;

%------------------------ plotting figures -------------------------------%
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
% legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt);
% 
% figure(9)
% plot(t,Vgd*ones(size(t)),'Color',blue,'linewidth',lw);hold on;grid on;
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
% 
% 
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

function out = fun_parabola_withwind(t,x)
global Vgd k2 kv wx wy k_kai k_kaidot a h  k

dydx = 2*a.*(x(1) - h);
kai_p = atan(dydx);
epsilon = -((x(2) - k) - a.*(x(1) - h)^2);
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2)) ;

psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
 
Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));

out(1,1) = Vg*cos(x(3)) ;
out(2,1) = Vg*sin(x(3)) ;

if epsilon< 0
    kaid = kai_p - kai_o;
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = (2*a*xdot)./(1 + (tan(kai_p))^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = -(ydot - 2*a*(x(1)-h)*xdot);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = kai_p + kai_o;
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = (2*a*xdot)./(1 + (tan(kai_p))^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = -(ydot - 2*a*(x(1)-h)*xdot);
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end
 Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%
%  psid = kai + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
%
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


% out(1,1) = Vgd*cos(x(3));
% out(2,1) = Vgd*sin(x(3)) ;
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
out(5,1) = kv*(Vad - x(5));
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

function [kai_p,kai_o,kai_des, kaidot_des] = fun_propkaid(x1,x2,x3,x5)
global  Vgd k2 wx wy  a h k



for i = 1:length(x1)


   dydx(i,:) = 2*a.*(x1(i,:) - h);
kai_p(i,:) = atan(dydx(i,:));
epsilon(i,:) = -((x2(i,:) - k) - a.*(x1(i,:) - h).^2);
kai_o(i,:) = pi/2 - asin(1./(1 + k2*(epsilon(i,:)).^2)) ;

psi = x3(i,:) + asin((wx*sin(x3(i,:))-wy*cos(x3(i,:)))/x5(i,:));
 
Vg = sqrt(x5(i,:)^2  + wx^2 +wy^2 + 2*x5(i,:)*(wx*cos(psi) + wy*sin(psi)));

    if epsilon(i,:)< 0
        kai_des(i,:)     = kai_p(i,:) - kai_o(i,:);
        xdot(i,:)        = Vg*cos(x3(i,:));
        ydot(i,:)        = Vg*sin(x3(i,:));
        kaip_dot(i,:)    = (2*a*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
        factor_1(i,:)    = (2*k2)/((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
        epsilon_dot(i,:) = -(ydot(i,:) - 2*a*(x1(i,:) - h)*xdot(i,:));
        kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
        kaidot_des(i,:)  = kaip_dot(i,:) - kaio_dot(i,:);
    else
        kai_des(i,:)     = kai_p(i,:) + kai_o(i,:);
        xdot(i,:)        = Vg*cos(x3(i,:));
        ydot(i,:)        = Vg*sin(x3(i,:));
        kaip_dot(i,:)    = (2*a*xdot(i,:))./(1 + (tan(kai_p(i,:)))^2);
        factor_1(i,:)    = (2*k2)/((1 + k2*(epsilon(i,:)).^2).*(sqrt(2*k2 + (k2*epsilon(i,:)).^2)));
        epsilon_dot(i,:) = -(ydot(i,:) - 2*a*(x1(i,:) - h)*xdot(i,:));
        kaio_dot(i,:)    = factor_1(i,:)*epsilon_dot(i,:);
        kaidot_des(i,:)  = kaip_dot(i,:) + kaio_dot(i,:);
    end
end

end

function [value,isterminal,direction] = stopping_parabolic(t,x)
global  omega A

value(1) = ( (A*sin(omega*x(1))) - x(2)) + 0.02; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1) = 0; % negative direction(as R decreases from positive to zero d=-1;If R increases from negative to zero d=+1;d=0 implies no need of direction )

end


function [kai0,kaidot0] = get_initialvalue(xa_0,ya_0)
global  Vgd kv wx wy k_kai k_kaidot k2 a h k

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
end