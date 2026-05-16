%-------------------------------------------------------------------------%
%------               12th Setember 2022                        ---------%
%--------           Autopilot design with no wind                 --------%
%-------              Outside initial position                    --------%
%-------     Proposed vector field for circular path following    --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f1 = figure;
ax1 = axes;
f2 = figure;
ax2 = axes;
f3 = figure;
ax3 = axes;
% f4 = figure;
% ax4 = axes;
f5 = figure;
ax5 = axes;
f6 = figure;
ax6 = axes;
f7 = figure;
ax7 = axes;
f8 = figure;
ax8 = axes;
f9 = figure;
ax9 = axes;
f10 = figure;
ax10 = axes;
f11 = figure;
ax11 = axes;
f12 = figure;
ax12 = axes;
f13 = figure;
ax13 = axes;
f14 = figure;
ax14 = axes;
f15 = figure;
ax15 = axes;
f16 = figure;
ax16 = axes;
f17 = figure;
ax17 = axes;
%--------------- declaration of constant parameters ----------------------%
global rd Vgd k_kai k_kaidot x_c y_c flag1 t_flag1 flag2 t_flag2 
rd = 50;
Vgd = 10;
k_kai = 150;
k_kaidot = sqrt(k_kai);
% flag1 = true;
% flag2 = true;
%--------------------- standoff circle    --------------------------------%
x_c = 0;
y_c = 0;
theta = 0:2*pi/1000:2*pi;
x_sc = x_c + rd*cos(theta);
y_sc = y_c + rd*sin(theta);
%-------------------------------------------------------------------------%

%------------------- initial conditions   --------------------------------%

tspan = [0 300];
r_0 = 80;
gamma_0 = -deg2rad(90);
x0 = r_0*cos(gamma_0) + x_c;
y0 = r_0*sin(gamma_0) + y_c;
% k_prop = [0.0007 0.0013 0.0018 0.0023 0.0029 0.0034 0.0039 0.0045 0.005 0.006];
% k_prop = 0.05;
% k_prop = 1*10^(-4):(0.01-2*10^(-5))/19 :0.01;
k_prop = 0.0006;
flag1 = true(length(k_prop));
t_flag1 = zeros(1,length(k_prop));
flag2 = true(length(k_prop));
t_flag2 = zeros(1,length(k_prop));
for i = 1:length(k_prop)
% k_prop = 0.004;
if r_0 < rd
    kai0_prop(i) = gamma_0 + asin(1./(1+k_prop(i)*(r_0-rd).^2)) ;
else
    kai0_prop(i) = gamma_0 + pi - asin(1./(1+k_prop(i)*(r_0-rd).^2)) ;
end
gammadot0_prop(i) = (Vgd./r_0).*sin(kai0_prop(i) - gamma_0) ;
kaidot0_prop(i) = gammadot0_prop(i) - (2*k_prop(i)*Vgd.*(r_0-rd))./((1+k_prop(i)*(r_0-rd).^2).^2);

kai0_deg(i) = rad2deg(wrapToPi(kai0_prop(i)));

k_Nelson(i) = sqrt(k_prop(i)^2 * (r_0-rd)^2 + 2*k_prop(i));
kai0_Nelson(i) = gamma_0 + pi/2 + atan(k_Nelson(i).*(r_0 - rd)) ;
gammadot0_Nelson(i) = (Vgd./r_0).*sin(kai0_Nelson(i) - gamma_0) ;
kaidot0_Nelson(i) = gammadot0_Nelson(i) - (k_Nelson(i)^2*Vgd.*(r_0-rd))./((1+k_Nelson(i)^2*(r_0-rd).^2).^(3/2));
x_initial = [r_0;gamma_0;x0;y0;kai0_prop(i);kaidot0_prop(i)];
y_initial = [r_0;gamma_0;x0;y0;kai0_Nelson(i);kaidot0_Nelson(i)];

%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
% range  = -100:4:100;
% [X,Y] = meshgrid(range);
% r = sqrt((X-x_c).^2 + (Y-y_c).^2);
% gama = atan2((Y-y_c),(X-x_c));
% [xdot, ydot] = fun_vf(r,gama);
% 
% for i = 1:length(X)
%     for j = 1:length(X)
%         if r(i,j) < rd
%             X_in(i,j) = X(i,j);
%             Y_in(i,j) = Y(i,j);
%            xdot_in(i,j) = xdot(i,j);
%            ydot_in(i,j) = ydot(i,j);
%         else
%             X_out(i,j) = X(i,j);
%             Y_out(i,j) = Y(i,j);
%            xdot_out(i,j) = xdot(i,j);
%            ydot_out(i,j) = ydot(i,j);
%         end
%     end
% end
%-------------------------------------------------------------------------%
%---------------- ode solver without stopping event ----------------------%
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% 
% [t,x] = ode45(@(t,x) fun_proposed_nowind(t,x,k_prop(i),i), tspan, x_initial, options);
% [t1,y] = ode45(@(t1,y) fun_Nelson_nowind(t1,y,k_Nelson(i),i), tspan, y_initial, options);

%----------------- ode solver with stopping event ------------------------%

option_x=odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(t,x) stopping_circular_prop(t,x));
option_y=odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(t,y) stopping_circular_Nelson(t,y));

%% ode solver equation  %%
[t,x, te, ze, ie] = ode45(@(t,x) fun_proposed_nowind(t,x,k_prop(i),i), tspan, x_initial, option_x);
[t1,y, te1, ze1, ie1] = ode45(@(t1,y) fun_Nelson_nowind(t1,y,k_Nelson(i),i), tspan, y_initial, option_y);
%----------------------  parameters to plot   ----------------------------%
if (~isempty(te))
x_stop(i) = x(i,1);%% ze contain all states to apply event in any of the desired %% 
t_stop_prop(i) = te; %% time when selected state achieves the required criteria  %%
end

if (~isempty(te1))
y_stop(i) = y(i,1);%% ze contain all states to apply event in any of the desired %% 
t_stop_Nelson(i) = te1; %% time when selected state achieves the required criteria  %%
end

r_prop = x(:,1) ;
error_r = (r_prop - rd) ;

x_ini = x(1,3) + x_c ;
y_ini = x(1,4) + y_c ;
x_end = x(end,3) + x_c ;
y_end = x(end,4) + y_c ;

% for i = 1:length(t)
%     if x(i,1) < rd
%         kaid_prop(i,1) = x(i,2) + asin(1./(1+k_prop*(x(i,1)-rd).^2)) ;
%         gamma_dot_prop(i,1) = (Vgd./x(i,1)).*sin(kaid_prop(i,1) - x(i,2)) ;
%         kaid_dot_prop(i,1) = gamma_dot_prop(i,1) - 2*k_prop*Vgd.*(x(i,1)-rd)./((1+k_prop*(x(i,1)-rd).^2).^2);
%     else
%         kaid_prop(i,1) = x(i,2) + pi - asin(1./(1+k_prop*(x(i,1)-rd).^2)) ;
%         gamma_dot_prop(i,1) = (Vgd./x(i,1)).*sin(kaid_prop(i,1) - x(i,2)) ;
%         kaid_dot_prop(i,1) = gamma_dot_prop(i,1) - 2*k_prop*Vgd.*(x(i,1)-rd)./((1+k_prop*(x(i,1)-rd).^2).^2);
%     end
% end
% kappa_des_prop = kaid_dot_prop./Vgd ;

kappa_actual_prop = x(:,6)./Vgd ;
control_effort_prop(i) = trapz(t,(Vgd*x(:,6)).^2);
% control_input_prop = k_kai*(kaid_prop - x(:,5)) + k_kaidot*(kaid_dot_prop - x(:,6))  ;
[val_max_prop(i) ,ind_prop(i)] = max(abs(kappa_actual_prop));
x_max_prop(i) = x(ind_prop(i),1);

%----------------------- Nelson et al ------------------------------------%
r_Nelson = y(:,1);
error_Nelson = (r_Nelson - rd) ;
gamma_Nelson = y(:,2);
% kaid_Nelson = gamma_Nelson + pi/2 + atan(k_Nelson.*(r_Nelson - rd)) ;
% gamma_dot_Nelson = (Vgd./r_Nelson).*sin(kaid_Nelson - gamma_Nelson) ;
% kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(r_Nelson-rd))./((1+k_Nelson^2*(r_Nelson-rd).^2).^(3/2));
% kappa_des_Nelson = kaid_dot_Nelson./Vgd ;

kappa_actual_Nelson = y(:,6)./Vgd ;
control_effort_Nelson(i) = trapz(t1,(Vgd*y(:,6)).^2);
% % control_input_Nelson = k_kai*(kaid_Nelson - y(:,5)) + k_kaidot*(kaid_dot_Nelson - y(:,6))  ;
[val_max_Nelson(i), ind_Nelson(i)] = max(abs(kappa_actual_Nelson));
x_max_Nelson(i) = y(ind_Nelson(i),1);

percent_red_curvature(i) = ((val_max_Nelson(i) - val_max_prop(i))./val_max_Nelson(i)) *100;
percent_red_controleffort(i) = ((control_effort_Nelson(i) - control_effort_prop(i))./control_effort_Nelson(i)) *100;
percent_red_settlingtime(i) = -((t_stop_Nelson(i) - t_stop_prop(i))./t_stop_Nelson(i)) *100;

ratio_curvature(i) = val_max_prop(i)./val_max_Nelson(i) ;
ratio_controleffort(i) = control_effort_prop(i)./control_effort_Nelson(i) ;
ratio_settlingtime(i) = t_stop_prop(i)./t_stop_Nelson(i);
%-------------------------------------------------------------------------%
f1;
% if r_0 < rd
%     quiver(X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% else
%     quiver(X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% end
plot(ax1,x_sc,y_sc,'k','linewidth',2);hold(ax1,'on');
plot(ax1,x(:,3),x(:,4),'r','linewidth',2);hold(ax1,'on');
plot(ax1,y(:,3),y(:,4),'b','linewidth',2);hold(ax1,'on');
h2 = plot(ax1,x_c,y_c,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','magenta');hold(ax1,'on');grid(ax1,'on');
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h3 = plot(ax1,x_ini,y_ini,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green');hold(ax1,'on');grid(ax1,'on');
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(ax1,x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold(ax1,'on');
ax1.FontSize = 13;
xlabel(ax1,' $ x, $ \ (m)','Fontsize',18);
ylabel(ax1,'$ y, $ \ (m)','Fontsize',18);
% legend('Vector field','Circular \ orbit','Proposed','Nelson \ et al. 2007','Centre of orbit','Initial \ position','Final \ position','Fontsize',14);
legend(ax1,'Circular \ orbit','Proposed','Nelson \ et al. 2007','Fontsize',14)
axis(ax1,'equal')

f2;
plot(ax2,t,r_prop,'r','LineWidth',2);hold(ax2,'on');grid(ax2,'on');
plot(ax2,t1,r_Nelson,'b','LineWidth',2);hold(ax2,'on');grid(ax2,'on');
ax2.FontSize = 13;
xlabel(ax2,'Time, \ $t$, \ (s)','Fontsize',18)
ylabel(ax2,'Radial \ distance, \ $r$, \ (m)','Fontsize',18);
legend(ax2,'Proposed','Nelson \ et al. 2007','Fontsize',14);

f3;
plot(ax3,t,error_r,'r','LineWidth',2);hold(ax3,'on');grid(ax3,'on');
plot(ax3,t1,error_Nelson,'b','LineWidth',2);hold(ax3,'on');grid(ax3,'on');
ax3.FontSize = 13;
xlabel(ax3,'Time, \ $t$, \ (s)','Fontsize',18)
ylabel(ax3,'Radial \ error, \ $e = (r - r_d)$, \ (m)','Fontsize',18)
legend(ax3,'Proposed','Nelson \ et al. 2007','Fontsize',14);

% f4;
% plot(ax4,t,wrapToPi(x(:,2))*(180/pi),'r','linewidth',2);hold(ax4,'on');grid(ax4,'on');
% plot(ax4,t1,wrapToPi(y(:,2))*(180/pi),'b','linewidth',2);hold(ax4,'on');grid(ax4,'on');
% ax4.FontSize = 13;
% xlabel(ax4,'Time, \ $ t,$ \ (s)','Fontsize',18)
% ylabel(ax4,'Bearing \ angle, \ $ \gamma$, \ $ (^\circ) $','Fontsize',18)
% legend(ax4,'Proposed','Nelson \ et al. 2007','Fontsize',14);

f5;
plot(ax5,t,wrapToPi(x(:,5))*(180/pi),'r','linewidth',2);hold(ax5,'on');grid(ax5,'on');
plot(ax5,t1,wrapToPi(y(:,5))*(180/pi),'b','linewidth',2);hold(ax5,'on');grid(ax5,'on');
ax5.FontSize = 13;
xlabel(ax5,'Time,  $ t,$ \ (s)','Fontsize',18)
ylabel(ax5,'Course  angle, \ $ \chi$, \ deg.','Fontsize',18)
legend(ax5,'Proposed','Nelson \ et al. 2007','Fontsize',14);
% 
f6;
plot(ax6,t,kappa_actual_prop,'r','LineWidth',2);hold(ax6,'on');grid(ax6,'on');
plot(ax6,t1,kappa_actual_Nelson,'b','LineWidth',2);hold(ax6,'on');grid(ax6,'on');
ax6.FontSize = 13;
xlabel(ax6,'Time, \ $ t,$ \ (s)','Fontsize',18)
ylabel(ax6,'Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
legend(ax6,'Proposed','Nelson \ et al. 2007','Fontsize',14);

f7;
plot(ax7,r_prop,kappa_actual_prop,'r','LineWidth',2);hold(ax7,'on');grid(ax7,'on');
plot(ax7,r_Nelson,kappa_actual_Nelson,'b','LineWidth',2);hold(ax7,'on');grid(ax7,'on');
% h7 = plot(ax7,x_max_prop,val_max_prop,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','cyan'); hold(ax7,'on');grid(ax7,'on');
% h7.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(ax7,x_max_Nelson,val_max_Nelson,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','green'); hold on;grid on;
% plot(ax7,[x_max_prop(i) x_max_prop(i)],[0 val_max_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax7,'on');grid(ax7,'on');
% plot(ax7,[x_max_Nelson x_max_Nelson],[0 val_max_Nelson],'--','color',[0.75 0.75 0.75],'linewidth',2);hold on;grid on;
ax7.FontSize = 13;
xlabel(ax7,'Radial \ distance, \ $r$, \ (m)','Fontsize',18);
ylabel(ax7,'Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
% legend('Proposed','Nelson \ et al. 2007','Fontsize',14);
legend(ax7,'Proposed','Nelson \ et al. 2007','Fontsize',14);

f8;
plot(ax8,t,x(:,6),'r','LineWidth',2);hold(ax8,'on');grid(ax8,'on');
plot(ax8,t1,y(:,6),'b','LineWidth',2);hold(ax8,'on');grid(ax8,'on');
ax8.FontSize = 13;
xlabel(ax8,'Time,  $ t,$  s','Fontsize',18)
ylabel(ax8,'Course  rate,   $\dot{\chi}$,   rad./s ','Fontsize',18)
legend(ax8,'Proposed','Nelson  et al. 2007','Fontsize',14);

f9;
plot(ax9,r_prop,x(:,6),'r','LineWidth',2);hold(ax9,'on');grid(ax9,'on');
plot(ax9,r_Nelson,y(:,6),'b','LineWidth',2);hold(ax9,'on');grid(ax9,'on');
ax9.FontSize = 13;
xlabel(ax9,'Radial \ distance, \ $r$, \ (m)','Fontsize',18);
ylabel(ax9,'Course \ rate,  \ $\dot{\chi}$, \  rad./s ','Fontsize',18)
legend(ax9,'Proposed','Nelson \ et al. 2007','Fontsize',14);

% f10;
% plot(ax10,error_r,wrapToPi(x(:,5))*(180/pi),'r','linewidth',2);hold(ax10,'on');grid(ax10,'on');
% plot(ax10,error_Nelson,wrapToPi(y(:,5))*(180/pi),'b','linewidth',2);hold(ax10,'on');grid(ax10,'on');
% ax10.FontSize = 13;
% xlabel(ax10,'Radial \ error, \ $e = (r - r_d)$, \ (m)','Fontsize',18)
% ylabel(ax10,'Course \ angle, \ $ \chi$, \ deg.','Fontsize',18)
% legend(ax10,'Proposed','Nelson \ et al. 2007','Fontsize',14);

end

f11;
plot(ax11,k_prop,ratio_controleffort,'o-','MarkerFaceColor','green',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax11,'on');grid(ax11,'on');
ax11.FontSize = 18;
xlabel(ax11,' Proposed \ gain, \ $k_{s}$ ','Fontsize',23);
ylabel(ax11,' Control effort \ ratio','Fontsize',23');
axis(ax11,[0 k_prop(end) 0 2.5])

f12;
plot(ax12,k_prop,ratio_curvature,'s-','MarkerFaceColor','red',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax12,'on');grid(ax12,'on');
ax12.FontSize = 18;
xlabel(ax12,' Proposed \ gain, \ $k_{s}$ ','Fontsize',23);
ylabel(ax12,' Maximum \ curvature \ ratio','Fontsize',23');
axis(ax12,[0 k_prop(end) 0 1.5])

f13;
plot(ax13,k_prop,ratio_settlingtime,'s-','MarkerFaceColor','red',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax13,'on');grid(ax13,'on');
ax13.FontSize = 18;
xlabel(ax13,' Proposed \ gain, \ $k_{s}$ ','Fontsize',23);
ylabel(ax13,' Settling \ time \ ratio','Fontsize',23');
axis(ax13,[0 k_prop(end) 0 2])


f14;
plot(ax14,k_prop,percent_red_settlingtime,'d-','MarkerFaceColor','cyan',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax14,'on');grid(ax14,'on');
ax14.FontSize = 13;
xlabel(ax14,' Proposed \ gain, \ $k_{s}$ ','Fontsize',23);
ylabel(ax14,' $ \% $ Reduction \ in \ settling \ time ','Fontsize',23');
% axis(ax14,[0 k_prop(end) 0 50])

f15;
plot(ax15,k_prop,percent_red_curvature,'s-','MarkerFaceColor','red','MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax15,'on');grid(ax15,'on');
ax15.FontSize = 13;
xlabel(ax15,' Proposed \ gain, \ $k_{s}$ ','Fontsize',18);
ylabel(ax15,' $ \% $ Reduction \ in \ max. \ curvature ','Fontsize',18');
axis(ax15,[0 k_prop(end) -20 100])

f16;
plot(ax16,k_prop,percent_red_controleffort,'o-','MarkerFaceColor','green',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax16,'on');grid(ax16,'on');
ax16.FontSize = 13;
xlabel(ax16,' Proposed \ gain, \ $k_{s}$ ','Fontsize',18);
ylabel(ax16,' $ \% $ Reduction \ in \ control effort','Fontsize',18');

f17;
yyaxis left
plot(ax17,k_prop,kai0_deg,'s-','MarkerFaceColor','red','MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax17,'on');grid(ax17,'on');
axis(ax17,[0 k_prop(end) 0 90])
ax17.FontSize = 13;
ylabel(ax17,'Initial \ course \ angle, \ $ \chi_{o}$,\ ( deg.)','Fontsize',18');hold(ax17,'on');grid(ax17,'on');
yyaxis right
ylabel(ax17,' Nelson et al. \ gain, \ $ k$ ','Fontsize',18');hold(ax17,'on');grid(ax17,'on');
plot(ax17,k_prop,k_Nelson,'o-','MarkerFaceColor','green','MarkerEdgeColor','black', 'MarkerIndices', 1:length(k_prop),'linewidth',2);hold(ax17,'on');grid(ax17,'on');
axis(ax17,[0 k_prop(end) 0 1.5]);
xlabel(ax17,' Proposed \ gain, \ $k_{s}$ ','Fontsize',18);
legend(ax17,'$\chi_{o}$','$k$','Fontsize',17)


% function [xdot, ydot] = fun_vf(r,gama)
% global Vgd k_prop rd
% for i = 1:length(r)
%     for j = 1:length(r)
%         if r(i,j) < rd
%             kaid_proposed(i,j) = gama(i,j) + asin(1./(1+k_prop*(r(i,j)-rd).^2)) ;
%             %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
%             %    kaid_dot = gamma_dot + 2*k_prop*Vgd.*(r-rd)./((1+k_prop*(r-rd).^2).^2);
%         else
%             kaid_proposed(i,j) = gama(i,j) + pi - asin(1./(1+k_prop*(r(i,j)-rd).^2)) ;
%             %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
%             %    kaid_dot = gamma_dot - 2*k_prop*Vgd.*(r-rd)./((1+k_prop*(r-rd).^2).^2);
%         end
%     end
% end
% xdot  = Vgd*cos(kaid_proposed);
% ydot  = Vgd*sin(kaid_proposed);
% 
% end
function out = fun_proposed_nowind(t,x,k_prop,i)
global k_kai k_kaidot rd Vgd flag1 t_flag1 

if (flag1(i) == true) && (x(1)>(rd + 0.5))
    t_flag1(i) = t;
    flag1(i) = false;
end

if x(1) < rd
    kaid = x(2) + asin(1./(1+k_prop*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1)-rd))./((1+k_prop*(x(1) -rd).^2).^2);
else
    kaid = x(2) + pi - asin(1./(1+k_prop*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1) -rd))./((1+k_prop*(x(1) -rd).^2).^2);
end


out(1,1) = Vgd.*cos(x(5) - x(2));
out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
out(3,1) = Vgd*cos(x(5))  ;
out(4,1) = Vgd*sin(x(5))  ;
out(5,1) = x(6) ;
out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;

end
function out = fun_Nelson_nowind(t1,y,k_Nelson,i)
global k_kai k_kaidot rd Vgd flag2 t_flag2


if (flag2(i) == true) && (y(1)>(rd + 0.5))
    t_flag2(i) = t1;
    flag2(i) = false;
end


kaid_Nelson = y(2) + pi/2 + atan(k_Nelson.*(y(1) - rd)) ;
gamma_dot_Nelson = (Vgd./y(1) ).*sin(kaid_Nelson - y(2)) ;
kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(y(1) -rd))./((1+k_Nelson^2*(y(1) -rd).^2).^(3/2));

out(1,1) = Vgd.*cos(y(5) - y(2));
out(2,1) = (Vgd./y(1)).*sin(y(5) - y(2));
out(3,1) =  Vgd*cos(y(5))  ;
out(4,1) = Vgd*sin(y(5))  ;
out(5,1) = y(6) ;
out(6,1) = k_kai*(kaid_Nelson - y(5)) + k_kaidot*(kaid_dot_Nelson - y(6) )  ;


end