%-------------------------------------------------------------------------%
%------               12th Setember 2022                        ---------%
%--------           Autopilot design with no wind                 --------%
%-------              Outside initial position                    --------%
%-------     Proposed vector field for circular path following    --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%%%%%%%%%%%%%%%% declaration of constant parameters %%%%%%%%%%%%%%%%%%%%%%%
global rd Vgd k_kai k_kaidot k_prop k_Nelson x_c y_c 
rd = 100;
Vgd = 5;
k_kai = 200;
k_kaidot = sqrt(k_kai);

%%%%%%%%%%%%%%%%%%%%%% standoff circle    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_c = 0;
y_c = 0;
theta = 0:2*pi/1000:2*pi;
x_sc = x_c + rd*cos(theta);
y_sc = y_c + rd*sin(theta);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 50];

r_0 = 5;
gamma_0 = -deg2rad(90);
x0 = r_0*cos(gamma_0) + x_c;
y0 = r_0*sin(gamma_0) + y_c;

k_prop = 0.005;
if r_0 < rd
    kai0_prop = gamma_0 + asin(1./(1+k_prop*(r_0-rd).^2)) ;
else
    kai0_prop = gamma_0 + pi - asin(1./(1+k_prop*(r_0-rd).^2)) ;
end
gammadot0_prop = (Vgd./r_0).*sin(kai0_prop - gamma_0) ;
kaidot0_prop = gammadot0_prop - (2*k_prop*Vgd.*(r_0-rd))./((1+k_prop*(r_0-rd).^2).^2);


k_Nelson = sqrt(k_prop.^2 *(r_0 - rd).^2 + 2*k_prop);
kai0_Nelson = gamma_0 + pi/2 + atan(k_Nelson.*(r_0 - rd)) ;
gammadot0_Nelson = (Vgd./r_0).*sin(kai0_Nelson - gamma_0) ;
kaidot0_Nelson = gammadot0_Nelson - (k_Nelson^2*Vgd.*(r_0-rd))./((1+k_Nelson^2*(r_0-rd).^2).^(3/2));
x_initial = [r_0;gamma_0;x0;y0;kai0_prop;kaidot0_prop;r_0;gamma_0;x0;y0;kai0_Nelson;kaidot0_Nelson];

%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -100:4:100;
[X,Y] = meshgrid(range);
r = sqrt((X-x_c).^2 + (Y-y_c).^2);
gama = atan2((Y-y_c),(X-x_c));
[xdot, ydot] = fun_vf(r,gama);

for i =1:length(X)
    for j = 1:length(X)
        if r(i,j) < rd
            X_in(i,j) = X(i,j);
            Y_in(i,j) = Y(i,j);
           xdot_in(i,j) = xdot(i,j);
           ydot_in(i,j) = ydot(i,j);
        else
            X_out(i,j) = X(i,j);
            Y_out(i,j) = Y(i,j);
           xdot_out(i,j) = xdot(i,j);
           ydot_out(i,j) = ydot(i,j);
        end
    end
end
%-------------------------------------------------------------------------%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@(t,x)fun_comp_autopilot_nowind(t,x) ,tspan, x_initial,options);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_prop = x(:,1) ;
error_r = (r_prop - rd) ;
% Vg_act = sqrt(x(:,7).^2 + wx^2 + wy^2 + 2*wx*x(:,7).*cos(x(:,5)) + 2*wy*x(:,7).*sin(x(:,5)));
x_ini = x(1,3) + x_c ;
y_ini = x(1,4) + y_c ;
x_end = x(end,3) + x_c ;
y_end = x(end,4) + y_c ;

for i = 1:length(t)
    if x(i,1) < rd
        kaid_prop(i,1) = x(i,2) + asin(1./(1+k_prop*(x(i,1)-rd).^2)) ;
        gamma_dot_prop(i,1) = (Vgd./x(i,1)).*sin(kaid_prop(i,1) - x(i,2)) ;
        kaid_dot_prop(i,1) = gamma_dot_prop(i,1) - 2*k_prop*Vgd.*(x(i,1)-rd)./((1+k_prop*(x(i,1)-rd).^2).^2);
    else
        kaid_prop(i,1) = x(i,2) + pi - asin(1./(1+k_prop*(x(i,1)-rd).^2)) ;
        gamma_dot_prop(i,1) = (Vgd./x(i,1)).*sin(kaid_prop(i,1) - x(i,2)) ;
        kaid_dot_prop(i,1) = gamma_dot_prop(i,1) - 2*k_prop*Vgd.*(x(i,1)-rd)./((1+k_prop*(x(i,1)-rd).^2).^2);
    end
end
kappa_des_prop = kaid_dot_prop./Vgd ;
%  Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%  psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% kappa_des = 1./(r_prop.*(1 + k2.*(r_prop - rd).^2)) - (2*k2.*(r_prop - rd))./((1 + k2.*(r_prop - rd).^2).^2);
kappa_actual_prop = x(:,6)./Vgd ;
control_input_prop = k_kai*(kaid_prop - x(:,5)) + k_kaidot*(kaid_dot_prop - x(:,6))  ;
[val_max_prop ,ind_prop] = max(abs(kappa_actual_prop));
x_max_prop = x(ind_prop,1);

%----------------------- Nelson et al ------------------------------------%
r_Nelson = x(:,7);
error_Nelson = (r_Nelson - rd) ;
gamma_Nelson = x(:,8);
kaid_Nelson = gamma_Nelson + pi/2 + atan(k_Nelson.*(r_Nelson - rd)) ;
gamma_dot_Nelson = (Vgd./r_Nelson).*sin(kaid_Nelson - gamma_Nelson) ;
kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(r_Nelson-rd))./((1+k_Nelson^2*(r_Nelson-rd).^2).^(3/2));
kappa_des_Nelson = kaid_dot_Nelson./Vgd ;
%  Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%  psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% kappa_des = 1./(r_prop.*(1 + k2.*(r_prop - rd).^2)) - (2*k2.*(r_prop - rd))./((1 + k2.*(r_prop - rd).^2).^2);
kappa_actual_Nelson = x(:,12)./Vgd ;
control_input_Nelson = k_kai*(kaid_Nelson - x(:,11)) + k_kaidot*(kaid_dot_Nelson - x(:,12))  ;
[val_max_Nelson, ind_Nelson] = max(abs(kappa_actual_Nelson));
x_max_Nelson = x(ind_Nelson,7);
%-------------------------------------------------------------------------%
figure(1)
% if r_0 < rd
%     quiver(X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% else
%     quiver(X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% end
plot(x_sc,y_sc,'k','linewidth',2);hold on;grid on;
plot(x(:,3),x(:,4),'r','linewidth',2);hold on;
plot(x(:,9),x(:,10),'b','linewidth',2);hold on;
h11 = plot(x_c,y_c,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','magenta'); hold on;
h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12 = plot(x_ini,y_ini,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green'); hold on;
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
h13 = plot(x_end,y_end,'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold on;
h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax1 = gca;
ax1.FontSize = 13;
xlabel(' $ x, $ \ (m)','Fontsize',18);
ylabel('$ y, $ \ (m)','Fontsize',18);
% legend('Vector field','Circular \ orbit','Proposed','Nelson \ et al. 2007','Centre of orbit','Initial \ position','Final \ position','Fontsize',14);
legend('Circular \ orbit','Proposed','Nelson \ et al. 2007','Fontsize',14)
axis equal


figure(2)
plot(t,r_prop,'r','LineWidth',2);hold on;grid on;
plot(t,r_Nelson,'b','LineWidth',2);hold on;grid on;
ax2 = gca;
ax2.FontSize = 13;
xlabel('Time, \ $t$, \ (s)','Fontsize',18)
ylabel('Radial \ distance, \ $r$, \ (m)','Fontsize',18);
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

figure(3)
plot(t,error_r,'r','LineWidth',2);hold on;grid on;
plot(t,error_Nelson,'b','LineWidth',2);hold on;grid on;
ax3 = gca;
ax3.FontSize = 13;
xlabel('Time, \ $t$, \ (s)','Fontsize',18)
ylabel('Radial \ error, \ $e = (r - r_d)$, \ (m)','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

figure(4)
plot(t,wrapToPi(x(:,2))*(180/pi),'r','linewidth',2);hold on;grid on;
plot(t,wrapToPi(x(:,8))*(180/pi),'b','linewidth',2);hold on;grid on;
ax4 = gca;
ax4.FontSize = 13;
xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
ylabel('Bearing \ angle, \ $ \gamma$, \ $ (^\circ) $','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

% figure(5)
% plot(r_prop,wrapToPi(kaid_dot_prop)*(180/pi),'r','LineWidth',2);hold on;grid on;
% plot(r_prop,wrapToPi(x(:,6))*(180/pi),'b','LineWidth',2);grid on;
% xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',18)
% ylabel('Course \ rate,  \ $\dot{\chi}$, \ ( $ ^\circ $/s) ','Fontsize',18)
% legend('Desired','Actual','Fontsize',14);

% figure(6)
% plot(r_prop,kappa_des_prop,'r','LineWidth',2);hold on;grid on;
% plot(r_prop,kappa_actual_prop,'b','LineWidth',2);hold on;grid on;
% xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',18)
% ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
% legend('Desired','Actual','Fontsize',14);

% figure(7)
% plot(t,wrapToPi(kaid_prop)*(180/pi),'r','linewidth',2);hold on;grid on;
% plot(t,wrapToPi(x(:,5))*(180/pi),'b','linewidth',2);grid on;
% xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
% ylabel('Course \ angle, \ $ \chi$, \ $ (^\circ) $','Fontsize',18)
% legend('Desired','Actual','Fontsize',14);

% figure(8)
% plot(t,kappa_des_prop,'r','LineWidth',2);hold on;grid on;
% plot(t,kappa_actual_prop,'b','LineWidth',2);hold on;grid on;
% xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
% ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
% legend('Desired','Actual','Fontsize',14);

figure(9)
plot(t,control_input_prop,'r','linewidth',2);hold on;grid on;
plot(t,control_input_Nelson,'b','linewidth',2);hold on;grid on;
ax9 = gca;
ax9.FontSize = 13;
xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
ylabel('Control \ input, \ $ u$, \  (units) ','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

figure(10)
plot(t,kappa_actual_prop,'r','LineWidth',2);hold on;grid on;
plot(t,kappa_actual_Nelson,'b','LineWidth',2);hold on;grid on;
ax10 = gca;
ax10.FontSize = 13;
xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

figure(11)
h = plot(r_prop,kappa_actual_prop,'r','LineWidth',2);hold on;grid on;
h1 = plot(r_Nelson,kappa_actual_Nelson,'b','LineWidth',2);hold on;grid on;
plot(x_max_prop,val_max_prop,'-s','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','cyan'); hold on;grid on;
plot(x_max_Nelson,val_max_Nelson,'-s','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold on;grid on;
% plot([x_max_prop x_max_prop],[0 val_max_prop],'--','color',[0.75 0.75 0.75],'linewidth',2);hold on;grid on;
% plot([x_max_Nelson x_max_Nelson],[0 val_max_Nelson],'--','color',[0.75 0.75 0.75],'linewidth',2);hold on;grid on;
dt_prop = datatip(h,x_max_prop,val_max_prop);hold(ax4,'on');grid(ax4,'on');
h.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h.DataTipTemplate.DataTipRows(2).Format = '%.3f';
dt_Nelson = datatip(h1,x_max_Nelson,val_max_Nelson);hold(ax4,'on');grid(ax4,'on');
h1.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h1.DataTipTemplate.DataTipRows(2).Format = '%.3f';
ax11 = gca;
ax11.FontSize = 13;
xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',18);
ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',18)
% legend('Proposed','Nelson \ et al. 2007','Fontsize',14);
legend('Proposed','Nelson \ et al. 2007','','','Fontsize',14);

figure(12)
plot(t,x(:,6),'r','LineWidth',2);hold on;grid on;
plot(t,x(:,12),'b','LineWidth',2);hold on;grid on;
ax12 = gca;
ax12.FontSize = 13;
xlabel('Time, \ $ t,$ \ (s)','Fontsize',18)
ylabel('Course \ rate,  \ $\dot{\chi}$, \ ( rad./s) ','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);

figure(13)
plot(r_prop,x(:,6),'r','LineWidth',2);hold on;grid on;
plot(r_Nelson,x(:,12),'b','LineWidth',2);hold on;grid on;
ax13 = gca;
ax13.FontSize = 13;
xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',18);
ylabel('Course \ rate,  \ $\dot{\chi}$, \ ( rad./s) ','Fontsize',18)
legend('Proposed','Nelson \ et al. 2007','Fontsize',14);



function out = fun_comp_autopilot_nowind(t,x)
global k_kai k_kaidot rd k_prop Vgd k_Nelson


if x(1) < rd
    kaid = x(2) + asin(1./(1+k_prop*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1)-rd))./((1+k_prop*(x(1) -rd).^2).^2);
else
    kaid = x(2) + pi - asin(1./(1+k_prop*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1) -rd))./((1+k_prop*(x(1) -rd).^2).^2);
end
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos(x(5)) + 2*wy*x(7).*sin(x(5)));

out(1,1) = Vgd.*cos(x(5) - x(2));
out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
out(3,1) = Vgd*cos(x(5))  ;
out(4,1) = Vgd*sin(x(5))  ;
out(5,1) = x(6) ;
out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;

kaid_Nelson = x(8) + pi/2 + atan(k_Nelson.*(x(7) - rd)) ;
gamma_dot_Nelson = (Vgd./x(7) ).*sin(kaid_Nelson - x(8)) ;
kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(x(7) -rd))./((1+k_Nelson^2*(x(7) -rd).^2).^(3/2));

out(7,1) = Vgd.*cos(x(11) - x(8));
out(8,1) = (Vgd./x(7)).*sin(x(11) - x(8));
out(9,1) =  Vgd*cos(x(11))  ;
out(10,1) = Vgd*sin(x(11))  ;
out(11,1) = x(12) ;
out(12,1) = k_kai*(kaid_Nelson - x(11)) + k_kaidot*(kaid_dot_Nelson - x(12) )  ;


end
function [xdot, ydot] = fun_vf(r,gama)
global Vgd k_prop rd
for i = 1:length(r)
    for j = 1:length(r)
        if r(i,j) < rd
            kaid_proposed(i,j) = gama(i,j) + asin(1./(1+k_prop*(r(i,j)-rd).^2)) ;
            %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
            %    kaid_dot = gamma_dot + 2*k2*Vgd.*(r-rd)./((1+k2*(r-rd).^2).^2);
        else
            kaid_proposed(i,j) = gama(i,j) + pi - asin(1./(1+k_prop*(r(i,j)-rd).^2)) ;
            %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
            %    kaid_dot = gamma_dot - 2*k2*Vgd.*(r-rd)./((1+k2*(r-rd).^2).^2);
        end
    end
end
xdot  = Vgd*cos(kaid_proposed);
ydot  = Vgd*sin(kaid_proposed);

end