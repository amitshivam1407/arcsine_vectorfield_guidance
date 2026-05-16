%-------------------------------------------------------------------------%
%------               12th Setember 2022                         ---------%
%-------- Circular path with Autopilot design with no wind        --------%
%-------   Proposed vector field for oustide initial position     --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%---------------- declaration of constant parameters ---------------------%
global rd Vgd k_kai k_kaidot k2 xc yc tf
rd = 50;
Vgd = 12;
k_kai = 200;
k_kaidot = sqrt(k_kai);
red = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
%%%%%%%%%%%%%%%%%%%%%% standoff circle    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc = 0;
yc = 0;
theta = 0:2*pi/1000:2*pi;
x_sc = xc + rd*cos(theta);
y_sc = yc + rd*sin(theta);
%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%
tf = 60;
tspan = [0 tf];

r_0 = 100;
gamma_0 = deg2rad(225);
x0 = r_0*cos(gamma_0) + xc;
y0 = r_0*sin(gamma_0) + yc;
% k2 = 14.9260./(r_0 - rd).^2;
k2 = 0.005;
if r_0 < rd
    kai0 = gamma_0 + asin(1./(1+k2*(r_0-rd).^2)) ;
else
    kai0 = gamma_0 + pi - asin(1./(1+k2*(r_0-rd).^2)) ;
end
gamma_dot0 = (Vgd./r_0).*sin(kai0 - gamma_0) ;
kaidot0 = gamma_dot0 - (2*k2*Vgd.*(r_0-rd))./((1+k2*(r_0-rd).^2).^2);

kai0_deg = rad2deg(wrapToPi(kai0));
% x_initial = [r_0;gamma_0;x0;y0;kai0;kaidot0];
x_initial = [r_0;gamma_0;x0;y0;kai0;0;0.01;0.01;0.01;0.01];
%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
range  = -100:7:100;
[X,Y] = meshgrid(range);
r = sqrt((X-xc).^2 + (Y-yc).^2);
gama = atan2((Y-yc),(X-xc));
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
[t,x] = ode45(@(t,x)fun_circ_autopilot_nowind(t,x) ,tspan, x_initial,options);

%----------------------- parameters to plot   ----------------------------%

r_prop = x(:,1) ;
error_r = (r_prop - rd) ;
% Vg_act = sqrt(x(:,7).^2 + wx^2 + wy^2 + 2*wx*x(:,7).*cos(x(:,5)) + 2*wy*x(:,7).*sin(x(:,5)));
x_ini = x(1,3) + xc ;
y_ini = x(1,4) + yc ;
x_end = x(end,3) + xc ;
y_end = x(end,4) + yc ;

for i = 1:length(t)
    if x(i,1) < rd
        kaid(i,1) = x(i,2) + asin(1./(1+k2*(x(i,1)-rd).^2)) ;
        gamma_dot(i,1) = (Vgd./x(i,1)).*sin(kaid(i,1) - x(i,2)) ;
        kaid_dot(i,1) = gamma_dot(i,1) - 2*k2*Vgd.*(x(i,1)-rd)./((1+k2*(x(i,1)-rd).^2).^2);
    else
        kaid(i,1) = x(i,2) + pi - asin(1./(1+k2*(x(i,1)-rd).^2)) ;
        gamma_dot(i,1) = (Vgd./x(i,1)).*sin(kaid(i,1) - x(i,2)) ;
        kaid_dot(i,1) = gamma_dot(i,1) - 2*k2*Vgd.*(x(i,1)-rd)./((1+k2*(x(i,1)-rd).^2).^2);
    end
end
% kappa_des = kaid_dot./Vgd ;
% %  Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% %  psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% % psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% % kappa_des = 1./(r_prop.*(1 + k2.*(r_prop - rd).^2)) - (2*k2.*(r_prop - rd))./((1 + k2.*(r_prop - rd).^2).^2);
% kappa_actual = x(:,6)./Vgd ;
% [val, ind] = min(kappa_actual);
% rmax = x(ind,1);
% control_input = k_kai*(kaid - x(:,5)) + k_kaidot*(kaid_dot - x(:,6))  ;

% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = ax_fnt+2; % Label font size
leg_fnt = ax_fnt-1; % Legend font size
ax_lw = 3;        % Axis line width

%-------------- figures --------------------------------------------------%
figure(1)
if r_0 < rd
    quiver(X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
else
    quiver(X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
end
plot(x_sc,y_sc,'k','linewidth',lw_);hold on;
plot(x(:,3),x(:,4),'color',red,'linewidth',lw_);hold on;
plot(xc,yc,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','magenta'); hold on;
plot(x_ini,y_ini,'-o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green'); hold on;
plot(x_end,y_end,'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold on;
ax1 = gca;
ax1.FontSize = ax_fnt;
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3)
xlabel(ax1,' $ x, $  m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $  m','Fontsize',lbl_fnt);
legend(ax1,'Vector field','Circular orbit','UAV trajectory','','','','Fontsize',leg_fnt)
axis(ax1,'equal')

figure(2)
plot(t,r_prop,'color',red,'LineWidth',lw_);grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3)
xlabel(ax2,' $t$,  s','Fontsize',lbl_fnt)
ylabel(ax2,'Radial distance $$ r $$, m','Fontsize',lbl_fnt)

figure(3)
plot(t,error_r,'color',red,'LineWidth',lw_);grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3)
xlabel(ax3,' $t$,  s','Fontsize',lbl_fnt)
ylabel(ax3,'Radial error $$(r - r_{\mathrm{d}})$$, m','Fontsize',lbl_fnt)

figure(4)
plot(t,wrapToPi(x(:,2))*(180/pi),'color',red,'linewidth',lw_);grid on;
ax4 = gca;
ax4.FontSize = ax_fnt;
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',3)
xlabel(ax4,' $ t,$ s','Fontsize',lbl_fnt)
ylabel(ax4,' $ \gamma$,  deg.','Fontsize',lbl_fnt)

% figure(5)
% plot(r_prop,wrapToPi(kaid_dot),'-','color',red,'LineWidth',2);hold on;grid on;
% plot(r_prop,wrapToPi(x(:,6)),'--','color',red,'LineWidth',2);grid on;
% ax5 = gca;
% ax5.FontSize = ax_fnt;
% ax5.XColor = 'black';         % Box horizontal lines' color
% ax5.YColor = 'black';         % Box vertical lines' color
% set(ax5,'linewidth',3)
% xlabel(ax5,'Radial  distance, m','Fontsize',lbl_fnt)
% ylabel(ax5,'$\dot{\chi}$, rad./s','Fontsize',lbl_fnt)
% legend('Commanded','Achieved','Fontsize',leg_fnt);

% figure(6)
% % plot(r_prop,kappa_des,'-','color',red,'LineWidth',2);hold on;grid on;
% plot(r_prop,kappa_actual,'-','color',red,'LineWidth',lw_);hold on;grid on;
% ax6 = gca;
% ax6.FontSize = ax_fnt;
% ax6.XColor = 'black';         % Box horizontal lines' color
% ax6.YColor = 'black';         % Box vertical lines' color
% set(ax6,'linewidth',3)
% xlabel(ax6,'Radial distance $$r $$, m','Fontsize',lbl_fnt)
% ylabel(ax6,'Curvature $$\kappa $$, m $ ^{-1} $ ','Fontsize',lbl_fnt)
% axis(ax6,[50 100 -0.04 0.04])
% % legend(ax6,'Commanded','Achieved','Fontsize',leg_fnt);

figure(7)
plot(t,wrapToPi(kaid)*(180/pi),'-','color',red,'linewidth',lw_);hold on;grid on;
plot(t,wrapToPi(x(:,5))*(180/pi),'--','color',blue,'linewidth',lw_);grid on;
ax7 = gca;
ax7.FontSize = ax_fnt;
ax7.XColor = 'black';         % Box horizontal lines' color
ax7.YColor = 'black';         % Box vertical lines' color
set(ax7,'linewidth',3);
% set(ax7,'XTick',[-180, -120, -60, 0, 60, 120, 180],'xticklabel',{'-180' ,'-120', '-60', '0', '60', '120', '180'});
% xticks(ax7,[-180 -120 -60 0 60 120 180]);
% xticklabels(ax7,{'-180' ,'-120', '-60', '0', '60', '120', '180'});
xlabel(ax7,' $ t,$ s','Fontsize',lbl_fnt)
ylabel(ax7,' $ \chi$, deg. ','Fontsize',lbl_fnt)
legend(ax7,'Commanded','Achieved','Fontsize',leg_fnt);

% figure(8)
% % plot(t,kappa_des,'-','color',red,'LineWidth',2);hold on;grid on;
% plot(t,kappa_actual,'--','color',red,'LineWidth',2);hold on;grid on;
% ax8 = gca;
% ax8.FontSize = ax_fnt;
% ax8.XColor = 'black';         % Box horizontal lines' color
% ax8.YColor = 'black';         % Box vertical lines' color
% set(ax8,'linewidth',1.5)
% xlabel(ax8,' $ t,$ s','Fontsize',17)
% ylabel(ax8,'Curvature $$\kappa $$, m $ ^{-1} $ ','Fontsize',17)
% % legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt);
%
% figure(9)
% plot(t,control_input,'color',red,'linewidth',2);grid on;
% box on
% ax9 = gca;
% ax9.FontSize = ax_fnt;
% % Switch on the box around the axis
% ax9.XColor = 'black';         % Box horizontal lines' color
% ax9.YColor = 'black';         % Box vertical lines' color
% set(ax9,'linewidth',1.5)
% xlabel(ax9,'$ t,$  s','Fontsize',17)
% ylabel(ax9,'Control \ input,  rad./s ','Fontsize',17)


function out = fun_circ_autopilot_nowind(t,x)
global  k_kai k_kaidot rd k2 Vgd xc yc tf


% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos(x(5)) + 2*wy*x(7).*sin(x(5)));

kr = 2;

if (t>=0) && (t<=2)
    
    if x(1) < rd
        kaid = x(2) + asin(1./(1+k2*(x(1) -rd).^2)) ;
        %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
        %    kaid_dot = gamma_dot + 2*k2*Vgd.*(x(1)-rd)./((1+k2*(x(1) -rd).^2).^2);
    else
        kaid = x(2) + pi - asin(1./(1+k2*(x(1) -rd).^2)) ;
        %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
        %    kaid_dot = gamma_dot - 2*k2*Vgd.*(x(1) -rd)./((1+k2*(x(1) -rd).^2).^2);
    end

    out(1,1) = Vgd.*cos(x(5) - x(2));
    out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
    out(3,1) = Vgd*cos(x(5))  ;
    out(4,1) = Vgd*sin(x(5))  ;
    out(5,1) = k_kai*wrapToPi(kaid - x(5));
    
    out(6,1) = 0 ;  
    out(7,1) = 0;
    out(8,1) = 0;
    out(9,1) = 0;
    out(10,1) = 0;
    
elseif (t>2) && (t<=4)
     if x(1) < rd
        kaid = x(2) + asin(1./(1+k2*(x(1) -rd).^2)) ;
        %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
        %    kaid_dot = gamma_dot + 2*k2*Vgd.*(x(1)-rd)./((1+k2*(x(1) -rd).^2).^2);
    else
        kaid = x(2) + pi - asin(1./(1+k2*(x(1) -rd).^2)) ;
        %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
        %    kaid_dot = gamma_dot - 2*k2*Vgd.*(x(1) -rd)./((1+k2*(x(1) -rd).^2).^2);
    end

    out(1,1) = Vgd.*cos(x(5) - x(2));
    out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
    out(3,1) = Vgd*cos(x(5))  ;
    out(4,1) = Vgd*sin(x(5))  ;
    out(5,1) = k_kai*wrapToPi(kaid - x(5));
    
%     out(1,1) = Vgd.*cos(x(5) - x(2));
%     out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
%     out(3,1) = Vgd*cos(x(5))  ;
%     out(4,1) = Vgd*sin(x(5))  ;
%     out(5,1) = k_kai*wrapToPi(kaid - x(5));
    
   if x(6) < rd
    kaid_hat = x(7) + asin(1./(1+k2*(x(6) -rd).^2)) ;
    %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    %    kaid_dot = gamma_dot + 2*k2*Vgd.*(x(1)-rd)./((1+k2*(x(1) -rd).^2).^2);
else
    kaid_hat = x(7) + pi - asin(1./(1+k2*(x(6) -rd).^2)) ;
    %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    %    kaid_dot = gamma_dot - 2*k2*Vgd.*(x(1) -rd)./((1+k2*(x(1) -rd).^2).^2);
   end 
 Lr = 20;   Ltheta = 20; Lx = 20; Ly = 20; Lkai = 20;
 zr = Lr*sign(x(1) - x(6));
 ztheta = Ltheta*sign(x(2) - x(7));
 zx = Lx*sign(x(3) - x(8));
 zy = Ly*sign(x(4) - x(9));
 zkai = Lkai*sign(x(5) - x(10));
% ;
% - kr *(abs(x(9)- x(4))^(1/3))*sign(x(9)- x(4));

    out(6,1) = -kr *(abs(x(6) - rd)^(1/3))*sign(x(6) - rd) ;
    % out(7,1) = (Vgd./x(6)).*sin(x(10) - x(7)) ;
    out(7,1) = -kr *(abs(x(7))^(1/3))*sign(x(7)) ;
    out(8,1) =  Vgd *cos(x(10)) ;
    out(9,1) =  Vgd *sin(x(10))  ;
    out(10,1) = k_kai*wrapToPi(kaid_hat - x(10)) ;

elseif (t>4) && (t<=tf)
    if x(1) < rd
    kaid = x(2) + asin(1./(1+k2*(x(1) -rd).^2)) ;
    %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    %    kaid_dot = gamma_dot + 2*k2*Vgd.*(x(1)-rd)./((1+k2*(x(1) -rd).^2).^2);
else
    kaid = x(2) + pi - asin(1./(1+k2*(x(1) -rd).^2)) ;
    %    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    %    kaid_dot = gamma_dot - 2*k2*Vgd.*(x(1) -rd)./((1+k2*(x(1) -rd).^2).^2);
end
    out(1,1) = Vgd.*cos(x(5) - x(2));
    out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
    out(3,1) = Vgd*cos(x(5))  ;
    out(4,1) = Vgd*sin(x(5))  ;
    out(5,1) = k_kai*wrapToPi(kaid - x(5));
    
    out(6,1) = 0;
    out(7,1) = 0;
    out(8,1) = 0 ;
    out(9,1) =  0;
    out(10,1) = 0;
% else
    
    end
% out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;
% out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;

end

function [xdot, ydot] = fun_vf(r,gama)
global Vgd k2 rd
for i = 1:length(r)
    for j = 1:length(r)
        if r(i,j) < rd
            kaid_proposed(i,j) = gama(i,j) + asin(1./(1+k2*(r(i,j)-rd).^2)) ;
            %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
            %    kaid_dot = gamma_dot + 2*k2*Vgd.*(r-rd)./((1+k2*(r-rd).^2).^2);
        else
            kaid_proposed(i,j) = gama(i,j) + pi - asin(1./(1+k2*(r(i,j)-rd).^2)) ;
            %    gamma_dot = (Vgd./r).*sin(kaid_proposed - gamma) ;
            %    kaid_dot = gamma_dot - 2*k2*Vgd.*(r-rd)./((1+k2*(r-rd).^2).^2);
        end
    end
end
xdot  = Vgd*cos(kaid_proposed);
ydot  = Vgd*sin(kaid_proposed);

end