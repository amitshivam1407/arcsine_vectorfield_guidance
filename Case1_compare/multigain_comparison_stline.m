%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
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
f4 = figure;
ax4 = axes;
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
% f18 = figure;
% ax18 = axes;
%--------------- declaration of constant parameters ----------------------%
global  Vgd k_kai  k_kaidot flag1 t_flag1 flag2 t_flag2 
Vgd = 10;
k_kai = 500;  
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%

%-------------------------------------------------------------------------%

%------------------- initial conditions   --------------------------------%

tspan = [0 50];

x0 = 100;
y0 = -100;

% k_prop = [0.0001 0.0003 0.0005 0.0007 0.0009 0.0015 0.0022 0.003 0.0037 0.0044 0.0051 0.0058 0.0065 0.0072 0.0079 0.0087 0.0094 0.01];
k_prop = [0.0001 0.00018 0.0003 0.0005 0.0009 0.0014 0.0020 0.0026 0.0033 0.004 0.0048 0.0058 0.0065 0.0072 0.0079 0.0087 0.0094 0.01];
% k_prop = 1*10^(-4):0.001/5 :0.001;
flag1 = true(10);
t_flag1 = zeros(1,10);
flag2 = true(10);
t_flag2 = zeros(1,10);

for i = 1:length(k_prop)

if x0 < 0
   kaid0_prop(i) =  asin(1./(1+ k_prop(i)*(x0).^2)) ;
   kaid_dot0_prop(i) = (2*k_prop(i)*Vgd.*(x0))./((1+k_prop(i)*(x0).^2).^2);
else
   kaid0_prop(i) = pi - asin(1./(1 + k_prop(i)*(x0).^2)) ;  
   kaid_dot0_prop(i) = - (2*k_prop(i)*Vgd.*(x0))./((1 + k_prop(i)*(x0).^2).^2);
end
kai0_prop(i) = kaid0_prop(i);
kaidot0_prop(i) = kaid_dot0_prop(i);

kai0_deg(i) = rad2deg(wrapToPi(kai0_prop(i)));


k_Nelson(i) = sqrt(k_prop(i)^2 * x0^2 + 2*k_prop(i));

kai0_Nelson(i) = pi/2 + atan(k_Nelson(i).*x0);
kaidot0_Nelson(i) = -(k_Nelson(i)^2*Vgd*x0)./(1+k_Nelson(i)^2.*x0.^2).^(3/2);



%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
% range  = -110:10:110;
% [X,Y] = meshgrid(range);
% [kaid,kaid_dot] = vf_proposed(X);
% xdot  = Vgd*cos(kaid);
% ydot  = Vgd*sin(kaid);


%----------------  ode solver   ------------------------------------------%
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

x_initial = [x0;y0;kai0_prop(i);kaidot0_prop(i);x0;y0;kai0_Nelson(i);kaidot0_Nelson(i)];
[t,x] = ode45(@(t,x)fun_stline_upd(t,x,k_prop(i),k_Nelson(i),i) ,tspan, x_initial,options);
% 
% k_prop_plot = 0.0002:(0.005-0.0002)/(length(t)-1):0.005;
% k_prop_plot = k_prop_plot';
% for i = 1 :length(t)
%     if x0 < 0
%         kaid0_prop_plot(i) =  asin(1./(1+ k_prop_plot(i)*(x0).^2)) ;
%     else
%         kaid0_prop_plot(i) = pi - asin(1./(1 + k_prop_plot(i)*(x0).^2)) ;
%     end
%     kaid0_prop_plot(i) = rad2deg(wrapToPi(kaid0_prop_plot(i)));
%     k_Nelson_plot(i) = sqrt(k_prop_plot(i)^2 * x0^2 + 2*k_prop_plot(i));    
% end
% kaid0_prop_plot = kaid0_prop_plot' ;
% k_Nelson_plot =  k_Nelson_plot';
% 
% %----------------------  parameters to plot   ----------------------------%
% 
% % x_ini = x(1,1) ;
% % y_ini = x(1,2) ;
% % x_end = x(end,1) ;
% % y_end = x(end,2) ;
% 
% 

[kai_des_prop,kaidot_des] = fun_prop_kaid(x(:,1),k_prop(i));
% for j = 1:length(x)
% if x(j,1) < 0
%    kai_des_prop(j,:) =  asin(1./(1+ k_prop(i)*(x(j,1)).^2)) ;
%    kaidot_des(j,:) = -(2*k_prop(i)*Vgd.*x(j,1))./((1+k_prop(i)*(x(j,1)).^2).^2);
% else
%    kai_des_prop(j,:) = pi - asin(1./(1 + k_prop(i)*(x(j,1)).^2)) ;  
%    kaidot_des(j,:) = - (2*k_prop(i)*Vgd.*x(j,1))./((1 + k_prop(i)*(x(j,1)).^2).^2);
% end
% end
% 
kaidot_actual_prop = x(:,4);
% kappa_des = kaidot_des./Vgd;
kappa_actual_prop = x(:,4)./Vgd ;
control_effort_prop(i) = trapz(t,(Vgd*x(:,4)).^2);
[val_max_prop(i), ind_prop(i)] = min(kappa_actual_prop);
x_max_prop(i) = x(ind_prop(i),1);
%  
% 
kai_des_Nelson = pi/2 + atan(k_Nelson(i).*x(:,5));
kaidot_actual_Nelson = x(:,8);
control_effort_Nelson(i) = trapz(t,(Vgd*x(:,8)).^2);
% kappa_des = kaidot_des./Vgd;
kappa_actual_Nelson = x(:,8)./Vgd ;
[val_max_Nelson(i), ind_Nelson(i)] = min(kappa_actual_Nelson);
x_max_Nelson(i) = x(ind_Nelson(i),5);
%if val_max_prop == val_max_Nelson

percent_red_curvature(i) = ((val_max_Nelson(i) - val_max_prop(i))./val_max_Nelson(i)) *100;
percent_red_controleffort(i) = ((control_effort_Nelson(i) - control_effort_prop(i))./control_effort_Nelson(i)) *100;
percent_red_settlingtime(i) = -((t_flag2(i) - t_flag1(i))./t_flag2(i)) *100;

ratio_curvature(i) = val_max_prop(i)./val_max_Nelson(i) ;
ratio_controleffort(i) = control_effort_prop(i)./control_effort_Nelson(i) ;
ratio_settlingtime(i) = t_flag1(i)./t_flag2(i);
% 
% % for i = 1:length(t)
% %  [val_prop(i,:),val_Nelson(i,:)] = fun_kappa(x(i,1),x(i,5));
% %  ratio_curvature(i,:) = abs(val_prop(i,:))./abs(val_Nelson(i,:)) ;
% % end
% % ratio_curvature = abs(val_max_prop)./abs(val_max_Nelson);
% 
% 
% 
%----------------------- plotting figures --------------------------------%
%%
% Plot format control variables
lw = 3;             % Line width
ms = 6;             % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20; % Label font size
leg_fnt = 16; % Legend font size
ax_lw = 3;        % Axis line width
%%
f1;
% quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
path = xline(ax1,0,'k','linewidth',lw);hold(ax1,'on'); 
plot(ax1,x(:,1),x(:,2),'r','linewidth',lw);hold(ax1,'on'); 
plot(ax1,x(:,5),x(:,6),'b','linewidth',lw);hold(ax1,'on'); 
% plot(x_ini,y_ini,'-o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','green'); hold on;
% plot(x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
ax1.FontSize = ax_fnt;
box(ax1,'on')                     % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax1,' $ x, $  m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $  m','Fontsize',lbl_fnt);
legend(ax1,'Desired path','Proposed','Nelson et al. ','Fontsize',leg_fnt);
axis(ax1, 'equal');
% 
% 
f2;
h2 = plot(ax2,t,(x(:,1)),'r','LineWidth',lw);hold(ax2,'on');grid(ax2,'on');
% dt_prop = datatip(h2,t_flag1(i),0.2);hold(ax2,'on');grid(ax2,'on');
h21 = plot(ax2,t,(x(:,5)),'b','LineWidth',lw);hold(ax2,'on');grid(ax2,'on');
% dt_Nelson = datatip(h21,t_flag2(i),0.2);hold(ax2,'on');grid(ax2,'on');
% yline(ax2,0.2,'k','LineWidth',2);hold(ax2,'on');grid(ax2,'on');
ax2.FontSize = ax_fnt;
box(ax2,'on')                     % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax2,' $t$,  s','Fontsize',lbl_fnt);
ylabel(ax2,'Cross-track error, m','Fontsize',lbl_fnt);
legend(ax2,'Proposed','Nelson  et al. ','Fontsize',leg_fnt);
% 
f3;
plot(ax3,x(:,1),kaidot_actual_prop,'r','LineWidth',lw);hold(ax3,'on');grid(ax3,'on');
plot(ax3,x(:,5),kaidot_actual_Nelson,'b','LineWidth',lw);
ax3.FontSize = ax_fnt;
box(ax3,'on')                     % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax3,'Cross-track  error,  m','Fontsize',lbl_fnt);
ylabel(ax3,'$\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt);
legend(ax3,'Proposed','Nelson  et al. ','Fontsize',leg_fnt);
% 
% 
f4;
h4 = plot(ax4,x(:,1),kappa_actual_prop,'r','LineWidth',lw);hold(ax4,'on');grid(ax4,'on');
h41 = plot(ax4,x(:,5),kappa_actual_Nelson,'b','LineWidth',lw);
plot(ax4,x_max_prop(i),val_max_prop(i),'-s','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','cyan'); hold(ax4,'on');grid(ax4,'on');
plot(ax4,x_max_Nelson(i),val_max_Nelson(i),'-s','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold(ax4,'on');grid(ax4,'on');
plot(ax4,[x_max_prop(i) x_max_prop(i)],[0 val_max_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
plot(ax4,[x_max_Nelson(i) x_max_Nelson(i)],[0 val_max_Nelson(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
% dt_prop = datatip(h4,x_max_prop(i),val_max_prop(i));hold(ax4,'on');grid(ax4,'on');
% dt_Nelson = datatip(h41,x_max_Nelson(i),val_max_Nelson(i));hold(ax4,'on');grid(ax4,'on');
ax4.FontSize = ax_fnt;
box(ax4,'on')                     % Switch on the box around the axis
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax4,'Cross-track error, m','Fontsize',lbl_fnt);
ylabel(ax4,'Curvature, m $ ^{-1} $ ','Fontsize',lbl_fnt);
legend(ax4,'Proposed','Nelson et al. ','$\kappa_{m} $ Proposed','$\kappa_{m} $ Nelson et al. ','Fontsize',leg_fnt);
% 
f5;
plot(ax5,t,kaidot_actual_prop,'r','LineWidth',lw);hold(ax5,'on');grid(ax5,'on');
plot(ax5,t,kaidot_actual_Nelson,'b','LineWidth',lw);
ax5.FontSize = ax_fnt;
box(ax5,'on')                     % Switch on the box around the axis
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax5,' $t$, s','Fontsize',lbl_fnt);
ylabel(ax5,' $\dot{\chi}$, rad./s','Fontsize',lbl_fnt);
legend(ax5,'Proposed','Nelson et al. ','Fontsize',leg_fnt);
% 
f6;
plot(ax6,t,kappa_actual_prop,'LineWidth',lw);hold(ax6,'on');grid(ax6,'on');
% plot(ax6,t,kappa_actual_Nelson,'b','LineWidth',lw);
ax6.FontSize = ax_fnt;
box(ax6,'on')                     % Switch on the box around the axis
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax6,' $t$, s','Fontsize',lbl_fnt);
ylabel(ax6,'Curvature, m $ ^{-1} $ ','Fontsize',lbl_fnt);
legend(ax6,'Proposed','Fontsize',leg_fnt);
% 
f7;
plot(ax7,t,x(:,3)*(180/pi),'r','linewidth',lw);hold(ax7,'on');grid(ax7,'on');
plot(ax7,t,x(:,7)*(180/pi),'b','linewidth',lw);
ax7.FontSize = ax_fnt;
box(ax7,'on')                     % Switch on the box around the axis
ax7.XColor = 'black';         % Box horizontal lines' color
ax7.YColor = 'black';         % Box vertical lines' color
set(ax7,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax7,' $ t, $ s','Fontsize',lbl_fnt);
ylabel(ax7,' $ \chi$, deg.','Fontsize',lbl_fnt);
legend(ax7,'Proposed','Nelson et al. ','Fontsize',leg_fnt);

% 
f8;
plot(ax8,x(:,1),x(:,3)*(180/pi),'r','linewidth',lw);hold(ax8,'on');grid(ax8,'on');
plot(ax8,x(:,5),x(:,7)*(180/pi),'b','linewidth',lw);
ax8.FontSize = ax_fnt;
box(ax8,'on')                     % Switch on the box around the axis
ax8.XColor = 'black';         % Box horizontal lines' color
ax8.YColor = 'black';         % Box vertical lines' color
set(ax8,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax8,'Cross-track error,  m','Fontsize',lbl_fnt);
ylabel(ax8,'$ \chi$, deg.','Fontsize',lbl_fnt);
legend(ax8,'Proposed','Nelson et al. ','Fontsize',leg_fnt);
 
% 
f9;
plot(ax9,x(:,1),kai_des_prop.*(180/pi),'r','linewidth',lw);hold(ax9,'on');grid(ax9,'on');
plot(ax9,x(:,1),x(:,3)*(180/pi),'b','linewidth',lw);
ax9.FontSize = ax_fnt;
box(ax9,'on')                     % Switch on the box around the axis
ax9.XColor = 'black';         % Box horizontal lines' color
ax9.YColor = 'black';         % Box vertical lines' color
set(ax9,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax9,'Cross-track  error, m','Fontsize',lbl_fnt);
ylabel(ax9,'$ \chi$, deg.','Fontsize',lbl_fnt);
legend(ax9,'Desired proposed','Actual proposed','Fontsize',leg_fnt);

% 
f10;
plot(ax10,x(:,5),kai_des_Nelson.*(180/pi),'r','linewidth',lw);hold(ax10,'on');grid(ax10,'on');
plot(ax10,x(:,5),x(:,7)*(180/pi),'b','linewidth',lw);
ax10.FontSize = ax_fnt;
box(ax10,'on')                     % Switch on the box around the axis
ax10.XColor = 'black';         % Box horizontal lines' color
ax10.YColor = 'black';         % Box vertical lines' color
set(ax10,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax10,'Cross-track  error, m','Fontsize',lbl_fnt);
ylabel(ax10,' $ \chi$, deg.','Fontsize',lbl_fnt);
legend(ax10,'Desired Nelson','Actual Nelson','Fontsize',leg_fnt);

% f18;
% plot(ax18,x(:,5),kappa_actual_Nelson,'LineWidth',lw);hold(ax18,'on');grid(ax18,'on');
% ax18.FontSize = ax_fnt;
% box(ax18,'on')                     % Switch on the box around the axis
% ax18.XColor = 'black';         % Box horizontal lines' color
% ax18.YColor = 'black';         % Box vertical lines' color
% set(ax18,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax18,' $t$, s','Fontsize',lbl_fnt);
% ylabel(ax18,'Curvature, m $ ^{-1} $ ','Fontsize',lbl_fnt);
% legend(ax18,'Nelson et al. ','Fontsize',leg_fnt);



end

f11;
plot(ax11,k_prop,ratio_controleffort,'o-','MarkerFaceColor','green',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax11,'on');grid(ax11,'on');
ax11.FontSize = ax_fnt;
box(ax11,'on')                     % Switch on the box around the axis
ax11.XColor = 'black';         % Box horizontal lines' color
ax11.YColor = 'black';         % Box vertical lines' color
set(ax11,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax11,' Guidance  gain $$k_{\mathrm{s}}$$','Fontsize',lbl_fnt);
ylabel(ax11,' Control effort  ratio','Fontsize',lbl_fnt);
axis(ax11,[0 k_prop(end) 0 1.5])

f12;
plot(ax12,k_prop,ratio_curvature,'s-','MarkerFaceColor','red',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax12,'on');grid(ax12,'on');
ax12.FontSize = ax_fnt;
box(ax12,'on')                     % Switch on the box around the axis
ax12.XColor = 'black';         % Box horizontal lines' color
ax12.YColor = 'black';         % Box vertical lines' color
set(ax12,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax12,' Guidance  gain $$k_{\mathrm{s}}$$ ','Fontsize',lbl_fnt);
ylabel(ax12,' Maximum  curvature  ratio','Fontsize',lbl_fnt);
axis(ax12,[0.0001 .01 0 1.5])

f13;
plot(ax13,k_prop,ratio_settlingtime,'d-','MarkerFaceColor','cyan',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax13,'on');grid(ax13,'on');
ax13.FontSize = ax_fnt;
box(ax13,'on')                     % Switch on the box around the axis
ax13.XColor = 'black';         % Box horizontal lines' color
ax13.YColor = 'black';         % Box vertical lines' color
set(ax13,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax13,' Guidance  gain $$k_{\mathrm{s}}$$','Fontsize',lbl_fnt);
ylabel(ax13,' Settling time  ratio','Fontsize',lbl_fnt);
axis(ax13,[0 k_prop(end) 0 2])

f14;
plot(ax14,k_prop,percent_red_settlingtime,'d-','MarkerFaceColor','cyan',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax14,'on');grid(ax14,'on');
ax14.FontSize = ax_fnt;
box(ax14,'on')                     % Switch on the box around the axis
ax14.XColor = 'black';         % Box horizontal lines' color
ax14.YColor = 'black';         % Box vertical lines' color
set(ax14,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax14,' Guidance  gain $$k_{\mathrm{s}}$$','Fontsize',lbl_fnt);
ylabel(ax14,' Percentage ','Fontsize',lbl_fnt);
axis(ax14,[0 k_prop(end) 0 100])

f15;
plot(ax15,k_prop,percent_red_curvature,'s-','MarkerFaceColor','red','MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax15,'on');grid(ax15,'on');
ax15.FontSize = ax_fnt;
box(ax15,'on')                     % Switch on the box around the axis
ax15.XColor = 'black';         % Box horizontal lines' color
ax15.YColor = 'black';         % Box vertical lines' color
set(ax15,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax15,' Guidance  gain $$k_{\mathrm{s}}$$','Fontsize',lbl_fnt);
ylabel(ax15,' Percentage ','Fontsize',lbl_fnt);
% axis(ax15,[0.0001 0.01 0 100])

f16;
plot(ax16,k_prop,percent_red_controleffort,'o-','MarkerFaceColor','green',...
    'MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop),'linewidth',lw);hold(ax16,'on');grid(ax16,'on');
ax16.FontSize = ax_fnt;
box(ax16,'on')                     % Switch on the box around the axis
ax16.XColor = 'black';         % Box horizontal lines' color
ax16.YColor = 'black';         % Box vertical lines' color
set(ax16,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax16,' Guidance  gain $$k_{\mathrm{s}}$$','Fontsize',lbl_fnt);
ylabel(ax16,' Percentage','Fontsize',lbl_fnt);

f17;
yyaxis left
plot(ax17,k_prop,kai0_deg,'s-','MarkerFaceColor','red','MarkerEdgeColor','black','MarkerIndices', 1:length(k_prop)/18:length(k_prop),'linewidth',lw);hold(ax17,'on');grid(ax17,'on');
axis(ax17,[0 k_prop(end) 90 190])
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax17.FontSize = ax_fnt;
box(ax17,'on')                     % Switch on the box around the axis
ax17.XColor = 'black';         % Box horizontal lines' color
ax17.YColor = 'black';         % Box vertical lines' color
set(ax17,'linewidth',ax_lw)   % Axis linewidth (box and grid)
ylabel(ax17,'Commanded course angle $$(t = 0)$$, deg.','Fontsize',lbl_fnt-2);hold(ax17,'on');grid(ax17,'on');
yyaxis right
ylabel(ax17,'  Guidance  gain $$\hat{k}_{\mathrm{s}}$$ for Ref. [11]','Fontsize',lbl_fnt-2);hold(ax17,'on');grid(ax17,'on');
plot(ax17,k_prop,k_Nelson,'o-','MarkerFaceColor','green','MarkerEdgeColor','black', 'MarkerIndices', 1:length(k_prop)/18:length(k_prop),'linewidth',2);hold(ax17,'on');grid(ax17,'on');
axis(ax17,[0 k_prop(end) 0 1.25]);
xlabel(ax17,' Guidance  gain $$k_{\mathrm{s}}$$ ','Fontsize',lbl_fnt-2);
legend(ax17,'$\chi_{d_{0}}$','$$\hat{k}_{\mathrm{s}}$$','Fontsize',leg_fnt)



function out = fun_stline_upd(t,x,k_prop,k_Nelson,i)
global Vgd k_kai k_kaidot flag1 t_flag1 flag2 t_flag2 

if (flag1(i) == true) && (x(1)<0.2)
    t_flag1(i) = t;
    flag1(i) = false;
end

if (flag2(i) == true) && (x(5)<0.2)
    t_flag2(i) = t;
    flag2(i) = false;
end

% 
% if (flag2 == true) && (x(5)<0.2)
%     t_flag2 = t;
%     flag2 = false;
% end

if x(1) < 0
   kaid_prop =  asin(1./(1+ k_prop*(x(1)).^2)) ;
   kaid_dot_prop = -(2*k_prop*Vgd.*(x(1)))./((1+k_prop*(x(1)).^2).^2);
else
   kaid_prop = pi - asin(1./(1 + k_prop*(x(1)).^2)) ;  
   kaid_dot_prop = - (2*k_prop*Vgd.*(x(1)))./((1 + k_prop*(x(1)).^2).^2);
end

kaid_Nelson = pi/2 + atan(k_Nelson.*x(5));
kaid_dot_Nelson = -(k_Nelson^2*Vgd*x(5))./(1+k_Nelson^2.*x(5).^2).^(3/2);
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid_prop - x(3)) + k_kaidot*(kaid_dot_prop - x(4) )  ;

out(5,1) = Vgd*cos(x(7));
out(6,1) = Vgd*sin(x(7));
out(7,1) = x(8) ;
out(8,1) = k_kai*(kaid_Nelson - x(7)) + k_kaidot*(kaid_dot_Nelson - x(8))  ;
end

% function [val_prop,val_Nelson]  = fun_kappa(x_prop,x_Nelson,k_prop,k_Nelson)
% global Vgd 
% % for i = 1:length(x_prop)
%     if x_prop < 0
%         kaidot_des = -(2*k_prop*Vgd.*x_prop)./((1+k_prop*(x_prop).^2).^2);
%         kappa_des_prop = kaidot_des./Vgd ;
%     else      
%         kaidot_des = - (2*k_prop*Vgd.*x_prop)./((1 + k_prop*(x_prop).^2).^2);
%         kappa_des_prop = kaidot_des./Vgd ;
%     end
% [val_prop, ind_prop] = max(abs(kappa_des_prop));
% x_max_prop = x(ind_prop,1);
% kaid_dot_Nelson = -(k_Nelson^2*Vgd*x_Nelson)./(1+k_Nelson^2.*x_Nelson.^2).^(3/2);
% kappa_des_Nelson = kaid_dot_Nelson./Vgd ;
% [val_Nelson, ind_Nelson] = max(abs(kappa_des_Nelson));
% x_max_Nelson = x(ind_Nelson,1);
% % ratio_curvature(i,1) = abs(val_max_prop(i,1))./abs(val_max_Nelson(i,1)) ;
% %end
% end

function [kaid,kaid_dot] = vf_proposed(x)
global Vgd  k_prop 
for i = 1:length(x)
    for j = 1:length(x)
if x(i,j) < 0
   kaid(i,j) = asin(1./(1+k_prop*(x(i,j)).^2)) ;
   kaid_dot(i,j) =  2*k_prop*Vgd.*(x(i,j))./((1+k_prop*(x(i,j)).^2).^2);
else
   kaid(i,j) = pi - asin(1./(1+k_prop*(x(i,j)).^2)) ;
   kaid_dot(i,j) = - 2*k_prop*Vgd.*(x(i,j))./((1+k_prop*(x(i,j)).^2).^2);
end
    end
end
end

function [kaid,kaid_dot] = fun_prop_kaid(x,k_prop)
global Vgd 
for i = 1:length(x) 
if x(i,:) < 0
   kaid(i,:) = asin(1./(1+k_prop*(x(i,:)).^2)) ;
   kaid_dot(i,:) =  2*k_prop*Vgd.*(x(i,:))./((1+k_prop*(x(i,:)).^2).^2);
else
   kaid(i,:) = pi - asin(1./(1+k_prop*(x(i,:)).^2)) ;
   kaid_dot(i,:) = - 2*k_prop*Vgd.*(x(i,:))./((1+k_prop*(x(i,:)).^2).^2);
end
end
end

