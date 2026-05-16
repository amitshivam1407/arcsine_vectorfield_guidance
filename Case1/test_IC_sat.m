%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%--------------- declaration of constant parameters ----------------------%
global  Vgd k_kai k_kaidot k2 figure_type_good kaidot_max te kaidot_ode kaidot kaid_sim  kai_sim x_sim
Vgd = 10;
k_kai = 2 ; 
% k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
k_kaidot = 2.95;  % 20 for better curvature profile  %%
figure_type_good = 1;
kaidot_max = 0.98;
%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%

tspan = [0 25];

x0 = 5;
y0 = -100;
k2 = 0.002;
if x0 < 0
   kaid0 =  asin(1./(1+ k2 *(x0).^2)) ;
%    kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
   kaid0 = pi - asin(1./(1 + k2 *(x0).^2)) ;  
%   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
kaid0_deg = rad2deg(kaid0);
% kai0 = pi + kaid0;
kai0 = deg2rad(270);
% kai0 = normalise(deg2rad(-140));
kai0_deg = rad2deg(kai0);
% kai0 = 3*pi/2;

% kaidot0 = -(2*k2*Vgd.*x0)./((1+k2*(x0).^2).^2);
kaidot0 = 0;
% k2 = abs(14.92597/x0^2);

% x_initial = [x0;y0;kai0;kaidot0];
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
for i = 1:length(t)
if x(i,1) < 0
   kai_des(i,1) =  asin(1./(1+ k2*(x(i,1)).^2)) ;
   kaidot_des(i,1) = -2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
   kai_des(i,1) = pi - asin(1./(1 + k2*(x(i,1)).^2)) ;  
   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
end
kaidot_actual = x(:,4);
kappa_des = kaidot_des./Vgd;
kappa_actual = x(:,4)./Vgd ;
% [kai_des,kaidot_des] = vf_proposed(x(:,1));

% Va_des = sqrt((Vgd*cos(kai_des) - wx).^2 + (Vgd*sin(kai_des) - wy).^2) ;
% psi_des = kai_des + asin((wx*sin(kai_des) - wy*cos(kai_des))./ Va_des );
% % psiddot_des = kaidot_des.*((Va_des.*cos(psi_des - kai_des) + (wx*cos(kai_des) + wy*sin(kai_des)))./(Va_des.*cos(psi_des - kai_des)));
% c1 = (1./x(:,5)).*kaidot_des.*(wx*cos(kai_des)+wy*sin(kai_des));
% c2 = (wx*sin(kai_des)-wy*cos(kai_des)).*((1./x(:,5)).^2).*(Vgd./Va_des).*kaidot_des.*(-wx*sin(kai_des)+wy*cos(kai_des));
% psiddot_des = kaidot_des + (1./sqrt(1 - ((1./x(:,5)).*(wx*sin(kai_des)-wy*cos(kai_des))).^2)).*(c1 + c2);
% psi_actual = x(:,3);
% psidot_actual = x(:,4);
% Vg_actual = sqrt((x(:,5).*cos(x(:,3)) + wx).^2 + (x(:,5).*sin(x(:,3)) + wy).^2) ;
% % Vg_actual = (wx.*cos(x(:,3)) + wy.*sin(x(:,3))) + sqrt(x(:,5).^2 - (wx.*sin(x(:,3)) - wy.*cos(x(:,3))).^2);
% kai_actual = x(:,3) + asin((-wx*sin(x(:,3)) + wy*cos(x(:,3)))./ Vg_actual );
% kaidot_actual = x(:,4).*((Vg_actual.*cos(kai_actual - x(:,3)) - (wx*cos(x(:,3)) + wy*sin(x(:,3))))./(Vg_actual.*cos(kai_actual - x(:,3))));
% kappa_actual = kaidot_actual./Vg_actual;
%---------------------- plotting figures ---------------------------------%

% figure_plot = fun_stlineplot(t,X,Y,xdot,ydot,x(:,1),x(:,2),x_ini,y_ini,x_end,y_end,...
%     kaidot_des,kaidot_actual,kappa_des,kappa_actual,kai_des,x(:,3));

f1 = figure;
ax1 = axes;
f2 = figure;
ax2 = axes;
% f3 = figure;
% ax3 = axes;
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
% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20; % Label font size
leg_fnt = 16; % Legend font size
ax_lw = 3;        % Axis line width
% Colors
% blue = [0 0.4470 0.7410];
% red = [0.8500 0.3250 0.0980];
blue = [0 0 1];
red = [1 0 0];
orange = [0.9290 0.6940 0.1250];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
cyan = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840];
black = [0 0 0];

f1;
h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
xl = xline(ax1,0,'k','linewidth',3);hold(ax1,'on'); 
h1 = plot(ax1,x(:,1),x(:,2),'r','linewidth',3);hold(ax1,'on');
plot(ax1,x_ini,y_ini,'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold(ax1,'on');
plot(ax1,x_end,y_end,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','cyan'); hold(ax1,'on');
% ax1.FontSize = 13;
l = legend(ax1,[h xl h1] ,'Vector field','Desired path','UAV trajectory');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size
% ax1 = gca;     
ax1.FontSize = ax_fnt;
axis(ax1,'equal')
% Outer box setup
box on                      % Switch on the box around the axis
ax1.XColor = black;         % Box horizontal lines' color
ax1.YColor = black;         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax1,'$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
axis(ax1,'equal')

f2;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
h3 = plot(ax2,t,x(:,1),'r','linewidth',3);hold(ax2,'on');grid(ax2,'on');  
ax2.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax2.XColor = black;         % Box horizontal lines' color
ax2.YColor = black;         % Box vertical lines' color
set(ax2,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax2,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax2,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)

% f3;
% % h31 = yline(ax4,kaidot_max,'--k','linewidth',3);hold(ax3,'on');grid(ax3,'on'); 
% % h34 = yline(ax4,-kaidot_max,'--k','linewidth',3);hold(ax3,'on');grid(ax3,'on');
% % h34.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h32 = plot(ax3,te,kaidot_ode,'r','linewidth',3);hold(ax3,'on');grid(ax3,'on');
% % h31 = plot(ax3,te,kaidot,'--k','linewidth',3);hold(ax3,'on');grid(ax3,'on');
% h32 = plot(ax3,t,kaidot_des,'r','linewidth',3);hold(ax3,'on');grid(ax3,'on');
% % h4.Color = blue;
% h33 = plot(ax3,t,x(:,4),'--b','linewidth',3);hold(ax3,'on');grid(ax3,'on');
% % ax4 = gca;
% ax3.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',3)
% xlabel(ax3,'  $$t$$, s','Fontsize',lbl_fnt);
% % yticks(ax3,[-1.25 -0.65 0 0.65 1.25])
% % yticklabels(ax3,{'-1.25','-0.65','0','0.65','1.25'})
% ylabel(ax3,' $$\dot{\chi} $$, rad/s','Fontsize',lbl_fnt);
% % legend(ax3,'UAV $$\dot{\chi}_{\mathrm{max}}$$','Commanded','Achieved','Fontsize',leg_fnt)
% legend(ax3,'Commanded','Achieved','Fontsize',leg_fnt)
% % axis(ax3,[0 t(end) -1.3 1.3])

% f4;
% % figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
% h4 = plot(ax4,t,kappa_des,'r','linewidth',3);hold(ax4,'on') ;grid(ax4,'on')
% h41 = plot(ax4,t,kappa_actual,'b','linewidth',3);hold(ax4,'on') ;grid(ax4,'on')
% l = legend(ax4,[h4 h41],'Commanded','Achieved');
% l.Orientation = 'Vertical';
% l.Location = 'best';
% l.FontSize = leg_fnt;
% l.Interpreter = 'Latex';
% %Setup axis font size
% % ax4 = gca;     
% ax4.FontSize = ax_fnt;
% % axis(ax4,[0 max(t(end),t(end)) -5 5]);
% %Outer box setup
% box on                      % Switch on the box around the axis
% ax4.XColor = black;         % Box horizontal lines' color
% ax4.YColor = black;         % Box vertical lines' color
% set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% %Label and limit setup
% xlabel(ax4,'  $$t$$, s','Fontsize',lbl_fnt);
% ylabel(ax4,'Curvature $$\kappa$$, m $^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)
% % ax4.XLim = [0 100];
% % ax4.YLim = [-0.05 0.01];
% % 

f5;
% h5 = plot(ax5,t,wrapToPi(kai_des)*(180/pi),'r','linewidth',3);hold(ax5,'on') ;grid(ax5,'on')
% h51 = plot(ax5,t,wrapToPi(x(:,3))*(180/pi),'b','linewidth',3);hold(ax5,'on') ;grid(ax5,'on');
% h5 = plot(ax5,t,kai_des*(180/pi),'r','linewidth',3);hold(ax5,'on') ;grid(ax5,'on')
% h51 = plot(ax5,t,x(:,3)*(180/pi),'b','linewidth',3);hold(ax5,'on') ;grid(ax5,'on');
h5 = plot(ax5,te,kaid_sim*(180/pi),'r','linewidth',3);hold(ax5,'on') ;grid(ax5,'on')
h51 = plot(ax5,te,kai_sim*(180/pi),'b','linewidth',3);hold(ax5,'on') ;grid(ax5,'on');
l = legend(ax5,[h5 h51],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax5.FontSize = ax_fnt;
yticks(ax5,[0 60 90 120 180 240 300 360])                           % Y-axis grid line positioning
yticklabels(ax5,{'0', '60', '90', '120', '180', '240', '300', '360'})
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax5.XColor = black;         % Box horizontal lines' color
ax5.YColor = black;         % Box vertical lines' color
set(ax5,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax5,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax5,' $$ \chi$$,  deg.','Fontsize',lbl_fnt);
% 
f6;
% h6 = plot(ax6,x(:,1),wrapToPi(kai_des)*(180/pi),'r','linewidth',3);hold(ax6,'on') ;grid(ax6,'on')
% h61 = plot(ax6,x(:,1),wrapToPi(x(:,3))*(180/pi),'b','linewidth',3);hold(ax6,'on') ;grid(ax6,'on');
% h6 = plot(ax6,x(:,1),kai_des*(180/pi),'r','linewidth',3);hold(ax6,'on') ;grid(ax6,'on');
% h61 = plot(ax6,x(:,1),x(:,3)*(180/pi),'b','linewidth',3);hold(ax6,'on') ;grid(ax6,'on');
h6 = plot(ax6,x_sim,kaid_sim*(180/pi),'r','linewidth',3);hold(ax6,'on') ;grid(ax6,'on');
h61 = plot(ax6,x_sim,kai_sim*(180/pi),'b','linewidth',3);hold(ax6,'on') ;grid(ax6,'on');
l = legend(ax6,[h6 h61],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
%Setup axis font size   
ax6.FontSize = ax_fnt;
% yticks(ax6,[-180 -120 -60 0 60 120 180])                           % Y-axis grid line positioning
% yticklabels(ax6,{'-180', '-120', '-60', '0', '60', '120', '180'})
% axis(ax6,[0 200 -40 60]);
%Outer box setup
box on                      % Switch on the box around the axis
ax6.XColor = black;         % Box horizontal lines' color
ax6.YColor = black;         % Box vertical lines' color
set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
%Label and limit setup
xlabel(ax6,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax6,'$$ \chi$$,  deg.','Fontsize',lbl_fnt);

f7;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
% h7 = plot(ax7,x(:,1),kappa_des,'r','linewidth',3);hold(ax7,'on') ;grid(ax7,'on')
% h71 = plot(ax7,x(:,1),kappa_actual,'b','linewidth',3);hold(ax7,'on') ;grid(ax7,'on')
h7 = plot(ax7,x_sim,kaidot_ode/Vgd,'r','linewidth',3);hold(ax7,'on') ;grid(ax7,'on')
h71 = plot(ax7,x_sim,kaidot/Vgd,'b','linewidth',3);hold(ax7,'on') ;grid(ax7,'on')
l = legend(ax7,[h7 h71],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
%Setup axis font size
% ax7 = gca;     
ax7.FontSize = ax_fnt;
% axis(ax7,[0 max(t(end),t(end)) -5 5]);
%Outer box setup
box on                      % Switch on the box around the axis
ax7.XColor = black;         % Box horizontal lines' color
ax7.YColor = black;         % Box vertical lines' color
set(ax7,'linewidth',ax_lw)   % Axis linewidth (box and grid)
%Label and limit setup
xlabel(ax7,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax7,'Curvature $$\kappa$$, m $^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)
% % ax7.XLim = [0 100];
% % ax7.YLim = [-0.05 0.01];

f8;
h81 = yline(ax8,kaidot_max,'--b','linewidth',3);hold(ax8,'on');grid(ax8,'on'); 
h82 = yline(ax8,-kaidot_max,'--b','linewidth',3);hold(ax8,'on');grid(ax8,'on');
h82.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h83 = plot(ax8,te,kaidot_ode,'r','linewidth',3);hold(ax8,'on');grid(ax8,'on');
h84 = plot(ax8,te,kaidot,'-r','linewidth',3);hold(ax8,'on');grid(ax8,'on');
% h85 = plot(ax8,t,kaidot_des,'r','linewidth',3);hold(ax8,'on');grid(ax8,'on');
% h4.Color = blue;
% h83 = plot(ax8,t,x(:,4),'--b','linewidth',3);hold(ax8,'on');grid(ax8,'on');
% ax4 = gca;
ax8.FontSize = ax_fnt;
box(ax8, 'on');                       % Switch on the box around the axis
ax8.XColor = 'black';         % Box horizontal lines' color
ax8.YColor = 'black';         % Box vertical lines' color
set(ax8,'linewidth',3)
xlabel(ax8,'  $$t$$, s','Fontsize',lbl_fnt);
% yticks(ax3,[-1.25 -0.65 0 0.65 1.25])
% yticklabels(ax3,{'-1.25','-0.65','0','0.65','1.25'})
ylabel(ax8,' $$\dot{\chi} $$, rad/s','Fontsize',lbl_fnt);
% legend(ax8,'UAV $$\dot{\chi}_{\mathrm{max}}$$','Commanded','Achieved','Fontsize',leg_fnt)
legend(ax8,'$$\dot{\chi}_{\mathrm{max}}$$','$$\dot{\chi}$$','Fontsize',leg_fnt)
% legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt)
axis(ax8,[0 t(end) -1.6 1.6])

f9;
h91 = yline(ax9,kaidot_max/Vgd,'--b','linewidth',3);hold(ax9,'on');grid(ax9,'on'); 
h92 = yline(ax9,-kaidot_max/Vgd,'--b','linewidth',3);hold(ax9,'on');grid(ax9,'on');
h92.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h93 = plot(ax9,te,kaidot_ode/Vgd,'r','linewidth',3);hold(ax9,'on');grid(ax9,'on');
h94 = plot(ax9,te,kaidot/Vgd,'-r','linewidth',3);hold(ax9,'on');grid(ax9,'on');
% h95 = plot(ax8,t,kaidot_des,'r','linewidth',3);hold(ax8,'on');grid(ax8,'on');
% h4.Color = blue;
% h83 = plot(ax8,t,x(:,4),'--b','linewidth',3);hold(ax8,'on');grid(ax8,'on');
% ax4 = gca;
ax9.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax9.XColor = 'black';         % Box horizontal lines' color
ax9.YColor = 'black';         % Box vertical lines' color
set(ax9,'linewidth',3)
xlabel(ax9,'  $$t$$, s','Fontsize',lbl_fnt);
% yticks(ax3,[-1.25 -0.65 0 0.65 1.25])
% yticklabels(ax3,{'-1.25','-0.65','0','0.65','1.25'})
ylabel(ax9,' $$\kappa $$, m$^{-1}$','Fontsize',lbl_fnt);
% legend(ax9,' $$\kappa_{\mathrm{max}}$$','Commanded','Achieved','Fontsize',leg_fnt);
legend(ax9,' $$\kappa_{\mathrm{max}}$$','$$\kappa$$','Fontsize',leg_fnt);
% legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt)
axis(ax9,[0 t(end) -0.16 0.16])



function out = fun_stline_upd(t,x)
global Vgd k2 k_kai k_kaidot kaidot_max te kaidot_ode kaidot kaid_sim  kai_sim x_sim

if t==0
    te = [];   
    kaidot_ode = [];
    kaidot = [];
    kaid_sim = [];
    kai_sim = [];
    x_sim = [];
end


if x(1) < 0
   kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
   kaid_dot = -(2*k2*Vgd.*(x(1)))./((1+k2*(x(1)).^2).^2);
else
   kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
   kaid_dot = - (2*k2*Vgd.*(x(1)))./((1 + k2*(x(1)).^2).^2);
end

% kaid = normalise(kaid);

% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));

% if (kaid_dot < -kaidot_max)
%     kaid_dot = -kaidot_max;
% 
% end
% if (kaid_dot > kaidot_max) 
%     kaid_dot = kaidot_max;    
% end


% out(3,1) = x(4) ;

% kai = normalise( x(3));

if abs(x(4))>kaidot_max
    p = sign(x(4))*kaidot_max;
else
    p = x(4);
end

out(3,1) = p ;

% if x(4) < -kaidot_max 
% 
%     x(4) = -kaidot_max ;
% end
% if x(4) > kaidot_max
%     x(4) = kaidot_max ;
% end
% x(3) = normalise(x(3));

out(4,1) = k_kai*wrapToPi(kaid - x(3)) + k_kaidot*(kaid_dot - p )  ;
% out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - p )  ;
% out(5,1) = kv*(Vad - x(5));

te = [te t];
kaidot_ode = [kaidot_ode kaid_dot];
kaidot = [kaidot p];
kaid_sim = [kaid_sim  kaid];
kai_sim = [kai_sim  x(3)];
x_sim = [x_sim x(1)];
end

% function out = fun_stline_upd(t,x)
% global Vgd k2 wx wy  kv k_kai k_kaidot
% 
% 
% if x(1) < 0
%    kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
%    kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
% else
%    kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
%    kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
% end
% 
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% % psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% c1 = (1/x(5))*kaid_dot*(wx*cos(kaid)+wy*sin(kaid));
% c2 = (wx*sin(kaid)-wy*cos(kaid))*((1/x(5))^2)*(Vgd/Vad)*kaid_dot*(-wx*sin(kaid)+wy*cos(kaid));
% psiddot = kaid_dot + (1/sqrt(1 - ((1/x(5))*(wx*sin(kaid)-wy*cos(kaid)))^2))*(c1 + c2);
% % psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% % Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));
% 
% 
% out(1,1) = x(5)*cos(x(3))  + wx ;
% out(2,1) = x(5)*sin(x(3))  + wy ;
% out(3,1) = x(4) ;
% out(4,1) = k_kai*(psid - x(3)) + k_kaidot*(psiddot - x(4))  ;
% out(5,1) = kv*(Vad - x(5));
% end

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

function y = normalise(x)
y = atan2(sin(x),cos(x));
% y = atan(x);
end