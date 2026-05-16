%-------------------------------------------------------------------------%
%------               12th May 2024                         --------------%
%-------     Proposed vector field for sinusoidal path following  --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% f1 = figure;
% ax1 = axes;
% f2 = figure;
% ax2 = axes;
% f3 = figure;
% ax3 = axes;
% f4 = figure;
% ax4 = axes;
% f5 = figure;
% ax5 = axes;
% f6 = figure;
% ax6 = axes;
% f7 = figure;
% ax7 = axes;
% f8 = figure;
% ax8 = axes;
% f9 = figure;
% ax9 = axes;
% f11 = figure;
% ax11 = axes;
% f12 = figure;
% ax12 = axes;
% f13 = figure;
% ax13 = axes;
% f14 = figure;
% ax14 = axes;

%%%%%%%%%%%%%%%% declaration of constant parameters %%%%%%%%%%%%%%%%%%%%%%%
global  Vgd k_kai k_kaidot k2 flag1 t_flag1 flag2 t_flag2 k_aseem omega A
omega = 0.06; 
A = 5;  
Vgd = 10;
k_kai = 80;
k_kaidot = 50;

%-------------------------------------------------------------------------%
%---------------- Desired sinusoidal path---------------------------------%
xd = -130:.005:130;
yd = A*sin(omega.*xd);

%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 30];

x0 = -90;
y0 = -80;
epsilon0 =  A*sin(omega*x0) - y0;
% k_prop = [2*10^(-5) 2.4*10^(-4) 4.6*10^(-4) 6.8*10^(-4) 9*10^(-4) 0.0011 0.0013 0.0016 0.0018 0.0020 0.003 0.004 0.005 0.006 0.007 0.008];
k_prop = 0.005;
k2 = k_prop ;
flag1 = true(length(k_prop));
t_flag1 = zeros(1,length(k_prop));
flag2 = true(length(k_prop));
t_flag2 = zeros(1,length(k_prop));

% for i = 1:length(k_prop)

[kaid0,kaid_dot0] = get_initialvalue(x0,y0);


kai0_deg = rad2deg(wrapToPi(kaid0));

k_aseem = sqrt(k_prop^2 * (epsilon0)^2 + 2*k_prop);

[kaid0_aseem,kaid_dot0_aseem] = get_initialvalue_aseem(x0,y0);

x_initial = [x0;y0;kaid0;kaid_dot0];
y_initial = [x0;y0;kaid0_aseem;kaid_dot0_aseem];

%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range = -130:5:130;
                             % semi minor axis
[X,Y] = meshgrid(range);
dy_ddx_des = A*omega*cos(omega.*X);
kai_p_des  = atan(dy_ddx_des);
epsilon    =  A*sin(omega.*X)-Y;
[kaid]     = vf_sinusoidal(X,Y);
Xdot       = Vgd*cos(kaid);
Ydot       = Vgd*sin(kaid);
%------
%-------------------------------------------------------------------------%

option_x=odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(t,x) stopping_prop(t,x));
option_y=odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(t,y) stopping_aseem(t,y));

% %% ode solver equation  %%
[t,x, te, ze, ie] = ode45(@(t,x) fun_prop(t,x,k_prop), tspan, x_initial, option_x);
[t1,y, te1, ze1, ie1] = ode45(@(t1,y) fun_aseem_et_al(t1,y), tspan, y_initial, option_y);
%----------------------  parameters to plot   ----------------------------%
% beta_xx = x(:,1).^2 /a^2 + x(:,2).^2 /b^2 ;
% error_prop =  beta_xx - 1;

if (~isempty(te))
% error_stop = error_prop;%% ze contain all states to apply event in any of the desired %% 
t_stop_prop = te; %% time when selected state achieves the required criteria  %%
end
% 
% 
% beta_yy = y(:,1).^2 /a^2 + y(:,2).^2 /b^2 ;
% error_aseem =  beta_yy - 1;
% 
if (~isempty(te1))
% error_yy_stop = error_aseem;%% ze contain all states to apply event in any of the desired %% 
t_stop_aseem = te1; %% time when selected state achieves the required criteria  %%
end

% f=@(t,x) fun_prop(t,x,k_prop);
% 
% t0 = 0; tf = tspan(end); dt = 0.01;
% [x,t]=rk45(f,x_initial,t0,tf,dt);
% t = t';
% x=x.';
% 
% fy=@(t,y) fun_aseem_et_al(t,y);
% 
% t10 = 0; tf = tspan(end); dt = 0.01;
% [y,t1]=rk45(fy,y_initial,t10,tf,dt);
% t1 = t1';
% y  = y.';


error_prop =  A*sin(omega.*x(:,1)) - x(:,2);
% Vg_act = sqrt(x(:,7).^2 + wx^2 + wy^2 + 2*wx*x(:,7).*cos(x(:,5)) + 2*wy*x(:,7).*sin(x(:,5)));
x_ini = x(1,1)  ;
y_ini = x(1,2)  ;
x_end = x(end,1)  ;
y_end = x(end,2)  ;


kappa_actual_prop = x(:,4)./Vgd ;
control_effort_prop = trapz(t,(Vgd*x(:,4)).^2);

% [val_max_prop ,ind_prop] = max(abs(kappa_actual_prop));
% exx_max_prop = error_prop(ind_prop,1);
% [kaidot_max_prop ,ind_prop] = max(abs(x(:,4)));
% t_max_prop = t(ind_prop,1);

[val_max_prop ,ind_prop] = min((kappa_actual_prop));
exx_max_prop = error_prop(ind_prop,1);
[kaidot_max_prop ,ind_prop] = min((x(:,4)));
t_max_prop = t(ind_prop,1);

%----------------------- Aseem et al ------------------------------------%



error_aseem =  A*sin(omega.*y(:,1)) - y(:,2);
kappa_actual_aseem = y(:,4)./Vgd ;
control_effort_aseem = trapz(t1,(Vgd*y(:,4)).^2);


% [val_max_aseem, ind_aseem_kappa] = max(abs(kappa_actual_aseem));
% eyy_max_aseem = error_aseem(ind_aseem_kappa,1);

[val_max_aseem, ind_aseem_kappa] = min((kappa_actual_aseem));
eyy_max_aseem = error_aseem(ind_aseem_kappa,1);

[kaidot_max_aseem ,ind_aseem] = min((y(:,4)));
t_max_aseem = t1(ind_aseem,1);

% [kaidot_max_aseem ,ind_aseem] = max(abs(y(:,4)));
% t_max_aseem = t1(ind_aseem,1);


percent_red_curvature = ((val_max_aseem - val_max_prop)./val_max_aseem) *100;
percent_red_controleffort = ((control_effort_aseem - control_effort_prop)./control_effort_aseem) *100;
percent_red_settlingtime = -((t_stop_aseem - t_stop_prop)./t_stop_aseem) *100;



%-------------------------------------------------------------------------%
%% plot format control variables
lw = 3;            % Line width
ms = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 19; % Label font size
leg_fnt = 15; % Legend font size
ax_lw = 3;        % Axis line width
% Colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
% blue = [0 0 1];
% red = [1 0 0];
orange = [0.9290 0.6940 0.1250];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
cyan = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840];
black = [0 0 0];

figure(1)
% fig = gcf;
% fig.WindowState = 'maximized';
% fig.Position = [10 50 1000 950];
% fig.Position = figure('Position', [20 20 1920 1080]);
% fig.Position =  [20 20 1920 1080];
xp = plot(xd,yd,'k','linewidth',lw);hold on;grid on;
% quiver(X,Y,Xdot,Ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
plot(x(:,1),x(:,2),'color',maroon,'linewidth',lw);hold on;
% plot(x(:,5),x(:,6),'color',blue,'linewidth',lw);hold on;
plot(y(:,1),y(:,2),'color',blue,'linewidth',lw);hold on;
h12 = plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10); hold on;
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h13 = plot(x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
% h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax1 = gca;
ax1.FontSize = ax_fnt;
% Outer box setup
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax1,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax1,'$$ y, $$  m','Fontsize',lbl_fnt);
% legend(ax1,'Vector field','Elliptic orbit','Proposed','Griffiths et al 2006. 2019','Fontsize',leg_fnt);
legend(ax1,'Sinusoidal path','Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)
axis(ax1,'equal')


figure(2)
% fig = gcf;
% fig.WindowState = 'maximized';
plot(t,error_prop,'r','linewidth',3);hold on ;grid on;
plot(t1,error_aseem,'b','linewidth',3);hold on ;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax2,'  $$t$$, s','Fontsize',lbl_fnt);
ylabel(ax2,'  $$\epsilon $$,  m','Fontsize',lbl_fnt);
legend(ax2,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)
% legend('Desired error ','Fontsize',14);
% axis(ax2,'tight')
% set(gca, 'XLim', [t(1,1) t(end,1)],'YLim', [error_prop(end)-0.01 error_prop(1)]);
% hold(ax2,'on'); grid(ax2,'on');


figure(3)
% fig = gcf;
% fig.WindowState = 'maximized';
h3 = plot(t,x(:,4),'r','linewidth',3);hold on ;grid on;
% plot(t,x(:,8)./Vgd,'b','linewidth',3);hold on ;grid on;
h31 = plot(t1,y(:,4),'b','linewidth',3);hold on ;grid on;
dt_prop3 = datatip(h3,t_max_prop,kaidot_max_prop);hold on ;grid on;
% h4.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% h4.DataTipTemplate.DataTipRows(2).Format = '%.3f';
dt_aseem3 = datatip(h31,t_max_aseem,kaidot_max_aseem);hold on ;grid on;
% h41.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% h41.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% plot([t_max_prop t_max_prop],[0 val_max_prop/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
% plot([t_max_aseem t_max_aseem],[0 val_max_aseem/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax3,'  $$t$$, s','Fontsize',lbl_fnt);
ylabel(ax3,'  $$\dot{\chi} $$,  rad/s','Fontsize',lbl_fnt);
legend(ax3,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)

figure(4)
% fig = gcf;
% fig.WindowState = 'maximized';
plot(t,wrapToPi(x(:,3)),'r','linewidth',3);hold on ;grid on;
% plot(t,x(:,8)./Vgd,'b','linewidth',3);hold on ;grid on;
plot(t1,wrapToPi(y(:,3)),'b','linewidth',3);hold on ;grid on;
ax4 = gca;
ax4.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax4,'  $$t$$, s','Fontsize',lbl_fnt);
ylabel(ax4,'  $$\chi $$,  rad','Fontsize',lbl_fnt);
legend(ax4,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)


figure(5)
% fig = gcf;
% fig.WindowState = 'maximized';
plot(error_prop,x(:,3),'r','linewidth',3);hold on ;grid on;
% plot(t,x(:,8)./Vgd,'b','linewidth',3);hold on ;grid on;
plot(error_aseem,y(:,3),'b','linewidth',3);hold on ;grid on;
ax5 = gca;
ax5.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax5,' $$\epsilon $$,  m','Fontsize',lbl_fnt);
ylabel(ax5,'  $$\chi $$,  rad','Fontsize',lbl_fnt);
legend(ax5,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)


figure(6)
% fig = gcf;
% fig.WindowState = 'maximized';
h6 = plot(t,x(:,4)./Vgd,'r','linewidth',3);hold on ;grid on;
% plot(t,x(:,8)./Vgd,'b','linewidth',3);hold on ;grid on;
h61 = plot(t1,y(:,4)./Vgd,'b','linewidth',3);hold on ;grid on;
dt_prop6 = datatip(h6,t_max_prop,kaidot_max_prop);hold on ;grid on;
h6.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h6.DataTipTemplate.DataTipRows(2).Format = '%.3f';
dt_aseem6 = datatip(h61,t_max_aseem,kaidot_max_aseem);hold on ;grid on;
h61.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h61.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% plot([t_max_prop t_max_prop],[0 val_max_prop/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
% plot([t_max_aseem t_max_aseem],[0 val_max_aseem/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
ax6 = gca;
ax6.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax6,'  $$t$$, s','Fontsize',lbl_fnt);
ylabel(ax6,'  $$\kappa $$,  m$^{-1}$','Fontsize',lbl_fnt);
legend(ax6,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)


figure(7)
% fig = gcf;
% fig.WindowState = 'maximized';
h7 = plot(error_prop,x(:,4)./Vgd,'r','linewidth',3);hold on ;grid on;
% plot(t,x(:,8)./Vgd,'b','linewidth',3);hold on ;grid on;
h71 = plot(error_aseem,y(:,4)./Vgd,'b','linewidth',3);hold on ;grid on;
dt_prop = datatip(h7,exx_max_prop,val_max_prop);hold on ;grid on;
h7.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h7.DataTipTemplate.DataTipRows(2).Format = '%.3f';
dt_aseem = datatip(h71,eyy_max_aseem,val_max_aseem);hold on ;grid on;
h71.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h71.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% plot([t_max_prop t_max_prop],[0 val_max_prop/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
% plot([t_max_aseem t_max_aseem],[0 val_max_aseem/Vgd],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
ax7 = gca;
ax7.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax7.XColor = 'black';         % Box horizontal lines' color
ax7.YColor = 'black';         % Box vertical lines' color
set(ax7,'linewidth',3)        % Axis linewidth (box and grid)
xlabel(ax7,'  $$\epsilon $$,  m','Fontsize',lbl_fnt);
ylabel(ax7,'  $$\kappa $$,  m$^{-1}$','Fontsize',lbl_fnt);
legend(ax7,'Proposed','Griffiths et al. 2006','Fontsize',leg_fnt)



% set(gca, 'XLim', [t(1,1) t(end,1)],'YLim', [-0.1 0.1]);
% hold(ax3,'on'); grid(ax3,'on');
% % set(gca, 'XLim', [x(end,1) x(1,1)],'YLim', [-0.25 0]);


% %% Animate Trajectory (used variables: x, y, si)
% isACmarker = true;
% markerFileName = 'plane.jpeg';
% delayTime = 0.5;
% 
% noLOS = min(200,length(x));
% jumpStep = floor(length(x)/noLOS);
% % jumpStep = 5;
% 
% %% Setup Video creation
% 
% vidFlnm = 'Plane_marker Elliptic_Trajectory';
% % vidFlnm = 'Tracking error_ellipse';
% % vidFlnm = 'Elliptic_Curvature';
% 
% if(ischar(vidFlnm))
% 
%     vid = VideoWriter(vidFlnm,'MPEG-4');
%     vid.Quality = 100;
%     vid.FrameRate = 15;
% 
%     open(vid)
% 
% end

% Plot the initial coordinates
% if(isACmarker)
%     % Find the initial coordinates of the marker WITHOUT orientation
%     [xm0, ym0] = setupPlotMarker(markerFileName, 1/4, 0.5);
% 
%     % Find the initial coordinates of the marker WITH orientation
%     [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(1,3));
% %     [xmT, ymT] = calcMarkerCoord(xm0,ym0,x(1,7));
% 
%     h1 = plot(ax1,xmM+x(1,1),ymM+x(1,2),'s','MarkerSize',0.5,'Color','red','MarkerFaceColor','red');
% %     h2 = plot(ax1,xmT+x(1,9),ymT+x(1,10),'o','MarkerSize',0.5,'Color','blue','MarkerFaceColor','red');
% else
%     h1 = plot(ax1,x(1,1),x(1,2),'ob','MarkerSize',8,'MarkerFaceColor','blue');
% %    h2 = plot(ax1,x(1,9),x(1,10),'or','MarkerSize',8,'MarkerFaceColor','red');
% end

%----- trajectory comparison----------
% h11 = plot(ax1,x(1,1),x(1,2),'-r','LineWidth',4);
% h12 =  plot(ax1,x(1,9),x(1,10),'-b','LineWidth',4);
% legend(ax1,[xp h11 h12],'Desired path','Proposed','Nelson et al.','Fontsize',17);
% legend(ax1,[xp hvf h11],'Desired path','Vector field','UAV trajectory','Fontsize',17);

%----------- Tracking error comparison--------
% h21 = plot(ax2,t(1,1),error_prop(1,1),'-r','LineWidth',4);
% h22 =  plot(ax2,x(1,5),x(1,7),'-b','LineWidth',4);
% l = legend(ax2,[h21 h22],'Proposed','Nelson et al.','Fontsize',17,'Location','southeast');

%----------- Curvature comparison--------
% h31 = plot(ax3,t(1,1),x(1,4)/Vgd,'-r','LineWidth',4);
%  h32 =  plot(ax3,x(1,7),x(1,12)/Vgd,'-b','LineWidth',4);
%  legend(ax3,[h31 h32],'Proposed','Nelson et al.','Fontsize',17,'Location','southeast');

%% Animate the marker

% for i = 1+jumpStep:jumpStep:length(x)
%     %   for i = 1:length(x)
%     pause(delayTime)
% 
%            set(gca, 'XLim', [-100 100],'YLim', [-100 100]);
% %            set(gca, 'XLim', [t(1,1) t(end,1)],'YLim', [error_prop(end,1)-0.01 error_prop(1,1)]);
% %       set(gca, 'XLim', [t(1,1) t(end,1)],'YLim', [-0.1 0.1]);
% 
%     if i==1
% 
%              h11(i) = plot(ax1,x(i-jumpStep:i,1),x(i-jumpStep:i,2),'-r','LineWidth',4);
% %              h12(i) =  plot(ax1,x(i-jumpStep:i,9),x(i-jumpStep:i,10),'-b','LineWidth',4);
% 
% %                h21(i) = plot(ax2,t(i-jumpStep:i,1),error_prop(i-jumpStep:i,1),'-r','LineWidth',4);
% %              h22(i) =  plot(ax2,x(i-jumpStep:i,5),x(i-jumpStep:i,7),'-b','LineWidth',4);
% 
% %          h31(i) = plot(ax3,t(i-jumpStep:i,1),x(i-jumpStep:i,4)./Vgd,'-r','LineWidth',4);
% %         h32(i) =  plot(ax3,x(i-jumpStep:i,7),x(i-jumpStep:i,12)./Vgd,'-b','LineWidth',4);
% 
%     else
%              h11(i) = plot(ax1,x(i-jumpStep:i,1),x(i-jumpStep:i,2),'-r','LineWidth',4);
%             h11(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
% %            h12(i) =  plot(ax1,x(i-jumpStep:i,9),x(i-jumpStep:i,10),'-b','LineWidth',4);
% %            h12(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% %             h21(i) = plot(ax2,t(i-jumpStep:i,1),error_prop(i-jumpStep:i,1),'-r','LineWidth',4);
% %            h21(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
% %            h22(i) =  plot(ax2,x(i-jumpStep:i,5),x(i-jumpStep:i,7),'-b','LineWidth',4);
% %            h22(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% %          h31(i) = plot(ax3,t(i-jumpStep:i,1),x(i-jumpStep:i,4)./Vgd,'-r','LineWidth',4);
% %          h31(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
% %         h32(i) =  plot(ax3,x(i-jumpStep:i,7),x(i-jumpStep:i,12)./Vgd,'-b','LineWidth',4);
% %         h32(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
%         %     legend(ax1,'off');
%     end
% %     pos=get(gca,'Position');
% %     if i>=777
% %         h33 = plot(ax3,x_max_prop,val_max_prop,'-s','MarkerSize',20,...
% %             'MarkerEdgeColor','black',...
% %             'MarkerFaceColor','cyan'); hold(ax3,'on');grid(ax3,'on');
% %         h33.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %         h33.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% %         h33.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% %         dt_prop = datatip(h33,x_max_prop,val_max_prop);hold(ax3,'on');grid(ax3,'on');
% % %         a1 = annotation('textarrow',[pos(1),pos(1)-2*pos(1)],[pos(2),pos(1)-2*pos(4)],'String','$$\kappa_{\mathrm{max}}$$ ');hold(ax3,'on');grid(ax3,'on');
% % %         a1.Color = 'red';
% % %         a1.FontSize = 17;
% % %         a1.Interpreter = 'latex';
% %     end
% % 
% %     if i >=835
% %         h34 = plot(ax3,x_max_Nelson,val_max_Nelson,'-s','MarkerSize',20,...
% %             'MarkerEdgeColor','black',...
% %             'MarkerFaceColor','green'); hold(ax3,'on');grid(ax3,'on');
% %         h34.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %         h34.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% %         h34.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% %         dt_Nelson = datatip(h34,x_max_Nelson,val_max_Nelson);hold(ax3,'on');grid(ax3,'on');
% % %         a2 = annotation('textarrow',[0.3,0.5],[0.6,0.5]);hold(ax3,'on');grid(ax3,'on');
% % %         a2.Color = 'blue';
% % %         a2.FontSize = 17;
% %     end
%     % h41.DataTipTemplate.DataTipRows(1).Format = '%.2f';
%     % h41.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% 
%     %     if i==9468
%     %         plot(ax3,'oc','MarkerSize',40,'MarkerFaceColor','red');hold on;
%     %         annotation('textarrow',x(i,1),x(i,4)./Vgd,'String','$$\kappa_{\mathrm{max}}$$ ');hold on;
%     %     end
%     %
%     %     if i==10013
%     %         plot(ax3,'sm','MarkerSize',40,'MarkerFaceColor','blue');hold on;
%     %         annotation('textarrow',x(i,5),x(i,8)./Vgd,'String','$$\kappa_{\mathrm{max}}$$ ');hold on;
%     %     end
% 
%         % Marker Position Update
%         if(isACmarker)
%             [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(i,3));
% %             [xmT, ymT] = calcMarkerCoord(xm0,ym0,x(i,11));
%     
%             set(h1,'XData',xmM+x(i,1),'YData',ymM+x(i,2))
% %             set(h2,'XData',xmT+x(i,9),'YData',ymT+x(i,10))
%         else
%             set(h1,'XData',x(i,1),'YData',x(i,2))
% %             set(h2,'XData',x(i,9),'YData',x(i,10))
%         end
%     %
%     drawnow
% 
%     if(ischar(vidFlnm))
%         
%         F = getframe(gcf);
%         writeVideo(vid,F)
% 
%     end
%      if ~ishghandle(fig)
%         break
%     end
%     ...
% 
% end

% plot(ax1,x(i-jumpStep:end,1),x(i-jumpStep:end,2),'-r','LineWidth',4)
% plot(ax1,x(i-jumpStep:end,5),x(i-jumpStep:end,6),'-b','LineWidth',4)
% l = legend(ax1,'Fontsize',17);
% l = legend(ax1,[xl h11 h12],'Desired path','Proposed','Nelson et al','Fontsize',17);
% if(isACmarker)
%     [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(end,3));
%     [xmT, ymT] = calcMarkerCoord(xm0,ym0,x(end,7));
%     set(h1,'XData',xmM+x(end,1),'YData',ymM+x(end,2))
%     set(h2,'XData',xmT+x(end,5),'YData',ymT+x(end,6))
% else
%     set(h1,'XData',x(end,1),'YData',x(end,2))
%     set(h2,'XData',x(end,5),'YData',x(end,6))
% end

% F = getframe(gcf);
% writeVideo(vid,F)
% 
% close(vid)

% figure(2)
% plot(t,r_prop,'color',maroon,'LineWidth',lw);hold on;grid on;
% plot(t,r_Nelson,'color',blue,'LineWidth',lw);hold on;grid on;
% ax2 = gca;
% ax2.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax2.XColor = 'black';         % Box horizontal lines' color
% ax2.YColor = 'black';         % Box vertical lines' color
% set(ax2,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax2,' $$t$$,  s','Fontsize',lbl_fnt)
% ylabel(ax2,'Radial distance $$r$$,  m','Fontsize',lbl_fnt);
% legend(ax2,'Proposed','Nelson et al. ','Fontsize',leg_fnt);
% 
% figure(3)
% plot(t,error_r,'color',maroon,'LineWidth',lw);hold on;grid on;
% plot(t,error_Nelson,'color',blue,'LineWidth',lw);hold on;grid on;
% ax3 = gca;
% ax3.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax3,' $$t$$, s','Fontsize',lbl_fnt)
% ylabel(ax3,'Radial  error $$(r - r_d)$$, m','Fontsize',lbl_fnt)
% legend(ax3,'Proposed','Nelson et al. ','Fontsize',leg_fnt);
% 
% figure(4)
% plot(t,wrapToPi(x(:,2))*(180/pi),'color',maroon,'linewidth',lw);hold on;grid on;
% plot(t,wrapToPi(x(:,8))*(180/pi),'color',blue,'linewidth',lw);hold on;grid on;
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax4,'$$ t,$$ s','Fontsize',lbl_fnt)
% ylabel(ax4,'Bearing  angle $ \gamma$,   deg.','Fontsize',lbl_fnt)
% legend(ax4,'Proposed','Nelson et al. ','Fontsize',leg_fnt);

% figure(5)
% plot(r_prop,wrapToPi(kaid_dot_prop)*(180/pi),'r','LineWidth',lw);hold on;grid on;
% plot(r_prop,wrapToPi(x(:,6))*(180/pi),'b','LineWidth',lw);grid on;
% xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',lbl_fnt)
% ylabel('Course \ rate,  \ $\dot{\chi}$, \ ( $ ^\circ $/s) ','Fontsize',lbl_fnt)
% legend('Desired','Actual','Fontsize',leg_fnt);

% figure(6)
% plot(r_prop,kappa_des_prop,'r','LineWidth',lw);hold on;grid on;
% plot(r_prop,kappa_actual_prop,'b','LineWidth',lw);hold on;grid on;
% xlabel('Radial \ distance, \ $r$, \ (m)','Fontsize',lbl_fnt)
% ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',lbl_fnt)
% legend('Desired','Actual','Fontsize',leg_fnt);

% figure(7)
% plot(t,wrapToPi(kaid_prop)*(180/pi),'r','linewidth',lw);hold on;grid on;
% plot(t,wrapToPi(x(:,5))*(180/pi),'b','linewidth',lw);grid on;
% xlabel('Time, \ $ t,$ \ (s)','Fontsize',lbl_fnt)
% ylabel('Course \ angle, \ $ \chi$, \ $ (^\circ) $','Fontsize',lbl_fnt)
% legend('Desired','Actual','Fontsize',leg_fnt);

% figure(8)
% plot(t,kappa_des_prop,'r','LineWidth',lw);hold on;grid on;
% plot(t,kappa_actual_prop,'b','LineWidth',lw);hold on;grid on;
% xlabel('Time, \ $ t,$ \ (s)','Fontsize',lbl_fnt)
% ylabel('Curvature, \ $\kappa$, \ (m $ ^-1 $) ','Fontsize',lbl_fnt)
% legend('Desired','Actual','Fontsize',leg_fnt);

% figure(9)
% plot(t,control_input_prop,'r','linewidth',lw);hold on;grid on;
% plot(t,control_input_Nelson,'b','linewidth',lw);hold on;grid on;
% ax9 = gca;
% ax9.FontSize = ax_fnt;
% xlabel('Time, \ $ t,$ \ (s)','Fontsize',lbl_fnt)
% ylabel('Control \ input, \ $ u$, \  (units) ','Fontsize',lbl_fnt)
% legend('Proposed','Nelson \ et al. 2007','Fontsize',leg_fnt);

% figure(10)
% plot(t,kappa_actual_prop,'color',maroon,'LineWidth',lw);hold on;grid on;
% plot(t,kappa_actual_Nelson,'color',blue,'LineWidth',lw);hold on;grid on;
% ax10 = gca;
% ax10.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax10.XColor = 'black';         % Box horizontal lines' color
% ax10.YColor = 'black';         % Box vertical lines' color
% set(ax10,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax10,' $$ t,$$ s','Fontsize',lbl_fnt)
% ylabel(ax10,'Curvature,  $\kappa$,  m $ ^{-1} $ ','Fontsize',lbl_fnt)
% legend(ax10,'Proposed','Nelson  et al. ','Fontsize',leg_fnt);
% 
% figure(11)
% h = plot(r_prop,kappa_actual_prop,'color',maroon,'LineWidth',lw);hold on;grid on;
% h1 = plot(r_Nelson,kappa_actual_Nelson,'color',blue,'LineWidth',lw);hold on;grid on;
% plot(x_max_prop,val_max_prop,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','cyan'); hold on;grid on;
% plot(x_max_Nelson,val_max_Nelson,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','green'); hold on;grid on;
% % plot([x_max_prop x_max_prop],[0 val_max_prop],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
% % plot([x_max_Nelson x_max_Nelson],[0 val_max_Nelson],'--','color',[0.75 0.75 0.75],'linewidth',lw);hold on;grid on;
% dt_prop = datatip(h,x_max_prop,val_max_prop);hold(ax4,'on');grid(ax4,'on');
% h.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% h.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% dt_Nelson = datatip(h1,x_max_Nelson,val_max_Nelson);hold(ax4,'on');grid(ax4,'on');
% h1.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% h1.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% ax11 = gca;
% ax11.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax11.XColor = 'black';         % Box horizontal lines' color
% ax11.YColor = 'black';         % Box vertical lines' color
% set(ax11,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax11,'Radial  distance,  $r$,  m','Fontsize',lbl_fnt);
% ylabel(ax11,'Curvature,  $\kappa$, m $ ^{-1} $ ','Fontsize',lbl_fnt)
% % legend('Proposed','Nelson \ et al. 2007','Fontsize',leg_fnt);
% legend(ax11,'Proposed','Nelson et al. ','','','Fontsize',leg_fnt);
% 
% figure(12)
% plot(t,x(:,6),'color',maroon,'LineWidth',lw);hold on;grid on;
% plot(t,x(:,12),'color',blue,'LineWidth',lw);hold on;grid on;
% ax12 = gca;
% ax12.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax12.XColor = 'black';         % Box horizontal lines' color
% ax12.YColor = 'black';         % Box vertical lines' color
% set(ax12,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax12,'Time,  $ t,$  s','Fontsize',lbl_fnt)
% ylabel(ax12,'Course  rate,   $\dot{\chi}$,  rad./s ','Fontsize',lbl_fnt)
% legend(ax12,'Proposed','Nelson et al. ','Fontsize',leg_fnt);
% 
% figure(13)
% plot(r_prop,x(:,6),'color',maroon,'LineWidth',lw);hold on;grid on;
% plot(r_Nelson,x(:,12),'color',blue,'LineWidth',lw);hold on;grid on;
% ax13 = gca;
% ax13.FontSize = ax_fnt;
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax13.XColor = 'black';         % Box horizontal lines' color
% ax13.YColor = 'black';         % Box vertical lines' color
% set(ax13,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% xlabel(ax13,'Radial  distance,  $r$,  m','Fontsize',lbl_fnt);
% ylabel(ax13,'Course  rate,   $\dot{\chi}$,   rad./s ','Fontsize',lbl_fnt)
% legend(ax13,'Proposed','Nelson  et al. ','Fontsize',leg_fnt);
% end

% function out = fun_proposed_nowind(t,x,k_prop,i)
% global k_kai k_kaidot rd Vgd flag1 t_flag1 
% 
% if (flag1(i) == true) && (x(1)>(rd + 0.2))
%     t_flag1(i) = t;
%     flag1(i) = false;
% end
% 
% if x(1) < rd
%     kaid = x(2) + asin(1./(1+k_prop*(x(1) -rd).^2)) ;
%     gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
%     kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1)-rd))./((1+k_prop*(x(1) -rd).^2).^2);
% else
%     kaid = x(2) + pi - asin(1./(1+k_prop*(x(1) -rd).^2)) ;
%     gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
%     kaid_dot = gamma_dot - (2*k_prop*Vgd.*(x(1) -rd))./((1+k_prop*(x(1) -rd).^2).^2);
% end
% 
% 
% out(1,1) = Vgd.*cos(x(5) - x(2));
% out(2,1) = (Vgd./x(1)).*sin(x(5) - x(2));
% out(3,1) = Vgd*cos(x(5))  ;
% out(4,1) = Vgd*sin(x(5))  ;
% out(5,1) = x(6) ;
% out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;
% 
% end

function out = fun_aseem_et_al(t1,y)
global k_kai k_kaidot  Vgd flag2 t_flag2 k_aseem A omega

out(1,1) =  Vgd*cos(y(3))  ;
out(2,1) =  Vgd*sin(y(3))  ;


dydx = A*omega*cos(omega.*y(1));
kai_p =  atan(dydx);
epsilon =  A*sin(omega.*y(1)) - y(2);
kaio1 = atan(k_aseem.*(epsilon )) ;
kait1 = atan(dydx);
kaid_aseem = kait1  +  kaio1  ;
xdot1 = out(1,1);
ydot1 = out(2,1);
kait_dot1 = -(A*omega^2 *sin(omega*y(1))*xdot1)./(1 + (tan(kait1)).^2);
epsilon_dot = A*omega*cos(omega*y(1))*xdot1 - ydot1;
kaio_dot1 = ((k_aseem)./((1 + (tan(kaio1)).^2))).*epsilon_dot;
kaid_dot_aseem = kait_dot1 + kaio_dot1 ;

% kaid_Nelson = y(8) + pi/2 + atan(k_Nelson.*(y(7) - rd)) ;
% gamma_dot_Nelson = (Vgd./y(7) ).*sin(kaid_Nelson - y(8)) ;
% kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(y(7) -rd))./((1+k_Nelson^2*(y(7) -rd).^2).^(3/2));
% 


out(3,1) = y(4) ;
out(4,1) = k_kai*wrapToPi(kaid_aseem - y(3)) + k_kaidot*(kaid_dot_aseem - y(4) )  ;


 end

 function out = fun_prop(t,x,k2)
global  k_kai k_kaidot Vgd flag1 t_flag1 flag2 t_flag2 k_aseem A omega 

% if (flag1 == true) && (x(1)>( 0.2))
%     t_flag1 = t;
%     flag1 = false;
% end

dydx = A*omega*cos(omega.*x(1));
kai_p =  atan(dydx);
epsilon =  A*sin(omega.*x(1)) - x(2);
kai_o = pi/2 - asin(1./(1 + k2*(epsilon).^2)) ;
out(1,1) = Vgd*cos(x(3)) ;
out(2,1) = Vgd*sin(x(3)) ;
if epsilon < 0
    kaid = (kai_p) - kai_o;
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = -(A*omega^2 *sin(omega*x(1))*xdot)./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = A*omega*cos(omega*x(1))*xdot - ydot;
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot - kaio_dot;
else
    kaid = kai_p + kai_o ;
    xdot = out(1,1);
    ydot = out(2,1);
    kaip_dot = -(A*omega^2 *sin(omega*x(1))*xdot)./(1 + (tan(kai_p)).^2);
    factor_1 = (2*k2)/((1 + k2*(epsilon).^2).*(sqrt(2*k2 + (k2*epsilon).^2)));
    epsilon_dot = A*omega*cos(omega*x(1))*xdot - ydot;
    kaio_dot = factor_1*epsilon_dot;
    kaid_dot = kaip_dot + kaio_dot;
end

out(3,1) = x(4);
out(4,1) = k_kai*wrapToPi(kaid - x(3)) + k_kaidot*(kaid_dot - x(4)) ;

% out(5,1) =  Vgd*cos(x(7))  ;
% out(6,1) = Vgd*sin(x(7))  ;
% 
% beta1 = x(5)^2/a^2 + x(6)^2/b^2 ;
% fx1 = 2*x(5)./a^2 ;
% fxx1 = 2/a^2 ;
% fy1 = 2*x(6)./b^2 ;
% fyy1 = 2./b^2 ;
% kaio1 = atan(k_aseem.*(beta1 - 1)) ;
% kait1 = atan2(fx1,-fy1);
% kaid_aseem = kait1  +  kaio1  ;
% xdot1 = out(5,1);
% ydot1 = out(6,1);
% kait_dot1 = -((fy1.*fxx1.*xdot1 - fx1.*fyy1.*ydot1)./fy1^2).*(cos(kait1)).^2 ;
% beta_dot1 = fx1.*xdot1 + fy1.*ydot1 ;
% kaio_dot1 = ((k_aseem)./((1 + (tan(kaio1)).^2))).*beta_dot1;
% kaid_dot_aseem = kait_dot1 + kaio_dot1 ;
% 
% % kaid_Nelson = x(8) + pi/2 + atan(k_Nelson.*(x(7) - rd)) ;
% % gamma_dot_Nelson = (Vgd./x(7) ).*sin(kaid_Nelson - x(8)) ;
% % kaid_dot_Nelson = gamma_dot_Nelson - (k_Nelson^2*Vgd.*(x(7) -rd))./((1+k_Nelson^2*(x(7) -rd).^2).^(3/2));
% % 
% 
% 
% out(7,1) = x(8) ;
% out(8,1) = k_kai*wrapToPi(kaid_aseem - x(7)) + k_kaidot*(kaid_dot_aseem - x(8) )  ;

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

function [value,isterminal,direction] = stopping_prop(t,x)
global  Vgd A omega
% value(1)      = ((x(1).^2 /a^2 + x(2).^2 /b^2) - 1) - 0.0001; 
value(1)      = (A*sin(omega.*x(1)) - x(2))-0.0001; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1)  = 0; % negative direction(as R decreases from positive to zero d=-1;
                  % If R increases from negative to zero d=+1;d=0 implies no need of direction )
end

function [value,isterminal,direction] = stopping_aseem(t,y)
global  Vgd A omega
% value(1)      = ((y(1).^2 /a^2 + y(2).^2 /b^2) - 1) - 0.0001; 
value(1)      = (A*sin(omega.*y(1)) - y(2))-0.0001; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1)  = 0; % negative direction(as R decreases from positive to zero d=-1;
                  % If R increases from negative to zero d=+1;d=0 implies no need of direction )
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
function [kaid0,kaid_dot0] = get_initialvalue(x0,y0)
global Vgd k2 omega A


epsilon0 =  A*sin(omega*x0) - y0;
dydx0 = A*omega*cos(omega.*x0);
kai_p0 = atan(dydx0);
kai_o0 =  pi/2 - asin(1./(1 + k2*(epsilon0).^2));

if epsilon0 < 0
    kaid0          = (kai_p0) - kai_o0;
    xdot0         = Vgd*cos(kaid0);
    ydot0         = Vgd*sin(kaid0);
    kaip_dot0     = -(A*omega^2*sin(omega*x0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k2)/((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0  =  A*omega*cos(omega*x0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaid_dot0       = kaip_dot0 - kaio_dot0;
else
    kaid0          = kai_p0 + kai_o0 ;
    xdot0         = Vgd*cos(kaid0);
    ydot0         = Vgd*sin(kaid0);
    kaip_dot0     = -(A*omega^2*sin(omega*x0).*xdot0)./(1 + (tan(kai_p0))^2);
    factor_10     = (2*k2)/((1 + k2*(epsilon0).^2).*(sqrt(2*k2 + (k2*epsilon0).^2)));
    epsilon_dot0  = A*omega*cos(omega*x0)*xdot0 - ydot0;
    kaio_dot0     = factor_10*epsilon_dot0;
    kaid_dot0       = kaip_dot0 + kaio_dot0;
end
end

function [kaid0_aseem,kaid_dot0_aseem] = get_initialvalue_aseem(x0,y0)
global Vgd k_aseem omega A


epsilon0 =  A*sin(omega*x0) - y0;
dydx0 = A*omega*cos(omega.*x0);
kait0 = atan(dydx0);
kaio0 =  atan(k_aseem*(epsilon0)) ;
kaid0_aseem = kait0  +  kaio0  ;
xdot0 = Vgd*cos(kaid0_aseem);
ydot0 = Vgd*sin(kaid0_aseem);
kaip_dot0     = -(A*omega^2*sin(omega*x0).*xdot0)./(1 + (tan(kait0))^2);
factor_10     = k_aseem./(1 + (tan(kaio0)).^2) ;
epsilon_dot0  = A*omega*cos(omega*x0)*xdot0 - ydot0;
kaio_dot0     = factor_10*epsilon_dot0;
kaid_dot0_aseem = kaip_dot0 + kaio_dot0 ;

end