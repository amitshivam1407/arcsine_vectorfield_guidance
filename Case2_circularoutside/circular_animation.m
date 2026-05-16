%-------------------------------------------------------------------------%
%------               12th Setember 2022                         ---------%
%-------- Circular path with Autopilot design with no wind        --------%
%-------   Proposed vector field for oustide initial position     --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%---------------- declaration of constant parameters ---------------------%
global rd Vgd k_kai k_kaidot k2 x_c y_c
rd = 50;                            % desired radius
Vgd = 10;                           % ground speed
k_kai = 100;                        % autopilot proportional gain
k_kaidot = sqrt(k_kai);             % autopilot derivative gain 

red = [0.8500 0.3250 0.0980];       % color for plotting
blue = [0 0.4470 0.7410];           % color for plotting

%%%%%%%%%%%%%%%%%%%%%% desired circle    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_c = 0;                        % x-coordinate of the centre of desired circle
y_c = 0;                        % y-coordinate of the centre of desired circle
theta = 0:2*pi/1000:2*pi;       % desired circular path angular paramter
x_sc = x_c + rd*cos(theta);     % x-coordinate of the desired circular path
y_sc = y_c + rd*sin(theta);     % y-coordinate of the desired circular path
%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%

tspan = [0 25];                         % simulation run time

r_0 = 100;                              % initial radial distance in m
gamma_0 = deg2rad(225);                 % initial LOS angle in rad
x0 = r_0*cos(gamma_0) + x_c;            % initial UAV position x-coordinate in m
y0 = r_0*sin(gamma_0) + y_c;            % initial UAV position y-coordinate in m
% k2 = 14.9260./(r_0 - rd).^2;
k2 = 0.005;                             % Guidance gain
if r_0 < rd
   kai0 = gamma_0 + asin(1./(1+k2*(r_0-rd).^2)) ;      % Initial heading in rad
else
   kai0 = gamma_0 + pi - asin(1./(1+k2*(r_0-rd).^2)) ; % Initial heading in rad
end
gamma_dot0 = (Vgd./r_0).*sin(kai0 - gamma_0) ; % LOS rate
kaidot0 = gamma_dot0 - (2*k2*Vgd.*(r_0-rd))./((1+k2*(r_0-rd).^2).^2);  % Initial heading rate in rad/s

kai0_deg = rad2deg(wrapToPi(kai0));            % Initial heading in deg
x_initial = [r_0;gamma_0;x0;y0;kai0;kaidot0];  % initial condition array

%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
range  = -100:7:100;
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
[t,x] = ode45(@(t,x)fun_circ_autopilot_nowind(t,x) ,tspan, x_initial,options);

%----------------------- parameters to plot   ----------------------------%

r_prop = x(:,1) ;
error_r = (r_prop - rd) ;
% Vg_act = sqrt(x(:,7).^2 + wx^2 + wy^2 + 2*wx*x(:,7).*cos(x(:,5)) + 2*wy*x(:,7).*sin(x(:,5)));
x_ini = x(1,3) + x_c ;
y_ini = x(1,4) + y_c ;
x_end = x(end,3) + x_c ;
y_end = x(end,4) + y_c ;

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
kappa_des = kaid_dot./Vgd ;
%  Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%  psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))./Vad);
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% kappa_des = 1./(r_prop.*(1 + k2.*(r_prop - rd).^2)) - (2*k2.*(r_prop - rd))./((1 + k2.*(r_prop - rd).^2).^2);
kappa_actual = x(:,6)./Vgd ;
[val, ind] = min(kappa_actual);
rmax = x(ind,1);
control_input = k_kai*(kaid - x(:,5)) + k_kaidot*(kaid_dot - x(:,6))  ;

% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20;       % Label font size
leg_fnt = 16;       % Legend font size
ax_lw = 3;          % Axis line width      
        
%-------------- figures --------------------------------------------------%

% figure(1)
% if r_0 < rd
%     quiver(X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% else
%     quiver(X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
% end
% plot(x_sc,y_sc,'k','linewidth',lw_);hold on;
% plot(x(:,3),x(:,4),'color',red,'linewidth',lw_);hold on;
% plot(x_c,y_c,'-o','linewidth',2,'MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','red'); hold on;
% plot(x_ini,y_ini,'-o','linewidth',2,'MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','green'); hold on;
% plot(x_end,y_end,'-s','linewidth',2,'MarkerSize',10,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','cyan'); hold on;
% ax1 = gca;
% ax1.FontSize = ax_fnt;
% ax1.XColor = 'black';         % Box horizontal lines' color
% ax1.YColor = 'black';         % Box vertical lines' color
% set(ax1,'linewidth',3)
% xlabel(ax1,' $ x, $  m','Fontsize',lbl_fnt);
% ylabel(ax1,'$ y, $  m','Fontsize',lbl_fnt);
% legend(ax1,'Vector field','Circular orbit','UAV trajectory','','','','Fontsize',leg_fnt)
% axis(ax1,'equal')
% 
% figure(2)
% plot(t,r_prop,'color',red,'LineWidth',lw_);grid on;
% ax2 = gca;
% ax2.FontSize = ax_fnt;
% ax2.XColor = 'black';         % Box horizontal lines' color
% ax2.YColor = 'black';         % Box vertical lines' color
% set(ax2,'linewidth',3)
% xlabel(ax2,' $t$,  s','Fontsize',lbl_fnt)
% ylabel(ax2,'Radial distance $$ r $$, m','Fontsize',lbl_fnt)
% 
% figure(3)
% plot(t,error_r,'color',red,'LineWidth',lw_);grid on;
% ax3 = gca;
% ax3.FontSize = ax_fnt;
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',3)
% xlabel(ax3,' $t$,  s','Fontsize',lbl_fnt)
% ylabel(ax3,'Radial error $$(r - r_{\mathrm{d}})$$, m','Fontsize',lbl_fnt)
% 
% figure(4)
% plot(error_r,wrapToPi(kaid)*(180/pi),'-','color',red,'linewidth',lw_);hold on;grid on;
% plot(error_r,wrapToPi(x(:,5))*(180/pi),'--','color',blue,'linewidth',lw_);grid on;
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',3)
% xlabel(ax4,'Radial error $$(r - r_{\mathrm{d}})$$, m','Fontsize',lbl_fnt)
% ylabel(ax4,' $ \chi$,  deg.','Fontsize',lbl_fnt)
% 
% figure(5)
% plot(r_prop,kaid_dot,'-','color',red,'LineWidth',lw_);hold on;grid on;
% plot(r_prop,x(:,6),'--','color',blue,'LineWidth',lw_);grid on;
% ax5 = gca;
% ax5.FontSize = ax_fnt;
% ax5.XColor = 'black';         % Box horizontal lines' color
% ax5.YColor = 'black';         % Box vertical lines' color
% set(ax5,'linewidth',3)
% xlabel(ax5,'Radial  distance, m','Fontsize',lbl_fnt)
% ylabel(ax5,'$\dot{\chi}$, rad/s','Fontsize',lbl_fnt)
% legend(ax5,'Commanded','Achieved','Fontsize',leg_fnt);
% 
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
% 
% figure(7)
% plot(t,wrapToPi(kaid)*(180/pi),'-','color',red,'linewidth',lw_);hold on;grid on;
% plot(t,wrapToPi(x(:,5))*(180/pi),'--','color',blue,'linewidth',lw_);grid on;
% ax7 = gca;
% ax7.FontSize = ax_fnt;
% ax7.XColor = 'black';         % Box horizontal lines' color
% ax7.YColor = 'black';         % Box vertical lines' color
% set(ax7,'linewidth',3);
% % set(ax7,'XTick',[-180, -120, -60, 0, 60, 120, 180],'xticklabel',{'-180' ,'-120', '-60', '0', '60', '120', '180'});
% % xticks(ax7,[-180 -120 -60 0 60 120 180]);
% % xticklabels(ax7,{'-180' ,'-120', '-60', '0', '60', '120', '180'});
% xlabel(ax7,' $ t,$ s','Fontsize',lbl_fnt)
% ylabel(ax7,' $ \chi$, deg. ','Fontsize',lbl_fnt)
% legend(ax7,'Commanded','Achieved','Fontsize',leg_fnt);
% 
% figure(8)
% % plot(t,kappa_des,'-','color',red,'LineWidth',2);hold on;grid on;
% plot(t,kappa_actual,'-','color',red,'LineWidth',lw_);hold on;grid on;
% ax8 = gca;
% ax8.FontSize = ax_fnt;
% ax8.XColor = 'black';         % Box horizontal lines' color
% ax8.YColor = 'black';         % Box vertical lines' color
% set(ax8,'linewidth',3) 
% xlabel(ax8,' $ t,$ s','Fontsize',lbl_fnt)
% ylabel(ax8,'Curvature $$\kappa $$, m $ ^{-1} $ ','Fontsize',lbl_fnt)
% % legend(ax8,'Commanded','Achieved','Fontsize',leg_fnt);
% 
% figure(9)
% plot(t,control_input,'color',red,'linewidth',lw_);grid on;
% box on  
% ax9 = gca; 
% ax9.FontSize = ax_fnt;
% % Switch on the box around the axis
% ax9.XColor = 'black';         % Box horizontal lines' color
% ax9.YColor = 'black';         % Box vertical lines' color
% set(ax9,'linewidth',3) 
% xlabel(ax9,'$$ t,$$  s','Fontsize',lbl_fnt)
% ylabel(ax9,'Control \ input,  rad./s ','Fontsize',lbl_fnt)
% 
% figure(10)
% plot(t,kaid_dot,'-','color',red,'LineWidth',lw_);hold on;grid on;
% plot(t,x(:,6),'--','color',blue,'LineWidth',lw_);grid on;
% ax10 = gca;
% ax10.FontSize = ax_fnt;
% ax10.XColor = 'black';         % Box horizontal lines' color
% ax10.YColor = 'black';         % Box vertical lines' color
% set(ax10,'linewidth',3)
% xlabel(ax10,'$$ t,$$  s','Fontsize',lbl_fnt)
% ylabel(ax10,'$\dot{\chi}$, rad/s','Fontsize',lbl_fnt)
% legend(ax10,'Commanded','Achieved','Fontsize',leg_fnt);

% figure(11)
fig = gcf;
fig.WindowState = 'maximized';
% fig.Position = [10 50 1000 950];
% fig.Position = figure('Position', [20 20 1920 1080]);
% fig.Position =  [20 20 1920 1080];
% xp = plot(x_sc,y_sc,'k','linewidth',lw);hold on;grid on;
xp = plot(x_sc,y_sc,'k','linewidth',lw_);hold on;
if r_0 < rd
  hvf =  quiver(X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
else
    hvf = quiver(X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold on;
end

% plot(x(:,3),x(:,4),'color',maroon,'linewidth',lw);hold on;
% plot(x(:,9),x(:,10),'color',blue,'linewidth',lw);hold on;
h1c = plot(x_c,y_c,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
     'MarkerFaceColor','magenta'); hold on;
h1c.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12 = plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10); hold on;
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h13 = plot(x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
% h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax11 = gca;
ax11.FontSize = ax_fnt;
% Outer box setup
box on                      % Switch on the box around the axis
ax11.XColor = 'black';         % Box horizontal lines' color
ax11.YColor = 'black';         % Box vertical lines' color
set(ax11,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax11,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax11,'$$ y, $$  m','Fontsize',lbl_fnt);
% legend('Vector field','Circular \ orbit','Proposed','Nelson \ et al. 2007','Centre of orbit','Initial \ position','Final \ position','Fontsize',leg_fnt);
% legend(ax11,'Desired circular orbit','Vector field','Proposed','Fontsize',leg_fnt);hold(ax1,'on');grid(ax1,'on');
axis(ax11,'equal')



%% Setup Video creation

delayTime = 0.02;

noLOS = min(200,length(x));
% jumpStep = floor(length(x)/noLOS);
jumpStep = 10;

vidFlnm = 'Circular_Trajectory1';


if(ischar(vidFlnm))

    vid = VideoWriter(vidFlnm,'MPEG-4');
    vid.Quality = 100;
    vid.FrameRate = 15;

    open(vid)

end

% % %% Plot the initial coordinates

h1 = plot(ax11,x(1,3),x(1,4),'-s','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
     'MarkerFaceColor','cyan');

% %----- trajectory animation ----------
h11 = plot(ax11,x(1,3),x(1,4),'-r','LineWidth',4);
legend(ax11,[xp hvf h11],'Desired circular orbit','Vector field','UAV trajectory','Fontsize',17);
axis(ax11,'equal')

% %% Animate the marker

for i = 1+jumpStep:jumpStep:length(x)
    %   for i = 1:length(x)
    pause(delayTime)

           set(gca, 'XLim', [-105 105],'YLim', [-105 105]);

     if i==1
             h11(i) = plot(ax11,x(i-jumpStep:i,3),x(i-jumpStep:i,4),'-r','LineWidth',4); 
     else
            h11(i) = plot(ax11,x(i-jumpStep:i,3),x(i-jumpStep:i,4),'-r','LineWidth',4);
           h11(i).Annotation.LegendInformation.IconDisplayStyle = 'off';

     end

            set(h1,'XData',x(i,3),'YData',x(i,4))

    drawnow

    if(ischar(vidFlnm))
        
        F = getframe(gcf);
        writeVideo(vid,F)

    end
     if ~ishghandle(fig)
        break
    end
    ...

end



F = getframe(gcf);
writeVideo(vid,F)

close(vid)


function out = fun_circ_autopilot_nowind(t,x)
global  k_kai k_kaidot rd k2 Vgd


if x(1) < rd
   kaid = x(2) + asin(1./(1+k2*(x(1) -rd).^2)) ;
   gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
   kaid_dot = gamma_dot + 2*k2*Vgd.*(x(1)-rd)./((1+k2*(x(1) -rd).^2).^2);
else
   kaid = x(2) + pi - asin(1./(1+k2*(x(1) -rd).^2)) ;
   gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
   kaid_dot = gamma_dot - 2*k2*Vgd.*(x(1) -rd)./((1+k2*(x(1) -rd).^2).^2);
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