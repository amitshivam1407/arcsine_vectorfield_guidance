%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
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
% f7 = figure;
% ax7 = axes;
% f8 = figure;
% ax8 = axes;
% Plot format control variables
lw = 3;            % Line width
ms = 6;            % Marker size
ax_fnt = 23;        % Axis font size
lbl_fnt = 25; % Label font size
leg_fnt = 17; % Legend font size
ax_lw = 4;        % Axis line width
% Colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
orange = [0.9290 0.6940 0.1250];
violet = [0.4940 0.1840 0.5560];
green = [0 1 0];
cyan = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840];
black = [0 0 0];
%--------------- declaration of constant parameters ----------------------%
global  kv wx wy k_kai k_kaidot k2 figure_type_good

wx = -1;
wy = 2;
kv = 5;
k_kai = 30;  
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
figure_type_good = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 25];
Vgd = 10;
for i = 1:length(Vgd)
k2 = 0.002;
x0 = 100;
y0 = -100;

if x0 < 0
   kai0 =  asin(1./(1+ k2 *(x0).^2)) ;
   kaidot0 = - (2*k2*Vgd.*(x0))./((1 + k2*(x0).^2).^2);
else
   kai0 = pi - asin(1./(1 + k2 *(x0).^2)) ;  
   kaidot0 = - (2*k2*Vgd.*(x0))./((1 + k2*(x0).^2).^2);
end
kai0_deg = rad2deg(wrapToPi(kai0));
Va0 = 8;

x_initial = [x0;y0;kai0;kaidot0;Va0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -110:10:110;
[X,Y] = meshgrid(range);
[kaid,kaid_dot] = vf_proposed(X,Vgd(1));
xdot  = Vgd(1)*cos(kaid);
ydot  = Vgd(1)*sin(kaid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@(t,x)fun_stline_withwind(t,x,Vgd(i)) ,tspan, x_initial,options);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vg_act = sqrt(x(:,5).^2 + wx^2 + wy^2 + 2*wx*x(:,5).*cos(psi) + 2*wy*x(:,5).*sin(psi));
x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;
% for i = 1:length(t)
% if x(i,1) < 0
%    kai_des(i,1) =  asin(1./(1+ k2*(x(i,1)).^2)) ;
%    kaidot_des(i,1) = -2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
% else
%    kai_des(i,1) = pi - asin(1./(1 + k2*(x(i,1)).^2)) ;  
%    kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
% end
% end
[kai_des,kaidot_des] = fun_kaidstline(x(:,1),Vgd(i));
kappa_des = kaidot_des./Vgd(i);
Vad = sqrt((Vgd(i).*cos(kai_des) - wx).^2 + (Vgd(i)*sin(kai_des) - wy).^2) ;

psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
Vg = sqrt((x(:,5).*cos(psi) + wx).^2 + (x(:,5).*sin(psi) + wy).^2) ;
kaidot_actual = x(:,4);
kappa_actual = kaidot_actual./Vg;
%---------------------- plotting figures ---------------------------------%

%% Figure 1: Trajectory with Vector Field
figure(1)
fig = gcf;
fig.WindowState = 'maximized';
xp = xline(0,'-k','linewidth',3);hold on;grid on;
quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
h_traj = plot(x(:,1),x(:,2),'color',maroon,'linewidth',lw);hold on;
h12 = plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold on;
h13 = plot(x_end,y_end,'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',cyan); hold on;
h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax1 = gca;
ax1.FontSize = ax_fnt;
box on
ax1.XColor = 'black';
ax1.YColor = 'black';
set(ax1,'linewidth',ax_lw)
xlabel(ax1,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax1,'$$ y, $$  m','Fontsize',lbl_fnt);
legend(ax1,'Desired straight line $$x = 0$$','Vector field','Trajectory','Initial position','Fontsize',leg_fnt);
axis(ax1,'equal')

%% Figure 2: Cross-track error vs time
figure(2)
fig = gcf;
fig.WindowState = 'maximized';
ax2 = gca;
plot(t,x(:,1),'color',maroon,'LineWidth',lw);
ax2.FontSize = ax_fnt;
box on
ax2.XColor = 'black';
ax2.YColor = 'black';
set(ax2,'linewidth',ax_lw)
xlabel(ax2,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax2,'$$ x, $$ m','Fontsize',lbl_fnt);
grid on;

%% Figure 3: Course angle rate vs cross-track error
figure(3)
fig = gcf;
fig.WindowState = 'maximized';
ax3 = gca;
plot(x(:,1),kaidot_des,'color',maroon,'LineWidth',lw);hold on;
plot(x(:,1),x(:,4),'--','color',blue,'LineWidth',lw);
ax3.FontSize = ax_fnt;
box on
ax3.XColor = 'black';
ax3.YColor = 'black';
set(ax3,'linewidth',ax_lw)
xlabel(ax3,'$$ x, $$ m','Fontsize',lbl_fnt);
ylabel(ax3,'$$ \dot{\chi}, $$ rad/s','Fontsize',lbl_fnt);
legend('Desired','Actual','Fontsize',leg_fnt);
grid on;

%% Figure 4: Curvature vs cross-track error
figure(4)
fig = gcf;
fig.WindowState = 'maximized';
ax4 = gca;
plot(x(:,1),kappa_des,'color',maroon,'LineWidth',lw);hold on;
plot(x(:,1),kappa_actual,'--','color',blue,'LineWidth',lw);
ax4.FontSize = ax_fnt;
box on
ax4.XColor = 'black';
ax4.YColor = 'black';
set(ax4,'linewidth',ax_lw)
xlabel(ax4,'$$ x, $$ m','Fontsize',lbl_fnt);
ylabel(ax4,'$$ \kappa, $$ m$$^{-1}$$','Fontsize',lbl_fnt);
legend('Desired','Actual','Fontsize',leg_fnt);
grid on;

%% Figure 5: Course angle vs time
figure(5)
fig = gcf;
fig.WindowState = 'maximized';
ax5 = gca;
plot(t,kai_des*180/pi,'color',maroon,'LineWidth',lw);hold on;
plot(t,x(:,3)*180/pi,'--','color',blue,'LineWidth',lw);
ax5.FontSize = ax_fnt;
box on
ax5.XColor = 'black';
ax5.YColor = 'black';
set(ax5,'linewidth',ax_lw)
xlabel(ax5,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax5,'$$ \chi, $$ deg','Fontsize',lbl_fnt);
legend('Desired','Actual','Fontsize',leg_fnt);
grid on;

%% Figure 6: Animation with Fixed-Wing UAV Marker
figure(6)
fig = gcf;
fig.WindowState = 'maximized';
ax6 = gca;
xline(0,'-k','linewidth',3);hold on;grid on;
quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold on;
ax6.FontSize = ax_fnt;
box on
ax6.XColor = 'black';
ax6.YColor = 'black';
set(ax6,'linewidth',ax_lw)
xlabel(ax6,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax6,'$$ y, $$  m','Fontsize',lbl_fnt);
axis(ax6,'equal')

% Initialize trajectory trail and UAV marker
h_trail = plot(nan, nan, 'color', maroon, 'linewidth', lw);

% Create fixed-wing UAV shape (triangle pointing in direction of motion)
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav = fill(nan, nan, green, 'EdgeColor', black, 'LineWidth', 2);

legend('Desired line $$x = 0$$','Vector field','Initial position','Trajectory','UAV','Fontsize',leg_fnt);

% Setup video writer
video_filename = sprintf('wind_casestudy_animation.mp4');
vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 10;
vidObj.Quality = 95;
open(vidObj);

% Animation loop
dt_anim = 0.1;
t_anim = 0:dt_anim:t(end);
for i = 1:length(t_anim)
    % Interpolate position and heading at current animation time
    x_curr = interp1(t, x(:,1), t_anim(i));
    y_curr = interp1(t, x(:,2), t_anim(i));
    chi_curr = interp1(t, x(:,3), t_anim(i));
    
    % Update trajectory trail
    idx = find(t <= t_anim(i));
    set(h_trail, 'XData', x(idx,1), 'YData', x(idx,2));
    
    % Rotate and translate UAV shape
    R = [cos(chi_curr), -sin(chi_curr); sin(chi_curr), cos(chi_curr)];
    uav_rotated = R * [uav_shape_x; uav_shape_y];
    uav_x = uav_rotated(1,:) + x_curr;
    uav_y = uav_rotated(2,:) + y_curr;
    
    set(h_uav, 'XData', uav_x, 'YData', uav_y);
    
    title(ax6, sprintf('Time: %.2f s', t_anim(i)), 'Fontsize', lbl_fnt);
    drawnow;
    
    % Capture frame for video
    frame = getframe(gcf);
    writeVideo(vidObj, frame);
end

% Close video writer
close(vidObj);
fprintf('Animation saved to: %s\n', video_filename);

disp('Animation complete!');

%% Save trajectory data to .mat file
trajectory_data.time = t;
trajectory_data.position_x = x(:,1);
trajectory_data.position_y = x(:,2);
trajectory_data.course_angle = x(:,3);
trajectory_data.course_angle_rate = x(:,4);
trajectory_data.air_speed = x(:,5);
trajectory_data.course_angle_desired = kai_des;
trajectory_data.course_angle_rate_desired = kaidot_des;
trajectory_data.curvature_actual = kappa_actual;
trajectory_data.curvature_desired = kappa_des;
trajectory_data.initial_conditions = x_initial;
trajectory_data.parameters.Vgd = Vgd;
trajectory_data.parameters.k2 = k2;
trajectory_data.parameters.k_kai = k_kai;
trajectory_data.parameters.k_kaidot = k_kaidot;
trajectory_data.parameters.wx = wx;
trajectory_data.parameters.wy = wy;
trajectory_data.parameters.kv = kv;
trajectory_data.vector_field.X = X;
trajectory_data.vector_field.Y = Y;
trajectory_data.vector_field.xdot = xdot;
trajectory_data.vector_field.ydot = ydot;

% Save with timestamp
filename = sprintf('wind_casestudy_trajectory.mat');
save(filename, 'trajectory_data');
fprintf('Trajectory data saved to: %s\n', filename);
end

function out = fun_stline_withwind(t,x,Vgd )
global k2 wx wy  kv k_kai k_kaidot


if x(1) < 0
   kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
   kaid_dot = -2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
else
   kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
   kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
end
 
Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
%  psid = kai + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
 Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vg.*cos(x(3));
out(2,1) = Vg.*sin(x(3));
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
out(5,1) = kv*(Vad - x(5));
end


function [kaid,kaid_dot] = vf_proposed(x,Vgd)
global k2 
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

function [kai_des,kaidot_des] = fun_kaidstline(x,Vgd)
global k2 
for i = 1:length(x)
if x(i,1) < 0
   kai_des(i,1) =  asin(1./(1+ k2*(x(i,1)).^2)) ;
   kaidot_des(i,1) = -(2*k2*Vgd.*(x(i,1)))./((1+k2*(x(i,1)).^2).^2);
else
   kai_des(i,1) = pi - asin(1./(1 + k2*(x(i,1)).^2)) ;  
   kaidot_des(i,1) = - (2*k2*Vgd.*(x(i,1)))./((1 + k2*(x(i,1)).^2).^2);
end
end
end
