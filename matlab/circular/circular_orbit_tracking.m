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

tspan = [0 30];                         % simulation run time

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


%% plot format control variables
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

% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20;       % Label font size
leg_fnt = 16;       % Legend font size
ax_lw = 3;          % Axis line width      
        
%-------------- figures --------------------------------------------------%

 figure(1)
% fig = gcf;
% fig.WindowState = 'maximized';
xp = plot(x_sc,y_sc,'k','linewidth',lw);hold on;grid on;
if r_0 < rd
  hvf = quiver(X_in,Y_in,xdot_in,ydot_in,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
else
    hvf = quiver(X_out,Y_out,xdot_out,ydot_out,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
end
h_traj = plot(x(:,3),x(:,4),'color',maroon,'linewidth',lw);hold on;
h1c = plot(x_c,y_c,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','magenta'); hold on;
h1c.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12 = plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10); hold on;
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax1 = gca;
ax1.FontSize = ax_fnt;
box on
ax1.XColor = 'black';
ax1.YColor = 'black';
set(ax1,'linewidth',ax_lw)
xlabel(ax1,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax1,'$$ y, $$  m','Fontsize',lbl_fnt);
legend(ax1,'Desired orbit','Vector field','UAV trajectory','','','Fontsize',leg_fnt);
axis(ax1,'equal')

%% Figure 2: Radial distance vs time
figure(2)
fig = gcf;
fig.WindowState = 'maximized';
ax2 = gca;
plot(t,r_prop,'color',maroon,'LineWidth',lw);
ax2.FontSize = ax_fnt;
box on
ax2.XColor = 'black';
ax2.YColor = 'black';
set(ax2,'linewidth',ax_lw)
xlabel(ax2,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax2,'$$ r, $$ m','Fontsize',lbl_fnt);
grid on;

%% Figure 3: Radial error vs time
figure(3)
fig = gcf;
fig.WindowState = 'maximized';
ax3 = gca;
plot(t,error_r,'color',maroon,'LineWidth',lw);
ax3.FontSize = ax_fnt;
box on
ax3.XColor = 'black';
ax3.YColor = 'black';
set(ax3,'linewidth',ax_lw)
xlabel(ax3,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax3,'$$ r - r_d, $$ m','Fontsize',lbl_fnt);
grid on;

%% Figure 4: Course angle vs time
figure(4)
fig = gcf;
fig.WindowState = 'maximized';
ax4 = gca;
plot(t,wrapToPi(kaid)*(180/pi),'color',maroon,'LineWidth',lw);hold on;
plot(t,wrapToPi(x(:,5))*(180/pi),'--','color',blue,'LineWidth',lw);
ax4.FontSize = ax_fnt;
box on
ax4.XColor = 'black';
ax4.YColor = 'black';
set(ax4,'linewidth',ax_lw)
xlabel(ax4,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax4,'$$ \\chi, $$ deg','Fontsize',lbl_fnt);
legend('Desired','Actual','Fontsize',leg_fnt);
grid on;

%% Figure 5: Curvature vs time
figure(5)
fig = gcf;
fig.WindowState = 'maximized';
ax5 = gca;
plot(t,kappa_actual,'color',maroon,'LineWidth',lw);
ax5.FontSize = ax_fnt;
box on
ax5.XColor = 'black';
ax5.YColor = 'black';
set(ax5,'linewidth',ax_lw)
xlabel(ax5,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax5,'$$ \\kappa, $$ m$$^{-1}$$','Fontsize',lbl_fnt);
grid on;

%% Figure 6: Animation with Fixed-Wing UAV Marker
figure(6)
fig = gcf;
fig.WindowState = 'maximized';
ax6 = gca;
xp = plot(x_sc,y_sc,'k','linewidth',lw);hold on;grid on;
if r_0 < rd
  quiver(X_in,Y_in,xdot_in,ydot_in,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
else
    quiver(X_out,Y_out,xdot_out,ydot_out,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
end
plot(x_c,y_c,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','magenta'); hold on;
plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
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

legend('Desired orbit','Vector field','','','Trajectory','','Fontsize',leg_fnt);

% Setup video writer
video_filename = sprintf('circular_orbit_animation_%s.mp4', datestr(now,''));
vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 10;
vidObj.Quality = 95;
open(vidObj);

% Animation loop
dt_anim = 0.1;
t_anim = 0:dt_anim:t(end);
for i = 1:length(t_anim)
    % Interpolate position and heading at current animation time
    x_curr = interp1(t, x(:,3), t_anim(i));
    y_curr = interp1(t, x(:,4), t_anim(i));
    chi_curr = interp1(t, x(:,5), t_anim(i));
    
    % Update trajectory trail
    idx = find(t <= t_anim(i));
    set(h_trail, 'XData', x(idx,3), 'YData', x(idx,4));
    
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
trajectory_data.position_x = x(:,3);
trajectory_data.position_y = x(:,4);
trajectory_data.radial_distance = x(:,1);
trajectory_data.los_angle = x(:,2);
trajectory_data.course_angle = x(:,5);
trajectory_data.course_angle_rate = x(:,6);
trajectory_data.course_angle_desired = kaid;
trajectory_data.course_angle_rate_desired = kaid_dot;
trajectory_data.curvature_actual = kappa_actual;
trajectory_data.curvature_desired = kappa_des;
trajectory_data.initial_conditions = x_initial;
trajectory_data.parameters.rd = rd;
trajectory_data.parameters.Vgd = Vgd;
trajectory_data.parameters.k2 = k2;
trajectory_data.parameters.k_kai = k_kai;
trajectory_data.parameters.k_kaidot = k_kaidot;
trajectory_data.vector_field.X = X_out;
trajectory_data.vector_field.Y = Y_out;
trajectory_data.vector_field.xdot = xdot_out;
trajectory_data.vector_field.ydot = ydot_out;
trajectory_data.vector_field.kappa = kappa_des;

% Save with timestamp
filename = sprintf('circular_orbit_trajectory.mat');
save(filename, 'trajectory_data');
fprintf('Trajectory data saved to: %s\n', filename);

%% Create GIF animation: UAV moving on vector field

gif_name = 'circular_path_following.gif';

figure(10); clf;
set(gcf,'Color','w','Position',[100 100 900 650]);
ax_gif = gca;
hold(ax_gif,'on');
grid(ax_gif,'on');
box(ax_gif,'on');

% Plot fixed background once
plot(ax_gif,x_sc,y_sc,'k','linewidth',lw); hold(ax_gif,'on');
quiver(ax_gif,X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5); hold(ax_gif,'on');
plot(ax_gif,x_c,y_c,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','magenta'); hold(ax_gif,'on');
plot(ax_gif,x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold(ax_gif,'on');

xlabel(ax_gif,' $ x, $  m','Fontsize',lbl_fnt);
ylabel(ax_gif,' $ y, $  m','Fontsize',lbl_fnt);
% title(ax_gif,'Circular Orbit Following','Fontsize',lbl_fnt);

ax_gif.FontSize = ax_fnt;
ax_gif.XColor = 'black';
ax_gif.YColor = 'black';
set(ax_gif,'linewidth',ax_lw)
axis(ax_gif,'equal');
xlim(ax_gif,[-100 100]);
ylim(ax_gif,[-100 100]);

% animated objects
traj_line = plot(ax_gif,nan,nan,'color',maroon,'linewidth',lw);

% Create fixed-wing UAV shape (triangle pointing in direction of motion)
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav_gif = fill(ax_gif,nan, nan, green, 'EdgeColor', black, 'LineWidth', 2);

legend(ax_gif,'Desired orbit','Vector field','','UAV trajectory','', ...
    'Fontsize',leg_fnt,'Location','northwest');

% choose fewer frames for smaller GIF
skip = 8;

for k = 1:skip:length(t)
    
    set(traj_line,'XData',x(1:k,3),'YData',x(1:k,4));
    
    % Get current heading
    chi_curr = x(k,5);
    
    % Rotate and translate UAV shape
    R = [cos(chi_curr), -sin(chi_curr); sin(chi_curr), cos(chi_curr)];
    uav_rotated = R * [uav_shape_x; uav_shape_y];
    uav_x = uav_rotated(1,:) + x(k,3);
    uav_y = uav_rotated(2,:) + x(k,4);
    
    set(h_uav_gif, 'XData', uav_x, 'YData', uav_y);
    
    drawnow;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    
    if k == 1
        imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.1);
    end
end

fprintf('GIF animation saved to: %s\n', gif_name);

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