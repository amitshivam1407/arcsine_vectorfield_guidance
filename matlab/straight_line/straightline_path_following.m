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
global  Vgd kv k_kai k_kaidot k2 figure_type_good
Vgd = 10;
kv = 5;
k_kai = 50;  
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
figure_type_good = 1;

%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%

tspan = [0 25];

x0 = 100;
y0 = -100;
k2 = 0.002;
if x0 < 0
   kaid0 =  asin(1./(1+ k2 *(x0).^2)) ;
%    kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
   kaid0 = pi - asin(1./(1 + k2 *(x0).^2)) ;  
%   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
kai0 = kaid0;
kaidot0 = -(2*k2*Vgd.*x0)./((1+k2*(x0).^2).^2);

% k2 = abs(14.92597/x0^2);

x_initial = [x0;y0;kai0;kaidot0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -110:10:110;
[X,Y] = meshgrid(range);
[kaid,kaid_dot] = vf_proposed(X);
xdot  = Vgd*cos(kaid);
ydot  = Vgd*sin(kaid);

kaid_mod = vf_proposed_mod(X);
xdot_mod  = Vgd*cos(kaid_mod);
ydot_mod  = Vgd*sin(kaid_mod);

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
   kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
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

%% Figure 1: Trajectory with Vector Field
figure(1)
fig = gcf;
% fig.WindowState = 'maximized';
xp = xline(0,'-k','linewidth',3);hold on;grid on;
if x0 < 0
    hvf = quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
else
    hvf = quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
end
h_traj = plot(x(:,1),x(:,2),'color',maroon,'linewidth',lw);hold on;
h12 = plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold on;
h13 = plot(x_end,y_end,'-s','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','k',...
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

%% Figure 2: Course angle vs time
figure(2)
fig = gcf;
% fig.WindowState = 'maximized';
ax2 = gca;
plot(t,x(:,3)*180/pi,'color',blue,'LineWidth',lw);hold on;
plot(t,kai_des*180/pi,'--','color',red,'LineWidth',lw);
ax2.FontSize = ax_fnt;
box on
ax2.XColor = 'black';
ax2.YColor = 'black';
set(ax2,'linewidth',ax_lw)
xlabel(ax2,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax2,'$$ \chi, $$ deg','Fontsize',lbl_fnt);
legend('Actual','Desired','Fontsize',leg_fnt);
grid on;

%% Figure 3: Course angle rate vs time
figure(3)
% fig = gcf;
fig.WindowState = 'maximized';
ax3 = gca;
plot(t,kaidot_actual*180/pi,'color',blue,'LineWidth',lw);hold on;
plot(t,kaidot_des*180/pi,'--','color',red,'LineWidth',lw);
ax3.FontSize = ax_fnt;
box on
ax3.XColor = 'black';
ax3.YColor = 'black';
set(ax3,'linewidth',ax_lw)
xlabel(ax3,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax3,'$$ \dot{\chi}, $$ deg/s','Fontsize',lbl_fnt);
legend('Actual','Desired','Fontsize',leg_fnt);
grid on;

%% Figure 4: Curvature vs time
figure(4)
fig = gcf;
% fig.WindowState = 'maximized';
ax4 = gca;
plot(t,kappa_actual,'color',blue,'LineWidth',lw);hold on;
plot(t,kappa_des,'--','color',red,'LineWidth',lw);
ax4.FontSize = ax_fnt;
box on
ax4.XColor = 'black';
ax4.YColor = 'black';
set(ax4,'linewidth',ax_lw)
xlabel(ax4,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax4,'$$ \kappa, $$ m$$^{-1}$$','Fontsize',lbl_fnt);
legend('Actual','Desired','Fontsize',leg_fnt);
grid on;

%% Figure 5: x position vs time
figure(5)
fig = gcf;
% fig.WindowState = 'maximized';
ax5 = gca;
plot(t,x(:,1),'color',blue,'LineWidth',lw);
ax5.FontSize = ax_fnt;
box on
ax5.XColor = 'black';
ax5.YColor = 'black';
set(ax5,'linewidth',ax_lw)
xlabel(ax5,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax5,'$$ x, $$ m','Fontsize',lbl_fnt);
grid on;

%% Figure 6: y position vs time
figure(6)
fig = gcf;
fig.WindowState = 'maximized';
ax6 = gca;
plot(t,x(:,2),'color',blue,'LineWidth',lw);
ax6.FontSize = ax_fnt;
box on
ax6.XColor = 'black';
ax6.YColor = 'black';
set(ax6,'linewidth',ax_lw)
xlabel(ax6,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax6,'$$ y, $$ m','Fontsize',lbl_fnt);
grid on;

%% Figure 7: Animation with Fixed-Wing UAV Marker
figure(7)
fig = gcf;
fig.WindowState = 'maximized';
ax7 = gca;
xline(0,'-k','linewidth',3);hold on;grid on;
if x0 < 0
    quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
else
    quiver(X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
end
plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold on;
ax7.FontSize = ax_fnt;
box on
ax7.XColor = 'black';
ax7.YColor = 'black';
set(ax7,'linewidth',ax_lw)
xlabel(ax7,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax7,'$$ y, $$  m','Fontsize',lbl_fnt);
axis(ax7,'equal')
% xlim([min(x(:,1))-10, max(x(:,1))+10]);
% ylim([min(x(:,2))-10, max(x(:,2))+10]);
axis(ax7,'equal')

% Initialize trajectory line and UAV marker
h_trail = plot(nan, nan, 'color', maroon, 'linewidth', lw);

% Create fixed-wing UAV shape (triangle pointing in direction of motion)
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav = fill(nan, nan, green, 'EdgeColor', black, 'LineWidth', 2);

legend('Desired line','Vector field','','Trajectory','','Fontsize',leg_fnt);

% Setup video writer
% Fix figure size in pixels so video frames are consistent across monitors
set(gcf, 'Units', 'pixels', 'Position', [100 100 1074 648]);  % fixed size
drawnow;  % force resize before first getframe
video_filename = sprintf('straightline_animation_%s.mp4', datestr(now,''));
vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 10; % 10 fps for smooth animation
vidObj.Quality = 95; % High quality
open(vidObj);
% Capture first frame to lock the expected frame size
initFrame = getframe(gcf);
expectedSize = size(initFrame.cdata);
writeVideo(vidObj, initFrame);

% Animation loop
dt_anim = 0.1; % Animation time step
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
    
    title(ax7, sprintf('Time: %.2f s', t_anim(i)), 'Fontsize', lbl_fnt);
    drawnow;
    
    % Capture frame for video and resize to match expected size
    % (prevents frame-size mismatch across different monitors/scaling)
    frame = getframe(gcf);
    if ~isequal(size(frame.cdata), expectedSize)
        frame.cdata = imresize(frame.cdata, [expectedSize(1) expectedSize(2)]);
    end
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
trajectory_data.course_angle_desired = kai_des;
trajectory_data.course_angle_rate_desired = kaidot_des;
trajectory_data.curvature_actual = kappa_actual;
trajectory_data.curvature_desired = kappa_des;
trajectory_data.initial_conditions = x_initial;
trajectory_data.parameters.Vgd = Vgd;
trajectory_data.parameters.k2 = k2;
trajectory_data.parameters.k_kai = k_kai;
trajectory_data.parameters.k_kaidot = k_kaidot;
trajectory_data.vector_field.X = X;
trajectory_data.vector_field.Y = Y;
trajectory_data.vector_field.xdot = xdot;
trajectory_data.vector_field.ydot = ydot;

% Save with timestamp
filename = sprintf('straightline_trajectory_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
save(filename, 'trajectory_data');
fprintf('Trajectory data saved to: %s\n', filename);

%% Create GIF animation: UAV moving on vector field

gif_name = 'straight_line_path_following.gif';

figure(10); clf;
set(gcf,'Color','w','Position',[100 100 900 650]);
ax_gif = gca;
hold(ax_gif,'on');
grid(ax_gif,'on');
box(ax_gif,'on');

% Plot fixed background once
xline(ax_gif,0,'-k','linewidth',3); hold(ax_gif,'on');
quiver(ax_gif,X,Y,xdot,ydot,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5); hold(ax_gif,'on');
plot(ax_gif,x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',green,...
    'MarkerSize',10); hold(ax_gif,'on');

xlabel(ax_gif,' $ x, $  m','Fontsize',lbl_fnt);
ylabel(ax_gif,' $ y, $  m','Fontsize',lbl_fnt);
title(ax_gif,'Arcsine Vector Field Guidance','Fontsize',lbl_fnt);

ax_gif.FontSize = ax_fnt;
ax_gif.XColor = 'black';
ax_gif.YColor = 'black';
set(ax_gif,'linewidth',ax_lw)
axis(ax_gif,'equal');
xlim(ax_gif,[-120 120]);
ylim(ax_gif,[-120 120]);

% animated objects
traj_line = plot(ax_gif,nan,nan,'color',maroon,'linewidth',lw);

% Create fixed-wing UAV shape (triangle pointing in direction of motion)
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav_gif = fill(ax_gif,nan, nan, green, 'EdgeColor', black, 'LineWidth', 2);

legend(ax_gif,'Desired line','Vector field','','UAV trajectory','', ...
    'Fontsize',leg_fnt,'Location','northwest');

% choose fewer frames for smaller GIF
skip = 8;

for k = 1:skip:length(t)
    
    set(traj_line,'XData',x(1:k,1),'YData',x(1:k,2));
    
    % Get current heading
    chi_curr = x(k,3);
    
    % Rotate and translate UAV shape
    R = [cos(chi_curr), -sin(chi_curr); sin(chi_curr), cos(chi_curr)];
    uav_rotated = R * [uav_shape_x; uav_shape_y];
    uav_x = uav_rotated(1,:) + x(k,1);
    uav_y = uav_rotated(2,:) + x(k,2);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fun_stline_upd(t,x)
global Vgd k2 wx wy  kv k_kai k_kaidot

if x(1) < 0
   kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
   kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
else
   kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
   kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
end
% 
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
% 
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
% 
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
% out(5,1) = kv*(Vad - x(5));
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



function [kaid_mod,kaid_dot] = vf_proposed_mod(x)
global Vgd  k2 
q = 1.5;
for i = 1:length(x)
    for j = 1:length(x)
if x(i,j) < 0
   % kaid(i,j) = asin(1./(1+k2*(x(i,j)).^2)) ;
   kaid_mod(i,j) = pi/2 - asin(500*k2.*abs(x(i,j)).^q./(1 + 500*k2.*abs(x(i,j)).^q));
   % kaid_dot(i,j) =  2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
else
   % kaid(i,j) = pi - asin(1./(1+k2*(x(i,j)).^2)) ;
   % kaid_dot(i,j) = - 2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
kaid_mod(i,j) = pi/2 + asin(500*k2.*abs(x(i,j)).^q./(1 + 500*k2.*abs(x(i,j)).^q));
end
    end
end
end