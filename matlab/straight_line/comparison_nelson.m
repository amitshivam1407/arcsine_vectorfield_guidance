%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
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
% f10 = figure;
% ax10 = axes;
% f11 = figure;
% ax11 = axes;
%--------------- declaration of constant parameters ----------------------%
global  Vgd k_kai k_prop k_Nelson k_kaidot flag1 t_flag1 flag2 t_flag2 figure_type_good
flag1 = true;
flag2 = true;
Vgd = 10;
k_kai = 200;  
 k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
 % k_kaidot = 10;  % 20 for better curvature profile  %%
figure_type_good = 1;
%-------------------------------------------------------------------------%

%------------------- initial conditions   --------------------------------%

tspan = [0 20];

x0 = 100;
y0 = -90;

k_prop = 0.005;
% flag1 = true(1);
% t_flag1 = zeros(1,2);
% flag2 = true(1);
% t_flag2 = zeros(1,2);

% for i = 1:length(k_prop)

if x0 < 0
   kaid0_prop =  asin(1./(1+ k_prop*(x0).^2)) ;
   kaid_dot0_prop = (2*k_prop*Vgd.*(x0))./((1+k_prop*(x0).^2).^2);
else
   kaid0_prop = pi - asin(1./(1 + k_prop*(x0).^2)) ;  
   kaid_dot0_prop = - (2*k_prop*Vgd.*(x0))./((1 + k_prop*(x0).^2).^2);
end
kai0_prop = kaid0_prop;
kaidot0_prop = kaid_dot0_prop;

kai0_deg = rad2deg(wrapToPi(kai0_prop));


k_Nelson = sqrt(k_prop^2 * x0^2 + 2*k_prop);

kai0_Nelson = pi/2 + atan(k_Nelson.*x0);
kaidot0_Nelson = -(k_Nelson^2*Vgd*x0)./(1+k_Nelson^2.*x0.^2).^(3/2);



%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
% range  = -110:10:110;
% [X,Y] = meshgrid(range);
% [kaid,kaid_dot] = vf_proposed(X);
% xdot  = Vgd*cos(kaid);
% ydot  = Vgd*sin(kaid);


%----------------  ode solver   ------------------------------------------%
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

x_initial = [x0;y0;kai0_prop;kaidot0_prop;x0;y0;kai0_Nelson;kaidot0_Nelson];
[t,x] = ode45(@(t,x)fun_stline_upd(t,x,k_prop,k_Nelson) ,tspan, x_initial,options);

k_prop_plot = 0.0002:(0.005-0.0002)/(length(t)-1):0.005;
k_prop_plot = k_prop_plot';
k_prop_indices = 0.0002:(0.005-0.0002)/9:0.005;
for i = 1 :length(t)
    if x0 < 0
        kaid0_prop_plot(i) =  asin(1./(1+ k_prop_plot(i)*(x0).^2)) ;
    else
        kaid0_prop_plot(i) = pi - asin(1./(1 + k_prop_plot(i)*(x0).^2)) ;
    end
    kaid0_prop_plot(i) = rad2deg(wrapToPi(kaid0_prop_plot(i)));
    k_Nelson_plot(i) = sqrt(k_prop_plot(i)^2 * x0^2 + 2*k_prop_plot(i));    
end
kaid0_prop_plot = kaid0_prop_plot' ;
k_Nelson_plot =  k_Nelson_plot';

%----------------------  parameters to plot   ----------------------------%

x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;


for j = 1:length(t)
  kai_des_Nelson(j,:) = pi/2 + atan(k_Nelson.*x(j,5));
if x(j,1) < 0
   kai_des_prop(j,:) =  asin(1./(1+ k_prop*(x(j,1)).^2)) ;
   kaidot_des(j,:) = -(2*k_prop*Vgd.*x(j,1))./((1+k_prop*(x(j,1)).^2).^2);
else
   kai_des_prop(j,:) = pi - asin(1./(1 + k_prop*(x(j,1)).^2)) ;  
   kaidot_des(j,:) = - (2*k_prop*Vgd.*x(j,1))./((1 + k_prop*(x(j,1)).^2).^2);
end
end

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

kaidot_actual_prop = x(:,4);
% kappa_des = kaidot_des./Vgd;
kappa_actual_prop = x(:,4)./Vgd ;
control_effort_prop = trapz(t,(Vgd*x(:,4)).^2);
[val_max_prop, ind_prop] = min(kappa_actual_prop);
x_max_prop = x(ind_prop,1);
 

kaidot_actual_Nelson = x(:,8);
control_effort_Nelson = trapz(t,(Vgd*x(:,8)).^2);
% kappa_des = kaidot_des./Vgd;
kappa_actual_Nelson = x(:,8)./Vgd ;
[val_max_Nelson, ind_Nelson] = min(kappa_actual_Nelson);
x_max_Nelson = x(ind_Nelson,5);

for i = 1:length(t)
 [val_prop(i,:),val_Nelson(i,:)] = fun_kappa(x(i,1),x(i,5));
 ratio_curvature(i,:) = abs(val_prop(i,:))./abs(val_Nelson(i,:)) ;
end
% ratio_curvature = abs(val_max_prop)./abs(val_max_Nelson);



%---------------------- plotting figures ---------------------------------%

%% Figure 1: Trajectory Comparison with Vector Field
figure(1)
fig = gcf;
fig.WindowState = 'maximized';
xp = xline(0,'-k','linewidth',3);hold on;grid on;
% Vector field for proposed method
range  = -110:10:110;
[X,Y] = meshgrid(range);
[kaid_prop_vf,kaid_dot_prop_vf] = vf_proposed(X);
xdot_prop_vf  = Vgd*cos(kaid_prop_vf);
ydot_prop_vf  = Vgd*sin(kaid_prop_vf);
quiver(X,Y,xdot_prop_vf,ydot_prop_vf,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
% Plot trajectories
h_traj_prop = plot(x(:,1),x(:,2),'color',maroon,'linewidth',lw);hold on;
h_traj_nelson = plot(x(:,5),x(:,6),'color',blue,'linewidth',lw);hold on;
% Initial positions
h_ini_prop = plot(x(1,1),x(1,2),'^k','MarkerSize',12,'MarkerFaceColor',cyan);hold on;
h_ini_nelson = plot(x(1,5),x(1,6),'^k','MarkerSize',12,'MarkerFaceColor', green);hold on;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on
ax1.XColor = 'black';
ax1.YColor = 'black';
set(ax1,'linewidth',ax_lw)
xlabel(ax1,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax1,'$$ y, $$  m','Fontsize',lbl_fnt);
legend(ax1,'Desired line $$x = 0$$','Vector field','Proposed','Nelson et al.','Initial (Proposed)','Initial (Nelson)','Fontsize',leg_fnt);
axis(ax1,'equal')

%% Figure 2: Cross-track error comparison
figure(2)
fig = gcf;
fig.WindowState = 'maximized';
ax2 = gca;
plot(t,x(:,1),'color',maroon,'LineWidth',lw);hold on;
plot(t,x(:,5),'color',blue,'LineWidth',lw);
ax2.FontSize = ax_fnt;
box on
ax2.XColor = 'black';
ax2.YColor = 'black';
set(ax2,'linewidth',ax_lw)
xlabel(ax2,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax2,'$$ x, $$ m','Fontsize',lbl_fnt);
legend('Proposed','Nelson et al.','Fontsize',leg_fnt);
grid on;

%% Figure 3: Course angle comparison
figure(3)
fig = gcf;
fig.WindowState = 'maximized';
ax3 = gca;
plot(t,x(:,3)*180/pi,'color',maroon,'LineWidth',lw);hold on;
plot(t,x(:,7)*180/pi,'color',blue,'LineWidth',lw);
ax3.FontSize = ax_fnt;
box on
ax3.XColor = 'black';
ax3.YColor = 'black';
set(ax3,'linewidth',ax_lw)
xlabel(ax3,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax3,'$$ \chi, $$ deg','Fontsize',lbl_fnt);
legend('Proposed','Nelson et al.','Fontsize',leg_fnt);
grid on;

%% Figure 4: Curvature comparison
figure(4)
fig = gcf;
fig.WindowState = 'maximized';
ax4 = gca;
plot(t,kappa_actual_prop,'color',maroon,'LineWidth',lw);hold on;
plot(t,kappa_actual_Nelson,'color',blue,'LineWidth',lw);
ax4.FontSize = ax_fnt;
box on
ax4.XColor = 'black';
ax4.YColor = 'black';
set(ax4,'linewidth',ax_lw)
xlabel(ax4,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax4,'$$ \kappa, $$ m$$^{-1}$$','Fontsize',lbl_fnt);
legend('Proposed','Nelson et al.','Fontsize',leg_fnt);
grid on;

%% Figure 4: Curvature crosstrack error comparison
figure(6)
fig = gcf;
fig.WindowState = 'maximized';
ax6 = gca;
plot(x(:,1),kappa_actual_prop,'color',maroon,'LineWidth',lw);hold on;
plot(x(:,5),kappa_actual_Nelson,'color',blue,'LineWidth',lw);
ax6.FontSize = ax_fnt;
box on
ax6.XColor = 'black';
ax6.YColor = 'black';
set(ax6,'linewidth',ax_lw)
xlabel(ax6,' $$ x, $$ m','Fontsize',lbl_fnt);
ylabel(ax6,'$$ \kappa, $$ m$$^{-1}$$','Fontsize',lbl_fnt);
legend('Proposed','Nelson et al.','Fontsize',leg_fnt);
grid on;
%% Figure 5: Animation with Fixed-Wing UAV Markers
figure(5)
fig = gcf;
fig.WindowState = 'maximized';
ax5 = gca;
xline(0,'-k','linewidth',3);hold on;grid on;
quiver(X,Y,xdot_prop_vf,ydot_prop_vf,0.7,'color',[0.2 0.6 0.8],'linewidth',1.5);hold on;
plot(x(1,1),x(1,2),'^k','MarkerSize',12,'MarkerFaceColor',cyan);hold on;
plot(x(1,5),x(1,6),'^k','MarkerSize',12,'MarkerFaceColor', green);hold on;
ax5.FontSize = ax_fnt;
box on
ax5.XColor = 'black';
ax5.YColor = 'black';
set(ax5,'linewidth',ax_lw)
xlabel(ax5,' $$ x, $$  m','Fontsize',lbl_fnt);
ylabel(ax5,'$$ y, $$  m','Fontsize',lbl_fnt);
axis(ax5,'equal')

% Initialize trajectory trails and UAV markers
h_trail_prop = plot(nan, nan, 'color', maroon, 'linewidth', lw);
h_trail_nelson = plot(nan, nan, 'color', green, 'linewidth', lw);

% Create fixed-wing UAV shapes (triangles pointing in direction of motion)
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav_prop = fill(nan, nan, maroon, 'EdgeColor', black, 'LineWidth', 2);
h_uav_nelson = fill(nan, nan, green, 'EdgeColor', black, 'LineWidth', 2);

legend('Desired line $$x = 0$$','Vector field','Proposed','Nelson et al.','Fontsize',leg_fnt);

% Setup video writer
video_filename = sprintf('comparison_nelson_stline.mp4');
vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 10;
vidObj.Quality = 95;
open(vidObj);

% Animation loop
dt_anim = 0.1;
t_anim = 0:dt_anim:t(end);
for i = 1:length(t_anim)
    % Interpolate positions and headings for both methods
    x_curr_prop = interp1(t, x(:,1), t_anim(i));
    y_curr_prop = interp1(t, x(:,2), t_anim(i));
    chi_curr_prop = interp1(t, x(:,3), t_anim(i));
    
    x_curr_nelson = interp1(t, x(:,5), t_anim(i));
    y_curr_nelson = interp1(t, x(:,6), t_anim(i));
    chi_curr_nelson = interp1(t, x(:,7), t_anim(i));
    
    % Update trajectory trails
    idx = find(t <= t_anim(i));
    set(h_trail_prop, 'XData', x(idx,1), 'YData', x(idx,2));
    set(h_trail_nelson, 'XData', x(idx,5), 'YData', x(idx,6));
    
    % Rotate and translate UAV shapes for proposed
    R_prop = [cos(chi_curr_prop), -sin(chi_curr_prop); sin(chi_curr_prop), cos(chi_curr_prop)];
    uav_rotated_prop = R_prop * [uav_shape_x; uav_shape_y];
    uav_x_prop = uav_rotated_prop(1,:) + x_curr_prop;
    uav_y_prop = uav_rotated_prop(2,:) + y_curr_prop;
    set(h_uav_prop, 'XData', uav_x_prop, 'YData', uav_y_prop);
    
    % Rotate and translate UAV shapes for Nelson
    R_nelson = [cos(chi_curr_nelson), -sin(chi_curr_nelson); sin(chi_curr_nelson), cos(chi_curr_nelson)];
    uav_rotated_nelson = R_nelson * [uav_shape_x; uav_shape_y];
    uav_x_nelson = uav_rotated_nelson(1,:) + x_curr_nelson;
    uav_y_nelson = uav_rotated_nelson(2,:) + y_curr_nelson;
    set(h_uav_nelson, 'XData', uav_x_nelson, 'YData', uav_y_nelson);
    
    title(ax5, sprintf('Time: %.2f s', t_anim(i)), 'Fontsize', lbl_fnt);
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
trajectory_data.proposed.position_x = x(:,1);
trajectory_data.proposed.position_y = x(:,2);
trajectory_data.proposed.course_angle = x(:,3);
trajectory_data.proposed.course_angle_rate = x(:,4);
trajectory_data.proposed.curvature = kappa_actual_prop;
trajectory_data.proposed.course_angle_desired = kai_des_prop;
trajectory_data.nelson.position_x = x(:,5);
trajectory_data.nelson.position_y = x(:,6);
trajectory_data.nelson.course_angle = x(:,7);
trajectory_data.nelson.course_angle_rate = x(:,8);
trajectory_data.nelson.curvature = kappa_actual_Nelson;
trajectory_data.nelson.course_angle_desired = kai_des_Nelson;
trajectory_data.parameters.Vgd = Vgd;
trajectory_data.parameters.k_prop = k_prop;
trajectory_data.parameters.k_Nelson = k_Nelson;
trajectory_data.parameters.k_kai = k_kai;
trajectory_data.parameters.k_kaidot = k_kaidot;
trajectory_data.control_effort.proposed = control_effort_prop;
trajectory_data.control_effort.nelson = control_effort_Nelson;
trajectory_data.max_curvature.proposed = val_max_prop;
trajectory_data.max_curvature.nelson = val_max_Nelson;

% Save with timestamp
filename = sprintf('comparison_nelson_trajectory.mat');
save(filename, 'trajectory_data');
fprintf('Trajectory data saved to: %s\n', filename);



function out = fun_stline_upd(t,x,k_prop,k_Nelson)
global Vgd k_kai k_kaidot flag1 t_flag1 flag2 t_flag2 

if (flag1 == true) && (x(1)<0.2)
    t_flag1 = t;
    flag1 = false;
end

if (flag2 == true) && (x(5)<0.2)
    t_flag2 = t;
    flag2 = false;
end

% 
% if (flag2 == true) && (x(5)<0.2)
%     t_flag2 = t;
%     flag2 = false;
% end

if x(1) < 0
   kaid_prop =  asin(1./(1+ k_prop*(x(1)).^2)) ;
   kaid_dot_prop = 2*k_prop*Vgd.*(x(1))./((1+k_prop*(x(1)).^2).^2);
else
   kaid_prop = pi - asin(1./(1 + k_prop*(x(1)).^2)) ;  
   kaid_dot_prop = - 2*k_prop*Vgd.*(x(1))./((1 + k_prop*(x(1)).^2).^2);
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

function [val_prop,val_Nelson]  = fun_kappa(x_prop,x_Nelson)
global Vgd k_prop k_Nelson
% for i = 1:length(x_prop)
    if x_prop < 0
        kaidot_des = -(2*k_prop*Vgd.*x_prop)./((1+k_prop*(x_prop).^2).^2);
        kappa_des_prop = kaidot_des./Vgd ;
    else      
        kaidot_des = - (2*k_prop*Vgd.*x_prop)./((1 + k_prop*(x_prop).^2).^2);
        kappa_des_prop = kaidot_des./Vgd ;
    end
[val_prop, ind_prop] = max(abs(kappa_des_prop));
% x_max_prop = x(ind_prop,1);
kaid_dot_Nelson = -(k_Nelson^2*Vgd*x_Nelson)./(1+k_Nelson^2.*x_Nelson.^2).^(3/2);
kappa_des_Nelson = kaid_dot_Nelson./Vgd ;
[val_Nelson, ind_Nelson] = max(abs(kappa_des_Nelson));
% x_max_Nelson = x(ind_Nelson,1);
% ratio_curvature(i,1) = abs(val_max_prop(i,1))./abs(val_max_Nelson(i,1)) ;
%end
end

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

