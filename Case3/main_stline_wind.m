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
f7 = figure;
ax7 = axes;
f8 = figure;
ax8 = axes;
% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20; % Label font size
leg_fnt = 16; % Legend font size
ax_lw = 3;        % Axis line width
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

f1;
xl = xline(ax1,0,'k','linewidth',3);hold(ax1,'on');
h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
h1 = plot(ax1,x(:,1),x(:,2),'r');hold(ax1,'on');
% h1.Color = red;                 % Color of plot
h1.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h1.LineWidth = lw_;             % Plot line width
h1.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h1.MarkerSize = ms_;            % Marker Size
h1.MarkerEdgeColor = 'auto';    % Marker edge color
h1.MarkerFaceColor = 'none';    % Marker face color
plot(ax1,x_ini,y_ini,'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold(ax1,'on');
plot(ax1,x_end,y_end,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold(ax1,'on');
% Setup axis font size
%ax1 = gca;
ax1.FontSize = ax_fnt;   
l = legend(ax1,[xl h h1] ,'Desired path','Vector field','UAV trajectory');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Outer box setup
box(ax1,"on")                      % Switch on the box around the axis
ax1.XColor = black;         % Box horizontal lines' color
ax1.YColor = black;         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax1,'$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
axis(ax1,'equal')

f2;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
h2 = plot(ax2,t,x(:,1),'r');hold(ax2,'on') ;grid(ax2,'on');
% h3.Color = red;                % Color of plot
h2.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h2.LineWidth = lw_;             % Plot line width
h2.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h2.MarkerSize = ms_;            % Marker Size
h2.MarkerEdgeColor = 'auto';    % Marker edge color
h2.MarkerFaceColor = 'none';    % Marker face color
% 
% l = legend(ax2,h3,'Thakkar et al.','Park et al.');
% l.Orientation = 'Vertical';
% l.Location = 'best';
% l.FontSize = leg_fnt;
% l.Interpreter = 'Latex';
% Setup axis font size   
ax2.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax2.XColor = black;         % Box horizontal lines' color
ax2.YColor = black;         % Box vertical lines' color
set(ax2,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax2,'Time $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax2,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)


f3;
h3 = plot(ax3,x(:,1),kaidot_des,'r');hold(ax3,'on');grid(ax3,'on');
% h3.Color = red;                % Color of plot
h3.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h3.LineWidth = lw_;             % Plot line width
h3.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h3.MarkerSize = ms_;            % Marker Size
h3.MarkerEdgeColor = 'auto';    % Marker edge color
h3.MarkerFaceColor = 'none';    % Marker face color

h31 = plot(ax3,x(:,1),x(:,4),'r');hold(ax3,'on');grid(ax3,'on');
% h31.Color = blue;                % Color of plot
h31.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h31.LineWidth = lw_;             % Plot line width
h31.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h31.MarkerSize = ms_;            % Marker Size
h31.MarkerEdgeColor = 'auto';    % Marker edge color
h31.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax3,[h3 h31],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
%Setup axis font size   
ax3.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax3.XColor = black;         % Box horizontal lines' color
ax3.YColor = black;         % Box vertical lines' color
set(ax3,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax3,'Cross-track error  $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax3,'Course rate $$\dot{\chi}$$, rad./s','Interpreter','Latex','FontSize',lbl_fnt)
% 
f4;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
h4 = plot(ax4,x(:,1),kappa_des,'r');hold(ax4,'on') ;grid(ax4,'on')
% h4.Color = red;                % Color of plot
h4.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h4.LineWidth = lw_;             % Plot line width
h4.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h4.MarkerSize = ms_;            % Marker Size
h4.MarkerEdgeColor = 'auto';    % Marker edge color
h4.MarkerFaceColor = 'none';    % Marker face color

h41 = plot(ax4,x(:,1),kappa_actual,'r');hold(ax4,'on') ;grid(ax4,'on')
% h41.Color = blue;                % Color of plot
h41.LineStyle = '--';            % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h41.LineWidth = lw_;             % Plot line width
h41.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h41.MarkerSize = ms_;            % Marker Size
h41.MarkerEdgeColor = 'auto';    % Marker edge color
h41.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax4,[h4 h41],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size
% ax4 = gca;     
ax4.FontSize = ax_fnt;
% axis(ax4,[0 max(t(end),t1(end)) -5 5]);
% Outer box setup
box on                      % Switch on the box around the axis
ax4.XColor = black;         % Box horizontal lines' color
ax4.YColor = black;         % Box vertical lines' color
set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax4,'Cross-track error  $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax4,'Curvature $$\kappa$$,  m $ ^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)

f5;
h5 = plot(ax5,t,kai_des*(180/pi),'r');hold(ax5,'on') ;grid(ax5,'on')
% h5.Color = 'red';                % Color of plot
h5.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h5.LineWidth = lw_;             % Plot line width
h5.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h5.MarkerSize = ms_;            % Marker Size
h5.MarkerEdgeColor = 'auto';    % Marker edge color
h5.MarkerFaceColor = 'none';    % Marker face color
h51 = plot(ax5,t,x(:,3)*(180/pi),'b');hold(ax5,'on') ;grid(ax5,'on');
% h51.Color = 'blue';                % Color of plot
h51.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h51.LineWidth = lw_;             % Plot line width
h51.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h51.MarkerSize = ms_;            % Marker Size
h51.MarkerEdgeColor = 'auto';    % Marker edge color
h51.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax5,[h5 h51],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax5.FontSize = ax_fnt;
yticks(ax5,[90 105 120 135 150 165 180])                           % Y-axis grid line positioning
yticklabels(ax5,{'90', '105', '120', '135', '150', '165','180'})
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax5.XColor = black;         % Box horizontal lines' color
ax5.YColor = black;         % Box vertical lines' color
set(ax5,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax5,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax5,'$$ \chi$$,  deg.','Fontsize',lbl_fnt);

f6;
h6 = plot(ax6,x(:,1),kai_des*(180/pi),'r');hold(ax6,'on') ;grid(ax6,'on')
% h6.Color = red;                % Color of plot
h6.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h6.LineWidth = lw_;             % Plot line width
h6.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h6.MarkerSize = ms_;            % Marker Size
h6.MarkerEdgeColor = 'auto';    % Marker edge color
h6.MarkerFaceColor = 'none';    % Marker face color
h61 = plot(ax6,x(:,1),x(:,3)*(180/pi),'r');hold(ax6,'on') ;grid(ax6,'on');
% h61.Color = blue;                % Color of plot
h61.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h61.LineWidth = lw_;             % Plot line width
h61.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h61.MarkerSize = ms_;            % Marker Size
h61.MarkerEdgeColor = 'auto';    % Marker edge color
h61.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax6,[h6 h61],'Desired','Actual');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax6.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax6.XColor = black;         % Box horizontal lines' color
ax6.YColor = black;         % Box vertical lines' color
set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax6,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax6,' $$ \chi$$,  deg.','Fontsize',lbl_fnt);

f7;
h7 = plot(ax7,t,Vad,'r');hold(ax7,'on') ;grid(ax7,'on')
% h7.Color = red;                % Color of plot
h7.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h7.LineWidth = lw_;             % Plot line width
% h7(4).Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h7(5).MarkerSize = ms_;            % Marker Size
% h7(6).MarkerEdgeColor = 'auto';    % Marker edge color
% h7(7).MarkerFaceColor = 'none';    % Marker face color
h71 = plot(ax7,t,x(:,5),'r');hold(ax7,'on') ;grid(ax7,'on');
% h71.Color = blue;                % Color of plot
h71.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h71.LineWidth = lw_;             % Plot line width
% h71(4).Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h71(5).MarkerSize = ms_;            % Marker Size
% h71(6).MarkerEdgeColor = 'auto';    % Marker edge color
% h71(7).MarkerFaceColor = 'none';    % Marker face color
l = legend(ax7,[h7 h71],'Desired','Actual');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax7.FontSize = ax_fnt;
axis(ax7,[0 t(end) 0 15]);
% % Outer box setup
 box on                      % Switch on the box around the axis
 ax7.XColor = black;         % Box horizontal lines' color
 ax7.YColor = black;         % Box vertical lines' color
set(ax7,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
 xlabel(ax7,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
 ylabel(ax7,' $$ V_{\mathrm{a}}$$,  m/s','Fontsize',lbl_fnt);
% 
f8;
h8 = plot(ax8,t,Vgd(i)*ones(length(t),1),'r','linewidth',3);hold(ax8,'on') ;grid(ax8,'on')
% h8.Color = 'red';                 % Color of plot
h8.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h8.LineWidth = lw_;             % Plot line width
% h8.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h8.MarkerSize = ms_;            % Marker Size
% h8.MarkerEdgeColor = 'auto';    % Marker edge color
% h8.MarkerFaceColor = 'none';    % Marker face color
h81 = plot(ax8,t,Vg,'linewidth',3);hold(ax8,'on') ;grid(ax8,'on');
h81.Color = 'b';                % Color of plot
h81.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h81.LineWidth = lw_;             % Plot line width
% h81.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h81.MarkerSize = ms_;            % Marker Size
% h81.MarkerEdgeColor = 'auto';    % Marker edge color
% h81.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax8,[h8 h81],'Commanded','Achieved');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax8.FontSize = ax_fnt;
axis(ax8,[0 t(end) 0 15]);
% % Outer box setup
box on                     % Switch on the box around the axis
ax8.XColor = black;         % Box horizontal lines' color
ax8.YColor = black;         % Box vertical lines' color
set(ax8,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
xlabel(ax8,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax8,' $$ V_{\mathrm{g}}$$,  m/s','Fontsize',lbl_fnt);
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
