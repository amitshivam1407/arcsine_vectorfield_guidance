%-------------------------------------------------------------------------%
%------------            11th September 2022            ------------------%
%------------         Autopilot design with no wind     ------------------%
%------------         Straight line path following      ------------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%--------------- declaration of constant parameters ----------------------%
global  Vgd k_kai k_kaidot figure_type_good kaidot_max te kaidot_ode kaidot
Vgd = 10;
k_kai = 0.98;  
% k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
k_kaidot = 1.5;  % 20 for better curvature profile  %%
figure_type_good = 1;
kaidot_max = 0.98;
%% figures
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
%% plot format control variables
lw = 3;            % Line width
ms = 6;            % Marker size
ax_fnt = 17;        % Axis font size
lbl_fnt = 20; % Label font size
leg_fnt = 16; % Legend font size
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
color = [red;blue;green;orange];
%-------------------------------------------------------------------------%

%-------------------- initial conditions  --------------------------------%

tspan = [0 20];

x0 = [5 10 20];
y0 = -100;
% k2 = [0.05 0.02 0.005 0.002];
k2 = 0.002;
for i = 1:length(x0)

if x0(i) < 0
   kaid0(i) =  asin(1./(1+ k2*(x0(i)).^2)) ;
   kaidot0(i) = -(2*k2*Vgd.*x0(i))./((1+k2*(x0(i)).^2).^2);
else
   kaid0(i) = pi - asin(1./(1 + k2*(x0(i)).^2)) ;  
   kaidot0(i) = - (2*k2*Vgd.*x0(i))./((1 + k2*(x0(i)).^2).^2);
end
kaid0_deg(i) = rad2deg(wrapToPi(kaid0(i)));
% x_initial = [x0;y0;kai0(i) ;kaidot0(i)];
kai0 = deg2rad(270);
kaidot0 = 0;
x_initial = [x0(i);y0;kai0 ;kaidot0];
%-------------------------------------------------------------------------%

%----------------  vector field construction  ----------------------------%
% range  = -100:10:100;
% [X,Y] = meshgrid(range);
% [kaid,kaid_dot] = vf_proposed(X);
% xdot  = Vgd*cos(kaid);
% ydot  = Vgd*sin(kaid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@(t,x)fun_stline_upd(t,x,k2) ,tspan, x_initial,options);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ini(i) = x(1,1) ;
y_ini(i) = x(1,2) ;
x_end(i) = x(end,1) ;
y_end(i) = x(end,2) ;

kaidot_actual = x(:,4);

kappa_actual = x(:,4)./Vgd ;

[val_prop(i) , ind_prop(i)] = min(kappa_actual);

x_prop(i) = x(ind_prop(i),1);
y_prop(i) = x(ind_prop(i),2);

[max_kaidot(i), ind_max(i)] = min(x(:,4));


%---------------------- plotting figures ---------------------------------%
f1;
% quiver(X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold on;
% path = xline(ax1,0,'k','linewidth',2);hold(ax1,'on');
if i == 1
h(i) = xline(ax1,0,'k','linewidth',3,'DisplayName','Desired path');hold(ax1,'on');grid(ax1,'on');
h(i).Annotation.LegendInformation.IconDisplayStyle = 'on';
h1(i) = plot(ax1,x(:,1),x(:,2),'linewidth',3,'DisplayName',['$$x_{0} = \ $$', num2str(x0(i))]);hold(ax1,'on'); 
h11(i) = plot(ax1,x_ini(i),y_ini(i),'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold on;
h11(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
else
  h1(i) = plot(ax1,x(:,1),x(:,2),'linewidth',3,'DisplayName',['$$x_{0} = \ $$', num2str(x0(i))]);hold(ax1,'on'); 
%   h(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
  h11(i) = plot(ax1,x_ini(i),y_ini(i),'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold on;
  h11(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
end
h1(i).Color = color(i,:);

% plot(ax1,x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold on;
ax1.FontSize = ax_fnt;
box(ax1,'on')
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
xlabel(ax1,'$ x, $  m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $  m','Fontsize',lbl_fnt);
% legend('Vector field','Desired \ path','UAV \ trajectory','Initial \ point','Final \ point','Fontsize',14);
axis(ax1,'equal')


f2;
h2 = plot(ax2,t,(x(:,1)),'LineWidth',3);hold(ax2,'on');grid(ax2,'on');
h2.Color = color(i,:);
ax2.FontSize = ax_fnt;
box(ax2,'on')
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3)   % Axis linewidth (box and grid)
xlabel(ax2,'$t$, s','Fontsize',lbl_fnt);
ylabel(ax2,'Cross-track error $$x$$, m','Fontsize',lbl_fnt);

f3;
h3 = plot(ax3,t,kaidot_actual,'LineWidth',3);hold(ax3,'on');grid(ax3,'on');
h3.Color = color(i,:);
ax3.FontSize = ax_fnt;
box(ax3,'on')
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3)   % Axis linewidth (box and grid)
xlabel(ax3,' $$t$$, m','Fontsize',lbl_fnt);
ylabel(ax3,' $\dot{\chi}$,  rad/s ','Fontsize',lbl_fnt);
% ax3.XLim = [t(1,1) t(end,1)];
% ax3.YLim = [-2 0.1];

f4;
h4 = plot(ax4,x(:,1),kappa_actual,'LineWidth',3);hold(ax4,'on');grid(ax4,'on');
h4.Color = color(i,:);
% h4 = plot(ax4,[x_prop(i) x_prop(i)],[0 val_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
% h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h2 = plot(ax4,[x_prop(i) x_prop(i)],[0 val_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
h41 = plot(ax4,x_prop(i),val_prop(i),'-s','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','cyan'); hold(ax4,'on');grid(ax4,'on');
% plot(ax4,[x_prop(i) x_prop(i)],[0 val_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
h41.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h42 = plot(ax4,[0 x_prop(i)],[val_prop(i) val_prop(i)],'--','color',[0.75 0.75 0.75],'linewidth',2);hold(ax4,'on');grid(ax4,'on');
% h42.Annotation.LegendInformation.IconDisplayStyle = 'off';
% yticks([-val_prop(1) -(val_prop(1)+0.02) -(val_prop(1)+0.04) -(val_prop(1)+0.06) ...
%     -(val_prop(1)+0.06) -(val_prop(1)+0.08) -(val_prop(1)+0.10)])
% yticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%   datatip(x_prop(i),val_prop(i));hold(ax4,'on');grid(ax4,'on');
% dt_prop = datatip(h4,x_prop(i),val_prop(i));hold(ax4,'on');grid(ax4,'on');
% h4.DataTipTemplate.DataTipRows(1).Format = '%.2f';
% h4.DataTipTemplate.DataTipRows(2).Format = '%.3f';
% l.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax4.FontSize = ax_fnt;
box(ax4,'on')
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',3)   % Axis linewidth (box and grid)
xlabel(ax4,'Cross-track error $$x$$, m','Fontsize',lbl_fnt);
ylabel(ax4,'Curvature $$\kappa$$, m $ ^{-1} $ ','Fontsize',lbl_fnt);

f5;
h5 = plot(ax5,t,x(:,3)*(180/pi),'linewidth',3);hold(ax5,'on');grid(ax5,'on');
h5.Color = color(i,:);
ax5.FontSize = ax_fnt;
box(ax5,'on')
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',3)   % Axis linewidth (box and grid)
% yticks(ax5,90:30:180);
yticks(ax5,[-120 -90 -60 -30 0 30 60 90 120 150 180])                           % Y-axis grid line positioning
yticklabels(ax5,{'-120', '-90', '-60', '-30' '0', '30', '60', '90', '120', '150', '180'})
xlabel(ax5,' $ t, $  s','Fontsize',lbl_fnt);
ylabel(ax5,'$ \chi$,  deg.','Fontsize',lbl_fnt);

f6;
h6 = plot(ax6,x(:,1),x(:,3)*(180/pi),'linewidth',3);hold(ax6,'on');grid(ax6,'on');
h6.Color = color(i,:);
ax6.FontSize = ax_fnt;
box(ax6,'on')
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',3)   % Axis linewidth (box and grid)
% yticks(ax6,90:30:180);
yticks(ax6,[-120 -90 -60 -30 0 30 60 90 120 150 180])                           % Y-axis grid line positioning
yticklabels(ax6,{'-120', '-90', '-60', '-30' '0', '30', '60', '90', '120', '150', '180'})
xlabel(ax6,'Cross-track error $$x$$, m','Fontsize',lbl_fnt);
ylabel(ax6,' $ \chi$,  deg.','Fontsize',lbl_fnt);
grid on;

f7;
h7 = plot(ax7,t,kappa_actual,'LineWidth',3);hold(ax7,'on');grid(ax7,'on');
h7.Color = color(i,:);
ax7.FontSize = ax_fnt;
box(ax7,'on')
ax7.XColor = 'black';         % Box horizontal lines' color
ax7.YColor = 'black';         % Box vertical lines' color
set(ax7,'linewidth',3)   % Axis linewidth (box and grid)
xlabel(ax7,' $ t, $ s','Fontsize',lbl_fnt);
ylabel(ax7,'Curvature $$\kappa$$,  m $ ^{-1} $ ','Fontsize',lbl_fnt);
grid on;
% figure_plot = fun_multistlineplot(t,k2(i),x(:,1),x(:,2),...
%   kaidot_actual,kappa_actual,...
%      x(:,3));
end
% Legend = cell(length(k2),1);
% for i = 1:length(k2)
% %   str1 = sprintf('$$k_{s} = %.3f $$',k2(i)) ; 
% Legend{i} =  ['$$k_{\mathrm{s}} = \ $$', num2str(k2(i))];
% end
% legend(ax1,Legend,'Fontsize',15 )
legend(ax1,'FontSize',leg_fnt)

Legend = cell(length(x0),1);
for i = 1:length(x0)
% str2 = sprintf('$ k_{s} =  %.3f $ ',k2(i)) ; 
Legend{i} =  ['$$x_{0} = \ $$', num2str(x0(i))];
end
legend(ax2,Legend,'Fontsize',leg_fnt )
Legend = cell(length(x0),1);
for i = 1:length(x0)
% str3 = sprintf('$k_{s} = %.3f $',k2(i)) ; 
Legend{i} =  ['$$x_{0} = \ $$', num2str(x0(i))];
end
legend(ax3,Legend,'Fontsize',leg_fnt )
Legend = cell(length(x0),1);
for i = 1:length(x0)
% str4 = sprintf('$k_{s} = %.3f $',k2(i)) ; 
Legend{i} =  ['$$x_{0} = \ $$', num2str(x0(i))];
end
legend(ax4,Legend,'Fontsize',leg_fnt )

Legend = cell(length(x0),1);
for i = 1:length(x0)
% str5 = sprintf('$k_{s} = %.3f $',k2(i)) ; 
Legend{i} = ['$$x_{0} = \ $$', num2str(x0(i))];
end
legend(ax5,Legend,'Fontsize',leg_fnt )

Legend = cell(length(x0),1);
for i = 1:length(x0)
% str6 = sprintf('$k_{s} = %.3f $',k2(i)) ; 
Legend{i} =  ['$$x_{0} = \ $$', num2str(x0(i))];
end
legend(ax6,Legend,'Fontsize',leg_fnt );

Legend = cell(length(x0),1);
for i = 1:length(x0)
% str7 = sprintf('$k_{s} = %.3f $',k2(i)) ; 
Legend{i} =  ['$$x_{0} =  \ $$', num2str(x0(i))];
end
legend(ax7,Legend,'Fontsize',leg_fnt )

function out = fun_stline_upd(t,x,k2)
global Vgd k_kai k_kaidot kaidot_max te kaidot_ode kaidot

if t==0
    te = [];   
    kaidot_ode = [];
    kaidot = [];
end

if x(1) < 0
   kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
   kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
else
   kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;  
   kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
end



out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));

if abs(x(4))>kaidot_max
    p = sign(x(4))*kaidot_max;
else
    p = x(4);
end

out(3,1) = p ;

% out(3,1) = x(4) ;

% out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;

out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - p )  ;

te = [te t];
kaidot_ode = [kaidot_ode kaid_dot];
kaidot = [kaidot p];

end

% function figure_plot = fun_multistlineplot(t,k2,x1,x2,...
%   kaidot_a,kappa_a,...
%      x3)
% global figure_type_good p_x0 p_y0 p_width p_height
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
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
% % Plot format control variables
% lw_ = 2;            % Line width
% ms_ = 6;            % Marker size
% ax_fnt = 16;        % Axis font size
% lbl_fnt = ax_fnt+1; % Label font size
% leg_fnt = ax_fnt-1; % Legend font size
% ax_lw = 1.5;        % Axis line width
% % Colors
% % blue = [0 0.4470 0.7410];
% % red = [0.8500 0.3250 0.0980];
% blue = [0 0 1];
% red = [1 0 0];
% orange = [0.9290 0.6940 0.1250];
% violet = [0.4940 0.1840 0.5560];
% green = [0.4660 0.6740 0.1880];
% cyan = [0.3010 0.7450 0.9330];
% maroon = [0.6350 0.0780 0.1840];
% black = [0 0 0];
% 
% if(figure_type_good == 1)
% % figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
% 
% f1;
% % h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
% xl = xline(ax1,0,'k','linewidth',2);hold(ax1,'on'); 
% xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h1 = plot(ax1,x1,x2);hold(ax1,'on');
% % h1.Color = red;                 % Color of plot
% h1.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h1.LineWidth = lw_;             % Plot line width
% h1.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h1.MarkerSize = ms_;            % Marker Size
% h1.MarkerEdgeColor = 'auto';    % Marker edge color
% h1.MarkerFaceColor = 'none';    % Marker face color
% % plot(ax1,x_ini,y_ini,'-o','MarkerSize',10,...
% %     'MarkerEdgeColor','blue',...
% %     'MarkerFaceColor','green'); hold(ax1,'on');
% % plot(ax1,x_end,y_end,'-s','MarkerSize',10,...
% %     'MarkerEdgeColor','blue',...
% %     'MarkerFaceColor','cyan'); hold(ax1,'on');
% ax1.FontSize = 13;
% % l = legend(ax1,h1);
% % l.Orientation = 'Vertical';
% % l.Location = 'best';
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% % Setup axis font size
% % ax1 = gca;     
% ax1.FontSize = ax_fnt;
% axis(ax1,'equal')
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax1.XColor = black;         % Box horizontal lines' color
% ax1.YColor = black;         % Box vertical lines' color
% set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax1,'$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
% axis(ax1,'equal')
% 
% f2;
% % figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
% h3 = plot(ax2,t,x1);hold(ax2,'on') ;grid(ax2,'on');
% h3.Color = red;                % Color of plot
% h3.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h3.LineWidth = lw_;             % Plot line width
% h3.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h3.MarkerSize = ms_;            % Marker Size
% h3.MarkerEdgeColor = 'auto';    % Marker edge color
% h3.MarkerFaceColor = 'none';    % Marker face color
% % 
% % l = legend(ax2,h3,'Thakkar et al.','Park et al.');
% % l.Orientation = 'Vertical';
% % l.Location = 'best';
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% % Setup axis font size   
% ax2.FontSize = ax_fnt;
% % axis(ax,[0 200 -40 60]);
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax2.XColor = black;         % Box horizontal lines' color
% ax2.YColor = black;         % Box vertical lines' color
% set(ax2,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax2,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax2,'Cross-track \ error, m','Interpreter','Latex','FontSize',lbl_fnt)
% 
% 
% f3;
% % h3 = plot(ax3,x1,kaidot_d,'DisplayName',['n = ',num2str(n(i))]);hold(ax3,'on');grid(ax3,'on');
% % h3.Color = red;                % Color of plot
% % h3.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% % h3.LineWidth = lw_;             % Plot line width
% % h3.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% % h3.MarkerSize = ms_;            % Marker Size
% % h3.MarkerEdgeColor = 'auto';    % Marker edge color
% % h3.MarkerFaceColor = 'none';    % Marker face color
% 
% h31 = plot(ax3,x1,kaidot_a);hold(ax3,'on');grid(ax3,'on');
% h31.Color = blue;                % Color of plot
% h31.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h31.LineWidth = lw_;             % Plot line width
% h31.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h31.MarkerSize = ms_;            % Marker Size
% h31.MarkerEdgeColor = 'auto';    % Marker edge color
% h31.MarkerFaceColor = 'none';    % Marker face color
% % l = legend(ax3,h31);
% % l.Orientation = 'Vertical';
% % l.Location = 'best'
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% %Setup axis font size   
% ax3.FontSize = ax_fnt;
% % axis(ax,[0 200 -40 60]);
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax3.XColor = black;         % Box horizontal lines' color
% ax3.YColor = black;         % Box vertical lines' color
% set(ax3,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax3,'Cross-track \ error, m','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax3,' $$\dot{\chi}$$, rad./s','Interpreter','Latex','FontSize',lbl_fnt)
% % 
% f4;
% % figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
% % h4 = plot(ax4,x1,kappa_d,'DisplayName',['n = ',num2str(n(i))]);hold(ax4,'on') ;grid(ax4,'on')
% % h4.Color = red;                % Color of plot
% % h4.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% % h4.LineWidth = lw_;             % Plot line width
% % h4.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% % h4.MarkerSize = ms_;            % Marker Size
% % h4.MarkerEdgeColor = 'auto';    % Marker edge color
% % h4.MarkerFaceColor = 'none';    % Marker face color
% 
% h41 = plot(ax4,x1,kappa_a);hold(ax4,'on') ;grid(ax4,'on')
% h41.Color = blue;                % Color of plot
% h41.LineStyle = '-';            % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h41.LineWidth = lw_;             % Plot line width
% h41.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h41.MarkerSize = ms_;            % Marker Size
% h41.MarkerEdgeColor = 'auto';    % Marker edge color
% h41.MarkerFaceColor = 'none';    % Marker face color
% % l = legend(ax4,[h4 h41],'Desired','Actual');
% % l.Orientation = 'Vertical';
% % l.Location = 'best';
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% % Setup axis font size
% % ax4 = gca;     
% ax4.FontSize = ax_fnt;
% % axis(ax4,[0 max(t(end),t1(end)) -5 5]);
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax4.XColor = black;         % Box horizontal lines' color
% ax4.YColor = black;         % Box vertical lines' color
% set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax4,'Cross-track \ error, m','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax4,'Curvature, m $^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)
% 
% f5;
% % h5 = plot(ax5,t,kai_d*(180/pi));hold(ax5,'on') ;grid(ax5,'on')
% % h5.Color = red;                % Color of plot
% % h5.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% % h5.LineWidth = lw_;             % Plot line width
% % h5.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% % h5.MarkerSize = ms_;            % Marker Size
% % h5.MarkerEdgeColor = 'auto';    % Marker edge color
% % h5.MarkerFaceColor = 'none';    % Marker face color
% h51 = plot(ax5,t,x3*(180/pi));hold(ax5,'on') ;grid(ax5,'on');
% h51.Color = red;                % Color of plot
% h51.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h51.LineWidth = lw_;             % Plot line width
% h51.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h51.MarkerSize = ms_;            % Marker Size
% h51.MarkerEdgeColor = 'auto';    % Marker edge color
% h51.MarkerFaceColor = 'none';    % Marker face color
% % l = legend(ax5,h51,'Desired','Actual');
% % l.Orientation = 'Vertical';
% % l.Location = 'best';
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% % Setup axis font size   
% ax5.FontSize = ax_fnt;
% ax5.FontSize = ax_fnt;
% yticks(ax5,[90 105 120 135 150 165 180])                           % Y-axis grid line positioning
% yticklabels(ax5,{'90', '105', '120', '135', '150', '165','180'})
% % axis(ax,[0 200 -40 60]);
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax5.XColor = black;         % Box horizontal lines' color
% ax5.YColor = black;         % Box vertical lines' color
% set(ax5,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax5,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax5,' $$ \chi$$,  deg.','Fontsize',lbl_fnt);
% 
% f6;
% % h6 = plot(ax6,x1,kai_d*(180/pi),'DisplayName',['n = ',num2str(n(i))]);hold(ax6,'on') ;grid(ax6,'on')
% % h6.Color = red;                % Color of plot
% % h6.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% % h6.LineWidth = lw_;             % Plot line width
% % h6.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% % h6.MarkerSize = ms_;            % Marker Size
% % h6.MarkerEdgeColor = 'auto';    % Marker edge color
% % h6.MarkerFaceColor = 'none';    % Marker face color
% h61 = plot(ax6,x1,x3*(180/pi));hold(ax6,'on') ;grid(ax6,'on');
% h61.Color = blue;                % Color of plot
% h61.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h61.LineWidth = lw_;             % Plot line width
% h61.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h61.MarkerSize = ms_;            % Marker Size
% h61.MarkerEdgeColor = 'auto';    % Marker edge color
% h61.MarkerFaceColor = 'none';    % Marker face color
% % l = legend(ax6,[h6 h61],'Desired','Actual');
% % l.Orientation = 'Vertical';
% % l.Location = 'best';
% % l.FontSize = leg_fnt;
% % l.Interpreter = 'Latex';
% % Setup axis font size   
% ax6.FontSize = ax_fnt;
% yticks(ax6,[90 105 120 135 150 165 180])                           % Y-axis grid line positioning
% yticklabels(ax6,{'90', '105', '120', '135', '150', '165','180'})     % Y-axis grid line labels
% % axis(ax,[0 200 -40 60]);
% % Outer box setup
% box on                      % Switch on the box around the axis
% ax6.XColor = black;         % Box horizontal lines' color
% ax6.YColor = black;         % Box vertical lines' color
% set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% % Label and limit setup
% xlabel(ax6,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
% ylabel(ax6,'$$ \chi$$,  deg.','Fontsize',lbl_fnt);
% 
% 
%     %-----------------------------------------------------%
% 
% 
% 
% %     
% %     % Setup axis font size
% %     ax = gca;
% %     ax.FontSize = ax_fnt;
% %     
% %     % Grid setup
% %     grid on
% %     ax.GridLineStyle = '-';
% %     ax.GridAlpha = 0.2;                                     % Grid transperancy
% %     xticks([-50 0 50 100 8 10])                                  % x-axis grid line positioning
% %     xticklabels({'0', '2', '4', '6', '8', '10'})            % x-axis grid line labels
% %     yticks([-1 -0.5 0 0.5 1 1.5])                           % Y-axis grid line positioning
% %     yticklabels({'-1', '-0.5', '0', '0.5', '1', '1.5'})     % Y-axis grid line labels
% %     
% %     % Outer box setup
% %     box on                      % Switch on the box around the axis
% %     ax.XColor = black_;         % Box horizontal lines' color
% %     ax.YColor = black_;         % Box vertical lines' color
% %     set(ax,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% %     
% %     % Label and limit setup
% %     xlabel('$$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
% %     ylabel('$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
% %     xlim([min(t) max(t)])
% %     ylim([1.2*min([x3 y4]) 1.5*max([x4 y5])])
% %     %-----------------------------------------------------%
% end
% figure_plot = [f1 f2 f3 f4 f5 f6];
% end