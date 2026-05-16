function figure_plot = fun_comp_stlineplot(t,x1,x2,x5,x6,...
    kaidot_a_prop,kaidot_a_Nelson,kappa_a_prop,kappa_a_Nelson,...
    kai_a_prop,kai_a_Nelson,x_max_prop,val_max_prop,x_max_Nelson,val_max_Nelson)
global figure_type_good p_x0 p_y0 p_width p_height

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
% Plot format control variables
lw_ = 3;            % Line width
ms_ = 6;            % Marker size
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

if(figure_type_good == 1)
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 

f1;
% h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
xl = xline(ax1,0,'-k','linewidth',3);hold(ax1,'on'); 
h1 = plot(ax1,x1,x2,'Color',red);hold(ax1,'on');grid(ax1,'on');
h1.Color = red;                 % Color of plot
h1.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h1.LineWidth = lw_;             % Plot line width
h1.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h1.MarkerSize = ms_;            % Marker Size
h1.MarkerEdgeColor = 'auto';    % Marker edge color
h1.MarkerFaceColor = 'none';    % Marker face color
% plot(ax1,x_ini,y_ini,'-o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','green'); hold(ax1,'on');
% plot(ax1,x_end,y_end,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','cyan'); hold(ax1,'on');
h11 = plot(ax1,x5,x6,'Color',blue);hold(ax1,'on');
h11.Color = blue;                 % Color of plot
h11.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h11.LineWidth = lw_;             % Plot line width
h11.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h11.MarkerSize = ms_;            % Marker Size
h11.MarkerEdgeColor = 'auto';    % Marker edge color
h11.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax1,[xl h1 h11] ,'Desired path','Proposed','Nelson et al.');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size
% ax1 = gca;     
ax1.FontSize = ax_fnt;
% axis(ax1,'equal')
% Outer box setup
box(ax1,'on')                     % Switch on the box around the axis
ax1.XColor = black;         % Box horizontal lines' color
ax1.YColor = black;         % Box vertical lines' color
set(ax1,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax1,'$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
axis(ax1,'equal')

f2;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
h2 = plot(ax2,t,x1);hold(ax2,'on') ;grid(ax2,'on');
h2.Color = red;                % Color of plot
h2.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h2.LineWidth = lw_;             % Plot line width
h2.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h2.MarkerSize = ms_;            % Marker Size
h2.MarkerEdgeColor = 'auto';    % Marker edge color
h2.MarkerFaceColor = 'none';    % Marker face color
h21 = plot(ax2,t,x5);hold(ax2,'on') ;grid(ax2,'on');
h21.Color = blue;                % Color of plot
h21.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h21.LineWidth = lw_;             % Plot line width
h21.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h21.MarkerSize = ms_;            % Marker Size
h21.MarkerEdgeColor = 'auto';    % Marker edge color
h21.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax2,[h2 h21],'Proposed','Nelson et al. ');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
%Setup axis font size   
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


f3;
h3 = plot(ax3,x1,kaidot_a_prop);hold(ax3,'on');grid(ax3,'on');
h3.Color = red;                % Color of plot
h3.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h3.LineWidth = lw_;             % Plot line width
h3.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h3.MarkerSize = ms_;            % Marker Size
h3.MarkerEdgeColor = 'auto';    % Marker edge color
h3.MarkerFaceColor = 'none';    % Marker face color

h31 = plot(ax3,x5,kaidot_a_Nelson);hold(ax3,'on');grid(ax3,'on');
h31.Color = blue;                % Color of plot
h31.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h31.LineWidth = lw_;             % Plot line width
h31.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h31.MarkerSize = ms_;            % Marker Size
h31.MarkerEdgeColor = 'auto';    % Marker edge color
h31.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax3,[h3 h31],'Proposed','Nelson et al. ');
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
xlabel(ax3,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax3,' $$\dot{\chi}$$, rad./s','Interpreter','Latex','FontSize',lbl_fnt)
% 
f4;
% figure_pos = figure('position', [p_x0, p_y0, p_width, p_height]); 
h4 = plot(ax4,x1,kappa_a_prop);hold(ax4,'on') ;grid(ax4,'on')
h4.Color = red;                % Color of plot
h4.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h4.LineWidth = lw_;             % Plot line width
h4.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h4.MarkerSize = ms_;            % Marker Size
h4.MarkerEdgeColor = 'auto';    % Marker edge color
h4.MarkerFaceColor = 'none';    % Marker face color

h41 = plot(ax4,x5,kappa_a_Nelson);hold(ax4,'on') ;grid(ax4,'on')
h41.Color = blue;                % Color of plot
h41.LineStyle = '-';            % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h41.LineWidth = lw_;             % Plot line width
h41.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h41.MarkerSize = ms_;            % Marker Size
h41.MarkerEdgeColor = 'auto';    % Marker edge color
h41.MarkerFaceColor = 'none';    % Marker face color
h42 = plot(ax4,x_max_prop,val_max_prop,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','cyan'); hold(ax4,'on');grid(ax4,'on');
h42.Annotation.LegendInformation.IconDisplayStyle = 'off';
h43 = plot(ax4,x_max_Nelson,val_max_Nelson,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green'); hold(ax4,'on');grid(ax4,'on');
h43.Annotation.LegendInformation.IconDisplayStyle = 'off';
dt_prop = datatip(h4,x_max_prop,val_max_prop);hold(ax4,'on');grid(ax4,'on');
h4.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h4.DataTipTemplate.DataTipRows(2).Format = '%.3f';
dt_Nelson = datatip(h41,x_max_Nelson,val_max_Nelson);hold(ax4,'on');grid(ax4,'on');
h41.DataTipTemplate.DataTipRows(1).Format = '%.2f';
h41.DataTipTemplate.DataTipRows(2).Format = '%.3f';

l = legend(ax4,[h4 h41],'Proposed','Nelson et al.');
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
xlabel(ax4,'Cross-track error $$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax4,'Curvature $$\kappa$$,  m $$^{-1}$$ ','Interpreter','Latex','FontSize',lbl_fnt)

f5;
h5 = plot(ax5,t,kai_a_prop*(180/pi));hold(ax5,'on') ;grid(ax5,'on')
h5.Color = red;                % Color of plot
h5.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h5.LineWidth = lw_;             % Plot line width
h5.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h5.MarkerSize = ms_;            % Marker Size
h5.MarkerEdgeColor = 'auto';    % Marker edge color
h5.MarkerFaceColor = 'none';    % Marker face color
h51 = plot(ax5,t,kai_a_Nelson*(180/pi));hold(ax5,'on') ;grid(ax5,'on');
h51.Color = blue;                % Color of plot
h51.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h51.LineWidth = lw_;             % Plot line width
h51.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h51.MarkerSize = ms_;            % Marker Size
h51.MarkerEdgeColor = 'auto';    % Marker edge color
h51.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax5,[h5 h51],'Proposed','Nelson et al.');
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
ylabel(ax5,' $$ \chi$$,  deg.','Fontsize',lbl_fnt);

f6;
h6 = plot(ax6,x1,kai_a_prop*(180/pi));hold(ax6,'on') ;grid(ax6,'on')
h6.Color = red;                % Color of plot
h6.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h6.LineWidth = lw_;             % Plot line width
h6.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h6.MarkerSize = ms_;            % Marker Size
h6.MarkerEdgeColor = 'auto';    % Marker edge color
h6.MarkerFaceColor = 'none';    % Marker face color
h61 = plot(ax6,x5,kai_a_Nelson*(180/pi));hold(ax6,'on') ;grid(ax6,'on');
h61.Color = blue;                % Color of plot
h61.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h61.LineWidth = lw_;             % Plot line width
h61.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h61.MarkerSize = ms_;            % Marker Size
h61.MarkerEdgeColor = 'auto';    % Marker edge color
h61.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax6,[h6 h61],'Proposed','Nelson et al. ');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax6.FontSize = ax_fnt;
yticks(ax6,[90 105 120 135 150 165 180])                           % Y-axis grid line positioning
yticklabels(ax6,{'90', '105', '120', '135', '150', '165','180'})
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax6.XColor = black;         % Box horizontal lines' color
ax6.YColor = black;         % Box vertical lines' color
set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax6,'Cross-track error $$x$$,  m','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax6,'$$ \chi$$,  deg.','Fontsize',lbl_fnt);


f7;
h7 = plot(ax7,t,kaidot_a_prop);hold(ax7,'on') ;grid(ax7,'on')
h7.Color = red;                % Color of plot
h7.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h7.LineWidth = lw_;             % Plot line width
h7.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h7.MarkerSize = ms_;            % Marker Size
h7.MarkerEdgeColor = 'auto';    % Marker edge color
h7.MarkerFaceColor = 'none';    % Marker face color
h71 = plot(ax7,t,kaidot_a_Nelson);hold(ax7,'on') ;grid(ax7,'on');
h71.Color = blue;                % Color of plot
h71.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h71.LineWidth = lw_;             % Plot line width
h71.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h71.MarkerSize = ms_;            % Marker Size
h71.MarkerEdgeColor = 'auto';    % Marker edge color
h71.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax7,[h7 h71],'Proposed','Nelson et al. ');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax7.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax7.XColor = black;         % Box horizontal lines' color
ax7.YColor = black;         % Box vertical lines' color
set(ax7,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax7,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax7,' $$\dot{\chi}$$, rad/s','Interpreter','Latex','FontSize',lbl_fnt)
% 

f8;
h8 = plot(ax8,t,kappa_a_prop);hold(ax8,'on') ;grid(ax8,'on')
h8.Color = red;                % Color of plot
h8.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h8.LineWidth = lw_;             % Plot line width
h8.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h8.MarkerSize = ms_;            % Marker Size
h8.MarkerEdgeColor = 'auto';    % Marker edge color
h8.MarkerFaceColor = 'none';    % Marker face color
h81 = plot(ax8,t,kappa_a_Nelson);hold(ax8,'on') ;grid(ax8,'on');
h81.Color = blue;                % Color of plot
h81.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
h81.LineWidth = lw_;             % Plot line width
h81.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
h81.MarkerSize = ms_;            % Marker Size
h81.MarkerEdgeColor = 'auto';    % Marker edge color
h81.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax8,[h8 h81],'Proposed','Nelson et al. ');
l.Orientation = 'Vertical';
l.Location = 'best';
l.FontSize = leg_fnt;
l.Interpreter = 'Latex';
% Setup axis font size   
ax8.FontSize = ax_fnt;
% axis(ax,[0 200 -40 60]);
% Outer box setup
box on                      % Switch on the box around the axis
ax8.XColor = black;         % Box horizontal lines' color
ax8.YColor = black;         % Box vertical lines' color
set(ax8,'linewidth',ax_lw)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax8,'$$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
ylabel(ax8,'Curvature $$\kappa$$,  m $$^{-1}$$ ','Interpreter','Latex','FontSize',lbl_fnt)
    %-----------------------------------------------------%



%     
%     % Setup axis font size
%     ax = gca;
%     ax.FontSize = ax_fnt;
%     
%     % Grid setup
%     grid on
%     ax.GridLineStyle = '-';
%     ax.GridAlpha = 0.2;                                     % Grid transperancy
%     xticks([-50 0 50 100 8 10])                                  % x-axis grid line positioning
%     xticklabels({'0', '2', '4', '6', '8', '10'})            % x-axis grid line labels
%     yticks([-1 -0.5 0 0.5 1 1.5])                           % Y-axis grid line positioning
%     yticklabels({'-1', '-0.5', '0', '0.5', '1', '1.5'})     % Y-axis grid line labels
%     
%     % Outer box setup
%     box on                      % Switch on the box around the axis
%     ax.XColor = black_;         % Box horizontal lines' color
%     ax.YColor = black_;         % Box vertical lines' color
%     set(ax,'linewidth',ax_lw)   % Axis linewidth (box and grid)
%     
%     % Label and limit setup
%     xlabel('$$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
%     ylabel('$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
%     xlim([min(t) max(t)])
%     ylim([1.2*min([x3 y4]) 1.5*max([x4 y5])])
%     %-----------------------------------------------------%
end
figure_plot = [f1 f2 f3 f4];
end