%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


f8 = figure;
ax8 = axes;
% Plot format control variables
lw_ = 2;            % Line width
ms_ = 6;            % Marker size
ax_fnt = 16;        % Axis font size
lbl_fnt = ax_fnt+1; % Label font size
leg_fnt = ax_fnt-1; % Legend font size
ax_lw = 1.5;        % Axis line width
black = [0 0 0];
%--------------- declaration of constant parameters ----------------------%
global  kv k_kai k_kaidot k2 figure_type_good
kv = 5;
k_kai = 30;
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
figure_type_good = 1;
Vgd = 10;
k2 = 0.002;

as = []; cte = [];vg = [];
for i = 1:500
    i
wx = -2+4*rand(1);
wy = -2+4*rand(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 30];

        
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

        options = odeset('RelTol',1e-11,'AbsTol',1e-11);  
        %--------- ode solver   ------------------------------------------%
        [t,x] = ode45(@(t,x)fun_stline_withwind(t,x,Vgd,wx,wy) ,tspan, x_initial,options);

        CTE = x(:,1);
        cte(i) = CTE(end);
        % cte(i) = mean(CTE);
        AIRSPEED = x(:,5);
        as(i) = AIRSPEED(end);

        psi = x(:,3) + asin((wx.*sin(x(:,3)) - wy.*cos(x(:,3)))./ x(:,5) );
        VG = sqrt((x(:,5).*cos(psi) + wx).^2 + (x(:,5).*sin(psi) + wy).^2) ;
        vg(i) = VG(end);
       
end
d1 = mean(as);
d2 = mean(cte);
d3 = mean(vg);


f8;
h8 = plot(ax8,t,Vgd*ones(length(t),1),'r','linewidth',2);hold(ax8,'on') ;grid(ax8,'on')
% h8.Color = 'red';                 % Color of plot
% h8.LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h8.LineWidth = lw_;             % Plot line width
% h8.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h8.MarkerSize = ms_;            % Marker Size
% h8.MarkerEdgeColor = 'auto';    % Marker edge color
% h8.MarkerFaceColor = 'none';    % Marker face color
h81 = plot(ax8,t,VG,'b','linewidth',2);hold(ax8,'on') ;grid(ax8,'on');
% h81.Color = 'blue';                % Color of plot
h81.LineStyle = '--';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
% h81.LineWidth = lw_;             % Plot line width
% h81.Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
% h81.MarkerSize = ms_;            % Marker Size
% h81.MarkerEdgeColor = 'auto';    % Marker edge color
% h81.MarkerFaceColor = 'none';    % Marker face color
l = legend(ax8,[h8 h81],'Desired','Actual');
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

function out = fun_stline_withwind(t,x,Vgd,wx,wy )
global k2  kv k_kai k_kaidot


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





