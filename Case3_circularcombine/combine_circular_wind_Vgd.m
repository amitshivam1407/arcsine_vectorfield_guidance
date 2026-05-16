%-------------------------------------------------------------------------%
%------               13th October 2022                          ---------%
%--------           Autopilot design with no wind                 --------%
%------- Circular path following using both initial conditions    --------%
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
f9 = figure;
ax9 = axes;

%%%%%%%%%%%%%%%% declaration of constant parameters %%%%%%%%%%%%%%%%%%%%%%%
global rd k_kai k_kaidot kv k2 wx wy x_c y_c
rd = 50;
wx = -1;
wy = 2;
kv = 10;
k_kai = 80;
k_kaidot = sqrt(k_kai);
k2 = 0.005;
%%%%%%%%%%%%%%%%%%%%%% standoff circle    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_c = 0;
y_c = 0;
theta = 0:2*pi/1000:2*pi;
x_sc = x_c + rd*cos(theta);
y_sc = y_c + rd*sin(theta);
%-------------------------------------------------------------------------%
%----------------  vector field construction  ----------------------------%
range  = -130:8:130;
[X,Y] = meshgrid(range);
r = sqrt((X-x_c).^2 + (Y-y_c).^2);
gama = atan2((Y-y_c),(X-x_c));
[xdot, ydot] = fun_vf(r,gama,10);
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
%%%%%%%%%%%%%%%%%%% initial conditions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 20];


Vgd = 12;
% for k = 1:length(Vgd)
r_0 = [100 10];
gamma_0 = deg2rad(225);
for i = 1:length(r_0)
    x0(i) = r_0(i)*cos(gamma_0) + x_c;
    y0(i) = r_0(i)*sin(gamma_0) + y_c;

    if r_0(i) < rd
%         Vgd = 12;
        kai0(i) = gamma_0 + asin(1./(1+k2*(r_0(i)-rd).^2)) ;
        gamma_dot0(i) = (Vgd./r_0(i)).*sin(kai0(i) - gamma_0) ;
        kaidot0(i) = gamma_dot0(i) - (2*k2*Vgd.*(r_0(i)-rd))./((1+k2*(r_0(i)-rd).^2).^2);
    else
%         Vgd = 14;
        kai0(i) = gamma_0 + pi - asin(1./(1+k2*(r_0(i)-rd).^2)) ;
        gamma_dot0(i) = (Vgd./r_0(i)).*sin(kai0(i) - gamma_0) ;
        kaidot0(i) = gamma_dot0(i) - (2*k2*Vgd.*(r_0(i)-rd))./((1+k2*(r_0(i)-rd).^2).^2);
    end


    kai0_deg(i) = rad2deg(wrapToPi(kai0(i)));
    Va0 = 10;

    x_initial(i,:) = [r_0(i);gamma_0;x0(i);y0(i);kai0(i);kaidot0(i);Va0];

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);

    %%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t,x] = ode45(@(t,x)fun_circular_wind(t,x,Vgd) ,tspan, x_initial(i,:),options);

    r_prop = x(:,1) ;
    error_r = (r_prop - rd) ;
    x_ini = x(1,3) + x_c ;
    y_ini = x(1,4) + y_c ;
    x_end = x(end,3) + x_c ;
    y_end = x(end,4) + y_c ;

    [kaid, kaidot,Vad] = fun_propkaid(x(:,1),x(:,2),Vgd);
    kappa_des = kaidot./Vgd ;
    kappa_actual = x(:,6)./Vgd ;
    [val(i), ind(i)] = max(kappa_actual);
    rmax(i) = x(ind(i),1);


    psi = x(:,5) + asin((wx.*sin(x(:,5))-wy.*cos(x(:,5)))./x(:,7));
    Vg = sqrt(x(:,7).^2  + wx^2 + wy^2 + 2*(x(:,7).*wx.*cos(psi ) + x(:,7).*wy.*sin(psi )));
    %-------------------------------------------------------------------------%
    % Plot format control variables
    lw = 3;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 17;        % Axis font size
    lbl_fnt = 20; % Label font size
    leg_fnt = 16; % Legend font size
    ax_lw = 3;        % Axis line width
    %% Colors
    blue = [0 0.4470 0.7410];
    red = [0.8500 0.3250 0.0980];
    orange = [0.9290 0.6940 0.1250];
    violet = [0.4940 0.1840 0.5560];
    green = [0.4660 0.6740 0.1880];
    cyan = [0.3010 0.7450 0.9330];
    maroon = [0.6350 0.0780 0.1840];
    black = [0 0 0];
    color = [green;maroon];
    %%
    f1;
    

    if i==1
         xd = plot(ax1,x_sc,y_sc,'k','linewidth',3,'DisplayName','Desired orbit');hold(ax1,'on');
          h11 = quiver(ax1,X,Y,xdot,ydot,'color',[0.75 0.75 0.75],'linewidth',1,'DisplayName','Vector field');hold(ax1,'on');
        % xd.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h1(i) = plot(ax1,x(:,3),x(:,4),'DisplayName','From outside the circle');hold on;
%         h1(i) = plot(ax1,x(:,3),x(:,4),'DisplayName','Outside, $$  r_{0} = 100 $$ m');hold on;
    else
        h1(i) = plot(ax1,x(:,3),x(:,4),'DisplayName','From inside the circle');hold on;
    end
    h1(i).Color = color(i,:);                 % Color of plot
    h1(i).LineStyle = '-';             % Line style ('--' -> Dash, '-.' -> Dot+Dash, ':' -> Dot)
    h1(i).LineWidth = lw;             % Plot line width
    h1(i).Marker = 'none';             % Marker type ('s' -> Square, 'o' -> Circle, 'd' -> Diamond)
    h1(i).MarkerSize = ms;            % Marker Size
    h1(i).MarkerEdgeColor = 'auto';    % Marker edge color
    h1(i).MarkerFaceColor = 'none';    % Marker face color
    hc = plot(ax1,x_c,y_c,'-o','linewidth',2,'MarkerSize',8,...
        'MarkerEdgeColor','blue',...
        'MarkerFaceColor','magenta'); hold(ax1,'on');
    hc.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_in = plot(ax1,x_ini,y_ini,'-o','linewidth',2,'MarkerSize',8,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','green'); hold(ax1,'on');
    h_in.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     if r_0(i) < rd
%         h11 = quiver(ax1,X_in,Y_in,xdot_in,ydot_in,'color',[0.75 0.75 0.75],'linewidth',1,'DisplayName','Vector field');hold(ax1,'on');
%     else
%         h22 =  quiver(ax1,X_out,Y_out,xdot_out,ydot_out,'color',[0.75 0.75 0.75],'linewidth',1);hold(ax1,'on');
%     end
%     % h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     h22.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % plot(ax1,x_end,y_end,'-s','MarkerSize',10,...
    %     'MarkerEdgeColor','blue',...
    %     'MarkerFaceColor','cyan'); hold(ax1,'on');

    % Setup axis font size
    ax1.FontSize = ax_fnt;
    axis(ax1,'equal')
    % Outer box setup
    box on                        % Switch on the box around the axis
    ax1.XColor = 'black';         % Box horizontal lines' color
    ax1.YColor = 'black';         % Box vertical lines' color
    set(ax1,'linewidth',ax_lw)    % Axis linewidth (box and grid)
    % Label and limit setup
    xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax1,'$$y$$, m','Interpreter','Latex','FontSize',lbl_fnt)
    axis(ax1,'equal')

    f2;
    h2 = plot(ax2,t,r_prop,'-','LineWidth',lw,'Color',color(i,:));hold(ax2,'on') ;grid(ax2,'on');

    % Setup axis font size
    ax2.FontSize = ax_fnt;

    % Outer box setup
    box on                        % Switch on the box around the axis
    ax2.XColor = 'black';         % Box horizontal lines' color
    ax2.YColor = 'black';         % Box vertical lines' color
    set(ax2,'linewidth',ax_lw)    % Axis linewidth (box and grid)
    % Label and limit setup
    xlabel(ax2,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax2,'Radial distance, m','Interpreter','Latex','FontSize',lbl_fnt)

    f3;
    h3 = plot(ax3,t,error_r,'-','LineWidth',lw,'Color',color(i,:));hold(ax3,'on');grid(ax3,'on');
    %Setup axis font size
    ax3.FontSize = ax_fnt;
    % axis(ax,[0 200 -40 60]);
    % Outer box setup
    box on                      % Switch on the box around the axis
    ax3.XColor = 'black';         % Box horizontal lines' color
    ax3.YColor = 'black';         % Box vertical lines' color
    set(ax3,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    % Label and limit setup
    xlabel(ax3,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax3,'Radial error, m','Interpreter','Latex','FontSize',lbl_fnt)


    f4;
    if i ==1
        h4(i) = plot(ax4,r_prop,kaidot,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 100 $$ m');hold(ax4,'on');grid(ax4,'on');
        h41(i) = plot(ax4,r_prop,x(:,6),'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 100 $$ m');grid(ax4,'on');
    else
        h4(i) = plot(ax4,r_prop,kaidot,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 10 $$ m');hold(ax4,'on');grid(ax4,'on');
        h41(i) = plot(ax4,r_prop,x(:,6),'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 10 $$ m');grid(ax4,'on');
    end

    ax4.FontSize = ax_fnt;
    % axis(ax4,[0 max(t(end),t1(end)) -5 5]);
    % Outer box setup
    box on                      % Switch on the box around the axis
    ax4.XColor = 'black';         % Box horizontal lines' color
    ax4.YColor = 'black';         % Box vertical lines' color
    set(ax4,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    xlabel(ax4,'Radial distance, m','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax4,'$\dot{\chi}$,  deg./s ','Interpreter','Latex','FontSize',lbl_fnt)

    f5;
    if i ==1
        h5(i) = plot(ax5,r_prop,kappa_des,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 100 $$ m');hold(ax5,'on');grid(ax5,'on');
        h51(i) = plot(ax5,r_prop,kappa_actual,'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 100 $$ m');grid(ax5,'on');
    else
        h5(i) = plot(ax5,r_prop,kappa_des,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 10 $$ m');hold(ax5,'on');grid(ax5,'on');
        h51(i) = plot(ax5,r_prop,kappa_actual,'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 10 $$ m');grid(ax5,'on');
    end

    %% Setup axis font size
    ax5.FontSize = ax_fnt;

    % Outer box setup
    box on                      % Switch on the box around the axis
    ax5.XColor = 'black';         % Box horizontal lines' color
    ax5.YColor = 'black';         % Box vertical lines' color
    set(ax5,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    xlabel(ax5,'Radial distance, m','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax5,'Curvature, m $ ^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)
    %
    %%
    f6;
    if i ==1
        h6(i) = plot(ax6,t,kappa_des,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 100 $$ m');hold(ax6,'on');grid(ax6,'on');
        h61(i) = plot(ax6,t,kappa_actual,'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 100 $$ m');grid(ax6,'on');
    else
        h6(i) = plot(ax6,t,kappa_des,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 10 $$ m');hold(ax6,'on');grid(ax6,'on');
        h61(i) = plot(ax6,t,kappa_actual,'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 10 $$ m');grid(ax6,'on');
    end

    %% Setup axis font size
    ax6.FontSize = ax_fnt;
    box on                      % Switch on the box around the axis
    ax6.XColor = 'black';         % Box horizontal lines' color
    ax6.YColor = 'black';         % Box vertical lines' color
    set(ax6,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    xlabel(ax6,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax6,'Curvature, m $ ^{-1} $ ','Interpreter','Latex','FontSize',lbl_fnt)

    f7;
    if i ==1
        h7(i) = plot(ax7,t,Vad,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 100 $$ m');hold(ax7,'on');grid(ax7,'on');
        h71(i) = plot(ax7,t,x(:,7),'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 100 $$ m');grid(ax7,'on');
    else
        h7(i) = plot(ax7,t,Vad,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$ r_{0} = 10 $$ m');hold(ax7,'on');grid(ax7,'on');
        h71(i) = plot(ax7,t,x(:,7),'--','LineWidth',lw,'Color',color(i,:),'DisplayName','Achieved, $$ r_{0} = 10 $$ m');grid(ax7,'on');
    end

    % Setup axis font size
    ax7.FontSize = ax_fnt;
    %     axis(ax7,[0 t(end) 0 15]);
    %% Outer box setup
    box on                      % Switch on the box around the axis
    ax7.XColor = 'black';         % Box horizontal lines' color
    ax7.YColor = 'black';         % Box vertical lines' color
    set(ax7,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    %% Label and limit setup
    xlabel(ax7,' $$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax7,' $$ V_{\mathrm{a}}$$,  m/s','Fontsize',lbl_fnt);

    f8;
    if r_0(i) < rd        
        h81(i) = plot(ax8,t,Vg,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','From inside the circle');hold(ax8,'on') ;grid(ax8,'on');
       
    else
         h8(i)  = plot(ax8,t,Vgd*ones(size(t)),'--','LineWidth',lw,'Color','b','DisplayName','Desired ground speed');hold(ax8,'on') ;grid(ax8,'on')
%         h8(i)  = plot(ax8,t,Vgd*ones(size(t)),'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$r_{0} = 100 $$ m');hold(ax8,'on') ;grid(ax8,'on')
        h81(i) = plot(ax8,t,Vg,'-','LineWidth',lw,'Color',color(i,:),'DisplayName','From outside the circle');hold(ax8,'on') ;grid(ax8,'on');

    end

    % Setup axis font size
    ax8.FontSize = ax_fnt;
    axis(ax8,[0 t(end) 0 20]);
    % % Outer box setup
    box on                     % Switch on the box around the axis
    ax8.XColor = 'black';         % Box horizontal lines' color
    ax8.YColor = 'black';         % Box vertical lines' color
    set(ax8,'linewidth',ax_lw)   % Axis linewidth (box and grid)
%     yticks(ax8,[0 5 10 12  15 20]);
%     yticklabels(ax8,{'0' '5' '10' '12' '15' '20'});
    % % Label and limit setup
    xlabel(ax8,'$$t$$, s','Interpreter','Latex','FontSize',lbl_fnt)
    ylabel(ax8,' $$ V_{\mathrm{g}}$$,  m/s','Fontsize',lbl_fnt);
    %axis(ax8,[])

    f9;
    if r_0(i) < rd
%         h9(i) = plot(ax9,t,wrapToPi(kaid)*(180/pi),'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$r_{0} = 10 $$ m');hold(ax9,'on') ;grid(ax9,'on')
         h91(i) = plot(ax9,t,wrapToPi(x(:,5))*(180/pi),'-','LineWidth',lw,'Color',color(i,:),'DisplayName','From inside the circle');hold(ax9,'on') ;grid(ax9,'on');

    else
        %h9(i) = plot(ax9,t,wrapToPi(kaid)*(180/pi),'-','LineWidth',lw,'Color',color(i,:),'DisplayName','Commanded, $$r_{0} = 100 $$ m');hold(ax9,'on') ;grid(ax9,'on')
        h91(i) = plot(ax9,t,wrapToPi(x(:,5))*(180/pi),'-','LineWidth',lw,'Color',color(i,:),'DisplayName','From outside the circle');hold(ax9,'on') ;grid(ax9,'on');
    end
    ax9.FontSize = ax_fnt;
    yticks(ax9,[-180 -120 -60 0 60 120 180 ]);
    yticklabels(ax9,{'-180' '-120' '-60' '0' '60' '120' '180'});
    % % Outer box setup
    box on                     % Switch on the box around the axis
    ax9.XColor = 'black';         % Box horizontal lines' color
    ax9.YColor = 'black';         % Box vertical lines' color
    set(ax9,'linewidth',ax_lw)   % Axis linewidth (box and grid)
    xlabel(ax9,' $ t,$  s','Fontsize',18)
    ylabel(ax9,' $ \chi$,  deg. ','Fontsize',18)

end

legend(ax1,'NumColumns',1);
Legend = cell(length(r_0),1);
for i = 1:length(r_0)
    str2 = sprintf('$$r_{0} = %d $$ m',r_0(i)) ;
    Legend{i} = str2;
end
legend(ax2,Legend,'FontSize',leg_fnt );

Legend = cell(length(r_0),1);
for i = 1:length(r_0)
    str3 = sprintf('$$r_{0} = %d $$ m',r_0(i)) ;
    Legend{i} = str3;
end
legend(ax3,Legend,'FontSize',leg_fnt );
legend(ax4)
legend(ax5)
legend(ax6)
legend(ax7)
legend(ax8)
legend(ax9)
% Legend9 = cell(length(r_0),1);
% for i = 1:length(r_0)
%     Legend9{i} =  ['$$r_{0} = %d $$ m', num2str(r_0(i))];
% end
% legend(kaid_plot,Legend9,'Fontsize',15 )

function out = fun_circular_wind(t,x,Vgd)
global  k_kai k_kaidot rd k2  kv wx wy


if x(1) < rd
    kaid = x(2) + asin(1./(1+k2*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k2*Vgd.*(x(1)-rd))./((1+k2*(x(1) -rd).^2).^2);
    Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
    psi = x(5) + asin((wx*sin(x(5))-wy*cos(x(5)))./x(7));
    Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos( psi) + 2*wy*x(7).*sin( psi));

else
    kaid = x(2) + pi - asin(1./(1+k2*(x(1) -rd).^2)) ;
    gamma_dot = (Vgd./x(1)).*sin(kaid - x(2)) ;
    kaid_dot = gamma_dot - (2*k2*Vgd.*(x(1) -rd))./((1+k2*(x(1) -rd).^2).^2);
    Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
    psi = x(5) + asin((wx*sin(x(5))-wy*cos(x(5)))./x(7));
    Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos( psi) + 2*wy*x(7).*sin( psi));

end


%     Vad = sqrt(Vgd_in^2  + wx^2 +wy^2 - 2*(Vgd_in*wx*cos(kaid) + Vgd_in*wy*sin(kaid)));
out(1,1) = Vg.*cos(x(5) - x(2));
out(2,1) = (Vg./x(1)).*sin(x(5) - x(2));
out(3,1) = x(7)*cos(psi) + wx ;
out(4,1) = x(7)*sin(psi) + wy ;
out(5,1) = x(6) ;
out(6,1) = k_kai*(kaid - x(5)) + k_kaidot*(kaid_dot - x(6))  ;
out(7,1) = kv*(Vad - x(7));
% else
%     Vad = sqrt(Vgd_out^2  + wx^2 +wy^2 - 2*(Vgd_out*wx*cos(kaid) + Vgd_out*wy*sin(kaid)));
% end
% psi = x(5) + asin((wx*sin(x(5))-wy*cos(x(5)))./x(7));
% Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos( psi) + 2*wy*x(7).*sin( psi));
% Vg = sqrt(x(7)^2  + wx^2 + wy^2 + 2*(x(7)*wx*cos(x(5)) + x(7)*wy*sin(x(5))));
% psi = x(5) + asin((wx*sin(x(5))-wy*cos(x(5)))./x(7));
% psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% Vg = sqrt(x(7).^2 + wx^2 + wy^2 + 2*wx*x(7).*cos(x(5)) + 2*wy*x(7).*sin(x(5)));




end

function [xdot, ydot] = fun_vf(r,gama,Vgd)
global k2 rd
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
function [kaid, kaidot,Vad] = fun_propkaid(r,gamma,Vgd)
global k2 rd wx wy
for i = 1:length(r)
    if r(i,:) < rd

        kaid(i,:) = gamma(i,:) + asin(1./(1+k2*(r(i,:)-rd).^2)) ;
        gamma_d(i,:)= (Vgd ./r(i,:)).*sin(kaid(i,:) - gamma(i,:)) ;
        kaidot(i,:) = gamma_d(i,:) - (2*k2*Vgd .*(r(i,:)-rd))./((1+k2*(r(i,:)-rd).^2).^2);
        Vad(i,:) = sqrt(Vgd ^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid(i,:)) + Vgd*wy*sin(kaid(i,:))));
    else

        kaid(i,:) = gamma(i,:) + pi - asin(1./(1+k2*(r(i,:)-rd).^2)) ;
        gamma_d(i,:) = (Vgd./r(i,:)).*sin(kaid(i,:) - gamma(i,:)) ;
        kaidot(i,:) = gamma_d(i,:) - (2*k2*Vgd.*(r(i,:)-rd))./((1+k2*(r(i,:)-rd).^2).^2);
        Vad(i,:) = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid(i,:)) + Vgd*wy*sin(kaid(i,:))));
    end

end
end