%-------------------------------------------------------------------------%
%------------                4th November 2024          ------------------%
%------------ Variable L0 based straight line following  -----------------%
%-------------------------------------------------------------------------%

close all;clear all; 

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%Parameters
v = 1;%m/s
Li = 0.1;lambda = 1;%L0 profile parameters

Lf_arr = [2 3 4];

for i = 1:length(Lf_arr)
    % for i = 1:2
    Lf = Lf_arr(i);
% Lf = 3 ; 


%path definition
yref = 0;%m

%constraints
cmax=0.5;
cmin=-cmax;

%Simulation Time and Time Step
TT=20;%Simulation Time
NN=1000*TT;%Number of Points
hh=TT/NN; %Time Step

%Figure time vector
tout=0:NN-1;
tout=tout*hh;

%Initial Conditions
x0=0;
y0=1;
theta0=0;

dx0=v*cos(theta0);
dy0=v*sin(theta0);
dtheta0=0;

 

%Simulation Initial State
x=x0;
y=y0;
theta=theta0;

dx=dx0;
dy=dy0;
dtheta=dtheta0;


 for k=1:NN %Starting Simulation

%Attributing values for trajectory vectors
X(:,k)=[x; y; theta];
dX(:,k)=[dx; dy; dtheta];

%Minimum distance and distance vector
 d=sqrt((y-yref)^2);
 d_v(k)=d;
 
 L0_var = (Lf - Li)*exp(-lambda*d) + Li;
L0var(k) = L0_var ;

 %Defining point Q (closest point) and point R (reference point)
 Q=[x; yref];

 
 R=Q+[L0_var; 0];
 

 %Defining the L1 vector
 L1=R-[x; y];

 %Calculating eta (heading angle deviation from L1)
n=[1;0;0];
dp=[0; dx; dy];
L1_L=[0; L1(1); L1(2)];  
aux1=cross(dp,L1_L);
aux2=aux1/(norm(dp)*norm(L1_L));
aux3=dot(dp,L1_L);
aux4=aux3/(norm(dp)*norm(L1_L));
eta=atan2(dot(aux2,n),aux4);
eta_v(k)=eta;

%computing curvature
c=2*sin(eta)/norm(L1);
%saturating c
if c>cmax
    c=cmax;
elseif c<cmin
    c=cmin;
end
c_v(k)=c;%c vector

%Dynamics
dx=v*cos(theta);
dy=v*sin(theta);
dtheta=v*c;
%Integration (Euler) 
x=x+hh*dx;
y=y+hh*dy;
theta=theta+hh*dtheta;

% if y<=0.2
%     t_stop = 
% else

 end

% plot parameter
ax_fnt = 20;
lbl_fnt = 24;
ax_wdth = 3;
lgd_fnt = 17;

%% Colors
        % blue = [0 0.4470 0.7410];
        % red = [0.8500 0.3250 0.0980];
        red = [1 0 0];
        orange = [0.9290 0.6940 0.1250];
        violet = [0.4940 0.1840 0.5560];
        % green = [0.4660 0.6740 0.1880];
        green = [0 1 0];
        blue = [0 0 1];
        cyan = [0.3010 0.7450 0.9330];
        maroon = [0.6350 0.0780 0.1840];
        black = [0 0 0];

        color = [blue;green;cyan];



% figure(1)
% h11 = yline(0,'color',maroon,'LineWidth',3,'DisplayName','Desired path'); hold on;grid on;
% h12 = plot(X(1,:),X(2,:),'-','color',red,'LineWidth',3,'DisplayName',['Varying $$ L_0 $$']) ;hold on;grid on;
% ax1 = gca;
% ax1.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax1.XColor = 'black';         % Box horizontal lines' color
% ax1.YColor = 'black';         % Box vertical lines' color
% set(ax1,'linewidth',3)
% xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
% ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
% % legend(ax1,'Desired path','UAV trajectory','Fontsize',14,'NumColumns',1)
% axis(ax1, 'tight');
% 
% figure(2)
% h2 = plot(tout,d_v,'-','color',red,'LineWidth',3,'DisplayName',['Varying $$ L_0 $$']);hold on;grid on;
% ax2 = gca;
% ax2.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax2.XColor = 'black';         % Box horizontal lines' color
% ax2.YColor = 'black';         % Box vertical lines' color
% set(ax2,'linewidth',3)
% xlabel(ax2,'  $$t$$, s','Fontsize',lbl_fnt);
% ylabel(ax2,' $$y$$, m','Fontsize',lbl_fnt);
% 
% figure(3)
% h32 = plot(tout,wrapToPi(eta_v)*(180/pi),'-','color',red,'LineWidth',3,'DisplayName',['Varying $$ L_0 $$']) ;hold on;grid on;
% ax3 = gca;
% ax3.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax3.XColor = 'black';         % Box horizontal lines' color
% ax3.YColor = 'black';         % Box vertical lines' color
% set(ax3,'linewidth',3)
% xlabel(ax3,'  $$t$$, s','Fontsize',lbl_fnt);
% yticks(ax3,[-120 -90 -45 0 45 90 120])
% yticklabels(ax3,{'-120','-90','-45','0','45','90','120'})
% ylabel(ax3,'  $$\eta $$, deg.','Fontsize',lbl_fnt);
% % legend(ax3,'Commanded','Achieved','Fontsize',lgd_fnt)
% 
% figure(4)
% h4 = plot(tout,c_v,'-','color',red,'LineWidth',3,'DisplayName',['Varying $$ L_0 $$']) ;hold on;grid on;
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',3)
% xlabel(ax4,'  $$t$$, s','Fontsize',lbl_fnt);
% ylabel(ax4,'Lateral acceleration $$a$$, m/s$^2$','Fontsize',lbl_fnt);
% axis(ax4,[0 20 -1 1])
% 
% 
% 
% legend(ax1, 'Fontsize',lgd_fnt ,'NumColumns',1)
% legend(ax2, 'Fontsize',lgd_fnt ,'NumColumns',1)
% legend(ax3, 'Fontsize',lgd_fnt ,'NumColumns',1)
% legend(ax4, 'Fontsize',lgd_fnt ,'NumColumns',1)

figure(1)
if i == 1
% plot(X(1,:), 0*ones(size(tout)),'k','LineWidth',3)
h11 = yline(0,'color',maroon,'LineWidth',3,'DisplayName','Desired path'); hold on;grid on;
h12 = plot(X(1,:),X(2,:),'-','color',color(i,:),'LineWidth',3,'DisplayName',['$$ L_\mathrm{f} = \ $$' num2str(Lf_arr(i)), '\ m']) ;hold on;grid on;
else
 h12 = plot(X(1,:),X(2,:),'-','color',color(i,:),'LineWidth',3,'DisplayName',['$$ L_\mathrm{f} = \ $$' num2str(Lf_arr(i)), '\ m']) ;hold on;grid on;   
end
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3)
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
% legend(ax1,'Desired path','UAV trajectory','Fontsize',14,'NumColumns',1)
axis(ax1, [0 20 -1.5 1.5]);

figure(2)
h2 = plot(tout,X(2,:),'-','color',color(i,:),'LineWidth',3,'DisplayName',['$$ L_\mathrm{f} = \ $$' num2str(Lf_arr(i)), '\ m']) ;hold on;grid on;
% h2.
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3)
xlabel(ax2,'  $$t$$, s','Fontsize',lbl_fnt);
ylabel(ax2,' $$y$$, m','Fontsize',lbl_fnt);

figure(3)
h32 = plot(tout,wrapToPi(eta_v)*(180/pi),'-','color',color(i,:),'linewidth',3,'DisplayName',['$$ L_\mathrm{f} = \ $$' num2str(Lf_arr(i)), '\ m']) ;hold on;grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3)
xlabel(ax3,'  $$t$$, s','Fontsize',lbl_fnt);
yticks(ax3,[-120 -90 -45 0 45 90 120])
yticklabels(ax3,{'-120','-90','-45','0','45','90','120'})
ylabel(ax3,'  $$\eta $$, deg.','Fontsize',lbl_fnt);
% legend(ax3,'Commanded','Achieved','Fontsize',lgd_fnt)
% 
% figure(4)
% h4 = plot(tout,a_v,'--','color',color(i,:),'LineWidth',3,'DisplayName',['$$ L_\mathrm{f} = \ $$' num2str(Lf_arr(i)), '\ m']) ;hold on;grid on;
% ax4 = gca;
% ax4.FontSize = ax_fnt;
% box on                        % Switch on the box around the axis
% ax4.XColor = 'black';         % Box horizontal lines' color
% ax4.YColor = 'black';         % Box vertical lines' color
% set(ax4,'linewidth',3)
% xlabel(ax4,'  $$t$$, s','Fontsize',lbl_fnt);
% ylabel(ax4,' $$\kappa$$, m/s$^{-1}$','Fontsize',lbl_fnt);
% axis(ax4,[0 20 -1 1])


% figure(3)
% subplot(2,1,1)
% plot(tout, 0*ones(size(tout)),'k')
% hold all
% plot(tout, d_v, '--')
% subplot(2,1,2)
% plot(tout, 0*ones(size(tout)),'k')
% hold all
% plot(tout, eta_v, '--')
end

legend(ax1, 'Fontsize',lgd_fnt ,'NumColumns',1)
legend(ax2, 'Fontsize',lgd_fnt ,'NumColumns',1)
legend(ax3, 'Fontsize',lgd_fnt ,'NumColumns',1)
% legend(ax4, 'Fontsize',lgd_fnt ,'NumColumns',1)

function [value,isterminal,direction] = stopping_L0(t,x)
global rd
value(1)      = x(2) - .2; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1)  = 0; % negative direction(as R decreases from positive to zero d=-1;
                  % If R increases from negative to zero d=+1;d=0 implies no need of direction )
end