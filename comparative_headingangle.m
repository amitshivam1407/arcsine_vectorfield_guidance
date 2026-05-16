%-------------------------------------------------------------------------%
%------------                    18th June 2023            ---------------%
%----------------- Arcsine versus ArcTan heading angle    ----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

x = -100:.1:100;
% k_s = [0.0002 0.002 0.02 0.2];
% k_s = [ 0.002 0.02 ];
k_s = 0.002 ;
% kappa = zeros(size(x));
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
% f6 = figure;
% ax6 = axes;
kappa = zeros(size(x));
gama = zeros(size(x));
kaid = zeros(size(x));
a = (pi/2).*ones(size(x));
b = pi.*ones(size(x));
c = 0.*ones(size(x));
q = 2;
for i = 1:length(k_s)
    k1 = k_s(:,i);
    k2(:,i) = sqrt(k_s(:,i)^2 *x(:,i)^2 + 2*k_s(:,i));
    xini_Nelson(:,i) = tan(0.01*pi)/k2(:,i) ;
    xini_prop(:,i) = sqrt((1/k1)*(1/sin(0.02*pi) - 1));
    for j = 1:length(x)
        kaid_Nelson(j,i) = pi/2 + atan(k2(:,i)*x(:,j));
        kappa_Nelson(j,i) = -((k2(:,i).*x(:,j))./(1 + k2(:,i).^2 *x(:,j).^2).^(3/2));
    if x(:,j) < 0
        kappa = -((2*k1.*x(:,j))./(1 + k1.*x(:,j).^2).^2);
        Curvature_prop(j,i) = kappa;
        gama = 1./(1 + k1.*x(:,j).^2);
        gamma_d(j,i) = gama ;
        kaid = asin(gama);
        kaid_prop(j,i) = kaid ;
        kaid_mod(j,i) = pi/2-asin(k1.*abs(x(:,j)).^q./(1 + k1.*abs(x(:,j)).^q));
%         epsilon =  asin(gama);
%         epsilon_prop(j,i) = epsilon;
        
     else
        kappa = -((2*k1.*x(:,j))./(1 + k1.*x(:,j).^2).^2);
        Curvature_prop(j,i) = kappa;
        gama = 1./(1 + k1.*x(:,j).^2);
        gamma_d(j,i) = gama ;
        kaid = (pi - asin(gama));
        kaid_prop(j,i) = kaid ;
        kaid_mod(j,i) = pi/2 + asin(k1.*abs(x(:,j)).^q./(1 + k1.*abs(x(:,j)).^q));
%         epsilon =  asin(gama);
%         epsilon_prop(j,i) = epsilon;
    end
    end
%     x0 = [-50 50]
%     %----Nelson et al kaid for x0
%     kai_cal_Nelson = rad2deg(pi/2 + atan(k2*x0)) ;
%     %----- Proposed kaid for x0
%     if x0 <0
%     kai_cal_prop = rad2deg(asin(1/(1 + k1*x0^2)));
%     else
%      kai_cal_prop = rad2deg(pi - asin(1/(1 + k1*x0^2)));   
%     end

    % Plot format control variables
    lw = 2;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 16;        % Axis font size
    lbl_fnt = ax_fnt+1; % Label font size
    leg_fnt = ax_fnt-1; % Legend font size
    ax_lw = 1.5;        % Axis line width
% figure
plot(ax1,x,Curvature_prop(:,i),'linewidth',2);hold(ax1, 'on');grid(ax1, 'on')
ax1.FontSize = 16;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3) ;
xlabel(ax1,'Cross-track error, m','interpreter','latex','Fontsize',17)
ylabel(ax1,'Curvature,  m $^{-1}$ ','interpreter','latex','Fontsize',17)
% title(ax1,'\boldmath $ Curvature \ Profile $','interpreter','latex')
axis(ax1,[-100 100 -0.5 0.5])
title(ax1, 'Proposed vector field guidance','Fontsize',20)

plot(ax2,x,kappa_Nelson(:,i),'linewidth',2);hold(ax2, 'on');grid(ax2, 'on')
ax2.FontSize = 16;
box on                      % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3) ;
xlabel(ax2,'Cross-track error, m','interpreter','latex','Fontsize',17)
ylabel(ax2,'Curvature,  m $^{-1}$ ','interpreter','latex','Fontsize',17)
title(ax2, 'Nelson et al. vector field guidance','Fontsize',20)
axis(ax2,[-100 100 -0.5 0.5])

plot(ax3,x,kaid_prop(:,i).*(180/pi),'linewidth',3);hold(ax3, 'on');grid(ax3, 'on')
h31 = xline(ax3,-50,'--b','linewidth',3);hold(ax3, 'on');grid(ax3, 'on');
h31.Annotation.LegendInformation.IconDisplayStyle = 'off';
h32 = xline(ax3,50,'--b','linewidth',3);hold(ax3, 'on');grid(ax3, 'on');
h32.Annotation.LegendInformation.IconDisplayStyle = 'off';
h33 = yline(ax3,1.123,'--g','linewidth',3);hold(ax3, 'on');grid(ax3, 'on');
h33.Annotation.LegendInformation.IconDisplayStyle = 'off';
h34 = yline(ax3,178.76,'--g','linewidth',3);hold(ax3, 'on');grid(ax3, 'on');
h34.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax3.FontSize = 16;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3) ;
yticks(ax3,[0 45 90 135 180]);
yticklabels(ax3,{'0','45','90','135','180'});
xlabel(ax3,'Cross-track error, m','interpreter','latex','Fontsize',20)
ylabel(ax3,' $ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',20)
title(ax3, 'Proposed vector field guidance','Fontsize',20)
axis(ax3,[-100 100 -10 190])
x03 = [-50 50 50 -50]; y03 = [1.123 1.123 178.76 178.76] ;c3 = 1;
h35 = fill(ax3,x03,y03,'yellow','FaceAlpha',0.3);hold(ax4, 'on');grid(ax4, 'on');
h35.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(ax4,x,kaid_Nelson(:,i).*(180/pi),'linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
h1 = xline(ax4,-50,'--b','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h11 = xline(ax4,50,'--b','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12 = yline(ax4,0.568,'--g','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';
h13 = yline(ax4,179.3,'--g','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax4.FontSize = 17;
box on                        % Switch on the box around the axis
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',3) ;
yticks(ax4,[0 45 90 135 180]);
yticklabels(ax4,{'0','45','90','135','180'});
xlabel(ax4,'Cross-track error, m','interpreter','latex','Fontsize',20)
ylabel(ax4,' $ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',20)
title(ax4, 'Nelson et al vector field guidance','Fontsize',20)
axis(ax4,[-100 100 -10 190])
x0 = [-50 50 50 -50]; y0 = [0.57 0.57 179.3 179.3] ;c = 1;
h14 = fill(ax4,x0,y0,'cyan','FaceAlpha',0.3);hold(ax4, 'on');grid(ax4, 'on');
h14.Annotation.LegendInformation.IconDisplayStyle = 'off';


plot(ax5,x,kaid_prop(:,i).*(180/pi),'linewidth',3);hold(ax5, 'on');grid(ax5, 'on')
h51 = xline(ax5,-50,'--b','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h51.Annotation.LegendInformation.IconDisplayStyle = 'off';
h52 = xline(ax5,50,'--b','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h52.Annotation.LegendInformation.IconDisplayStyle = 'off';
h53 = yline(ax5,1.123,'--g','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h53.Annotation.LegendInformation.IconDisplayStyle = 'off';
h54 = yline(ax5,178.76,'--g','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h54.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(ax5,x,kaid_Nelson(:,i).*(180/pi),'linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
plot(ax5,x,kaid_mod(:,i).*(180/pi),'linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
% h1 = xline(ax4,-50,'--b','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
% h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h11 = xline(ax4,50,'--b','linewidth',3);hold(ax4, 'on');grid(ax4, 'on');
% h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
h55 = yline(ax5,0.568,'--g','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h55.Annotation.LegendInformation.IconDisplayStyle = 'off';
h56 = yline(ax5,179.3,'--g','linewidth',3);hold(ax5, 'on');grid(ax5, 'on');
h56.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax5.FontSize = 17;
box on                        % Switch on the box around the axis
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',3) ;
yticks(ax5,[0 45 90 135 180]);
yticklabels(ax5,{'0','45','90','135','180'});
xlabel(ax5,'Cross-track error, m','interpreter','latex','Fontsize',19)
ylabel(ax5,' $ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',19)
legend(ax5,'Proposed','Nelson et al.','Fontsize',16)
x0 = [-50 50 50 -50]; y0 = [0.57 0.57 179.3 179.3] ;c = 1;
h57 = fill(ax5,x0,y0,'cyan','FaceAlpha',0.3);hold(ax5, 'on');grid(ax5, 'on');
h57.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
legendStr = [] ;
for i = 1:length(k_s)
    legendStr = [legendStr; {['$ k_{\mathrm{s}} $ = ', num2str(k_s(i))]}];
end
legend(ax1 , legendStr ,'interpreter','latex','location','SouthEast','Fontsize',15);hold on;
legendStr = [] ;
for i = 1:length(k_s)
    legendStr = [legendStr; {['$ \hat{k}_{\mathrm{s}} $ = ', num2str(k2(i))]}];  
end
legend(ax2 , legendStr ,'interpreter','latex','location','SouthEast','Fontsize',15);hold on;
legendStr = [] ;
for i = 1:length(k_s)
    legendStr = [legendStr; {['$ k_{\mathrm{s}} $ = ', num2str(k_s(i))]}];
end
legend(ax3 , legendStr ,'interpreter','latex','location','SouthEast','Fontsize',15);hold on;
legendStr = [] ;
for i = 1:length(k_s)
    legendStr = [legendStr; {['$ \hat{k}_{\mathrm{s}} $ = ', num2str(k2(i))]}];
end
legend(ax4 , legendStr ,'interpreter','latex','location','SouthEast','Fontsize',15);hold on;
