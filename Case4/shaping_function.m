%-------------------------------------------------------------------------%
%------------------------     30th March 2024      -----------------------%
%-----------  Sinusoidal path following using autopilot model  -----------%
%-------------------------------------------------------------------------%

close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global kg

% kk = [0.005 0.05 0.5];

kk = [ 0.05 ];
for i = 1: length(kk)
kg = kk(i);

epsilon = -100:100;
kai_o = pi/2 -  asin(1./(1 + kg*epsilon.^2))  ;

% plot parameter
ax_fnt = 17;
lbl_fnt = 19;
ax_wdth = 3;
lgd_fnt = 15;

figure(2)
plot(epsilon,kai_o,'LineWidth',3,'DisplayName',['$$k_{\mathrm{g}} = $$' num2str(kk(i))]);hold on;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',ax_wdth) ;
xlabel(ax2,' Tracking error $\epsilon $, m','Fontsize',lbl_fnt);
ylabel(ax2,' Shaping function $$\chi_{\mathrm{o}}$$, rad.','Fontsize',lbl_fnt);
axis(ax2,[epsilon(1) epsilon(end) 0 2])
end
legend(ax2)
