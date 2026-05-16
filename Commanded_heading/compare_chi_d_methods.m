%-------------------------------------------------------------------------%
%-------------------------   Comparison of chi_d methods   ---------------%
%-------------------------------------------------------------------------%
close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

x = -100:1:100;
Cs_values = [0.1, 0.2, 0.3 0.5 0.7];

lw = 3;
ax_fnt = 17;
lbl_fnt = 20;
leg_fnt = 16;

f1 = figure('Position', [100, 100, 800, 600]);
%% -------------------- Color definitions --------------------
blue      = [0 0.4470 0.7410];
red       = [0.8500 0.3250 0.0980];
green     = [0.4660 0.6740 0.1880];
orange    = [0.9290 0.6940 0.1250];
black     = [0 0 0];
colors = {blue, red, green ,orange,black};

for idx = 1:length(Cs_values)
    Cs = Cs_values(idx);
    
    ks = zeros(size(x));
    chi_d_arcsine = zeros(size(x));
    chi_d_other = zeros(size(x));
    
    for j = 1:length(x)
        x0 = x(j);
           
        
        if x0 < 0
            % ks(j) = (1/x0^2) * (1/sqrt(1 - (1 - exp(Cs*x0))^2) - 1);
            ks(j) = Cs^2/2 ;
            gama = 1./(1 + ks(j).*x0.^2);
            chi_d_arcsine(j) = asin(gama);
            chi_d_other(j) = acos(1- exp(Cs * x0));
        else
            % ks(j) = (1/x0^2) * (1/sqrt(1 - (1 - exp(-Cs*x0))^2) - 1);
            ks(j) = Cs^2/2 ;
            gama = 1./(1 + ks(j).*x0.^2);
            chi_d_arcsine(j) = pi - asin(gama);
            chi_d_other(j) = pi - acos(1- exp(-Cs * x0));
        end
        
    end
    
    plot(x, chi_d_arcsine.*(180/pi), '-', 'Color', colors{idx}, 'linewidth', lw); hold on;
    plot(x, chi_d_other.*(180/pi), '--', 'Color', colors{idx}, 'linewidth', lw);
end

grid on;
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax,'linewidth',3);

yticks([0 45 90 135 180]);
yticklabels({'0','45','90','135','180'});

xlabel('Cross-track error, m','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',lbl_fnt);

legendStr = {};
for i = 1:length(Cs_values)
    legendStr{end+1} = ['Arcsine, $ C_{\mathrm{s}} $ = ', num2str(Cs_values(i))];
    legendStr{end+1} = ['Arccos, $ C_{\mathrm{s}} $ = ', num2str(Cs_values(i))];
end
legend(legendStr, 'interpreter','latex','location','best','Fontsize',leg_fnt);
