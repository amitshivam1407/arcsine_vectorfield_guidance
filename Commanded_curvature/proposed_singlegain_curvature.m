%-------------------------------------------------------------------------%
%-------------------------         23rd June 2022       ------------------%
%---------------------   Proposed curvature for straight line   ----------%
%-------------------------------------------------------------------------%
close all;clear all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
f1 = figure;
ax = axes;
x = -100:0.01:100;
% k_s = [0.0002 0.002 0.02 0.2];
k_s = 0.002;


for i = 1:length(k_s)
    %     k1 = k_s(:,i);
    for j = 1:length(x)
        if x(:,j) < 0
            kappa = -((2*k_s(:,i).*x(:,j))./(1 + k_s(:,i).*x(:,j).^2).^2);
            Curvature_prop(j,i) = kappa;
            [val_max(i) ,ind_max(i) ] = max(Curvature_prop(:,i));
            x_max(i)  = x(ind_max(i) );
        else
            kappa = -((2*k_s(:,i).*x(:,j))./(1 + k_s(:,i).*x(:,j).^2).^2);
            Curvature_prop(j,i) = kappa;
            [val_min(i) , ind_min(i) ] = min(Curvature_prop(:,i));
            x_min(i)  = x(ind_min(i) );
        end



    end

    str_min = '$ -\frac{1}{\sqrt {3 k_{\mathrm{s}}}}$';
    str_max = '$ \frac{1}{\sqrt {3 k_{\mathrm{s}}}}$';
    kappa_max = '$ \frac{9}{8}\sqrt{\frac{k_{\mathrm{s}}}{3}}$';
    kappa_min = '$ -\frac{9}{8}\sqrt{\frac{k_{\mathrm{s}}}{3}}$';

    f1;
    
    if i ==1
        h2(i) = plot(ax,x_max(i),val_max(i),'-s','linewidth',2,'MarkerSize',10,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','cyan','DisplayName','Maxima '); hold(ax, 'on');

        h3(i) = plot(ax,x_min(i),val_min(i),'-s','linewidth',2,'MarkerSize',10,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','green','DisplayName','Minima '); hold(ax, 'on');
    else
        h2(i) = plot(ax,x_max(i),val_max(i),'-s','linewidth',2,'MarkerSize',10,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','cyan'); hold(ax, 'on');
        h2(i).Annotation.LegendInformation.IconDisplayStyle = 'off';

        h3(i) = plot(ax,x_min(i),val_min(i),'-s','linewidth',2,'MarkerSize',10,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','green'); hold(ax, 'on');
        h3(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    h1 = plot(ax,x',Curvature_prop(:,i),'r','linewidth',3,'DisplayName',['$$k_{\mathrm{s}} = \ $$', num2str(k_s(i))]);hold(ax, 'on');grid(ax, 'on');
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % plot(x,Curvature_prop(:,i),'linewidth',2);hold 'on';grid 'on';
    % h2 = plot(x_max,0,'-o','MarkerSize',8,...
    %     'MarkerEdgeColor','black',...
    %     'MarkerFaceColor','red'); hold on;
    % h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % plot(x_max,0,'-o','MarkerSize',8,...
    %     'MarkerEdgeColor','black',...
    %     'MarkerFaceColor','red'); hold on;
    % plot(x_min,0,'-o','MarkerSize',8,...
    %     'MarkerEdgeColor','black',...
    %     'MarkerFaceColor','green'); hold on;
    % yline(0,'--k','linewidth',2);hold 'on';grid 'on';
    h4 = plot(ax,[x_max(i)  x_max(i) ],[0 val_max(i) ],'--','color',[0.75 0.75 0.75],'linewidth',2);hold 'on';grid 'on';
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    text(x_max-19.09,-0.005,str_min,'interpreter','latex','Fontsize',18); hold on;
    text(x_max-11.09,val_max+0.008,kappa_max,'interpreter','latex','Fontsize',18); hold on;
    h5 = plot(ax,[x_min(i)  x_min(i) ],[0 val_min(i) ],'--','color',[0.75 0.75 0.75],'linewidth',2);hold 'on';grid 'on';
    h5.Annotation.LegendInformation.IconDisplayStyle = 'off';
    text(x_min-10.09,-0.005,str_max,'interpreter','latex','Fontsize',18); hold on;
    text(x_min-14.09,val_min-0.007,kappa_min,'interpreter','latex','Fontsize',18); hold on;
    ax.FontSize = 17;
    box on                      % Switch on the box around the axis
    ax.XColor = 'black';         % Box horizontal lines' color
    ax.YColor = 'black';         % Box vertical lines' color
    set(ax,'linewidth',3) ;
    xlabel(ax,'Cross-track error $$x$$, m','Fontsize',19)
    ylabel(ax,'Curvature $$ \kappa, $$  m$^{-1} $','Fontsize',19);
    %     legend(ax,'Maximum','Minimum','Fontsize',15)
    axis(ax,[-100 100 -0.045 0.045])
end
legend(ax)
