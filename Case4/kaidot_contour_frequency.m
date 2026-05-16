%-------------------------------------------------------------------------%
%------               12th May 2024                         --------------%
%-------     Proposed vector field for sinusoidal path following  --------%
%-------------------------------------------------------------------------%

close all;clear all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global  Vgd  omega A


Vgd = 13;
% omega_arr = [0.0008 0.002 0.008 0.02];
% A = 5*ones(length(omega_arr)); 
A = linspace(0,100,13);
omega_arr = linspace(0,0.02,14);
[X,Y] = meshgrid(A,omega_arr);
Z = X.*Y.^2.*Vgd ;

% Va0 = 10;
% Va0 = 8;
% for i = 1:length(omega_arr)
% omega = omega_arr(i); 
% kaip_dot = A(i)*omega^2*Vgd ;

% Plot format control variables
    lw = 3;            % Line width
    ms = 6;            % Marker size
    ax_fnt = 16;        % Axis font size
    lbl_fnt = 18; % Label font size
    leg_fnt = 15; % Legend font size
    ax_lw = 3;        % Axis line width
%---------    
figure(1)
% surf(X,Y,Z); hold on;grid on;
contour(X,Y,Z,'LineWidth',lw); hold on;grid on;
% colorbar;
% contour(A,omega_arr,kaip_dot); hold on;
% plot(t,path_error,'Color',color(i,:),'LineWidth',lw,'DisplayName',['$$k_{\mathrm{g}} = \ $$', num2str(k(i))]);hold on;grid on;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',ax_lw) ;
xlabel(ax1,' Amplitude m','Fontsize',lbl_fnt);
ylabel(ax1, ' Frequency $$f $$, cycles/m','Fontsize',lbl_fnt);
% 
% end