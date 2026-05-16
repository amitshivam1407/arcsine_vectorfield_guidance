%-------------------------------------------------------------------------%
%------------                    14th August 2024            -------------%
%-------------------------------------------------------------------------%

close all;clear all; clc

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

x = linspace(-10,10,1000);
y = linspace(-10,10,1000);
p = atan2(y,1);
p1 = atan2(sin(y),cos(y));
% p = 2*atan(x);

f1 = figure;
ax1 = axes;

f1; 
h1 = plot(ax1,y,p,'r','linewidth',3);hold(ax1,'on');grid(ax1,'on');
h2 = plot(ax1,y,p1,'--b','linewidth',3);hold(ax1,'on');
% Setup axis font size
% ax1 = gca;     
ax1.FontSize = 17;
axis(ax1,'equal')
% Outer box setup
box on                      % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3)   % Axis linewidth (box and grid)
% Label and limit setup
xlabel(ax1,'$$x$$, m','Interpreter','Latex','FontSize',19)
ylabel(ax1,'$$p$$, rad','Interpreter','Latex','FontSize',19)
% axis(ax1,'equal')