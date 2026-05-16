%% demo_straight_line_following.m
% Straight line following using cos^-1 course angle law with unicycle kinematics
% Based on equation (4): chi_d = cos^-1(1 - e^(Cs*x)) for x <= 0
%                               = pi - cos^-1(1 - e^(Cs*x)) otherwise

clear; clc; close all;

% Set LaTeX interpreter for better plotting
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Plot parameters
ax_fnt = 25;
lbl_fnt = 29;
ax_wdth = 3;
lgd_fnt = 19;

%% -------------------- Color definitions --------------------
blue      = [0 0.4470 0.7410];
red       = [0.8500 0.3250 0.0980];
green     = [0.4660 0.6740 0.1880];
orange    = [0.9290 0.6940 0.1250];
black     = [0 0 0];

%% -------------------- Simulation parameters --------------------
Tend = 30;      % Simulation time
dt   = 0.01;    % Time step
t    = 0:dt:Tend;
N    = numel(t);

%% -------------------- Unicycle parameters --------------------
v = 2.0;        % Constant forward velocity (m/s)

%% -------------------- Course angle law parameters --------------------
% Test different Cs values as shown in the figure
Cs_values = [0.1, 0.3, 0.5, 0.7];
n_Cs = length(Cs_values);
Cs_colors = [blue; red; green; orange];

%% -------------------- Initial conditions --------------------
% Different initial lateral errors
x0_values = [-100, -50, 50];  % Initial lateral errors (m)
n_IC = length(x0_values);
IC_colors = [blue; red; green];

%% -------------------- Desired course angle function (Equation 4) --------------------
chi_d = @(x, Cs) (x <= 0) .* acos(1 - exp(Cs*x)) + ...
                 (x > 0)  .* (pi - acos(1 - exp(-Cs*x)));

%% ============================================================
%% Simulation 1: Different Cs values with fixed initial condition
%% ============================================================
fprintf('\n=== Simulation 1: Effect of Cs parameter ===\n');

x0_fixed = -100;  % Fixed initial lateral error
y0 = 0;           % Initial y position
chi0 = pi/2;      % Initial course angle (pointing along y-axis)

% Storage for different Cs values
X_Cs = zeros(n_Cs, N);
Y_Cs = zeros(n_Cs, N);
Chi_Cs = zeros(n_Cs, N);

for i_Cs = 1:n_Cs
    Cs = Cs_values(i_Cs);
    fprintf('Simulating with Cs = %.1f\n', Cs);
    
    % Initialize state
    x = zeros(1, N);
    y = zeros(1, N);
    chi = zeros(1, N);
    
    x(1) = x0_fixed;
    y(1) = y0;
    chi(1) = chi0;
    
    % Simulation loop
    for k = 1:N-1
        % Desired course angle from equation (4)
        chi_desired = chi_d(x(k), Cs);
        
        % Control law: simple proportional control for course angle
        k_chi = 2.0;  % Course angle gain
        omega = -k_chi * (chi(k) - chi_desired);
        
        % Unicycle kinematics
        % dx/dt = v*sin(chi)
        % dy/dt = v*cos(chi)
        % dchi/dt = omega
        x(k+1) = x(k) + dt * v * sin(chi(k));
        y(k+1) = y(k) + dt * v * cos(chi(k));
        chi(k+1) = chi(k) + dt * omega;
        
        % Wrap angle to [-pi, pi]
        chi(k+1) = atan2(sin(chi(k+1)), cos(chi(k+1)));
    end
    
    % Store results
    X_Cs(i_Cs, :) = x;
    Y_Cs(i_Cs, :) = y;
    Chi_Cs(i_Cs, :) = chi;
    
    fprintf('  Final lateral error: %.4f m\n', x(end));
    fprintf('  Final course angle: %.4f rad (%.2f deg)\n', chi(end), rad2deg(chi(end)));
end

%% ============================================================
%% Simulation 2: Different initial conditions with fixed Cs
%% ============================================================
fprintf('\n=== Simulation 2: Different initial conditions ===\n');

Cs_fixed = 0.3;  % Fixed Cs value

% Storage for different initial conditions
X_IC = zeros(n_IC, N);
Y_IC = zeros(n_IC, N);
Chi_IC = zeros(n_IC, N);

for i_IC = 1:n_IC
    x0 = x0_values(i_IC);
    fprintf('Simulating with x0 = %.1f m\n', x0);
    
    % Initialize state
    x = zeros(1, N);
    y = zeros(1, N);
    chi = zeros(1, N);
    
    x(1) = x0;
    y(1) = y0;
    chi(1) = chi0;
    
    % Simulation loop
    for k = 1:N-1
        % Desired course angle from equation (4)
        chi_desired = chi_d(x(k), Cs_fixed);
        
        % Control law
        k_chi = 2.0;
        omega = -k_chi * (chi(k) - chi_desired);
        
        % Unicycle kinematics
        x(k+1) = x(k) + dt * v * sin(chi(k));
        y(k+1) = y(k) + dt * v * cos(chi(k));
        chi(k+1) = chi(k) + dt * omega;
        
        % Wrap angle
        chi(k+1) = atan2(sin(chi(k+1)), cos(chi(k+1)));
    end
    
    % Store results
    X_IC(i_IC, :) = x;
    Y_IC(i_IC, :) = y;
    Chi_IC(i_IC, :) = chi;
    
    fprintf('  Final lateral error: %.4f m\n', x(end));
    fprintf('  Final course angle: %.4f rad (%.2f deg)\n', chi(end), rad2deg(chi(end)));
end

%% ============================================================
%% Plotting: Effect of Cs parameter
%% ============================================================

% Figure 1: Lateral error vs y position for different Cs
figure(1);
hold on; grid on;
for i_Cs = 1:n_Cs
    plot(Y_Cs(i_Cs, :), X_Cs(i_Cs, :), '-', 'color', Cs_colors(i_Cs, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$C_s = %.1f$', Cs_values(i_Cs)));
end
plot([0, max(Y_Cs(:))], [0, 0], 'k--', 'linewidth', 2, 'DisplayName', 'Desired path');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ y, $ m', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);
axis equal;

% Figure 2: Course angle vs lateral error for different Cs
figure(2);
hold on; grid on;
for i_Cs = 1:n_Cs
    plot(X_Cs(i_Cs, :), rad2deg(Chi_Cs(i_Cs, :)), '-', 'color', Cs_colors(i_Cs, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$C_s = %.1f$', Cs_values(i_Cs)));
end
yline(90, 'k--', 'linewidth', 2, 'DisplayName', '$\chi_d = \pi/2$');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ \chi, $ deg', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);

% Figure 3: Lateral error vs time for different Cs
figure(3);
hold on; grid on;
for i_Cs = 1:n_Cs
    plot(t, X_Cs(i_Cs, :), '-', 'color', Cs_colors(i_Cs, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$C_s = %.1f$', Cs_values(i_Cs)));
end
plot([0, Tend], [0, 0], 'k--', 'linewidth', 2, 'HandleVisibility', 'off');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ t, $ s', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);

%% ============================================================
%% Plotting: Different initial conditions
%% ============================================================

% Figure 4: Lateral error vs y position for different ICs
figure(4);
hold on; grid on;
for i_IC = 1:n_IC
    plot(Y_IC(i_IC, :), X_IC(i_IC, :), '-', 'color', IC_colors(i_IC, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$x_0 = %.0f$ m', x0_values(i_IC)));
end
plot([0, max(Y_IC(:))], [0, 0], 'k--', 'linewidth', 2, 'DisplayName', 'Desired path');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ y, $ m', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);
axis equal;

% Figure 5: Course angle vs lateral error for different ICs
figure(5);
hold on; grid on;
for i_IC = 1:n_IC
    plot(X_IC(i_IC, :), rad2deg(Chi_IC(i_IC, :)), '-', 'color', IC_colors(i_IC, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$x_0 = %.0f$ m', x0_values(i_IC)));
end
yline(90, 'k--', 'linewidth', 2, 'DisplayName', '$\chi_d = \pi/2$');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ \chi, $ deg', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);

% Figure 6: Desired course angle function visualization
figure(6);
hold on; grid on;
x_range = linspace(-100, 100, 1000);
for i_Cs = 1:n_Cs
    chi_desired_plot = arrayfun(@(x) chi_d(x, Cs_values(i_Cs)), x_range);
    plot(x_range, rad2deg(chi_desired_plot), '-', 'color', Cs_colors(i_Cs, :), ...
        'linewidth', 3, 'DisplayName', sprintf('$C_s = %.1f$', Cs_values(i_Cs)));
end
yline(90, 'k--', 'linewidth', 2, 'DisplayName', '$\chi_d = \pi/2$');
ax = gca;
ax.FontSize = ax_fnt;
box on;
ax.XColor = 'black';
ax.YColor = 'black';
set(ax, 'linewidth', ax_wdth);
xlabel(ax, ' $ x, $ m', 'Fontsize', lbl_fnt);
ylabel(ax, ' $ \chi_d, $ deg', 'Fontsize', lbl_fnt);
legend(ax, 'Location', 'best', 'Fontsize', lgd_fnt);
grid on;

fprintf('\n=== Simulation Complete ===\n');
