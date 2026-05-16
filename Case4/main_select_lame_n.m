%% main_select_lame_n.m
% Curvature analysis and automatic selection of Lamé exponent n
% for geometry-actuation compatible safety design.

clear; clc; close all;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Parameters
P.a    = 10;
P.b    = 6;
P.v    = 1.0;
P.umax = 0.25;

% Vehicle curvature limit
P.kappa_max = P.umax / P.v;

% Candidate n values
n_values = [2 2.5 3 4 5 6 8 10];

% Boundary sampling
Ntheta = 2000;
theta = linspace(0, 2*pi, Ntheta);

% Critical region selection:
% full boundary, or side-to-corner region only
% Example: first quadrant side-to-corner region
% theta_crit = linspace(0.02, pi/2 - 0.02, 700);
theta_crit = linspace(0.02, pi/2 - 0.02, 700);

%% Analyze candidates
results = struct([]);
for i = 1:length(n_values)
    n = n_values(i);

    R = analyze_lame_curvature(n, theta, theta_crit, P);

    results(i).n                    = n;
    results(i).kappa_all_max        = R.kappa_all_max;
    results(i).kappa_crit_max       = R.kappa_crit_max;
    results(i).kappa_crit_mean      = R.kappa_crit_mean;
    results(i).kappa_side_mean      = R.kappa_side_mean;
    results(i).kappa_corner_max     = R.kappa_corner_max;
    results(i).feasible_all         = R.feasible_all;
    results(i).feasible_critical    = R.feasible_critical;
    results(i).score                = R.score;
    results(i).theta                = R.theta;
    results(i).x                    = R.x;
    results(i).y                    = R.y;
    results(i).kappa                = R.kappa;
    results(i).theta_crit           = R.theta_crit;
    results(i).kappa_crit           = R.kappa_crit;
end

%% Display results table
fprintf('\n================ CURVATURE ANALYSIS ================\n');
fprintf('a = %.3f, b = %.3f, v = %.3f, umax = %.3f, kappa_max = %.5f\n\n', ...
    P.a, P.b, P.v, P.umax, P.kappa_max);

fprintf('%8s %15s %15s %15s %15s %12s %12s %12s\n', ...
    'n', 'kappa_all_max', 'kappa_crit_max', 'kappa_crit_mean', ...
    'kappa_corner_max', 'feas_all', 'feas_crit', 'score');

for i = 1:length(results)
    fprintf('%8.2f %15.6f %15.6f %15.6f %15.6f %12d %12d %12.6f\n', ...
        results(i).n, ...
        results(i).kappa_all_max, ...
        results(i).kappa_crit_max, ...
        results(i).kappa_crit_mean, ...
        results(i).kappa_corner_max, ...
        results(i).feasible_all, ...
        results(i).feasible_critical, ...
        results(i).score);
end

%% Select best n
best_idx = select_best_n(results, P);
best_n   = results(best_idx).n;

fprintf('\nRecommended n = %.2f\n', best_n);

%% Plot curvature vs boundary parameter
plot_curvature_profiles(results, P);

%% Plot boundary shapes
plot_boundary_shapes(results, P);

%% Plot summary metrics vs n
plot_n_selection_metrics(results, P);


function R = analyze_lame_curvature(n, theta, theta_crit, P)
% Analyze Lamé boundary curvature for exponent n

% Full boundary
[x, y] = lame_boundary_param(theta, n, P.a, P.b);
kappa  = numerical_curvature(theta, x, y);

% Critical region boundary
[xc, yc] = lame_boundary_param(theta_crit, n, P.a, P.b);
kappa_crit = numerical_curvature(theta_crit, xc, yc);

% Define subregions inside first quadrant critical set
% Side region: near top side, smaller theta (close to pi/2) or right side (close to 0)
% Here we focus on upper side-to-corner transition in first quadrant.
%
% Let us define:
% side-like region   = theta in [0.65*pi/2, 0.95*pi/2]
% corner-like region = theta in [0.20*pi/2, 0.65*pi/2]
%
% You can change these depending on how you define "critical".
th_side   = theta_crit(theta_crit >= 0.65*(pi/2) & theta_crit <= 0.95*(pi/2));
th_corner = theta_crit(theta_crit >= 0.20*(pi/2) & theta_crit <= 0.65*(pi/2));

if isempty(th_side)
    kappa_side_mean = NaN;
else
    [xs, ys] = lame_boundary_param(th_side, n, P.a, P.b);
    ks = numerical_curvature(th_side, xs, ys);
    kappa_side_mean = mean(ks,'omitnan');
end

if isempty(th_corner)
    kappa_corner_max = NaN;
else
    [xco, yco] = lame_boundary_param(th_corner, n, P.a, P.b);
    kc = numerical_curvature(th_corner, xco, yco);
    kappa_corner_max = max(kc, [], 'omitnan');
end

kappa_all_max   = max(kappa, [], 'omitnan');
kappa_crit_max  = max(kappa_crit, [], 'omitnan');
kappa_crit_mean = mean(kappa_crit, 'omitnan');

% Feasibility with respect to available curvature
feasible_all      = (kappa_all_max  <= P.kappa_max);
feasible_critical = (kappa_crit_max <= P.kappa_max);

% Score design:
% penalize infeasibility heavily, then prefer
% - low critical mean curvature
% - low corner peak
% - larger n only if still feasible (mild reward)
penalty = 0;
if ~feasible_critical
    penalty = penalty + 1e3*(kappa_crit_max - P.kappa_max + 1);
end
if ~feasible_all
    penalty = penalty + 5e2*(kappa_all_max - P.kappa_max + 1);
end

% Mild reward for being more geometry-aligned than ellipse
% But do not let very large n dominate if curvature spikes
reward_n = -0.05*n;

score = penalty + 3.0*kappa_crit_mean + 2.0*kappa_corner_max + reward_n;

R.n                 = n;
R.theta             = theta;
R.x                 = x;
R.y                 = y;
R.kappa             = kappa;
R.theta_crit        = theta_crit;
R.kappa_crit        = kappa_crit;
R.kappa_all_max     = kappa_all_max;
R.kappa_crit_max    = kappa_crit_max;
R.kappa_crit_mean   = kappa_crit_mean;
R.kappa_side_mean   = kappa_side_mean;
R.kappa_corner_max  = kappa_corner_max;
R.feasible_all      = feasible_all;
R.feasible_critical = feasible_critical;
R.score             = score;
end


function [x, y] = lame_boundary_param(theta, n, a, b)

c = cos(theta);
s = sin(theta);

x = a * sign_nonzero(c) .* abs(c).^(2/n);
y = b * sign_nonzero(s) .* abs(s).^(2/n);

end

function z = sign_nonzero(v)
z = zeros(size(v));
z(v > 0)  = 1;
z(v < 0)  = -1;
z(v == 0) = 0;
end

function kappa = numerical_curvature(theta, x, y)

dx  = gradient(x, theta);
dy  = gradient(y, theta);
ddx = gradient(dx, theta);
ddy = gradient(dy, theta);

num = abs(dx .* ddy - dy .* ddx);
den = (dx.^2 + dy.^2).^(3/2);

kappa = num ./ den;

% Clean up numerical issues
kappa(~isfinite(kappa)) = NaN;
end


function best_idx = select_best_n(results, P)

scores = [results.score];
feas_crit = [results.feasible_critical];

% Prefer feasible candidates first
idx_feas = find(feas_crit);

if ~isempty(idx_feas)
    [~, loc] = min(scores(idx_feas));
    best_idx = idx_feas(loc);
else
    % if none feasible in critical region, choose least bad
    [~, best_idx] = min(scores);
end
end

function plot_curvature_profiles(results, P)

figure('Color','w');

subplot(2,1,1); hold on; grid on;
for i = 1:length(results)
    plot(results(i).theta, results(i).kappa, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('$n=%.1f$', results(i).n));
end
yline(P.kappa_max, 'k--', 'LineWidth', 2.0, ...
    'DisplayName', '$\kappa_{\max}=u_{\max}/v$');
xlabel('$\theta$ [rad]');
ylabel('$\kappa(\theta)$');
title('Boundary Curvature Along Full Lam\''e / Elliptic Boundary');
legend('Location','eastoutside');
box on;

subplot(2,1,2); hold on; grid on;
for i = 1:length(results)
    plot(results(i).theta_crit, results(i).kappa_crit, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('$n=%.1f$', results(i).n));
end
yline(P.kappa_max, 'k--', 'LineWidth', 2.0, ...
    'DisplayName', '$\kappa_{\max}$');
xlabel('$\theta$ [rad]');
ylabel('$\kappa(\theta)$');
title('Curvature in Critical Side-to-Corner Region');
legend('Location','eastoutside');
box on;
end

function plot_boundary_shapes(results, P)

figure('Color','w'); hold on; grid on; axis equal;

for i = 1:length(results)
    plot(results(i).x, results(i).y, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('$n=%.1f$', results(i).n));
end

xlabel('$x$');
ylabel('$y$');
title('Lam\''e / Elliptic Boundaries for Candidate $n$');
legend('Location','eastoutside');
box on;
end

function plot_n_selection_metrics(results, P)

nvals            = [results.n];
kappa_all_max    = [results.kappa_all_max];
kappa_crit_max   = [results.kappa_crit_max];
kappa_crit_mean  = [results.kappa_crit_mean];
kappa_corner_max = [results.kappa_corner_max];
score            = [results.score];

figure('Color','w');

subplot(2,2,1);
plot(nvals, kappa_all_max, '-o', 'LineWidth', 2.0, 'DisplayName','Full-boundary max');
hold on; grid on;
plot(nvals, kappa_crit_max, '--s', 'LineWidth', 2.0, 'DisplayName','Critical-region max');
yline(P.kappa_max, 'k--', 'LineWidth', 2.0, 'DisplayName','$\kappa_{\max}$');
xlabel('$n$');
ylabel('Curvature');
title('Maximum Curvature vs $n$');
legend('Location','best');
box on;

subplot(2,2,2);
plot(nvals, kappa_crit_mean, '-o', 'LineWidth', 2.0);
grid on;
xlabel('$n$');
ylabel('Curvature');
title('Mean Critical-Region Curvature vs $n$');
box on;

subplot(2,2,3);
plot(nvals, kappa_corner_max, '-o', 'LineWidth', 2.0);
hold on; grid on;
yline(P.kappa_max, 'k--', 'LineWidth', 2.0);
xlabel('$n$');
ylabel('Curvature');
title('Corner-Region Peak Curvature vs $n$');
box on;

subplot(2,2,4);
plot(nvals, score, '-o', 'LineWidth', 2.0);
grid on;
xlabel('$n$');
ylabel('Score');
title('Selection Score vs $n$');
box on;
end

