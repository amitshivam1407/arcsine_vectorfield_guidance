%% Compare arcsin vs arccos guidance laws from same initial lateral error x0
% - Uses your kinematics: xdot = V*cos(chi_d)
% - Integrates x(t) forward and computes chi_d(t), chi_d_dot(t), curvature kappa(t)
% - Two ways to "match" the laws:
%    (A) Slope-match near x=0:     ks = Cs^2/2
%    (B) Match at x0 exactly:      ks chosen so chi_sin(x0)=chi_cos(x0)
%
% Author: (generated)
clear; clc; close all;

%% User settings
V   = 15;          % speed (m/s)
x0  = 40;          % initial lateral error (m)  (choose + or -)
Cs  = 0.3;         % arccos shaping parameter
T   = 30;          % simulation time (s)
dt  = 0.01;        % time step (s)

matchMode = "slope";  % "slope" or "x0"

%% Compute ks from Cs
switch lower(matchMode)
    case "slope"
        ks = Cs^2/2;
    case "x0"
        ks = ks_match_at_x0(Cs, x0);
    otherwise
        error('matchMode must be "slope" or "x0".');
end

fprintf('Match mode: %s\nCs = %.4f, ks = %.6f, x0 = %.3f\n', matchMode, Cs, ks, x0);

%% Time vector
t = 0:dt:T;
N = numel(t);

%% Allocate
x_sin  = zeros(1,N);  x_cos  = zeros(1,N);
chi_sin = zeros(1,N); chi_cos = zeros(1,N);
chidot_sin = zeros(1,N); chidot_cos = zeros(1,N);
kappa_sin = zeros(1,N);  kappa_cos = zeros(1,N);

x_sin(1) = x0;
x_cos(1) = x0;

%% Simulate
for k = 1:N-1
    % --- ARCSIN law ---
    chi_sin(k) = chi_arcsin(x_sin(k), ks);
    dchidx_sin = dchi_dx_arcsin(x_sin(k), ks);
    xdot_sin   = V*cos(chi_sin(k));
    x_sin(k+1) = x_sin(k) + dt*xdot_sin;

    chidot_sin(k) = dchidx_sin * xdot_sin;
    kappa_sin(k)  = chidot_sin(k)/V;

    % --- ARCCOS law ---
    chi_cos(k) = chi_arccos(x_cos(k), Cs);
    dchidx_cos = dchi_dx_arccos(x_cos(k), Cs);
    xdot_cos   = V*cos(chi_cos(k));
    x_cos(k+1) = x_cos(k) + dt*xdot_cos;

    chidot_cos(k) = dchidx_cos * xdot_cos;
    kappa_cos(k)  = chidot_cos(k)/V;
end

% fill last samples
chi_sin(end) = chi_arcsin(x_sin(end), ks);
chi_cos(end) = chi_arccos(x_cos(end), Cs);

dchidx_sin_end = dchi_dx_arcsin(x_sin(end), ks);
dchidx_cos_end = dchi_dx_arccos(x_cos(end), Cs);

xdot_sin_end = V*cos(chi_sin(end));
xdot_cos_end = V*cos(chi_cos(end));

chidot_sin(end) = dchidx_sin_end * xdot_sin_end;
chidot_cos(end) = dchidx_cos_end * xdot_cos_end;

kappa_sin(end) = chidot_sin(end)/V;
kappa_cos(end) = chidot_cos(end)/V;

%% Report peak curvature (magnitude)
[ksin_max, idxS] = max(abs(kappa_sin));
[kcos_max, idxC] = max(abs(kappa_cos));

fprintf('\nPeak |kappa| arcsin  = %.6f 1/m at t=%.3f s, x=%.3f m\n', ksin_max, t(idxS), x_sin(idxS));
fprintf('Peak |kappa| arccos  = %.6f 1/m at t=%.3f s, x=%.3f m\n', kcos_max, t(idxC), x_cos(idxC));

%% Plots
figure; 
plot(t, x_sin, 'LineWidth', 1.5); hold on;
plot(t, x_cos, 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('x (m)');
legend('arcsin','arccos','Location','best');
title(sprintf('Lateral error x(t), x0=%.2f, Cs=%.3f, ks=%.4f (%s-match)', x0, Cs, ks, matchMode));

figure;
plot(t, rad2deg(chi_sin), 'LineWidth', 1.5); hold on;
plot(t, rad2deg(chi_cos), 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('\chi_d (deg)');
legend('arcsin','arccos','Location','best');
title('Commanded course angle \chi_d(t)');

figure;
plot(t, chidot_sin, 'LineWidth', 1.5); hold on;
plot(t, chidot_cos, 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('\dot{\chi}_d (rad/s)');
legend('arcsin','arccos','Location','best');
title('Commanded course rate \dot{\chi}_d(t)');

figure;
plot(t, kappa_sin, 'LineWidth', 1.5); hold on;
plot(t, kappa_cos, 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('\kappa = \dot{\chi}_d/V (1/m)');
legend('arcsin','arccos','Location','best');
title('Curvature \kappa(t)');

figure;
plot(x_sin, rad2deg(chi_sin), 'LineWidth', 1.5); hold on;
plot(x_cos, rad2deg(chi_cos), 'LineWidth', 1.5);
grid on; xlabel('x (m)'); ylabel('\chi_d (deg)');
legend('arcsin','arccos','Location','best');
title('Guidance map: \chi_d vs x along simulated trajectories');

%% ===== Local functions =====
function chi = chi_arcsin(x, ks)
    s = 1/(1 + ks*x^2);
    s = min(max(s,0),1); % safety
    if x <= 0
        chi = asin(s);
    else
        chi = pi - asin(s);
    end
end

function dchidx = dchi_dx_arcsin(x, ks)
    % derivative of piecewise arcsin law w.r.t x
    s = 1/(1 + ks*x^2);
    s = min(max(s,0),1);
    dsdx = -(2*ks*x)/(1 + ks*x^2)^2;
    % chi = asin(s) or pi-asin(s) => dchi/dx = (+/-) (1/sqrt(1-s^2))*ds/dx
    denom = sqrt(max(1 - s^2, 1e-12));
    if x <= 0
        dchidx = (1/denom)*dsdx;
    else
        dchidx = -(1/denom)*dsdx;
    end
end

function chi = chi_arccos(x, Cs)
    if x <= 0
        u = 1 - exp(Cs*x);
        u = min(max(u, -1), 1);
        chi = acos(u);
    else
        u = 1 - exp(-Cs*x);
        u = min(max(u, -1), 1);
        chi = pi - acos(u);
    end
end

function dchidx = dchi_dx_arccos(x, Cs)
    % derivative of piecewise arccos law w.r.t x
    if x <= 0
        u = 1 - exp(Cs*x);
        du = -Cs*exp(Cs*x);
        denom = sqrt(max(1 - u^2, 1e-12));
        dchidx = -(1/denom)*du;  % d/dx acos(u) = -u'/sqrt(1-u^2)
    else
        u = 1 - exp(-Cs*x);
        du = Cs*exp(-Cs*x);
        denom = sqrt(max(1 - u^2, 1e-12));
        dchidx = +(1/denom)*du;  % d/dx [pi - acos(u)] = +u'/sqrt(1-u^2)
    end
end

function ks = ks_match_at_x0(Cs, x0)
    % Solve for ks such that chi_arcsin(x0,ks) == chi_arccos(x0,Cs)
    % Closed-form from equating the principal angles at x0:
    % 1/(1+ks*x0^2) = sqrt(1-(1-exp(-Cs*|x0|))^2)
    xm = abs(x0);
    a  = exp(-Cs*xm);
    rhs = sqrt(max(1 - (1 - a)^2, 1e-12)); % sqrt(1-u^2) with u=1-a
    ks = (1/(rhs) - 1)/(xm^2);
end
