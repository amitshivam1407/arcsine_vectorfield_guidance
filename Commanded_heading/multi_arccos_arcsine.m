%% Multi-Cs comparison: arcsin (bold) vs arccos (dashed)
clear; clc; close all;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% User parameters
V   = 15;                  % speed (m/s)
x0  = 40;                  % initial lateral error (m)
Cs_vec = [0.08 0.1 0.2 0.3]; % values of Cs to compare
T   = 30;                  % simulation time (s)
dt  = 0.01;

matchMode = "x0";  % "slope" or "x0"

%% Time
t = 0:dt:T;
N = numel(t);

%% Colors and formatting
blue      = [0 0.4470 0.7410];
red       = [0.8500 0.3250 0.0980];
green     = [0.4660 0.6740 0.1880];
orange    = [0.9290 0.6940 0.1250];
cols = {blue, red, green, orange};

lw = 3;
ax_fnt = 17;
lbl_fnt = 20;
leg_fnt = 16;

%% Storage for legends
legText = {};

%% ---------- FIGURES ----------
fig1 = figure('Position', [100, 100, 800, 600]); hold on; grid on;
ax1 = gca;
ax1.FontSize = ax_fnt;
box on;
ax1.XColor = 'black';
ax1.YColor = 'black';
set(ax1,'linewidth',3);
xlabel('$ t, $ s','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ x, $ m','interpreter','latex','Fontsize',lbl_fnt);

fig2 = figure('Position', [150, 150, 800, 600]); hold on; grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on;
ax2.XColor = 'black';
ax2.YColor = 'black';
set(ax2,'linewidth',3);
xlabel('$ t, $ s','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',lbl_fnt);

fig3 = figure('Position', [200, 200, 800, 600]); hold on; grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
box on;
ax3.XColor = 'black';
ax3.YColor = 'black';
set(ax3,'linewidth',3);
xlabel('$ t, $ s','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ \kappa, $ m$^{-1}$','interpreter','latex','Fontsize',lbl_fnt);


fig4 = figure('Position', [200, 200, 800, 600]); hold on; grid on;
ax4 = gca;
ax4.FontSize = ax_fnt;
box on;
ax4.XColor = 'black';
ax4.YColor = 'black';
set(ax4,'linewidth',3);
xlabel('$ x, $ m','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ \kappa, $ m$^{-1}$','interpreter','latex','Fontsize',lbl_fnt);

fig5 = figure('Position', [200, 200, 800, 600]); hold on; grid on;
ax5 = gca;
ax5.FontSize = ax_fnt;
box on;
ax5.XColor = 'black';
ax5.YColor = 'black';
set(ax5,'linewidth',3);
xlabel('$ x, $ m','interpreter','latex','Fontsize',lbl_fnt);
ylabel('$ \chi_{\mathrm{d}}, $ deg.','interpreter','latex','Fontsize',lbl_fnt);
%% ---------- LOOP OVER Cs ----------
for i = 1:numel(Cs_vec)

    Cs = Cs_vec(i);

    % ks selection
    switch lower(matchMode)
        case "slope"
            ks = Cs^2/2;
        case "x0"
            ks = ks_match_at_x0(Cs, x0);
    end

    % Allocate
    x_sin = zeros(1,N);  x_cos = zeros(1,N);
    chi_sin = zeros(1,N); chi_cos = zeros(1,N);
    kappa_sin = zeros(1,N); kappa_cos = zeros(1,N);

    x_sin(1) = x0;
    x_cos(1) = x0;

    %% Time integration
    for k = 1:N-1

        % --- arcsin ---
        chi_sin(k) = chi_arcsin(x_sin(k), ks);
        dchidx = dchi_dx_arcsin(x_sin(k), ks);
        xdot = V*cos(chi_sin(k));
        x_sin(k+1) = x_sin(k) + dt*xdot;
        kappa_sin(k) = dchidx * cos(chi_sin(k));

        % --- arccos ---
        chi_cos(k) = chi_arccos(x_cos(k), Cs);
        dchidx = dchi_dx_arccos(x_cos(k), Cs);
        xdot = V*cos(chi_cos(k));
        x_cos(k+1) = x_cos(k) + dt*xdot;
        kappa_cos(k) = dchidx * cos(chi_cos(k));
    end

    % last samples
    chi_sin(end) = chi_arcsin(x_sin(end), ks);
    chi_cos(end) = chi_arccos(x_cos(end), Cs);

    %% ---------- PLOTS ----------
    figure(fig1);
    plot(t, x_sin, '-', 'LineWidth', lw, 'Color', cols{i});     % arcsin
    plot(t, x_cos, '--', 'LineWidth', lw, 'Color', cols{i}); % arccos

    figure(fig2);
    plot(t, rad2deg(chi_sin), '-', 'LineWidth', lw, 'Color', cols{i});
    plot(t, rad2deg(chi_cos), '--', 'LineWidth', lw, 'Color', cols{i});

    figure(fig3);
    plot(t, kappa_sin, '-', 'LineWidth', lw, 'Color', cols{i});
    plot(t, kappa_cos, '--', 'LineWidth', lw, 'Color', cols{i});

    figure(fig4);
    plot(x_sin, kappa_sin, '-', 'LineWidth', lw, 'Color', cols{i});
    plot(x_cos, kappa_cos, '--', 'LineWidth', lw, 'Color', cols{i});

    figure(fig5);
    plot(x_sin, rad2deg(chi_sin), '-', 'LineWidth', lw, 'Color', cols{i});
    plot(x_cos, rad2deg(chi_cos), '--', 'LineWidth', lw, 'Color', cols{i});


    legText{end+1} = sprintf('Arcsine, $ C_s $ = %.2f', Cs);
    legText{end+1} = sprintf('Arccos, $ C_s $ = %.2f', Cs);
end

%% Legends
figure(fig1); legend(legText, 'interpreter','latex','Location','best','Fontsize',leg_fnt);
figure(fig2); legend(legText, 'interpreter','latex','Location','best','Fontsize',leg_fnt);
figure(fig3); legend(legText, 'interpreter','latex','Location','best','Fontsize',leg_fnt);
figure(fig4); legend(legText, 'interpreter','latex','Location','best','Fontsize',leg_fnt);
figure(fig5); legend(legText, 'interpreter','latex','Location','best','Fontsize',leg_fnt);
%% ---------- LOCAL FUNCTIONS ----------
function chi = chi_arcsin(x, ks)
    s = 1/(1 + ks*x^2);
    s = min(max(s,0),1);
    if x <= 0
        chi = asin(s);
    else
        chi = pi - asin(s);
    end
end

function dchidx = dchi_dx_arcsin(x, ks)
    s = 1/(1 + ks*x^2);
    dsdx = -(2*ks*x)/(1 + ks*x^2)^2;
    denom = sqrt(max(1 - s^2, 1e-12));
    if x <= 0
        dchidx = (1/denom)*dsdx;
    else
        dchidx = -(1/denom)*dsdx;
    end
end

function chi = chi_arccos(x, Cs)
    if x <= 0
        chi = acos(1 - exp(Cs*x));
    else
        chi = pi - acos(1 - exp(-Cs*x));
    end
end

function dchidx = dchi_dx_arccos(x, Cs)
    if x <= 0
        u = 1 - exp(Cs*x);
        du = -Cs*exp(Cs*x);
        dchidx = -(1/sqrt(max(1-u^2,1e-12)))*du;
    else
        u = 1 - exp(-Cs*x);
        du = Cs*exp(-Cs*x);
        dchidx = +(1/sqrt(max(1-u^2,1e-12)))*du;
    end
end

function ks = ks_match_at_x0(Cs, x0)
    xm = abs(x0);
    a = exp(-Cs*xm);
    rhs = sqrt(max(1 - (1 - a)^2, 1e-12));
    ks = (1/rhs - 1)/(xm^2);
end
