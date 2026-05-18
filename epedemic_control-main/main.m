clear; close all; clc;

project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end

%% Closed-loop epidemic control
%
% Purpose:
%   Couples the two-city ABM with an MPC controller and compares the
%   closed-loop intervention policy against an uncontrolled baseline.
%
% Main outputs:
%   Figures comparing infection, hospitalisation, deaths, intervention
%   signal, and city-level hospital load.
%   integration/results_closed_loop.mat
%
% Model components:
%   Micro model: ABM with contact heterogeneity, hazard-based infection,
%   cross-city mixing, and independent city-level random streams.
%
%   Macro model: discrete plant from a linearised SIRVHDB model, used by
%   the MPC to regulate hospital load.
%
% Parameters:
%   macro/SIRVDB_model_parameters.m

addpath(fullfile(project_root, 'macro'));
addpath(fullfile(project_root, 'micro'));
addpath(fullfile(project_root, 'integration'));

%% Shared parameters
run(fullfile(project_root, 'macro', 'SIRVDB_model_parameters.m'));
project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end
gamma_val  = gamma;        % recovery rate; avoids MATLAB built-in conflict
H_ref_load = H_ref;        % hospital-load reference used by the controller

fprintf('Parameters loaded.\n');
fprintf('  beta1=%.4f  beta2=%.4f  k=%.2f\n', beta1, beta2, k);
fprintf('  H_cap=%.4f  H_ref=%.4f\n', H_cap1, H_ref_load);

%% ABM parameter structure
params.nu1      = nu1;
params.nu2      = nu2;
params.eta      = eta;
params.gamma    = gamma_val;
params.gamma_v  = gamma_v;
params.gamma_h  = gamma_h;
params.mu       = mu;
params.mu_v     = mu_v;
params.mu_h     = mu_h;
params.h        = h;
params.h_v      = h_v;
params.omega_V  = omega_V;
params.omega_R  = omega_R;
params.k        = k;
params.M11      = 0.85;    % city 1 contacts within city 1
params.M12      = 0.15;    % city 1 contacts with city 2
params.M21      = 0.10;    % city 2 contacts with city 1
params.M22      = 0.90;    % city 2 contacts within city 2

%% MPC setup
load(fullfile(project_root, 'macro', 'plant_ga.mat'));   % loads plant_ga

clear mpcobj;
mpcobj = mpc(plant_ga, 1);

mpcobj.PredictionHorizon                = 20;
mpcobj.ControlHorizon                   = 10;
mpcobj.Weights.ManipulatedVariables     = 0.001;
mpcobj.Weights.ManipulatedVariablesRate = 0.01;
mpcobj.Weights.OutputVariables          = 10;
mpcobj.MV.Min = 0;
mpcobj.MV.Max = 1;

xmpc = mpcstate(mpcobj);

fprintf('\nMPC ready.\n');
fprintf('  States=%d  max_eig=%.4f  DC_gain=%.4f\n', ...
        order(plant_ga), ...
        max(abs(eig(plant_ga.A))), ...
        dcgain(plant_ga));
fprintf('  PH=%d  CH=%d  H_ref=%.4f\n', ...
        mpcobj.PredictionHorizon, ...
        mpcobj.ControlHorizon, H_ref_load);

%% Simulation settings
T_total = 150;
N1      = 50000;
N2      = 50000;
I0_1    = 50;     % 0.1% initially infected city 1
I0_2    = 40;     % 0.08% initially infected city 2

fprintf('\nSimulation: T=%d days  N=%d+%d agents\n', T_total, N1, N2);
fprintf('  Hospital capacity : %.0f agents/city\n', H_cap1*N1);
fprintf('  H_ref agents      : %.0f total\n', ...
        H_ref_load * (H_cap1*N1 + H_cap2*N2));

%% Storage
all_days   = (1:T_total)';
all_u      = zeros(T_total, 1);
all_H_load = zeros(T_total, 1);
all_beta1  = zeros(T_total, 1);
all_beta2  = zeros(T_total, 1);

% Closed-loop states
cl_S  = zeros(T_total,1);  cl_I  = zeros(T_total,1);
cl_Iv = zeros(T_total,1);  cl_R  = zeros(T_total,1);
cl_V  = zeros(T_total,1);  cl_H  = zeros(T_total,1);
cl_D  = zeros(T_total,1);

% Per-city (needed for H_load feedback)
cl_c1_I = zeros(T_total,1);  cl_c1_H = zeros(T_total,1);
cl_c2_I = zeros(T_total,1);  cl_c2_H = zeros(T_total,1);

% Open-loop baseline (no MPC)
ol_I = zeros(T_total,1);
ol_H = zeros(T_total,1);
ol_D = zeros(T_total,1);

%% Initialise agents
% Identical seeds give a fair closed-loop/open-loop comparison.
rng(1);
agents_cl = init_agents(N1, N2, I0_1, I0_2, 0, 0, params);
rng(1);
agents_ol = init_agents(N1, N2, I0_1, I0_2, 0, 0, params);

fprintf('\n--- Running closed-loop MPC + ABM v2 (%d days) ---\n', T_total);

%% Main simulation loop
for t = 1:T_total

    %% Hospital-load feedback
    % Simulink-equivalent signal:
    % H_load = (H1_frac / H_cap1) + (H2_frac / H_cap2).
    if t == 1
        y_meas = 0;
    else
        y_meas = (cl_c1_H(t-1) / N1) / H_cap1 + ...
                 (cl_c2_H(t-1) / N2) / H_cap2;
    end

    all_H_load(t) = y_meas;

    %% MPC intervention
    % The controller predicts hospital load over the configured horizon
    % and regulates toward H_ref_load.
    u_opt = mpcmove(mpcobj, xmpc, y_meas, H_ref_load, []);
    u_now = max(0, min(1, u_opt));
    all_u(t) = u_now;

    %% Transmission modulation
    % beta(t) = beta0 * (1 - k*u) * seasonal_factor.
    sf      = 1 + season_amp * sin(2 * pi * t / season_period);
    beta1_t = beta1 * (1 - k * u_now) * sf;
    beta2_t = beta2 * (1 - k * u_now) * sf;

    all_beta1(t) = beta1_t;
    all_beta2(t) = beta2_t;

    %% Closed-loop ABM step
    % update_agents applies hazard-based infection, intervention-dependent
    % contact reduction, cross-city mixing, and independent city streams.
    agents_cl = update_agents(agents_cl, beta1_t, beta2_t, params, u_now);
    stats_cl  = aggregate_stats(agents_cl);

    cl_S(t)    = stats_cl.S;
    cl_I(t)    = stats_cl.I;
    cl_Iv(t)   = stats_cl.Iv;
    cl_R(t)    = stats_cl.R;
    cl_V(t)    = stats_cl.V;
    cl_H(t)    = stats_cl.H;
    cl_D(t)    = stats_cl.D;
    cl_c1_I(t) = stats_cl.city1.I_total;
    cl_c1_H(t) = stats_cl.city1.H;
    cl_c2_I(t) = stats_cl.city2.I_total;
    cl_c2_H(t) = stats_cl.city2.H;

    %% Open-loop ABM baseline
    beta1_ol = beta1 * sf;
    beta2_ol = beta2 * sf;

    agents_ol = update_agents(agents_ol, beta1_ol, beta2_ol, params, 0);
    stats_ol  = aggregate_stats(agents_ol);

    ol_I(t) = stats_ol.I + stats_ol.Iv;
    ol_H(t) = stats_ol.H;
    ol_D(t) = stats_ol.D;

end

%% Summary
N_total = N1 + N2;

fprintf('\n=== Simulation complete ===\n');
fprintf('  Closed-loop deaths  : %d (%.1f%%)\n', ...
        cl_D(end), 100*cl_D(end)/N_total);
fprintf('  Open-loop deaths    : %d (%.1f%%)\n', ...
        ol_D(end), 100*ol_D(end)/N_total);
fprintf('  Deaths prevented    : %d\n', ol_D(end) - cl_D(end));
fprintf('  Peak H_load (CL)    : %.4f\n', max(all_H_load));
fprintf('  H_ref_load target   : %.4f\n', H_ref_load);
fprintf('  Peak u(t)           : %.4f\n', max(all_u));
fprintf('  Mean u(t)           : %.4f\n', mean(all_u));
fprintf('  Days u > 0          : %d / %d\n', sum(all_u > 0.01), T_total);

%% Plots

% --- Plot 1: Infected ---
figure('Name','Infected: CL vs OL','Color','w');
plot(all_days, cl_I+cl_Iv, 'b-',  'LineWidth', 1.8); hold on;
plot(all_days, ol_I,       'b--', 'LineWidth', 1.5);
xlabel('Days'); ylabel('Total infected (agents)');
legend('Closed-loop (MPC)','Open-loop (no MPC)','Location','best');
title('Infected: closed-loop vs open-loop'); grid on;

% --- Plot 2: Hospitalised ---
figure('Name','Hospitalised: CL vs OL','Color','w');
plot(all_days, cl_H, 'r-',  'LineWidth', 1.8); hold on;
plot(all_days, ol_H, 'r--', 'LineWidth', 1.5);
yline(H_cap1*N1 + H_cap2*N2, '--k', 'Hospital capacity', 'LineWidth', 1.2);
xlabel('Days'); ylabel('Hospitalised (agents)');
legend('Closed-loop (MPC)','Open-loop (no MPC)','Location','best');
title('Hospitalised: closed-loop vs open-loop'); grid on;

% --- Plot 3: Deaths ---
figure('Name','Deaths: CL vs OL','Color','w');
plot(all_days, cl_D, 'k-',  'LineWidth', 1.8); hold on;
plot(all_days, ol_D, 'k--', 'LineWidth', 1.5);
xlabel('Days'); ylabel('Cumulative deaths (agents)');
legend('Closed-loop (MPC)','Open-loop (no MPC)','Location','best');
title('Deaths: closed-loop vs open-loop'); grid on;

% Hospital-load tracking is the primary control performance metric.
figure('Name','H_load vs MPC reference','Color','w');
plot(all_days, all_H_load, 'b-', 'LineWidth', 1.8); hold on;
yline(H_ref_load, '--r', ...
      sprintf('H_{ref} = %.2f', H_ref_load), 'LineWidth', 1.5);
yline(1/H_cap1 + 1/H_cap2, '--k', 'Max capacity', 'LineWidth', 1.0);
xlabel('Days'); ylabel('H\_load (ODE fraction units)');
title('Hospital load vs MPC reference — ABM feedback');
legend('H\_load (ABM)','Reference','Max capacity','Location','best');
grid on;

% --- Plot 5: MPC intervention signal ---
figure('Name','MPC intervention u(t)','Color','w');
stairs(all_days, all_u, 'b-', 'LineWidth', 1.8);
xlabel('Days');
ylabel('u(t)   [0 = no restriction,   1 = full lockdown]');
title('MPC intervention signal u(t)');
ylim([0 1.05]); grid on;

% --- Plot 6: Effective beta per city ---
figure('Name','Effective beta','Color','w');
plot(all_days, all_beta1, 'b-', 'LineWidth', 1.8); hold on;
plot(all_days, all_beta2, 'r-', 'LineWidth', 1.8);
yline(beta1, '--b', '\beta_0 city 1', 'LineWidth', 1.0);
yline(beta2, '--r', '\beta_0 city 2', 'LineWidth', 1.0);
xlabel('Days'); ylabel('\beta(t)');
title('Effective transmission rate under MPC control');
legend('City 1 \beta(t)','City 2 \beta(t)','Location','best');
grid on;

% --- Plot 7: Per-city H_load ---
figure('Name','Per-city hospital load','Color','w');
H_load1 = (cl_c1_H / N1) / H_cap1;
H_load2 = (cl_c2_H / N2) / H_cap2;

subplot(2,1,1);
plot(all_days, H_load1, 'b-', 'LineWidth', 1.8); hold on;
yline(1.0,        '--k', 'Capacity', 'LineWidth', 1.0);
yline(H_ref_load, ':r',  'H_{ref}',  'LineWidth', 1.2);
ylabel('H load city 1');
title('Per-city hospital load — closed-loop'); grid on;

subplot(2,1,2);
plot(all_days, H_load2, 'b-', 'LineWidth', 1.8); hold on;
yline(1.0,        '--k', 'Capacity', 'LineWidth', 1.0);
yline(H_ref_load, ':r',  'H_{ref}',  'LineWidth', 1.2);
xlabel('Days'); ylabel('H load city 2'); grid on;

% --- Plot 8: Full population state stack ---
figure('Name','Population state stack','Color','w');
area(all_days, ...
     [cl_S, cl_V, cl_I+cl_Iv, cl_R, cl_H, cl_D], ...
     'LineWidth', 0.5);
colororder([0.20 0.45 1.00;
            0.62 0.25 0.90;
            1.00 0.10 0.10;
            0.10 0.72 0.20;
            1.00 0.58 0.10;
            0.20 0.20 0.20]);
xlabel('Days'); ylabel('Agents');
legend('S','V','I+Iv','R','H','D','Location','best');
title('Full population state — closed-loop'); grid on;

% --- Plot 9: City comparison (stochastic divergence) ---
figure('Name','City comparison','Color','w');
subplot(2,1,1);
plot(all_days, cl_c1_I, 'b-', 'LineWidth', 1.8); hold on;
plot(all_days, cl_c2_I, 'r-', 'LineWidth', 1.8);
ylabel('Infected (agents)');
legend('City 1','City 2','Location','best');
title('City 1 vs City 2 — stochastic divergence (ABM v2)'); grid on;

subplot(2,1,2);
plot(all_days, cl_c1_H, 'b-', 'LineWidth', 1.8); hold on;
plot(all_days, cl_c2_H, 'r-', 'LineWidth', 1.8);
yline(H_cap1*N1, '--k', 'Capacity', 'LineWidth', 1.0);
xlabel('Days'); ylabel('Hospitalised (agents)');
legend('City 1','City 2','Location','best'); grid on;

%% Save results
results_cl.time    = all_days;
results_cl.S       = cl_S;   results_cl.I   = cl_I;
results_cl.Iv      = cl_Iv;  results_cl.R   = cl_R;
results_cl.V       = cl_V;   results_cl.H   = cl_H;
results_cl.D       = cl_D;
results_cl.city1_I = cl_c1_I; results_cl.city1_H = cl_c1_H;
results_cl.city2_I = cl_c2_I; results_cl.city2_H = cl_c2_H;

save(fullfile(project_root, 'integration', 'results_closed_loop.mat'), ...
    'results_cl', ...
    'all_u',  'all_H_load', 'all_beta1', 'all_beta2', 'all_days', ...
    'ol_I',   'ol_H',       'ol_D', ...
    'H_cap1', 'H_cap2',     'H_ref', 'H_ref_load', 'N1', 'N2');

fprintf('\nResults saved to integration/results_closed_loop.mat\n');
fprintf('Done.\n');
