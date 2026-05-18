%% Control-oriented plant calibration
%
% Purpose:
%   Calibrates the SIRVHDB ODE macro model to ABM ensemble behaviour using
%   a Genetic Algorithm, then linearises the fitted ODE for MPC design.
%
% Outputs:
%   macro/plant_ga.mat containing plant_ga, calibrated parameters, and
%   calibration diagnostics.
%
% Calibration rationale:
%   The ODE is fitted across multiple intervention scenarios so the MPC
%   plant remains representative over the policy range u in [0,1].
%
% USAGE:
%   run('integration/build_plant_ga.m')

clear; close all; clc;

% Resolve project root when the script is launched from another folder.
project_root = fileparts(fileparts(mfilename('fullpath')));
if isempty(project_root)
    project_root = pwd;
end
cd(project_root);

addpath(fullfile(project_root, 'macro'));
addpath(fullfile(project_root, 'micro'));
addpath(fullfile(project_root, 'integration'));

%% Parameters
run('macro/SIRVDB_model_parameters.m');
gamma_val = gamma;   % avoid MATLAB built-in conflict

params.nu1     = nu1;       params.nu2     = nu2;
params.eta     = eta;       params.gamma   = gamma_val;
params.gamma_v = gamma_v;   params.gamma_h = gamma_h;
params.mu      = mu;        params.mu_v    = mu_v;
params.mu_h    = mu_h;      params.h       = h;
params.h_v     = h_v;       params.omega_V = omega_V;
params.omega_R = omega_R;   params.k       = k;
params.M11     = 0.85;      params.M12     = 0.15;
params.M21     = 0.10;      params.M22     = 0.90;

N1   = 10000;   % agents per city
N2   = 10000;
T    = 150;
I0_1 = 10;
I0_2 = 8;
n_seeds = 5;

fprintf('Parameters loaded. N=%d per city, T=%d days.\n', N1, T);

%% Calibration scenarios
% Each scenario combines an intervention trajectory and transmission
% multiplier to cover relevant MPC operating conditions.
scenarios = struct();

% Uncontrolled epidemic baseline.
scenarios(1).name         = 'No control';
scenarios(1).u            = zeros(T,1);
scenarios(1).beta_mult    = 1.0;

% Sustained moderate national control.
scenarios(2).name         = 'Moderate national';
scenarios(2).u            = 0.3 * ones(T,1);
scenarios(2).beta_mult    = 1.0;

% Sustained strong national control.
scenarios(3).name         = 'Strong national';
scenarios(3).u            = 0.7 * ones(T,1);
scenarios(3).beta_mult    = 1.0;

% Alternating local-control proxy.
scenarios(4).name         = 'Local city control';
u_local = zeros(T,1);
u_local(1:60)   = 0.5;   % city 1 early lockdown
u_local(61:120) = 0.1;   % ease off
u_local(121:T)  = 0.4;   % reimpose partially
scenarios(4).u            = u_local;
scenarios(4).beta_mult    = 1.0;

% Delayed policy response.
scenarios(5).name         = 'Delayed response';
u_delayed = zeros(T,1);
u_delayed(30:end) = 0.6;
scenarios(5).u            = u_delayed;
scenarios(5).beta_mult    = 1.0;

% Increased seasonal transmission.
scenarios(6).name         = 'Seasonal increase';
scenarios(6).u            = zeros(T,1);
scenarios(6).beta_mult    = 1.4;

% Increased-transmissibility shock with delayed response.
scenarios(7).name         = 'Shock variant';
u_shock = zeros(T,1);
u_shock(50:end) = 0.5;   % policy response delayed by 10 days
scenarios(7).u            = u_shock;
scenarios(7).beta_mult    = 1.6;   % 60% more transmissible

n_scenarios = length(scenarios);
fprintf('\n%d calibration scenarios defined.\n', n_scenarios);

%% ODE initial conditions
x0 = zeros(14,1);
x0(1)  = 1 - I0_1/N1;   % S1
x0(2)  = I0_1/N1;        % I1
x0(8)  = 1 - I0_2/N2;   % S2
x0(9)  = I0_2/N2;        % I2
%% ABM ensemble for calibration targets
fprintf('\n--- Running ABM ensemble (%d scenarios x %d seeds) ---\n', ...
        n_scenarios, n_seeds);

H_abm = zeros(T, n_scenarios);   % hospital load (ODE units)
I_abm = zeros(T, n_scenarios);   % infected fraction
S_abm = zeros(T, n_scenarios);   % susceptible fraction

for sc = 1:n_scenarios

    u_sc      = scenarios(sc).u;
    b_mult    = scenarios(sc).beta_mult;
  
    H_runs = zeros(T, n_seeds);
    I_runs = zeros(T, n_seeds);
    S_runs = zeros(T, n_seeds);

    for s = 1:n_seeds
        rng(s * 100 + sc);
        agents = init_agents(N1, N2, I0_1, I0_2, 0, 0, params);

        for t = 1:T
            u_t  = u_sc(min(t, end));
            sf   = 1 + season_amp * sin(2*pi*t/season_period);
            b1_t = beta1 * b_mult * (1 - k*u_t) * sf;
            b2_t = beta2 * b_mult * (1 - k*u_t) * sf;

            agents = update_agents(agents, b1_t, b2_t, params, u_t);
            stats  = aggregate_stats(agents);

            H_runs(t,s) = (stats.city1.H/N1)/H_cap1 + ...
                          (stats.city2.H/N2)/H_cap2;
            I_runs(t,s) = (stats.I + stats.Iv) / (N1+N2);
            S_runs(t,s) = stats.S / (N1+N2);
        end
    end

    H_abm(:,sc) = mean(H_runs, 2);
    I_abm(:,sc) = mean(I_runs, 2);
    S_abm(:,sc) = mean(S_runs, 2);

    fprintf('  Scenario %d (%s): peak H_load=%.3f  peak I=%.1f%%\n', ...
            sc, scenarios(sc).name, max(H_abm(:,sc)), 100*max(I_abm(:,sc)));
end

fprintf('ABM ensemble complete.\n');

%% GA parameter optimisation
% Fits ODE parameters to match ABM across all calibration scenarios.
% Parameters optimised: beta1, beta2, gamma, gamma_v, h, h_v,
% mu, mu_h, omega_R, omega_V. Less sensitive parameters are fixed.
% The objective weights peak hospital load because this is the main MPC
% constraint-relevant quantity.
fprintf('\n--- Running GA optimisation ---\n');
fprintf('Fitting ODE to ABM across %d scenarios simultaneously.\n', n_scenarios);

% Fixed parameters.
p_fix.eta          = eta;
p_fix.nu1          = nu1;
p_fix.nu2          = nu2;
p_fix.gamma_h      = gamma_h;
p_fix.mu_v         = mu_v;
p_fix.H_cap1       = H_cap1;
p_fix.H_cap2       = H_cap2;
p_fix.k            = k;
p_fix.season_amp   = season_amp;
p_fix.season_period= season_period;
p_fix.M11          = 0.85;  p_fix.M12 = 0.15;
p_fix.M21          = 0.10;  p_fix.M22 = 0.90;

%          [beta1  beta2  gamma  gammav  h      hv     mu      mu_h   omgR   omgV ]
lb       = [0.15,  0.12,  1/20,  1/15,  0.001, 0.001, 0.0001, 0.001, 1/200, 1/200];
ub       = [0.80,  0.70,  1/7,   1/6,   0.020, 0.010, 0.003,  0.020, 1/50,  1/60 ];
n_params = length(lb);

obj_fn = @(theta) ga_objective(theta, p_fix, x0, ...
                                H_abm, I_abm, scenarios, T);

ga_opts = optimoptions('ga', ...
    'PopulationSize',   150, ...
    'MaxGenerations',   100, ...
    'FunctionTolerance', 1e-7, ...
    'CrossoverFraction', 0.8, ...
    'EliteCount',        10, ...
    'Display',          'iter', ...
    'UseParallel',      false);

fprintf('Starting GA (PopSize=150, MaxGen=200)...\n');
tic;
[theta_best, fval] = ga(obj_fn, n_params, [], [], [], [], ...
                         lb, ub, [], ga_opts);
elapsed = toc;

fprintf('\nGA complete in %.1f seconds. Final cost: %.6f\n', elapsed, fval);

%% Optimised parameter struct
p_opt = p_fix;
p_opt.beta1   = theta_best(1);
p_opt.beta2   = theta_best(2);
p_opt.gamma   = theta_best(3);
p_opt.gamma_v = theta_best(4);
p_opt.h       = theta_best(5);
p_opt.h_v     = theta_best(6);
p_opt.mu      = theta_best(7);
p_opt.mu_h    = theta_best(8);
p_opt.omega_R = theta_best(9);
p_opt.omega_V = theta_best(10);

fprintf('\nOptimised ODE parameters:\n');
fprintf('  beta1   = %.4f  (original: %.4f)\n', p_opt.beta1,   beta1);
fprintf('  beta2   = %.4f  (original: %.4f)\n', p_opt.beta2,   beta2);
fprintf('  gamma   = %.4f  (original: %.4f)\n', p_opt.gamma,   gamma_val);
fprintf('  gamma_v = %.4f  (original: %.4f)\n', p_opt.gamma_v, gamma_v);
fprintf('  h       = %.5f  (original: %.5f)\n', p_opt.h,       h);
fprintf('  h_v     = %.5f  (original: %.5f)\n', p_opt.h_v,     h_v);
fprintf('  mu      = %.6f  (original: %.6f)\n', p_opt.mu,      mu);
fprintf('  mu_h    = %.4f  (original: %.4f)\n', p_opt.mu_h,    mu_h);
fprintf('  omega_R = %.5f  (original: %.5f)\n', p_opt.omega_R, omega_R);
fprintf('  omega_V = %.5f  (original: %.5f)\n', p_opt.omega_V, omega_V);

%% ODE-ABM validation
fprintf('\n--- Validation plots ---\n');

figure('Name', 'GA validation — ODE vs ABM (all scenarios)', 'Color', 'w');

for sc = 1:n_scenarios
    [H_ode, t_ode] = run_ode_scenario(p_opt, scenarios(sc), T, x0);

    subplot(3, 3, sc);
    plot(1:T, H_abm(:,sc), 'b-',  'LineWidth', 2.0, ...
         'DisplayName', 'ABM mean');
    hold on;
    plot(t_ode+1, H_ode, 'r--', 'LineWidth', 2.0, ...
         'DisplayName', 'ODE fit');
    yline(H_cap1, ':k', 'Capacity', 'LineWidth', 1.0, 'FontSize', 7);

    err = abs(H_ode(1:T) - H_abm(:,sc));
    mape = 100 * mean(err(H_abm(:,sc) > 0.01)) / max(H_abm(:,sc) + 1e-6);

    title(sprintf('%s (MAPE=%.1f%%)', scenarios(sc).name, mape), ...
          'FontSize', 9);
    xlabel('Days'); ylabel('H\_load');
    if sc == 1, legend('Location','best','FontSize',7); end
    grid on;
end

subplot(3,3,8 ); axis off;
text(0.1, 0.7, sprintf('GA cost: %.5f', fval), 'FontSize', 11);
text(0.1, 0.5, sprintf('\\beta_1 = %.4f', p_opt.beta1), 'FontSize', 10);
text(0.1, 0.3, sprintf('\\gamma  = %.4f', p_opt.gamma), 'FontSize', 10);
text(0.1, 0.1, sprintf('h      = %.5f', p_opt.h), 'FontSize', 10);

sgtitle('Control-oriented calibration: ODE fit vs ABM (7 scenarios)', ...
        'FontSize', 12, 'FontWeight', 'bold');

%% Linearise fitted ODE
% Linearisation uses a controlled operating point, representative of the
% MPC regime rather than the uncontrolled initial condition.
fprintf('\n--- Linearising fitted ODE ---\n');

[~, t_op, x_op] = run_ode_scenario(p_opt, ...
    struct('u', 0.3*ones(60,1), 'beta_mult', 1.0, 'name', 'op'), 60, x0);
x_op_pt = x_op(31, :)';   % state at t=30 days
u_op    = 0.3;

% Numerical Jacobian by finite differences.
n_states = 14;
A_lin    = zeros(n_states, n_states);
B_lin    = zeros(n_states, 1);

f0  = sirvhdb_ode(30, x_op_pt, p_opt, u_op);
eps = 1e-7;

for i = 1:n_states
    xp    = x_op_pt;
    xp(i) = xp(i) + eps;
    A_lin(:,i) = (sirvhdb_ode(30, xp, p_opt, u_op) - f0) / eps;
end

fp    = sirvhdb_ode(30, x_op_pt, p_opt, u_op + eps);
B_lin = (fp - f0) / eps;

% Output matrix: H_load = H1/H_cap1 + H2/H_cap2.
C_lin       = zeros(1, n_states);
C_lin(6)    = 1 / p_opt.H_cap1;
C_lin(13)   = 1 / p_opt.H_cap2;
D_lin       = 0;

% Discretise with one-day sample time.
sys_cont  = ss(A_lin, B_lin, C_lin, D_lin);
plant_ga  = c2d(sys_cont, 1, 'zoh');

% Deflate marginal population-conservation eigenvalue for MPC stability.
eig_mags = abs(eig(plant_ga.A));
if max(eig_mags) >= 1.0
    plant_ga = ss(plant_ga.A * 0.9999, plant_ga.B, ...
                  plant_ga.C, plant_ga.D, 1);
    fprintf('Applied stability deflation (eig was %.6f)\n', max(eig_mags));
end

fprintf('Linearised plant:\n');
fprintf('  States      : %d\n',    order(plant_ga));
fprintf('  Max eig     : %.6f\n',  max(abs(eig(plant_ga.A))));
fprintf('  DC gain     : %.6f\n',  dcgain(plant_ga));
fprintf('  Ts          : %.1f day\n', plant_ga.Ts);

%% Save calibrated plant
save(fullfile('macro', 'plant_ga.mat'), ...
    'plant_ga', 'p_opt', 'theta_best', 'fval', ...
    'H_abm', 'I_abm', 'S_abm', 'scenarios');

fprintf('\nSaved to macro/plant_ga.mat\n');
fprintf('\nTo use in main_v4.m, replace the plant loading section:\n');
fprintf('  load(fullfile(''macro'',''plant_ga.mat''));\n');
fprintf('  clear mpcobj;\n');
fprintf('  mpcobj = mpc(plant_ga, 1);\n');
fprintf('  %% then set weights as before\n');

%% Local functions

function cost = ga_objective(theta, p_fix, x0, H_abm, I_abm, scenarios, T)
% Weighted MSE between ODE and ABM across calibration scenarios.

    p = p_fix;
    p.beta1   = theta(1);
    p.beta2   = theta(2);
    p.gamma   = theta(3);
    p.gamma_v = theta(4);
    p.h       = theta(5);
    p.h_v     = theta(6);
    p.mu      = theta(7);
    p.mu_h    = theta(8);
    p.omega_R = theta(9);
    p.omega_V = theta(10);

    cost = 0;
    n_sc = length(scenarios);

    for sc = 1:n_sc
        try
            [H_ode, t_ode] = run_ode_scenario(p, scenarios(sc), T, x0);

            H_daily = interp1(t_ode, H_ode, 0:T-1, 'pchip', 'extrap')';
            H_daily = max(0, H_daily);

            % Higher weight around peak hospitalisation.
            H_ref_sc = H_abm(:, sc);
            peak_val = max(H_ref_sc) + 1e-6;
            weights  = 1 + 4 * (H_ref_sc / peak_val);

            residual = H_daily(1:T) - H_ref_sc;
            cost     = cost + mean(weights .* residual.^2);

            % Secondary penalty on infected fraction.
            I_ref_sc = I_abm(:, sc);
            I_ode_sc = interp1(t_ode, ...
                get_I_frac(t_ode, p, scenarios(sc), T, x0), ...
                0:T-1, 'pchip', 'extrap')';
            cost = cost + 0.3 * mean((I_ode_sc(1:T) - I_ref_sc).^2);

        catch
            cost = cost + 1e6;
        end
    end
end


function [H_load, t_out, x_out] = run_ode_scenario(p, scenario, T, x0)
% Runs the SIRVHDB ODE for one calibration scenario.

    u_series = scenario.u;
    b_mult   = scenario.beta_mult;

    % Piecewise daily integration allows day-specific intervention values.
    t_out = [];
    x_out = [];
    x_cur = x0;

    for day = 1:T
        u_t     = u_series(min(day, end));
        t_span  = [day-1, day];

        p_sc        = p;
        p_sc.beta1  = p.beta1 * b_mult;
        p_sc.beta2  = p.beta2 * b_mult;

        opts = odeset('RelTol',1e-5,'AbsTol',1e-7,'NonNegative',1:14);
        [t_seg, x_seg] = ode45(@(t,x) sirvhdb_ode(t, x, p_sc, u_t), ...
                                t_span, x_cur, opts);

        t_out = [t_out; t_seg];
        x_out = [x_out; x_seg];
        x_cur = x_seg(end, :)';
    end

    [t_out, idx] = unique(t_out);
    x_out = x_out(idx, :);

    H1     = max(x_out(:,6),  0);
    H2     = max(x_out(:,13), 0);
    H_load = H1/p.H_cap1 + H2/p.H_cap2;
end


function I_frac = get_I_frac(t_ode, p, scenario, T, x0)
% Returns infected fraction from the ODE solution.

    [~, t_out, x_out] = run_ode_scenario(p, scenario, T, x0);
    I_frac = (max(x_out(:,2),0) + max(x_out(:,3),0) + ...
              max(x_out(:,9),0) + max(x_out(:,10),0)) / 2;
    I_frac = interp1(t_out, I_frac, t_ode, 'pchip', 'extrap');
end
