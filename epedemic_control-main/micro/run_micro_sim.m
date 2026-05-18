function results = run_micro_sim(beta_t1_series, beta_t2_series, H_cap1, H_cap2)
% RUN_MICRO_SIM
% Runs a compact ABM simulation for prescribed transmission-rate series.
%
% Inputs:
%   beta_t1_series,beta_t2_series - city transmission rates over time
%   H_cap1,H_cap2                 - optional hospital capacity fractions
%
% Output:
%   results - struct with compartment counts and feedback signals.

rng(1);

T = length(beta_t1_series);

%% Parameters
params.nu1 = 0.0000;
params.nu2 = 0.0000;

params.eta = 0.6;
params.gamma = 1/12;
params.gamma_v = 1/9;
params.gamma_h = 1/25;

params.mu = 0.001;
params.mu_v = 0.0005;
params.mu_h = 0.02;

params.h = 0.05;
params.h_v = 0.02;

params.omega_V = 1/140;
params.omega_R = 1/110;

%% Population
N1 = 500;
N2 = 500;

I0_1 = 10;
I0_2 = 8;
V0_1 = 0;
V0_2 = 0;

agents = init_agents(N1, N2, I0_1, I0_2, V0_1, V0_2);

%% Storage
results.time = (0:T-1)';

results.S  = zeros(T,1);
results.I  = zeros(T,1);
results.Iv = zeros(T,1);
results.R  = zeros(T,1);
results.V  = zeros(T,1);
results.H  = zeros(T,1);
results.D  = zeros(T,1);

results.I_total = zeros(T,1);
results.H_total = zeros(T,1);
results.H_load  = zeros(T,1);

results.city1_I = zeros(T,1);
results.city2_I = zeros(T,1);
results.city1_H = zeros(T,1);
results.city2_H = zeros(T,1);

results.city1_H_load = zeros(T,1);
results.city2_H_load = zeros(T,1);

%% Simulation loop
for t = 1:T
    beta_t1 = beta_t1_series(t);
    beta_t2 = beta_t2_series(t);

    agents = update_agents(agents, beta_t1, beta_t2, params);

    if nargin >= 3
        stats = aggregate_stats(agents, H_cap1, H_cap2);
    else
        stats = aggregate_stats(agents);
    end

    results.S(t)  = stats.S;
    results.I(t)  = stats.I;
    results.Iv(t) = stats.Iv;
    results.R(t)  = stats.R;
    results.V(t)  = stats.V;
    results.H(t)  = stats.H;
    results.D(t)  = stats.D;

    % Feedback signals used by macro-level control.
    results.I_total(t) = stats.I_total;
    results.H_total(t) = stats.H;

    if ~isempty(stats.H_load)
        results.H_load(t) = stats.H_load;
        results.city1_H_load(t) = stats.city1.H_load;
        results.city2_H_load(t) = stats.city2.H_load;
    end

    results.city1_I(t) = stats.city1.I_total;
    results.city2_I(t) = stats.city2.I_total;

    results.city1_H(t) = stats.city1.H;
    results.city2_H(t) = stats.city2.H;
end

%% Plot
figure;
plot(results.time, results.I_total, 'LineWidth', 1.8); hold on;
plot(results.time, results.H_total, 'LineWidth', 1.8);
plot(results.time, results.D, 'LineWidth', 1.8);
xlabel('Days');
ylabel('Agents');
legend('Total infected', 'Hospitalized', 'Dead');
title('Micro-model epidemic results');
grid on;

end
