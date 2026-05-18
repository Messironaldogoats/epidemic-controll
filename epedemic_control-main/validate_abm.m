clear; close all; clc;

%% ABM validation against ODE
%
% Purpose:
%   Compares the ABM ensemble mean with the two-city SIRVHDB ODE under
%   no intervention. This validation supports later GA calibration.
%
% Outputs:
%   Diagnostic figures for global and per-city dynamics.
%   integration/abm_v2_ensemble.mat

addpath('macro');
addpath('micro');

%% Parameters
run('macro/SIRVDB_model_parameters.m');

% Travel matrix used by update_agents for cross-city mixing.
params.nu1     = nu1;     params.nu2     = nu2;
params.eta     = eta;     params.gamma   = gamma;
params.gamma_v = gamma_v; params.gamma_h = gamma_h;
params.mu      = mu;      params.mu_v    = mu_v;
params.mu_h    = mu_h;    params.h       = h;
params.h_v     = h_v;     params.omega_V = omega_V;
params.omega_R = omega_R; params.k       = k;
params.M11 = 0.85;        params.M12 = 0.15;
params.M21 = 0.10;        params.M22 = 0.90;

N1 = 10000;
N2 = 10000;
T  = 150;
I0_1 = 10;
I0_2 = 8;
n_runs = 10;     % ensemble size
u_fixed = 0.0;   % no intervention

fprintf('Validating ABM v2 (N=%d per city, %d runs)...\n', N1, n_runs);

%% ABM ensemble
S_runs  = zeros(T, n_runs);
I_runs  = zeros(T, n_runs);
H_runs  = zeros(T, n_runs);
D_runs  = zeros(T, n_runs);

S1_runs = zeros(T, n_runs);
I1_runs = zeros(T, n_runs);
H1_runs = zeros(T, n_runs);

S2_runs = zeros(T, n_runs);
I2_runs = zeros(T, n_runs);
H2_runs = zeros(T, n_runs);

for r = 1:n_runs
    rng(r);
    agents = init_agents(N1, N2, I0_1, I0_2, 0, 0, params);

    for t = 1:T
        sf = 1 + season_amp * sin(2*pi*t/season_period);
        b1 = beta1 * (1 - k*u_fixed) * sf;
        b2 = beta2 * (1 - k*u_fixed) * sf;

        agents = update_agents(agents, b1, b2, params, u_fixed);
        stats  = aggregate_stats(agents);

        S_runs(t,r)  = stats.S  / (N1+N2);
        I_runs(t,r)  = stats.I_total / (N1+N2);
        H_runs(t,r)  = stats.H  / (N1+N2);
        D_runs(t,r)  = stats.D  / (N1+N2);

        S1_runs(t,r) = stats.city1.S / N1;
        I1_runs(t,r) = stats.city1.I_total / N1;
        H1_runs(t,r) = stats.city1.H / N1;

        S2_runs(t,r) = stats.city2.S / N2;
        I2_runs(t,r) = stats.city2.I_total / N2;
        H2_runs(t,r) = stats.city2.H / N2;
    end

    fprintf('  Run %d complete — peak I=%.1f%%\n', r, 100*max(I_runs(:,r)));
end

%% Ensemble statistics
days = (1:T)';

S_mean = mean(S_runs, 2);   S_std = std(S_runs, 0, 2);
I_mean = mean(I_runs, 2);   I_std = std(I_runs, 0, 2);
H_mean = mean(H_runs, 2);   H_std = std(H_runs, 0, 2);
D_mean = mean(D_runs, 2);

I1_mean = mean(I1_runs, 2);
I2_mean = mean(I2_runs, 2);
H1_mean = mean(H1_runs, 2);
H2_mean = mean(H2_runs, 2);

fprintf('\n=== ABM Ensemble Statistics ===\n');
fprintf('  Peak I (mean)  : %.2f%%\n', 100*max(I_mean));
fprintf('  Peak H (mean)  : %.4f\n',   max(H_mean));
fprintf('  Final D (mean) : %.2f%%\n', 100*D_mean(end));
fprintf('  I std at peak  : %.4f (%.1f%% of mean)\n', ...
        I_std(I_mean==max(I_mean)), ...
        100*I_std(I_mean==max(I_mean))/max(I_mean));

%% ODE comparison
fprintf('\nRunning ODE for comparison...\n');

% Initial conditions as population fractions.
x0 = [1 - I0_1/N1;   % S1
       I0_1/N1;       % I1
       0;             % Iv1
       0;             % R1
       0;             % V1
       0;             % H1
       0;             % D1
       1 - I0_2/N2;   % S2
       I0_2/N2;       % I2
       0; 0; 0; 0; 0];

p_ode.beta1 = beta1; p_ode.beta2 = beta2;
p_ode.gamma = gamma; p_ode.gamma_v = gamma_v;
p_ode.gamma_h = gamma_h;
p_ode.mu = mu; p_ode.mu_v = mu_v; p_ode.mu_h = mu_h;
p_ode.h = h; p_ode.h_v = h_v;
p_ode.eta = eta; p_ode.nu1 = nu1; p_ode.nu2 = nu2;
p_ode.omega_V = omega_V; p_ode.omega_R = omega_R;
p_ode.k = k;
p_ode.season_amp = season_amp;
p_ode.season_period = season_period;
p_ode.M11 = 0.85; p_ode.M12 = 0.15;
p_ode.M21 = 0.10; p_ode.M22 = 0.90;
p_ode.H_cap1 = H_cap1; p_ode.H_cap2 = H_cap2;

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:14);
[t_ode, x_ode] = ode45(@(t,x) sirvhdb_ode_val(t, x, p_ode, 0.0), ...
                         [0 T], x0, opts);

t_daily = 0:T-1;
x_daily = interp1(t_ode, x_ode, t_daily, 'pchip');

S_ode = (x_daily(:,1) + x_daily(:,8)) / 2;
I_ode = (x_daily(:,2) + x_daily(:,3) + x_daily(:,9) + x_daily(:,10)) / 2;
H_ode = (x_daily(:,6) + x_daily(:,13)) / 2;
D_ode = (x_daily(:,7) + x_daily(:,14)) / 2;

fprintf('  ODE peak I: %.2f%%\n', 100*max(I_ode));
fprintf('  ODE peak H: %.4f\n', max(H_ode));

%% Global ABM-ODE comparison
figure('Name', 'ABM v2 vs ODE — global', 'Color', 'w');

subplot(3,1,1);
for r = 1:n_runs
    plot(days, S_runs(:,r), 'b-', 'LineWidth', 0.5, ...
         'Color', [0.2 0.45 1.0 0.2], 'HandleVisibility', 'off');
    hold on;
end
plot(days, S_mean, 'b-',  'LineWidth', 2.0, 'DisplayName', 'ABM mean');
plot(t_daily+1, S_ode, 'r--', 'LineWidth', 2.0, 'DisplayName', 'ODE');
ylabel('S fraction'); legend('Location','best'); grid on;
title('ABM v2 validation: ensemble mean vs ODE (no intervention)');

subplot(3,1,2);
for r = 1:n_runs
    plot(days, I_runs(:,r), 'r-', 'LineWidth', 0.5, ...
         'Color', [1.0 0.1 0.1 0.15], 'HandleVisibility', 'off');
    hold on;
end
plot(days, I_mean, 'r-',  'LineWidth', 2.0, 'DisplayName', 'ABM mean');
plot(t_daily+1, I_ode, 'k--', 'LineWidth', 2.0, 'DisplayName', 'ODE');
ylabel('I fraction'); legend('Location','best'); grid on;

subplot(3,1,3);
for r = 1:n_runs
    plot(days, H_runs(:,r), 'g-', 'LineWidth', 0.5, ...
         'Color', [0.1 0.72 0.2 0.2], 'HandleVisibility', 'off');
    hold on;
end
plot(days, H_mean, 'g-',  'LineWidth', 2.0, 'DisplayName', 'ABM mean');
plot(t_daily+1, H_ode, 'm--', 'LineWidth', 2.0, 'DisplayName', 'ODE');
xlabel('Days'); ylabel('H fraction'); legend('Location','best'); grid on;

%% City comparison
figure('Name', 'ABM v2 — city comparison', 'Color', 'w');

subplot(2,1,1);
for r = 1:n_runs
    plot(days, I1_runs(:,r), 'b-', 'LineWidth', 0.5, ...
         'Color', [0.2 0.45 1.0 0.2], 'HandleVisibility', 'off');
    hold on;
    plot(days, I2_runs(:,r), 'r-', 'LineWidth', 0.5, ...
         'Color', [1.0 0.1 0.1 0.15], 'HandleVisibility', 'off');
end
plot(days, I1_mean, 'b-', 'LineWidth', 2.0, 'DisplayName', 'City 1 mean');
plot(days, I2_mean, 'r-', 'LineWidth', 2.0, 'DisplayName', 'City 2 mean');
ylabel('I fraction'); legend('Location','best'); grid on;
title('Per-city infected — independent stochastic dynamics');

subplot(2,1,2);
for r = 1:n_runs
    plot(days, H1_runs(:,r), 'b-', 'LineWidth', 0.5, ...
         'Color', [0.2 0.45 1.0 0.2], 'HandleVisibility', 'off');
    hold on;
    plot(days, H2_runs(:,r), 'r-', 'LineWidth', 0.5, ...
         'Color', [1.0 0.1 0.1 0.15], 'HandleVisibility', 'off');
end
plot(days, H1_mean, 'b-', 'LineWidth', 2.0, 'DisplayName', 'City 1 mean');
plot(days, H2_mean, 'r-', 'LineWidth', 2.0, 'DisplayName', 'City 2 mean');
xlabel('Days'); ylabel('H fraction'); legend('Location','best'); grid on;
title('Per-city hospitalised');

%% Stochastic variation
figure('Name', 'ABM v2 — stochastic variation', 'Color', 'w');

I_lo = I_mean - 2*I_std;
I_hi = I_mean + 2*I_std;

fill([days; flipud(days)], [I_lo; flipud(I_hi)], ...
     [0.2 0.45 1.0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
     'DisplayName', '\pm2\sigma band');
hold on;
plot(days, I_mean, 'b-',  'LineWidth', 2.0, 'DisplayName', 'ABM mean');
plot(t_daily+1, I_ode, 'r--', 'LineWidth', 2.0, 'DisplayName', 'ODE');
xlabel('Days'); ylabel('Infected fraction');
legend('Location','best'); grid on;
title('ABM stochastic variation vs ODE mean-field');

%% Error metrics
len = min(length(I_mean), length(I_ode));
err_I = abs(I_mean(1:len) - I_ode(1:len));
err_H = abs(H_mean(1:len) - H_ode(1:len));

fprintf('\n=== ABM vs ODE Agreement ===\n');
fprintf('  MAE infected  : %.6f (%.1f%% of peak I)\n', ...
        mean(err_I), 100*mean(err_I)/max(I_ode));
fprintf('  MAE hosp      : %.6f (%.1f%% of peak H)\n', ...
        mean(err_H), 100*mean(err_H)/max(H_ode));

if max(err_I)/max(I_ode) < 0.15
    fprintf('  ✓ ABM tracks ODE within 15%% — ACCEPTABLE\n');
else
    fprintf('  ✗ ABM deviates >15%% from ODE — GA tuning needed\n');
end

%% Save ensemble data for GA
save('integration/abm_v2_ensemble.mat', ...
     'S_runs', 'I_runs', 'H_runs', 'D_runs', ...
     'S1_runs', 'I1_runs', 'H1_runs', ...
     'S2_runs', 'I2_runs', 'H2_runs', ...
     'S_mean', 'I_mean', 'H_mean', 'D_mean', ...
     'I1_mean', 'I2_mean', 'H1_mean', 'H2_mean', ...
     'S_ode', 'I_ode', 'H_ode', 'D_ode', ...
     'N1', 'N2', 'T', 'days');

fprintf('\nEnsemble data saved to integration/abm_v2_ensemble.mat\n');
fprintf('Use this for GA parameter tuning.\n');

%% Local ODE function
function dxdt = sirvhdb_ode_val(t, x, p, u_val)
% Two-city SIRVHDB ODE used for validation.

S1=max(x(1),0);  I1=max(x(2),0);  Iv1=max(x(3),0);
R1=max(x(4),0);  V1=max(x(5),0);  H1=max(x(6),0); D1=max(x(7),0);
S2=max(x(8),0);  I2=max(x(9),0);  Iv2=max(x(10),0);
R2=max(x(11),0); V2=max(x(12),0); H2=max(x(13),0); D2=max(x(14),0);

N1 = max(S1+I1+Iv1+R1+V1+H1, 1e-10);
N2 = max(S2+I2+Iv2+R2+V2+H2, 1e-10);

sf    = 1 + p.season_amp * sin(2*pi*t/p.season_period);
b1    = p.beta1 * (1 - p.k*u_val) * sf;
b2    = p.beta2 * (1 - p.k*u_val) * sf;

f1 = (I1+Iv1)/N1;
f2 = (I2+Iv2)/N2;

lam1  = b1 * (p.M11*f1 + p.M12*f2);
lam2  = b2 * (p.M21*f1 + p.M22*f2);
lam1v = p.eta * lam1;
lam2v = p.eta * lam2;

dS1  = -lam1*S1  - p.nu1*S1  + p.omega_R*R1 + p.omega_V*V1;
dI1  =  lam1*S1  - (p.gamma + p.h + p.mu)*I1;
dIv1 =  lam1v*V1 - (p.gamma_v + p.h_v + p.mu_v)*Iv1;
dR1  =  p.gamma*I1 + p.gamma_v*Iv1 + p.gamma_h*H1 - p.omega_R*R1;
dV1  =  p.nu1*S1 - lam1v*V1 - p.omega_V*V1;
dH1  =  p.h*I1 + p.h_v*Iv1 - (p.gamma_h + p.mu_h)*H1;
dD1  =  p.mu*I1 + p.mu_v*Iv1 + p.mu_h*H1;

dS2  = -lam2*S2  - p.nu2*S2  + p.omega_R*R2 + p.omega_V*V2;
dI2  =  lam2*S2  - (p.gamma + p.h + p.mu)*I2;
dIv2 =  lam2v*V2 - (p.gamma_v + p.h_v + p.mu_v)*Iv2;
dR2  =  p.gamma*I2 + p.gamma_v*Iv2 + p.gamma_h*H2 - p.omega_R*R2;
dV2  =  p.nu2*S2 - lam2v*V2 - p.omega_V*V2;
dH2  =  p.h*I2 + p.h_v*Iv2 - (p.gamma_h + p.mu_h)*H2;
dD2  =  p.mu*I2 + p.mu_v*Iv2 + p.mu_h*H2;

dxdt = [dS1;dI1;dIv1;dR1;dV1;dH1;dD1;
        dS2;dI2;dIv2;dR2;dV2;dH2;dD2];
end
