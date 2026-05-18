clear; close all; clc;

%% Beta reconstruction validation
%
% Purpose:
%   Validates beta reconstruction on an open-loop ABM trajectory with a
%   known transmission schedule.
%
% Outputs:
%   Diagnostic reconstruction figures and
%   integration/results_validation.mat.
%
% Hospital-load units:
%   H_load = (H1_agents/N1)/H_cap1 + (H2_agents/N2)/H_cap2.

addpath('macro');
addpath('micro');
addpath('integration');

SIRVDB_model_parameters;

%% Settings
T     = 150;
N1    = 500;
N2    = 500;
N_pop = N1 + N2;

rec_window      = 14;
trust_threshold = 5;

%% Known beta schedule
beta_true1 = zeros(T, 1);
beta_true2 = zeros(T, 1);

for t = 1:T
    sf = 1 + season_amp * sin(2 * pi * t / season_period);
    if t <= 40
        b1 = beta1 * sf;
        b2 = beta2 * sf;
    elseif t <= 90
        b1 = beta1 * (1 - 0.6*k) * sf;
        b2 = beta2 * (1 - 0.6*k) * sf;
    else
        b1 = beta1 * (1 - 0.3*k) * sf;
        b2 = beta2 * (1 - 0.3*k) * sf;
    end
    beta_true1(t) = max(0, b1);
    beta_true2(t) = max(0, b2);
end

%% Storage
val_S   = zeros(T, 1); val_I   = zeros(T, 1);
val_Iv  = zeros(T, 1); val_R   = zeros(T, 1);
val_V   = zeros(T, 1); val_H   = zeros(T, 1);
val_D   = zeros(T, 1);

val_c1_S  = zeros(T,1); val_c1_I  = zeros(T,1);
val_c1_Iv = zeros(T,1); val_c1_R  = zeros(T,1);
val_c1_V  = zeros(T,1); val_c1_H  = zeros(T,1);
val_c1_D  = zeros(T,1);

val_c2_S  = zeros(T,1); val_c2_I  = zeros(T,1);
val_c2_Iv = zeros(T,1); val_c2_R  = zeros(T,1);
val_c2_V  = zeros(T,1); val_c2_H  = zeros(T,1);
val_c2_D  = zeros(T,1);

%% Open-loop ABM with known beta
agents = init_agents(N1, N2, 10, 8, 0, 0);
fprintf('Running open-loop ABM with known beta schedule...\n');

for t = 1:T
    agents = update_agents(agents, beta_true1(t), beta_true2(t), params);
    stats  = aggregate_stats(agents);

    val_S(t)  = stats.S;   val_I(t)  = stats.I;
    val_Iv(t) = stats.Iv;  val_R(t)  = stats.R;
    val_V(t)  = stats.V;   val_H(t)  = stats.H;
    val_D(t)  = stats.D;

    val_c1_S(t)  = stats.city1.S;    val_c1_I(t)  = stats.city1.I_total;
    val_c1_Iv(t) = stats.city1.Iv;   val_c1_R(t)  = stats.city1.R;
    val_c1_V(t)  = stats.city1.V;    val_c1_H(t)  = stats.city1.H;
    val_c1_D(t)  = stats.city1.D;

    val_c2_S(t)  = stats.city2.S;    val_c2_I(t)  = stats.city2.I_total;
    val_c2_Iv(t) = stats.city2.Iv;   val_c2_R(t)  = stats.city2.R;
    val_c2_V(t)  = stats.city2.V;    val_c2_H(t)  = stats.city2.H;
    val_c2_D(t)  = stats.city2.D;
end

fprintf('ABM run complete.\n');

%% Hospital-load feedback signal
H_load1 = (val_c1_H / N1) / H_cap1;
H_load2 = (val_c2_H / N2) / H_cap2;
H_load  = H_load1 + H_load2;

fprintf('\nH_load peak:  %.2f\n', max(H_load));
fprintf('H_ref_load:   %.2f\n',   H_ref_load);
fprintf('Capacity:     %.2f\n',   1/H_cap1 + 1/H_cap2);

%% Beta reconstruction
partial.S        = val_S;    partial.I   = val_I;
partial.Iv       = val_Iv;   partial.R   = val_R;
partial.V        = val_V;    partial.H   = val_H;
partial.D        = val_D;

partial.city1_S  = val_c1_S;  partial.city1_I  = val_c1_I;
partial.city1_Iv = val_c1_Iv; partial.city1_R  = val_c1_R;
partial.city1_V  = val_c1_V;  partial.city1_H  = val_c1_H;
partial.city1_D  = val_c1_D;

partial.city2_S  = val_c2_S;  partial.city2_I  = val_c2_I;
partial.city2_Iv = val_c2_Iv; partial.city2_R  = val_c2_R;
partial.city2_V  = val_c2_V;  partial.city2_H  = val_c2_H;
partial.city2_D  = val_c2_D;

[b1_raw, b2_raw, b1_smooth, b2_smooth] = reconstruct_beta(partial, params, rec_window);
fprintf('Reconstruction complete.\n');

%% Error metrics
days = (1:T)';
trusted1 = val_c1_I >= trust_threshold;
trusted2 = val_c2_I >= trust_threshold;

err1_raw    = abs(b1_raw    - beta_true1);
err2_raw    = abs(b2_raw    - beta_true2);
err1_smooth = abs(b1_smooth - beta_true1);
err2_smooth = abs(b2_smooth - beta_true2);

rel_err1 = err1_smooth(trusted1) ./ max(beta_true1(trusted1), 1e-6);
rel_err2 = err2_smooth(trusted2) ./ max(beta_true2(trusted2), 1e-6);

fprintf('\n--- Reconstruction accuracy (trusted days only) ---\n');
fprintf('City 1  MAE raw:    %.4f\n', mean(err1_raw(trusted1),    'omitnan'));
fprintf('City 1  MAE smooth: %.4f\n', mean(err1_smooth(trusted1), 'omitnan'));
fprintf('City 1  MAPE:       %.1f%%\n', 100*mean(rel_err1,         'omitnan'));
fprintf('City 2  MAE raw:    %.4f\n', mean(err2_raw(trusted2),    'omitnan'));
fprintf('City 2  MAE smooth: %.4f\n', mean(err2_smooth(trusted2), 'omitnan'));
fprintf('City 2  MAPE:       %.1f%%\n', 100*mean(rel_err2,         'omitnan'));
fprintf('--------------------------------------------------\n\n');

%% Plots

% Beta reconstruction.
figure('Name', 'Beta reconstruction validation', 'Color', 'w');
subplot(2,1,1);
plot(days, beta_true1, 'k-',  'LineWidth', 2.0, 'DisplayName', 'True'); hold on;
plot(days, b1_raw,     'b:',  'LineWidth', 1.2, 'DisplayName', 'Raw');
plot(days, b1_smooth,  'b-',  'LineWidth', 1.8, 'DisplayName', 'Smoothed');
fill_untrusted(days, ~trusted1);
xline(40, '--', 'Color', [0.6 0.6 0.6]);
xline(90, '--', 'Color', [0.6 0.6 0.6]);
ylabel('\beta City 1'); legend('Location','best'); grid on;
title(sprintf('Beta reconstruction (window=%d days)', rec_window));

subplot(2,1,2);
plot(days, beta_true2, 'k-',  'LineWidth', 2.0, 'DisplayName', 'True'); hold on;
plot(days, b2_raw,     'r:',  'LineWidth', 1.2, 'DisplayName', 'Raw');
plot(days, b2_smooth,  'r-',  'LineWidth', 1.8, 'DisplayName', 'Smoothed');
fill_untrusted(days, ~trusted2);
xline(40, '--', 'Color', [0.6 0.6 0.6]);
xline(90, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Days'); ylabel('\beta City 2');
legend('Location','best'); grid on;

% Reconstruction error.
figure('Name', 'Reconstruction error', 'Color', 'w');
subplot(2,1,1);
plot(days, err1_raw,    'b:', 'LineWidth', 1.2, 'DisplayName', 'Raw'); hold on;
plot(days, err1_smooth, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Smoothed');
fill_untrusted(days, ~trusted1);
ylabel('|\Delta\beta| City 1'); legend('Location','best'); grid on;
title('Absolute reconstruction error');

subplot(2,1,2);
plot(days, err2_raw,    'r:', 'LineWidth', 1.2, 'DisplayName', 'Raw'); hold on;
plot(days, err2_smooth, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Smoothed');
fill_untrusted(days, ~trusted2);
xlabel('Days'); ylabel('|\Delta\beta| City 2');
legend('Location','best'); grid on;

% Epidemic curves.
figure('Name', 'Epidemic curves (validation)', 'Color', 'w');
plot(days, val_I + val_Iv, 'LineWidth', 1.8); hold on;
plot(days, val_H,          'LineWidth', 1.8);
plot(days, val_D,          'LineWidth', 1.8);
xline(40, '--', 'Phase 2', 'Color', [0.6 0.6 0.6], 'LabelVerticalAlignment','bottom');
xline(90, '--', 'Phase 3', 'Color', [0.6 0.6 0.6], 'LabelVerticalAlignment','bottom');
xlabel('Days'); ylabel('Agents');
legend('Total infected','Hospitalized','Dead','Location','best');
title('Epidemic curves — validation (open-loop, known \beta)'); grid on;

% Hospital load in Simulink-equivalent units.
figure('Name', 'H_load validation', 'Color', 'w');
plot(days, H_load, 'LineWidth', 1.8); hold on;
yline(H_ref_load, '--r', sprintf('H_{ref\\_load} = %.1f', H_ref_load), 'LineWidth', 1.5);
yline(1/H_cap1 + 1/H_cap2, '--k', 'Max capacity', 'LineWidth', 1.0);
xlabel('Days'); ylabel('H\_load (ODE units)');
title('Hospital load — validation run (matches Simulink scale)');
legend('H\_load','Reference','Max capacity','Location','best'); grid on;

subplot_city_loads(days, H_load1, H_load2, H_ref_load/2);

% Trust threshold.
figure('Name', 'Trust threshold', 'Color', 'w');
subplot(2,1,1);
plot(days, val_c1_I, 'LineWidth', 1.8); hold on;
yline(trust_threshold, '--r', sprintf('Threshold=%d', trust_threshold), 'LineWidth',1.2);
ylabel('I_{total} City 1'); title('Infected count vs trust threshold'); grid on;
subplot(2,1,2);
plot(days, val_c2_I, 'LineWidth', 1.8); hold on;
yline(trust_threshold, '--r', sprintf('Threshold=%d', trust_threshold), 'LineWidth',1.2);
xlabel('Days'); ylabel('I_{total} City 2'); grid on;

%% Save
save('integration/results_validation.mat', ...
    'beta_true1','beta_true2','b1_raw','b2_raw','b1_smooth','b2_smooth', ...
    'err1_raw','err2_raw','err1_smooth','err2_smooth', ...
    'val_S','val_I','val_Iv','val_R','val_V','val_H','val_D', ...
    'val_c1_I','val_c2_I','val_c1_H','val_c2_H', ...
    'H_load','H_load1','H_load2','H_ref_load', ...
    'trusted1','trusted2','days');

fprintf('Saved to integration/results_validation.mat\n');


%% Local helpers
function fill_untrusted(days, mask)
    if ~any(mask), return; end
    yl     = ylim;
    idx    = find(mask);
    breaks = [0; find(diff(idx)>1); length(idx)];
    for b = 1:length(breaks)-1
        block = idx(breaks(b)+1:breaks(b+1));
        patch([days(block(1)) days(block(end)) days(block(end)) days(block(1))], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
              'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility','off');
    end
end

function subplot_city_loads(days, H_load1, H_load2, H_ref_city)
    figure('Name', 'City H_load validation', 'Color', 'w');
    subplot(2,1,1);
    plot(days, H_load1, 'LineWidth', 1.8); hold on;
    yline(1.0,        '--k', 'Capacity',   'LineWidth', 1.0);
    yline(H_ref_city, ':r',  'H_{ref}',    'LineWidth', 1.2);
    ylabel('H load City 1'); title('Per-city hospital load'); grid on;
    subplot(2,1,2);
    plot(days, H_load2, 'LineWidth', 1.8); hold on;
    yline(1.0,        '--k', 'Capacity',   'LineWidth', 1.0);
    yline(H_ref_city, ':r',  'H_{ref}',    'LineWidth', 1.2);
    xlabel('Days'); ylabel('H load City 2'); grid on;
end
