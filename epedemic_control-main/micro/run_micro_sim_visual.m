function results =









series, beta_t2_series)
% Visual agent-based epidemic simulation with inter-city movement and a
% dashboard for compartment counts.
%
% States:
% 1 = S (susceptible)
% 2 = I (infected)
% 3 = R (recovered)
% 4 = V (vaccinated)
% 5 = H (hospitalized)
% 6 = D (dead)

rng(1);

%% Settings
T = min(length(beta_t1_series), length(beta_t2_series));

N1 = 90;
N2 = 90;
N = N1 + N2;

move_step_normal = 1.0;
move_step_hospital = 0.08;
infection_radius = 2.8;

% Disease progression probabilities.
p_recover = 0.05;
p_hospital = 0.015;
p_die_from_I = 0.002;
p_die_from_H = 0.01;
p_leave_hospital = 0.06;

vacc_fraction = 0.15;

% Inter-city movement probability per day.
p_commute = 0.015;

contact_scale = 0.35;

vacc_protection = 0.35;

% City geometry for visualisation.
city1_box = [0, 0, 45, 40];
city2_box = [55, 0, 45, 40];

city1_center = [22.5, 20];
city2_center = [77.5, 20];

hospital1 = [8, 33];
hospital2 = [92, 33];

road_y = [17, 23];

%% Colours
cS = [0.20 0.45 1.00];
cI = [1.00 0.10 0.10];
cR = [0.10 0.72 0.20];
cV = [0.62 0.25 0.90];
cH = [1.00 0.58 0.10];
cD = [0.00 0.00 0.00];

%% Initialise agents
x = zeros(N,1);
y = zeros(N,1);
state = ones(N,1);
city = zeros(N,1);

home_city = zeros(N,1);
travel_days_left = zeros(N,1);

x(1:N1) = city1_box(1) + city1_box(3)*rand(N1,1);
y(1:N1) = city1_box(2) + city1_box(4)*rand(N1,1);
city(1:N1) = 1;
home_city(1:N1) = 1;

x(N1+1:end) = city2_box(1) + city2_box(3)*rand(N2,1);
y(N1+1:end) = city2_box(2) + city2_box(4)*rand(N2,1);
city(N1+1:end) = 2;
home_city(N1+1:end) = 2;

vacc_idx = randperm(N, round(vacc_fraction*N));
state(vacc_idx) = 4;

available = find(state == 1);
init_inf = available(randperm(length(available), 8));
state(init_inf) = 2;

%% Storage
results.time = (1:T)';
results.S = zeros(T,1);
results.I = zeros(T,1);
results.Iv = zeros(T,1);
results.R = zeros(T,1);
results.V = zeros(T,1);
results.H = zeros(T,1);
results.D = zeros(T,1);

results.city1_I = zeros(T,1);
results.city2_I = zeros(T,1);
results.city1_H = zeros(T,1);
results.city2_H = zeros(T,1);

%% Figure
fig = figure('Name','Micro epidemic simulation','Color','w');
set(fig, 'Position', [80 80 1450 650]);

for t = 1:T
    %% Inter-city movement
    mobile = (state ~= 5 & state ~= 6);

    for i = 1:N
        if ~mobile(i)
            continue;
        end

        if travel_days_left(i) > 0
            travel_days_left(i) = travel_days_left(i) - 1;

            if travel_days_left(i) == 0
                city(i) = home_city(i);
            end
        else
            if rand < p_commute
                if city(i) == 1
                    city(i) = 2;
                else
                    city(i) = 1;
                end
                travel_days_left(i) = randi([2 5]);
            end
        end
    end

    %% Move agents
    dx = zeros(N,1);
    dy = zeros(N,1);

    for i = 1:N
        if state(i) == 6
            continue;
        elseif state(i) == 5
            step = move_step_hospital;
            if city(i) == 1
                target = hospital1;
            else
                target = hospital2;
            end
        else
            step = move_step_normal;
            if city(i) == 1
                target = city1_center;
            else
                target = city2_center;
            end
        end

        drift = 0.08 * ([target(1) - x(i), target(2) - y(i)]);
        random_move = step * randn(1,2);
        move_vec = random_move + drift;

        dx(i) = move_vec(1);
        dy(i) = move_vec(2);
    end

    x = x + dx;
    y = y + dy;

    idx1 = city == 1;
    idx2 = city == 2;

    x(idx1) = max(city1_box(1), min(city1_box(1)+city1_box(3), x(idx1)));
    y(idx1) = max(city1_box(2), min(city1_box(2)+city1_box(4), y(idx1)));

    x(idx2) = max(city2_box(1), min(city2_box(1)+city2_box(3), x(idx2)));
    y(idx2) = max(city2_box(2), min(city2_box(2)+city2_box(4), y(idx2)));

    %% Infection
    new_infected = false(N,1);
    infected_idx = find(state == 2);

    for k = 1:length(infected_idx)
        i = infected_idx(k);

        if city(i) == 1
            beta_now = beta_t1_series(min(t, length(beta_t1_series)));
        else
            beta_now = beta_t2_series(min(t, length(beta_t2_series)));
        end

        beta_now = max(0, beta_now);

        susceptible_candidates = find((state == 1 | state == 4) & city == city(i));
        if isempty(susceptible_candidates)
            continue;
        end

        dist = sqrt((x(susceptible_candidates)-x(i)).^2 + ...
                    (y(susceptible_candidates)-y(i)).^2);

        near_mask = dist < infection_radius;
        neighbors = susceptible_candidates(near_mask);
        neighbor_dist = dist(near_mask);

        for j = 1:length(neighbors)
            n = neighbors(j);
            d = neighbor_dist(j);

            proximity_factor = max(0, 1 - d/infection_radius);
            p_infect = contact_scale * beta_now * proximity_factor;

            if state(n) == 4
                p_infect = p_infect * vacc_protection;
            end

            p_infect = min(max(p_infect, 0), 1);

            if rand < p_infect
                new_infected(n) = true;
            end
        end
    end

    state(new_infected) = 2;

    %% Disease progression
    infected_idx = find(state == 2);
    for k = 1:length(infected_idx)
        i = infected_idx(k);
        r = rand;

        if r < p_die_from_I
            state(i) = 6;
        elseif r < p_die_from_I + p_hospital
            state(i) = 5;
        elseif r < p_die_from_I + p_hospital + p_recover
            state(i) = 3;
        end
    end

    hosp_idx = find(state == 5);
    for k = 1:length(hosp_idx)
        i = hosp_idx(k);
        r = rand;

        if r < p_die_from_H
            state(i) = 6;
        elseif r < p_die_from_H + p_leave_hospital
            state(i) = 3;
        end
    end

    %% Store results
    results.S(t) = sum(state == 1);
    results.I(t) = sum(state == 2);
    results.Iv(t) = 0;
    results.R(t) = sum(state == 3);
    results.V(t) = sum(state == 4);
    results.H(t) = sum(state == 5);
    results.D(t) = sum(state == 6);

    results.city1_I(t) = sum(state == 2 & city == 1);
    results.city2_I(t) = sum(state == 2 & city == 2);
    results.city1_H(t) = sum(state == 5 & city == 1);
    results.city2_H(t) = sum(state == 5 & city == 2);

    %% Visualisation
    clf(fig);

    ax = axes('Parent', fig, 'Position', [0.05 0.10 0.72 0.82]);
    hold(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [-3 103]);
    ylim(ax, [-3 43]);
    box(ax, 'on');

    patch(ax, [-3 103 103 -3], [-3 -3 43 43], [0.95 0.97 0.95], 'EdgeColor', 'none');

    rectangle(ax, 'Position', city1_box, ...
        'FaceColor', [0.88 0.93 1.00], 'EdgeColor', [0.15 0.25 0.5], 'LineWidth', 2);

    rectangle(ax, 'Position', city2_box, ...
        'FaceColor', [1.00 0.92 0.88], 'EdgeColor', [0.5 0.2 0.15], 'LineWidth', 2);

    patch(ax, [45 55 55 45], [road_y(1) road_y(1) road_y(2) road_y(2)], ...
        [0.75 0.75 0.75], 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.5);
    plot(ax, [45 55], [20 20], '--w', 'LineWidth', 1.2);

    plot(ax, hospital1(1), hospital1(2), 'ks', 'MarkerSize', 12, ...
        'MarkerFaceColor', [1 1 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, hospital2(1), hospital2(2), 'ks', 'MarkerSize', 12, ...
        'MarkerFaceColor', [1 1 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');

    text(ax, 16, 41, sprintf('City 1   \\beta=%.3f', beta_t1_series(min(t,end))), ...
        'FontSize', 11, 'FontWeight', 'bold');
    text(ax, 70, 41, sprintf('City 2   \\beta=%.3f', beta_t2_series(min(t,end))), ...
        'FontSize', 11, 'FontWeight', 'bold');

    text(ax, hospital1(1)-2.2, hospital1(2)+2.1, 'Hospital', 'FontSize', 9, 'FontWeight', 'bold');
    text(ax, hospital2(1)-2.2, hospital2(2)+2.1, 'Hospital', 'FontSize', 9, 'FontWeight', 'bold');

    plot_agents(ax, x, y, state, 1, 'o', cS);
    plot_agents(ax, x, y, state, 2, 'o', cI);
    plot_agents(ax, x, y, state, 3, 'o', cR);
    plot_agents(ax, x, y, state, 4, 'o', cV);
    plot_agents(ax, x, y, state, 5, 'o', cH);
    plot_agents(ax, x, y, state, 6, 'x', cD);

    title(ax, sprintf('Micro epidemic simulation - Day %d', t), ...
        'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax, 'X position');
    ylabel(ax, 'Y position');

    %% Dashboard panel
    axInfo = axes('Parent', fig, 'Position', [0.80 0.10 0.18 0.82]);
    axis(axInfo, [0 1 0 1]);
    axis(axInfo, 'off');
    hold(axInfo, 'on');

    rectangle(axInfo, 'Position', [0.02 0.02 0.96 0.96], ...
        'FaceColor', [0.99 0.99 0.99], ...
        'EdgeColor', [0.82 0.82 0.82], ...
        'LineWidth', 1.2, ...
        'Curvature', 0.02);

    text(axInfo, 0.08, 0.94, 'Simulation Dashboard', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15]);

    text(axInfo, 0.08, 0.89, sprintf('Day %d', t), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.25 0.25 0.25]);

    plot(axInfo, [0.08 0.92], [0.86 0.86], 'Color', [0.85 0.85 0.85], 'LineWidth', 1);

    text(axInfo, 0.08, 0.82, 'States', ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    y0 = 0.77;
    dy = 0.06;
    draw_status_item(axInfo, 0.10, y0 - 0*dy, cS, 'o', 'Susceptible',   results.S(t));
    draw_status_item(axInfo, 0.10, y0 - 1*dy, cI, 'o', 'Infected',      results.I(t));
    draw_status_item(axInfo, 0.10, y0 - 2*dy, cR, 'o', 'Recovered',     results.R(t));
    draw_status_item(axInfo, 0.10, y0 - 3*dy, cV, 'o', 'Vaccinated',    results.V(t));
    draw_status_item(axInfo, 0.10, y0 - 4*dy, cH, 'o', 'Hospitalized',  results.H(t));
    draw_status_item(axInfo, 0.10, y0 - 5*dy, cD, 'x', 'Dead',          results.D(t));

    plot(axInfo, [0.08 0.92], [0.39 0.39], 'Color', [0.85 0.85 0.85], 'LineWidth', 1);

    num_travelers = sum(travel_days_left > 0 & state ~= 6);

    text(axInfo, 0.08, 0.35, 'City breakdown', ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    text(axInfo, 0.10, 0.30, sprintf('City 1: I = %d, H = %d', ...
        results.city1_I(t), results.city1_H(t)), ...
        'FontSize', 10, 'Color', [0.15 0.15 0.15]);

    text(axInfo, 0.10, 0.25, sprintf('City 2: I = %d, H = %d', ...
        results.city2_I(t), results.city2_H(t)), ...
        'FontSize', 10, 'Color', [0.15 0.15 0.15]);

    text(axInfo, 0.10, 0.18, sprintf('Travelers: %d', num_travelers), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.25 0.25 0.25]);

    rectangle(axInfo, 'Position', [0.08 0.05 0.84 0.08], ...
        'FaceColor', [0.96 0.97 1.00], ...
        'EdgeColor', [0.82 0.86 0.95], ...
        'LineWidth', 1.0, ...
        'Curvature', 0.02);

    text(axInfo, 0.12, 0.095, sprintf('\\beta_1 = %.3f    \\beta_2 = %.3f', ...
        beta_t1_series(min(t,end)), beta_t2_series(min(t,end))), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.15 0.25 0.45]);

    drawnow;
    pause(0.06);
end
end

function plot_agents(ax, x, y, state, state_id, marker_style, color)
idx = state == state_id;

if ~any(idx)
    return;
end

if marker_style == "x" || strcmp(marker_style, 'x')
    scatter(ax, x(idx), y(idx), 70, color, marker_style, 'LineWidth', 1.5);
else
    scatter(ax, x(idx), y(idx), 45 + 15*(state_id==2) + 25*(state_id==5), ...
        color, marker_style, 'filled');
end
end

function draw_status_item(ax, x, y, color, marker_style, label, value)
if marker_style == "x" || strcmp(marker_style, 'x')
    plot(ax, x, y, marker_style, 'Color', color, 'MarkerSize', 10, 'LineWidth', 1.8);
else
    scatter(ax, x, y, 55, color, marker_style, 'filled');
end

text(ax, x + 0.07, y, sprintf('%s: %d', label, value), ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 10, ...
    'Color', [0.15 0.15 0.15]);
end
