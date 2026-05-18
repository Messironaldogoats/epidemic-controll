function run_visual_abm(beta_t1_series, beta_t2_series, u_series)
% RUN_VISUAL_ABM
% Real-time visualisation using the same ABM functions as the main
% closed-loop simulation.
%
% Inputs:
%   beta_t1_series,beta_t2_series - optional city transmission series
%   u_series                      - optional MPC intervention series
%
% Control:
%   Close the figure window to stop.

addpath('micro');
addpath('macro');

%% Parameters
run('macro/SIRVDB_model_parameters.m');

gamma_rec = evalin('caller','gamma');   % recovery rate from parameter file

params.nu1     = nu1;
params.nu2     = nu2;
params.eta     = eta;
params.gamma   = gamma_rec;
params.gamma_v = gamma_v;
params.gamma_h = gamma_h;
params.mu      = mu;
params.mu_v    = mu_v;
params.mu_h    = mu_h;
params.h       = h;
params.h_v     = h_v;
params.omega_V = omega_V;
params.omega_R = omega_R;
params.k       = k;
params.M11     = 0.85;
params.M12     = 0.15;
params.M21     = 0.10;
params.M22     = 0.90;

% Small visual population for real-time rendering.
N1 = 150;
N2 = 150;
T  = 150;

%% Default inputs
if nargin < 1 || isempty(beta_t1_series)
    beta_t1_series = zeros(T,1);
    beta_t2_series = zeros(T,1);
    for t = 1:T
        sf = 1 + season_amp * sin(2*pi*t/season_period);
        beta_t1_series(t) = beta1 * sf;
        beta_t2_series(t) = beta2 * sf;
    end
end

if nargin < 3 || isempty(u_series)
    u_series = zeros(T,1);
end

T = min([T, length(beta_t1_series), length(beta_t2_series)]);

%% Initialise agents
rng(42);
agents = init_agents(N1, N2, 3, 2, 0, 0, params);
N = N1 + N2;

fprintf('Agents initialised via init_agents_v2.m\n');
fprintf('  N=%d  (City1=%d, City2=%d)\n', N, N1, N2);

%% Visual positions
% Display positions do not affect ABM state transitions.
city1_box = [0.03, 0.15, 0.37, 0.70];
city2_box = [0.60, 0.15, 0.37, 0.70];
bridge_x  = [0.40, 0.60];
bridge_y  = [0.43, 0.57];
hosp1_pos = [0.09, 0.76];
hosp2_pos = [0.88, 0.76];

px = zeros(N,1);
py = zeros(N,1);
vx = zeros(N,1);
vy = zeros(N,1);

c1 = (agents.city == 1);
c2 = (agents.city == 2);

px(c1) = city1_box(1) + city1_box(3)*(0.1 + 0.8*rand(sum(c1),1));
py(c1) = city1_box(2) + city1_box(4)*(0.1 + 0.8*rand(sum(c1),1));
px(c2) = city2_box(1) + city2_box(3)*(0.1 + 0.8*rand(sum(c2),1));
py(c2) = city2_box(2) + city2_box(4)*(0.1 + 0.8*rand(sum(c2),1));

vx = 0.004*(rand(N,1)-0.5);
vy = 0.004*(rand(N,1)-0.5);

traveling     = false(N,1);
travel_target = zeros(N,1);
travel_timer  = zeros(N,1);

%% Colours
COLS = [
    0.20 0.45 1.00;
    0.95 0.12 0.12;
    0.90 0.35 0.80;
    0.10 0.72 0.20;
    0.60 0.25 0.90;
    1.00 0.55 0.10;
    0.35 0.35 0.35;
];

LABELS = {'S susceptible','I infected','Iv vacc.inf.', ...
          'R recovered','V vaccinated','H hospitalised','D dead'};

%% Storage
hS  = nan(T,1); hI  = nan(T,1); hIv = nan(T,1);
hR  = nan(T,1); hV  = nan(T,1); hH  = nan(T,1); hD  = nan(T,1);
hI1 = nan(T,1); hI2 = nan(T,1);
hH1 = nan(T,1); hH2 = nan(T,1);

%% Figure
fig = figure('Name','Live ABM — micro model v2', ...
             'Color','w','Position',[30 30 1440 740],'NumberTitle','off');

ax_map   = axes('Parent',fig,'Position',[0.01 0.28 0.57 0.69]);
ax_curve = axes('Parent',fig,'Position',[0.63 0.54 0.35 0.42]);
ax_city  = axes('Parent',fig,'Position',[0.63 0.07 0.35 0.39]);
ax_dash  = axes('Parent',fig,'Position',[0.01 0.01 0.57 0.24]);

axis(ax_map,  [0 1 0 1]); axis(ax_map,  'off'); hold(ax_map,  'on');
axis(ax_dash, [0 1 0 1]); axis(ax_dash, 'off'); hold(ax_dash, 'on');
hold(ax_curve,'on'); grid(ax_curve,'on');
hold(ax_city, 'on'); grid(ax_city, 'on');

%% Static map elements
patch(ax_map,[0 1 1 0],[0 0 1 1],[0.93 0.95 0.92],'EdgeColor','none');

patch(ax_map, ...
    city1_box(1)+[0 city1_box(3) city1_box(3) 0], ...
    city1_box(2)+[0 0 city1_box(4) city1_box(4)], ...
    [0.87 0.92 1.00],'EdgeColor',[0.10 0.20 0.55],'LineWidth',1.8);

patch(ax_map, ...
    city2_box(1)+[0 city2_box(3) city2_box(3) 0], ...
    city2_box(2)+[0 0 city2_box(4) city2_box(4)], ...
    [1.00 0.91 0.87],'EdgeColor',[0.55 0.16 0.10],'LineWidth',1.8);

patch(ax_map, ...
    [bridge_x(1) bridge_x(2) bridge_x(2) bridge_x(1)], ...
    [bridge_y(1) bridge_y(1) bridge_y(2) bridge_y(2)], ...
    [0.79 0.79 0.79],'EdgeColor',[0.50 0.50 0.50],'LineWidth',1.0);

for bx = linspace(bridge_x(1)+0.022, bridge_x(2)-0.022, 5)
    plot(ax_map,[bx bx],bridge_y,'w--','LineWidth',0.8);
end

text(ax_map, mean(bridge_x), mean(bridge_y)+0.04, 'Bridge', ...
     'FontSize',9,'HorizontalAlignment','center','Color',[0.3 0.3 0.3],'FontWeight','bold');

draw_hosp(ax_map, hosp1_pos, 'Hospital 1');
draw_hosp(ax_map, hosp2_pos, 'Hospital 2');

text(ax_map, city1_box(1)+city1_box(3)/2, city1_box(2)+city1_box(4)+0.025, ...
     'City 1','FontSize',13,'FontWeight','bold','HorizontalAlignment','center','Color',[0.10 0.18 0.50]);
text(ax_map, city2_box(1)+city2_box(3)/2, city2_box(2)+city2_box(4)+0.025, ...
     'City 2','FontSize',13,'FontWeight','bold','HorizontalAlignment','center','Color',[0.50 0.14 0.10]);

for li = 1:7
    scatter(ax_map, 0.013, 0.097-(li-1)*0.013, 36, COLS(li,:), ...
            'filled','MarkerEdgeColor','none');
    text(ax_map, 0.027, 0.097-(li-1)*0.013, LABELS{li}, ...
         'FontSize',8,'Color',[0.2 0.2 0.2],'VerticalAlignment','middle');
end

h_agents = scatter(ax_map, px, py, 22, COLS(agents.state,:), ...
                   'filled','MarkerEdgeColor','none');

h_title = title(ax_map, 'Day 0  —  init_agents_v2 + update_agents_v2 + aggregate_stats_v2', ...
                'FontSize',10,'FontWeight','bold');

%% Dashboard
bw = 0.125;
h_counts = gobjects(7,1);
snames   = {'S','I','Iv','R','V','H','D'};

for si = 1:7
    bx = (si-1)*bw + 0.008;
    patch(ax_dash, bx+[0 bw-0.01 bw-0.01 0],[0.06 0.06 0.94 0.94], ...
          min(COLS(si,:)*0.28+0.72,1), 'EdgeColor',COLS(si,:)*0.55,'LineWidth',0.8);
    text(ax_dash, bx+bw/2-0.005, 0.74, snames{si}, ...
         'FontSize',12,'FontWeight','bold','HorizontalAlignment','center', ...
         'Color',COLS(si,:)*0.50);
    h_counts(si) = text(ax_dash, bx+bw/2-0.005, 0.34, '0', ...
         'FontSize',14,'FontWeight','bold','HorizontalAlignment','center', ...
         'Color',min(COLS(si,:)*0.65,1));
end

set(h_counts(1),'String',num2str(N));

h_info = text(ax_dash, 0.50, 0.04, '', ...
    'FontSize',10,'HorizontalAlignment','center','Color',[0.20 0.20 0.20]);

%% Main loop
fprintf('Running... close window to stop.\n');

for t = 1:T

    if ~isvalid(fig), break; end

    b1 = beta_t1_series(min(t,end));
    b2 = beta_t2_series(min(t,end));
    u  = u_series(min(t,end));

    %% ABM update
    agents = update_agents(agents, b1, b2, params, u);
    stats  = aggregate_stats(agents);

    %% Store
    hS(t)=stats.S; hI(t)=stats.I; hIv(t)=stats.Iv;
    hR(t)=stats.R; hV(t)=stats.V; hH(t)=stats.H; hD(t)=stats.D;
    hI1(t)=stats.city1.I_total; hI2(t)=stats.city2.I_total;
    hH1(t)=stats.city1.H;       hH2(t)=stats.city2.H;

    %% Update visual positions
    move_sc = 0.005 * max(0.15, 1 - 0.80*u);

    for i = 1:N
        s = agents.state(i);
        if s == 7, continue; end

        if s == 6
            tgt = hosp1_pos;
            if agents.city(i)==2, tgt = hosp2_pos; end
            px(i) = px(i) + 0.12*(tgt(1)-px(i)) + 0.002*randn;
            py(i) = py(i) + 0.12*(tgt(2)-py(i)) + 0.002*randn;
            continue;
        end

        if traveling(i)
            travel_timer(i) = travel_timer(i) - 1;
            if travel_target(i)==1
                tx = city1_box(1)+city1_box(3)/2;
            else
                tx = city2_box(1)+city2_box(3)/2;
            end
            px(i) = px(i) + 0.07*(tx-px(i)) + 0.004*randn;
            py(i) = py(i) + 0.07*(mean(bridge_y)-py(i)) + 0.003*randn;
            if travel_timer(i) <= 0
                traveling(i) = false;
                if travel_target(i)==1, b=city1_box; else, b=city2_box; end
                px(i) = b(1)+b(3)*(0.1+0.8*rand);
                py(i) = b(2)+b(4)*(0.1+0.8*rand);
            end
            continue;
        end

        vx(i) = 0.65*vx(i) + move_sc*randn;
        vy(i) = 0.65*vy(i) + move_sc*randn;
        px(i) = px(i) + vx(i);
        py(i) = py(i) + vy(i);

        if agents.city(i)==1, b=city1_box; else, b=city2_box; end
        px(i) = max(b(1)+0.01, min(b(1)+b(3)-0.01, px(i)));
        py(i) = max(b(2)+0.01, min(b(2)+b(4)-0.01, py(i)));
        if px(i)<=b(1)+0.015, vx(i)=abs(vx(i)); end
        if px(i)>=b(1)+b(3)-0.015, vx(i)=-abs(vx(i)); end
        if py(i)<=b(2)+0.015, vy(i)=abs(vy(i)); end
        if py(i)>=b(2)+b(4)-0.015, vy(i)=-abs(vy(i)); end

        p_tr = 0.009 * max(0.05, 1-0.90*u);
        if rand < p_tr
            traveling(i)     = true;
            travel_target(i) = 3 - agents.city(i);
            travel_timer(i)  = 7 + randi(7);
        end
    end

    %% Update map
    ag_colors = COLS(agents.state,:);
    ag_sizes  = 20*ones(N,1);
    ag_sizes(agents.state==2) = 44;
    ag_sizes(agents.state==3) = 34;
    ag_sizes(agents.state==6) = 38;
    ag_sizes(agents.state==7) = 7;

    set(h_agents,'XData',px,'YData',py,'CData',ag_colors,'SizeData',ag_sizes);
    set(h_title,'String', ...
        sprintf('Day %d / %d   |   update\\_agents\\_v2()  +  aggregate\\_stats\\_v2()', t, T));

    %% Update dashboard
    cnts = [stats.S, stats.I, stats.Iv, stats.R, stats.V, stats.H, stats.D];
    for si = 1:7
        set(h_counts(si),'String',num2str(cnts(si)));
    end

    if u<0.05, us='no restriction';
    elseif u<0.40, us='mild';
    elseif u<0.70, us='moderate lockdown';
    else, us='FULL LOCKDOWN'; end

    set(h_info,'String', ...
        sprintf('\\beta_1=%.3f   \\beta_2=%.3f   |   u=%.3f (%s)   |   travelers on bridge=%d', ...
        b1, b2, u, us, sum(traveling)));

    %% Update epidemic curve
if t >= 2
d = (1:t)';
cla(ax_curve); hold(ax_curve,'on');
area(ax_curve, d, [hS(1:t), hI(1:t)+hIv(1:t), hH(1:t), hR(1:t), hD(1:t)], ...
     'LineWidth',0.4);
colororder(ax_curve, [COLS(1,:);COLS(2,:);COLS(6,:);COLS(4,:);COLS(7,:)]);
yline(ax_curve, H_cap1*N1+H_cap2*N2,'--k','Hosp. cap','LineWidth',1.2,'FontSize',8);
xlim(ax_curve,[1 T]); ylim(ax_curve,[0 N]);
xlabel(ax_curve,'Days','FontSize',9);
ylabel(ax_curve,'Agents','FontSize',9);
title(ax_curve,sprintf('Day %d — population state',t),'FontSize',10);
legend(ax_curve,'S','I+Iv','H','R','D','Hosp.cap', ...
       'Location','eastoutside','FontSize',8);
grid(ax_curve,'on');
end

  %% Update city curves
if t >= 2
cla(ax_city); hold(ax_city,'on');
plot(ax_city,1:t,hI1(1:t),'b-','LineWidth',2.0, ...
     'DisplayName',sprintf('City 1  I=%d',hI1(t)));
plot(ax_city,1:t,hI2(1:t),'r-','LineWidth',2.0, ...
     'DisplayName',sprintf('City 2  I=%d',hI2(t)));
plot(ax_city,1:t,hH1(1:t),'b--','LineWidth',1.2, ...
     'DisplayName',sprintf('City 1  H=%d',hH1(t)));
plot(ax_city,1:t,hH2(1:t),'r--','LineWidth',1.2, ...
     'DisplayName',sprintf('City 2  H=%d',hH2(t)));
yline(ax_city,H_cap1*N1,':k','LineWidth',1.0,'HandleVisibility','off');
xlim(ax_city,[1 T]);
xlabel(ax_city,'Days','FontSize',9);
ylabel(ax_city,'Agents','FontSize',9);
title(ax_city,'City 1 vs City 2 — stochastic divergence','FontSize',10);
legend(ax_city,'Location','best','FontSize',8);
grid(ax_city,'on');
end

    drawnow limitrate;
    pause(0.03);
end

fprintf('\nDone.  D=%d  H=%d\n', stats.D, stats.H);
end


function draw_hosp(ax, pos, label)
    cx=pos(1); cy=pos(2); w=0.040; h=0.055;
    patch(ax,cx+[-w w w -w],cy+[0 0 h h],[1.0 1.0 0.84], ...
          'EdgeColor',[0.4 0.4 0.4],'LineWidth',1.0);
    plot(ax,[cx cx],[cy+h*0.25 cy+h*0.85],'r-','LineWidth',2.5);
    plot(ax,[cx-w*0.5 cx+w*0.5],[cy+h*0.55 cy+h*0.55],'r-','LineWidth',2.5);
    text(ax,cx,cy-0.024,label,'FontSize',8,'HorizontalAlignment','center', ...
         'FontWeight','bold','Color',[0.15 0.15 0.55]);
end
