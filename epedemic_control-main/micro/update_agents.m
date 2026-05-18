function agents = update_agents(agents, beta_t1, beta_t2, params, u_now)
% UPDATE_AGENTS
% Advances the two-city ABM by one day.
%
% Purpose:
%   Applies infection, vaccination, recovery, hospitalisation, mortality,
%   and waning-immunity transitions using heterogeneous contacts and
%   city-specific stochastic streams.
%
% Inputs:
%   agents         - ABM state struct from init_agents
%   beta_t1,beta_t2 - city transmission rates for the current day
%   params         - epidemiological and travel parameters
%   u_now          - intervention level in [0,1]
%
% Output:
%   agents         - updated ABM state struct
%
% Infection model:
%   p = 1 - exp(-lambda), with lambda adjusted by susceptibility,
%   contact count, intervention strength, and cross-city mixing.
%
% State codes: 1=S, 2=I, 3=Iv, 4=R, 5=V, 6=H, 7=D

%% Unpack parameters
gamma   = params.gamma;
gamma_v = params.gamma_v;
gamma_h = params.gamma_h;
mu      = params.mu;
mu_v    = params.mu_v;
mu_h    = params.mu_h;
h       = params.h;
h_v     = params.h_v;
omega_V = params.omega_V;
omega_R = params.omega_R;
eta     = params.eta;
nu1     = params.nu1;
nu2     = params.nu2;
k       = params.k;

M11 = params.M11;   % fraction of city1 contacts within city1
M12 = params.M12;   % fraction of city1 contacts with city2
M21 = params.M21;   % fraction of city2 contacts with city1
M22 = params.M22;   % fraction of city2 contacts within city2

%% Extract state vectors
state          = agents.state;
city           = agents.city;
susceptibility = agents.susceptibility;
n_contacts     = agents.n_contacts;
N              = length(state);

c1 = (city == 1);
c2 = (city == 2);

%% Infectious prevalence by city
I1_infectious = sum((state == 2 | state == 3) & c1);
I2_infectious = sum((state == 2 | state == 3) & c2);

N1_alive = sum(c1 & state ~= 7);
N2_alive = sum(c2 & state ~= 7);
N1_alive = max(N1_alive, 1);
N2_alive = max(N2_alive, 1);

f1 = I1_infectious / N1_alive;
f2 = I2_infectious / N2_alive;

%% Intervention-dependent contact reduction
contact_reduction = max(0, 1 - k * u_now);

%% City force of infection with travel mixing
lambda1_base = beta_t1 * (M11 * f1 + M12 * f2);
lambda2_base = beta_t2 * (M21 * f1 + M22 * f2);

%% State duration
agents.days_in_state(state ~= 7) = ...
    agents.days_in_state(state ~= 7) + 1;

%% Infection step
% lambda_i includes city force of infection, individual susceptibility,
% and relative contact intensity.

mean_c = mean(n_contacts);
n_eff  = max(1, round(n_contacts * contact_reduction));

lambda_agent = zeros(N, 1);
lambda_agent(c1) = lambda1_base * susceptibility(c1) .* ...
                   n_eff(c1) / mean_c;
lambda_agent(c2) = lambda2_base * susceptibility(c2) .* ...
                   n_eff(c2) / mean_c;

p_inf  = 1 - exp(-lambda_agent);
p_inf  = max(0, min(1, p_inf));

p_infv = 1 - exp(-eta * lambda_agent);
p_infv = max(0, min(1, p_infv));

%% State transitions
% Separate random streams preserve independent stochastic city dynamics.

N1_count = sum(c1);
N2_count = sum(c2);

r_c1 = rand(agents.stream1, N1_count, 6);
r_c2 = rand(agents.stream2, N2_count, 6);

r = zeros(N, 6);
r(c1, :) = r_c1;
r(c2, :) = r_c2;

new_state = state;

%% S to I or V
s_mask = (state == 1);

inf_mask = s_mask & (r(:,1) < p_inf);
new_state(inf_mask) = 2;

vacc_rate = zeros(N,1);
vacc_rate(c1 & ~inf_mask) = nu1;
vacc_rate(c2 & ~inf_mask) = nu2;
vacc_mask = s_mask & ~inf_mask & (r(:,2) < vacc_rate);
new_state(vacc_mask) = 5;

%% I to H, D, or R
i_mask = (state == 2);

hosp_mask = i_mask & (r(:,1) < h);
new_state(hosp_mask) = 6;

dead_mask = i_mask & ~hosp_mask & (r(:,2) < mu);
new_state(dead_mask) = 7;

rec_mask = i_mask & ~hosp_mask & ~dead_mask & (r(:,3) < gamma);
new_state(rec_mask) = 4;

%% Iv to H, D, or R
iv_mask = (state == 3);

hosp_v_mask = iv_mask & (r(:,1) < h_v);
new_state(hosp_v_mask) = 6;

dead_v_mask = iv_mask & ~hosp_v_mask & (r(:,2) < mu_v);
new_state(dead_v_mask) = 7;

rec_v_mask = iv_mask & ~hosp_v_mask & ~dead_v_mask & (r(:,3) < gamma_v);
new_state(rec_v_mask) = 4;

%% V to Iv or S
v_mask = (state == 5);

bk_mask = v_mask & (r(:,1) < p_infv);
new_state(bk_mask) = 3;

wane_mask = v_mask & ~bk_mask & (r(:,2) < omega_V);
new_state(wane_mask) = 1;

%% H to D or R
h_mask = (state == 6);

dead_h_mask = h_mask & (r(:,1) < mu_h);
new_state(dead_h_mask) = 7;

rec_h_mask = h_mask & ~dead_h_mask & (r(:,2) < gamma_h);
new_state(rec_h_mask) = 4;

%% R to S
r_mask = (state == 4);
wane_r_mask = r_mask & (r(:,1) < omega_R);
new_state(wane_r_mask) = 1;

%% Reset duration for transitioned agents
transitioned = (new_state ~= state);
agents.days_in_state(transitioned) = 0;

agents.state = new_state;

end
