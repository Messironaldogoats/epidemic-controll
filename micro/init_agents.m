function agents = init_agents(N1, N2, I0_1, I0_2, V0_1, V0_2, params)
% INIT_AGENTS
% Initialises a two-city ABM population with heterogeneous susceptibility
% and contact rates.
%
% Epidemiological assumptions:
%   Contact rates follow an overdispersed negative binomial distribution.
%   Susceptibility is gamma distributed to represent heterogeneity in
%   secondary transmission. The commute matrix favours within-city mixing.
%
% State codes:
%   1=S, 2=I, 3=Iv, 4=R, 5=V, 6=H, 7=D
%
% Inputs:
%   N1, N2      - population size per city
%   I0_1, I0_2  - initial infected per city
%   V0_1, V0_2  - initial vaccinated per city
%   params      - parameter struct (must contain M11, M12, M21, M22)
%
% Output:
%   agents - struct containing city assignment, compartment state,
%   susceptibility, daily contacts, state duration, and city random streams.

N = N1 + N2;

%% City assignment
agents.city  = [ones(N1,1); 2*ones(N2,1)];
agents.state = ones(N, 1);
agents.days_in_state = zeros(N, 1);

%% Independent city random streams
agents.stream1 = RandStream('mt19937ar', 'Seed', 42);
agents.stream2 = RandStream('mt19937ar', 'Seed', 137);

%% Susceptibility heterogeneity
% Gamma distributed with mean 1; lower shape values increase overdispersion.
k_susc = 2.0;   % moderate heterogeneity
agents.susceptibility = gamrnd(k_susc, 1/k_susc, N, 1);
agents.susceptibility = max(0.1, min(5.0, agents.susceptibility));

%% Contact-rate distribution
% Negative binomial contacts create overdispersed daily exposure.
mean_contacts = 10;
k_contacts    = 0.5;

p_nb = k_contacts / (k_contacts + mean_contacts);
agents.n_contacts = nbinrnd(k_contacts, p_nb, N, 1);
agents.n_contacts = max(1, min(50, agents.n_contacts));  % cap at 50

%% Initial infected, city 1
city1_idx = find(agents.city == 1);
perm1     = city1_idx(randperm(length(city1_idx)));
I0_1      = min(I0_1, length(city1_idx));
agents.state(perm1(1:I0_1)) = 2;

%% Initial infected, city 2
city2_idx = find(agents.city == 2);
perm2     = city2_idx(randperm(length(city2_idx)));
I0_2      = min(I0_2, length(city2_idx));
agents.state(perm2(1:I0_2)) = 2;

%% Initial vaccinated, city 1
if V0_1 > 0
    avail1    = city1_idx(agents.state(city1_idx) == 1);
    V0_1      = min(V0_1, length(avail1));
    v1_idx    = avail1(randperm(length(avail1), V0_1));
    agents.state(v1_idx) = 5;
end

%% Initial vaccinated, city 2
if V0_2 > 0
    avail2    = city2_idx(agents.state(city2_idx) == 1);
    V0_2      = min(V0_2, length(avail2));
    v2_idx    = avail2(randperm(length(avail2), V0_2));
    agents.state(v2_idx) = 5;
end

end
