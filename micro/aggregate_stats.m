function stats = aggregate_stats(agents, H_cap1, H_cap2)
% AGGREGATE_STATS
% Aggregates ABM compartment counts globally and per city.
%
% Inputs:
%   agents       - ABM state struct
%   H_cap1,H_cap2 - optional hospital capacity fractions
%
% Output:
%   stats - struct with global counts, city counts, fractions, and optional
%   hospital-load signals.
%
% State codes:
% 1=S, 2=I, 3=Iv, 4=R, 5=V, 6=H, 7=D

if nargin < 2
    H_cap1 = [];
    H_cap2 = [];
end

state = agents.state(:);
city  = agents.city(:);

c1 = (city == 1);
c2 = (city == 2);

alive_mask = (state ~= 7);

%% Global counts
stats.S  = sum(state == 1);
stats.I  = sum(state == 2);
stats.Iv = sum(state == 3);
stats.R  = sum(state == 4);
stats.V  = sum(state == 5);
stats.H  = sum(state == 6);
stats.D  = sum(state == 7);

stats.I_total = stats.I + stats.Iv;
stats.N_total = numel(state);
stats.alive   = sum(alive_mask);

if stats.alive > 0
    stats.I_frac = stats.I_total / stats.alive;
    stats.H_frac = stats.H / stats.alive;
else
    stats.I_frac = 0;
    stats.H_frac = 0;
end

%% City 1 counts
stats.city1.S  = sum(state(c1) == 1);
stats.city1.I  = sum(state(c1) == 2);
stats.city1.Iv = sum(state(c1) == 3);
stats.city1.R  = sum(state(c1) == 4);
stats.city1.V  = sum(state(c1) == 5);
stats.city1.H  = sum(state(c1) == 6);
stats.city1.D  = sum(state(c1) == 7);

stats.city1.I_total = stats.city1.I + stats.city1.Iv;
stats.city1.N_total = sum(c1);
stats.city1.alive   = sum(c1 & alive_mask);

if stats.city1.alive > 0
    stats.city1.I_frac = stats.city1.I_total / stats.city1.alive;
    stats.city1.H_frac = stats.city1.H / stats.city1.alive;
else
    stats.city1.I_frac = 0;
    stats.city1.H_frac = 0;
end

%% City 2 counts
stats.city2.S  = sum(state(c2) == 1);
stats.city2.I  = sum(state(c2) == 2);
stats.city2.Iv = sum(state(c2) == 3);
stats.city2.R  = sum(state(c2) == 4);
stats.city2.V  = sum(state(c2) == 5);
stats.city2.H  = sum(state(c2) == 6);
stats.city2.D  = sum(state(c2) == 7);

stats.city2.I_total = stats.city2.I + stats.city2.Iv;
stats.city2.N_total = sum(c2);
stats.city2.alive   = sum(c2 & alive_mask);

if stats.city2.alive > 0
    stats.city2.I_frac = stats.city2.I_total / stats.city2.alive;
    stats.city2.H_frac = stats.city2.H / stats.city2.alive;
else
    stats.city2.I_frac = 0;
    stats.city2.H_frac = 0;
end

%% Hospital-load signals
if ~isempty(H_cap1) && ~isempty(H_cap2)
    stats.city1.H_load = stats.city1.H / (H_cap1 * stats.city1.N_total);
    stats.city2.H_load = stats.city2.H / (H_cap2 * stats.city2.N_total);
    stats.H_load       = stats.H / ((H_cap1 + H_cap2) * stats.N_total / 2);
else
    stats.city1.H_load = [];
    stats.city2.H_load = [];
    stats.H_load       = [];
end

end
