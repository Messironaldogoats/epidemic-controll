%% SIRVDB_model_parameters.m
% Run this script to initialise all parameters for macro and micro models.
% Does NOT contain clear/close/clc so it is safe to call from other scripts.

%% -------------------------------------------------------
% City parameters
% -------------------------------------------------------
m       = 2;          % number of cities
nu1     = 0;          % vaccination rate city 1
nu2     = 0;          % vaccination rate city 2

travel  = [-0.02  0.02;
            0.01 -0.01];

% Commute contact matrix
M11 = 0.85;   % city 1 contacts within city 1
M12 = 0.15;   % city 1 contacts within city 2
M21 = 0.10;   % city 2 contacts within city 1
M22 = 0.90;   % city 2 contacts within city 2

%% -------------------------------------------------------
% Transmission
% -------------------------------------------------------
beta1         = 0.45;
beta2         = 0.45;
season_amp    = 0.4;
season_period = 365;
eta           = 0.6;    % infection breakthrough (vaccinated)
epsilon       = 0.6;    % vaccine efficiency

%% -------------------------------------------------------
% Recovery and mortality
% -------------------------------------------------------
gamma   = 1/12;    % recovery rate (unvaccinated infected)
gamma_v = 1/9;     % recovery rate (vaccinated infected)
gamma_h = 1/25;    % recovery rate (hospitalised)
mu      = 0.001;   % mortality rate (infected)
mu_v    = 0.0005;  % mortality rate (vaccinated infected)
mu_h    = 0.02;    % mortality rate (hospitalised)
h       = 0.05;    % hospitalisation rate (infected)
h_v     = 0.02;    % hospitalisation rate (vaccinated infected)

%% -------------------------------------------------------
% Waning immunity
% -------------------------------------------------------
omega_V = 1/140;   % immunity loss after vaccination
omega_R = 1/110;   % immunity loss after recovery

%% -------------------------------------------------------
% Healthcare capacity
% -------------------------------------------------------
H_cap1 = 0.005;                          % fraction of city 1 population
H_cap2 = 0.005;                          % fraction of city 2 population
H_ref  = 0.8 * (H_cap1 + H_cap2);       % 80% of total national capacity

%% -------------------------------------------------------
% MPC intervention
% -------------------------------------------------------
u1 = 0;    % 0 = no restrictions, 1 = full lockdown
u2 = 0;
k  = 0.7;  % intervention effectiveness

Ts = 1;    % sampling time (1 day)

%% -------------------------------------------------------
% Package params struct for micro model (update_agents)
% -------------------------------------------------------
params.nu1     = nu1;
params.nu2     = nu2;
params.eta     = eta;
params.gamma   = gamma;
params.gamma_v = gamma_v;
params.gamma_h = gamma_h;
params.mu      = mu;
params.mu_v    = mu_v;
params.mu_h    = mu_h;
params.h       = h;
params.h_v     = h_v;
params.omega_V = omega_V;
params.omega_R = omega_R;

%% -------------------------------------------------------
% MPC setup
% plant.mat lives in macro/ — addpath('macro') must be
% called before this script runs.
% -------------------------------------------------------
load(fullfile(fileparts(mfilename('fullpath')), 'plant.mat'));

mpcobj = mpc(plant, Ts);

mpcobj.PredictionHorizon                = 40;
mpcobj.ControlHorizon                   = 10;
mpcobj.Weights.ManipulatedVariables     = 0.05;
mpcobj.Weights.ManipulatedVariablesRate = 0.2;
mpcobj.Weights.OutputVariables          = 1;

mpcobj.MV.Min = 0;
mpcobj.MV.Max = 1;
mpcobj.OV.Max = 1;

xmpc = mpcstate(mpcobj);