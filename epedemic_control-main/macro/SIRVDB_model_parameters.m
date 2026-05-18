clear; close all; clc

%% SIRVDB model and MPC parameters
% Purpose:
%   Defines shared epidemiological, travel, healthcare, and MPC settings
%   for the two-city macro and micro simulations.
%
% Outputs:
%   Workspace parameters used by Simulink, MPC setup, and ABM scripts.

%% City and vaccination
m = 2;
nu1 = 0.001;                % vaccination rate, city 1
nu2 = 0.002;
travel = [-0.02 0.02; 0.01 -0.01];

%% Contact mixing matrix
M11 = 0.85;                 % city 1 contacts within city 1
M12 = 0.15;                 % city 1 contacts with city 2
M21 = 0.1;                  % city 2 contacts with city 1
M22 = 0.9;                  % city 2 contacts within city 2

%% Transmission
beta1 = 0.45;               % baseline transmission, city 1
beta2 = 0.45;
season_amp = 0.4;           % seasonal forcing amplitude
season_period = 365;        % seasonal period in days
eta = 0.6;                  % vaccinated breakthrough susceptibility
epsilon = 0.6;              % vaccine efficacy

%% Recovery and mortality
gamma = 1/12;               % recovery rate, infected
gamma_v = 1/9;              % recovery rate, vaccinated infected
gamma_h = 1/25;             % recovery rate, hospitalised
mu = 0.001;                 % mortality rate, infected
mu_v = 0.0005;              % mortality rate, vaccinated infected
mu_h = 0.02;                % mortality rate, hospitalised
h = 0.05;                   % hospitalisation rate, infected
h_v = 0.02;                 % hospitalisation rate, vaccinated infected

%% Waning immunity
omega_V = 1/140;            % immunity loss after vaccination
omega_R = 1/110;            % immunity loss after infection

%% Healthcare and control
H_cap1 = 0.005;             % hospital capacity fraction, city 1
H_cap2 = 0.005;


%% MPC
u1 = 0;
u2 = 0;
H_ref = 0.8;                % hospital-load reference
k = 0.7;                    % intervention effectiveness

Ts = 1;                     % sampling time in days


load(fullfile(fileparts(mfilename('fullpath')), 'plant.mat'));
mpcobj = mpc(plant, Ts);

mpcobj.PredictionHorizon = 20;
mpcobj.ControlHorizon = 10;
mpcobj.Weights.ManipulatedVariables = 1;     % intervention cost
mpcobj.Weights.ManipulatedVariablesRate = 0.5;
mpcobj.Weights.OutputVariables = 0.8;
mpcobj.Weights.ECR = 1;
mpcobj.MV.Min = 0;
mpcobj.MV.Max = 1;
mpcobj.OV.Max = 1;
xmpc = mpcstate(mpcobj);
