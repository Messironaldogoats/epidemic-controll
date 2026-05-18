function dxdt = sirvhdb_ode(t, x, p, u_val)
% SIRVHDB_ODE
% Two-city SIRVHDB compartmental ODE used for calibration and plant
% linearisation.
%
% States x (14 total, 7 per city):
%   City 1: x(1)=S1  x(2)=I1  x(3)=Iv1 x(4)=R1
%           x(5)=V1  x(6)=H1  x(7)=D1
%   City 2: x(8)=S2  x(9)=I2  x(10)=Iv2 x(11)=R2
%           x(12)=V2 x(13)=H2 x(14)=D2
%
% Inputs:
%   t     - current time (days)
%   x     - state vector (14x1), all compartments as fractions
%   p     - parameter struct (from SIRVDB_model_parameters.m)
%   u_val - intervention level [0,1]
%
% Output:
%   dxdt  - state derivative (14x1)

% Enforce non-negative compartment values during numerical integration.
S1  = max(x(1),  0);  I1  = max(x(2),  0);  Iv1 = max(x(3),  0);
R1  = max(x(4),  0);  V1  = max(x(5),  0);  H1  = max(x(6),  0);
D1  = max(x(7),  0);

S2  = max(x(8),  0);  I2  = max(x(9),  0);  Iv2 = max(x(10), 0);
R2  = max(x(11), 0);  V2  = max(x(12), 0);  H2  = max(x(13), 0);
D2  = max(x(14), 0);

N1 = max(S1+I1+Iv1+R1+V1+H1, 1e-10);
N2 = max(S2+I2+Iv2+R2+V2+H2, 1e-10);

% Seasonal transmission forcing.
sf    = 1 + p.season_amp * sin(2*pi*t / p.season_period);

% Effective beta under intervention.
b1    = p.beta1 * (1 - p.k * u_val) * sf;
b2    = p.beta2 * (1 - p.k * u_val) * sf;

f1 = (I1 + Iv1) / N1;
f2 = (I2 + Iv2) / N2;

% Force of infection with travel mixing.
lam1  = b1 * (p.M11*f1 + p.M12*f2);
lam2  = b2 * (p.M21*f1 + p.M22*f2);

% Breakthrough infection among vaccinated individuals.
lam1v = p.eta * lam1;
lam2v = p.eta * lam2;

%% City 1 ODEs
dS1  = -lam1*S1  - p.nu1*S1  + p.omega_R*R1 + p.omega_V*V1;
dI1  =  lam1*S1  - (p.gamma + p.h + p.mu)*I1;
dIv1 =  lam1v*V1 - (p.gamma_v + p.h_v + p.mu_v)*Iv1;
dR1  =  p.gamma*I1 + p.gamma_v*Iv1 + p.gamma_h*H1 - p.omega_R*R1;
dV1  =  p.nu1*S1 - lam1v*V1 - p.omega_V*V1;
dH1  =  p.h*I1  + p.h_v*Iv1  - (p.gamma_h + p.mu_h)*H1;
dD1  =  p.mu*I1 + p.mu_v*Iv1 + p.mu_h*H1;

%% City 2 ODEs
dS2  = -lam2*S2  - p.nu2*S2  + p.omega_R*R2 + p.omega_V*V2;
dI2  =  lam2*S2  - (p.gamma + p.h + p.mu)*I2;
dIv2 =  lam2v*V2 - (p.gamma_v + p.h_v + p.mu_v)*Iv2;
dR2  =  p.gamma*I2 + p.gamma_v*Iv2 + p.gamma_h*H2 - p.omega_R*R2;
dV2  =  p.nu2*S2 - lam2v*V2 - p.omega_V*V2;
dH2  =  p.h*I2  + p.h_v*Iv2  - (p.gamma_h + p.mu_h)*H2;
dD2  =  p.mu*I2 + p.mu_v*Iv2 + p.mu_h*H2;

dxdt = [dS1;dI1;dIv1;dR1;dV1;dH1;dD1;
        dS2;dI2;dIv2;dR2;dV2;dH2;dD2];
end
