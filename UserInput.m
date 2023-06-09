% Model name
model = 'dipole';
show_mesh = 1;
% Operating frequency
fmin = 1e9;
fmax = 9e9;
NoF  = 101;
% Requested number of eigenmodes
EigNum = 7;
show_eigspectrum = 1;
orth_opt = 0;
track = 0;
% Media definitions
er = 12.0*(1 - 1i*0.02);
ur = 1.0*(1 - 1i*0.00);
% Skin effects and IBC
zs = 0 + 1i*0;
% Quadrature rules
nP = 6;
nQ = 7;
% Excitations
E0 = 1;
% Domain decomposition
use_sc = 0;
% Field point definitions
fp = [0, 0, 0];