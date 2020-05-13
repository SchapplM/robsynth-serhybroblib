% Return Structural Kinematic Parameters of the Robot 
% palh3m2DE2
%
% Output:
% v_mdh [12x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [12x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [12x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x18) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = palh3m2DE2_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4; 1; 2; 7; 7; 4; 6; 8;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 0; 1; 0; 0; 0; 0; 0; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 9;
% Anzahl der Kinematikparameter
NKP = 18;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 4;
% Namen der Kinematikparameter
pkin_names = {'AB', 'BC', 'BE', 'BG', 'DC', 'DT2', 'EP', 'GH', 'GP', 'HW', 'OT1', 'T1A', 'T1T2', 'phi1', 'phi2', 'phi410', 'phi78', 'phi79'};
