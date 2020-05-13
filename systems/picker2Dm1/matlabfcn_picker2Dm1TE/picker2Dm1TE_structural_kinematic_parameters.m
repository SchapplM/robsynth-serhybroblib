% Return Structural Kinematic Parameters of the Robot 
% picker2Dm1TE
%
% Output:
% v_mdh [15x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [15x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [15x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x9) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = picker2Dm1TE_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 2; 0; 2; 0; 1; 3; 4; 8; 6; 7; 5; 9;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2; 2; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 11;
% Anzahl der Kinematikparameter
NKP = 9;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 2;
% Namen der Kinematikparameter
pkin_names = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'e', 'phi05', 'phi1'};
