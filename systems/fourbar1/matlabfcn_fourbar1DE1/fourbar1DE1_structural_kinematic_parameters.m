% Return Structural Kinematic Parameters of the Robot 
% fourbar1DE1
%
% Output:
% v_mdh [5x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [5x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [5x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x4) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = fourbar1DE1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 0; 2; 3;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 0; 0; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 4;
% Anzahl der Kinematikparameter
NKP = 4;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 1;
% Namen der Kinematikparameter
pkin_names = {'l1', 'l2', 'l3', 'l4'};
