% Return Structural Kinematic Parameters of the Robot 
% palh4m1DE2
%
% Output:
% v_mdh [9x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [9x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [9x1]
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
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = palh4m1DE2_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4; 5; 1; 4; 7;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 1; 0; 0; 0; 0; 0; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 0; 1; 0; 1; 1; 1; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 8;
% Anzahl der Kinematikparameter
NKP = 9;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 5;
% Namen der Kinematikparameter
pkin_names = {'AB', 'AD', 'CB', 'CE', 'EP', 'HC', 'OT', 'TA', 'TD'};
