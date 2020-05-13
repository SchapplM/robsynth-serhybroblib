% Return Structural Kinematic Parameters of the Robot 
% fourbar1turnDE2
%
% Output:
% v_mdh [6x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [6x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [6x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x5) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = fourbar1turnDE2_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 1; 3; 4;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 0; 0; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 5;
% Anzahl der Kinematikparameter
NKP = 5;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 2;
% Namen der Kinematikparameter
pkin_names = {'l1', 'l2', 'l3', 'l4', 'l5'};
