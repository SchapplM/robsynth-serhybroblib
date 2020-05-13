% Return Structural Kinematic Parameters of the Robot 
% palh2m2DE
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
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = palh2m2DE_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4; 5;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 0; 1; 0; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 7;
% Anzahl der Kinematikparameter
NKP = 5;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 4;
% Namen der Kinematikparameter
pkin_names = {'A2Off', 'A3Off', 'A4Off', 'L1', 'L2'};
