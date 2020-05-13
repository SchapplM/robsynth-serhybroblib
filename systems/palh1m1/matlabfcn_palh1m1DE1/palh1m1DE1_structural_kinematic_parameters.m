% Return Structural Kinematic Parameters of the Robot 
% palh1m1DE1
%
% Output:
% v_mdh [16x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [16x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [16x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x23) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = palh1m1DE1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4; 1; 2; 2; 8; 7; 7; 3; 4; 6; 9; 10;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2; 2; 2;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 11;
% Anzahl der Kinematikparameter
NKP = 23;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 4;
% Namen der Kinematikparameter
pkin_names = {'AB', 'AM', 'BC', 'BE', 'BG', 'BL', 'DA', 'DC', 'EP', 'GH', 'GP', 'HW', 'ML', 'OT2', 'T1D', 'T2A', 'T2T1', 'phi1', 'phi2', 'phi312', 'phi413', 'phi710', 'phi711'};
