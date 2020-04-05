% Berechne die Kinematikparameter aller möglicher Palettierer-Modelle aus
% gegebenen Parametern
% 
% TODO: Prüfen, ob mit A3Off funktioniert
% TODO: Dokumentation ausweiten

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function pkin_struct = palletizer_exampleparam(Par_Struct)

%% Eingabe verarbeiten
if nargin == 0
  % Neudefinition
  A2Off   = 300e-3;
  L1      = 1100e-3;
  A3Off   = 0e-3;
  L2      = 1200e-3;
  A4Off   = 200e-3;
else
  A2Off   = Par_Struct.A2Off;
  L1      = Par_Struct.L1;
  A3Off   = Par_Struct.A3Off;
  L2      = Par_Struct.L2;  
  A4Off   = Par_Struct.A4Off;
end
%% Ausgabe vorbereiten
pkin_struct = struct('pkin_palh1m1', [], 'pkin_palh2m1', [], 'pkin_palh2m2', [], 'pkin_S5RRRRR1', [], ...
  'theta_offset', []);

%% palh1m1: Kinematik-Parameter des 3-Parallelogramm-Roboters berechnen
% Siehe: palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar
% Anpassung der Kinematikparameter. Siehe [Bejaoui2018_S749]
% systems/palh1m1/codegen/palh1m1_kinematic_parameter_values.m

AB = sqrt(L1^2+A3Off^2); % gesetzt
% AM = A2Off;
% BC = A2Off;
BE = sqrt(2)*A4Off; % gesetzt // Warum hier A4?
BG = L2; % gesetzt
% BL = BC;% gesetzt
% DA = A2Off;
DC = AB;% gesetzt
EP = L2; % gesetzt
GH = A4Off; % gesetzt
GP = BE; % gesetzt aus Annahme
HW = 0; % % gesetzt
ML = AB; % gesetzt
OT2 = 0;% gesetzt
T1D = A2Off;% gesetzt - Vorzeichen eventuell negativ
T2A = A2Off;% gesetzt
T2T1 = -A2Off; % gesetzt

phi312=2*pi/3;
% konstaner Winkel fuer VGK Blau , fuer Rechnung der ZW erforderlich
T1S=((T2T1*T1D)/(T1D+T2A));
T2S=T2T1-T1S;
phi1=atan2(T1S,T1D); phi2=atan2(T2S,T2A);
AS=sqrt(T2A^2+T2S^2);
DS=sqrt(T1S^2+T1D^2);
DA=AS+DS;
AM = DA; % gesetzt
BC = DA; % gesetzt
BL = BC;% gesetzt
% konstanter Winkel fuer VGK GRUEN
BT3=A4Off; % gesetzt
phi711=acos(BT3/BC); 
phi413=pi/2; 
phi710=acos(BT3/BE);

% Definition der Roboterklasse für hybriden Palettierroboer
% TODO: Anderes Modell aus Modellliste nehmen, sobald verfügbar
% Siehe serhybrob-mdl/systems/palh1m1/models.csv
RS = hybroblib_create_robot_class('palh1m1', 'TE', 'palh1m1Bsp1');

pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
RS.update_mdh(pkin);

% Für den Winkel vorm Endeffektor kommen "krumme" Werte heraus, die durch
% Berechnung der Kinematik bestimmt werden können
% Winkelfehler mit alten Parametern wird mit direkter Kinematik bestimmt
q_test = zeros(4,1);
T_EE = RS.fkineEE(q_test);
phi_diff = r2eulxyz(t2r(T_EE));
phi413_neu = phi413+phi_diff(2);
pkin_it2=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413_neu,phi710,phi711]';
RS.pkin= pkin_it2;
pkin_struct.pkin_palh1m1 = pkin_it2;


%% palh3m1: Kinematik-Parameter des 2-Parallelogramm-Roboters berechnen
% Siehe: palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar
% Anpassung der Kinematikparameter. Siehe [Shan2019_S828]

AB = sqrt(L1^2+A3Off^2); % 1. Parallelogramm, Seite 1; gesetzt
% BC = 0.6; % 1. Parallelogramm, Seite 2
BE = sqrt(2)*A4Off; % gesetzt; Warum hier A4? 2. Parallelogramm, Seite 2
BG = L2; % 2. Parallelogramm; gesetzt
% DA = 0.6; % 1. Parallelogramm, Seite 2
DC = sqrt(L1^2+A3Off^2); % 1. Parallelogramm, Seite 1; gesetzt

EP = L2; % 2. Parallelogramm; gesetzt
GH = A4Off; % gesetzt
GP = BE; % gesetzt aus Annahme; 2. Parallelogramm, Seite 2
HW = 0; % gesetzt
OT1 = 0; % gesetzt
T1A = A2Off; % 0.42;
T1T2 = A2Off; % 0.48; % TODO
DT2 = sqrt(DA^2-T1T2^2)-T1A;%0.22;
phi1 = atan2(T1T2,T1A+DT2);
phi2 = phi1;
phi410 = pi/3;
phi78 = phi710; % Übernehmen von palh1m1-Beispiel. TODO: Eigene Berechnung hier
phi79 = phi711; % von oben. TODO: Nachvollziehen.

% Definition der Roboterklasse für hybriden Palettierroboer
% TODO: Anderes Modell aus Modellliste nehmen, sobald verfügbar
% Siehe serhybrob-mdl/systems/palh3m1/models.csv
RS = hybroblib_create_robot_class('palh3m1', 'TE', 'palh3m1Bsp1');

pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
RS.update_mdh(pkin);

% Für den Winkel vorm Endeffektor kommen "krumme" Werte heraus, die durch
% Berechnung der Kinematik bestimmt werden können
% Winkelfehler mit alten Parametern wird mit direkter Kinematik bestimmt
q_test = [0;pi/2;pi/2;0]; % Null-Stellung
T_EE = RS.fkineEE(q_test);
phi_diff = r2eulxyz(t2r(T_EE));
phi410_neu = phi410+phi_diff(2);
pkin_it2=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410_neu,phi78,phi79]';
RS.pkin= pkin_it2;
pkin_struct.pkin_palh3m1 = pkin_it2;

%% palh2m1:Parameter
% Parameter für serielles Modell, siehe 
% S5RRRRR1_fkine_fixb_rotmat_mdh_sym_varpar und
% palh2m1DE_fkine_fixb_rotmat_mdh_sym_varpar 
% pkin=[a2,a3,a4,a5,d1,d5]';
pkin = NaN(6,1);
pkin(1) = A2Off; % a2
pkin(2) = sqrt(L1^2+A3Off^2); % a3
pkin(3) = L2; % a4
pkin(4) = A4Off; % a5
pkin(5) = 0; % d1
pkin(6) = 0; % d5
pkin_struct.pkin_palh2m1 = pkin;
pkin_struct.pkin_S5RRRRR1 = pkin;
%% palh2m2: Parameter
% Wie palh2m1, aber ohne d1,d5; mit offset
% zusätzlich Anpassung des a3-Parameters an Parallel verschobenen A3off
% Siehe palh2m2DE_fkine_fixb_rotmat_mdh_sym_varpar 
% pkin=[A2Off,A3Off,A4Off,L1,L2]';
theta_offset = atan2(A3Off, L1);
pkin = NaN(5,1);
pkin(1) = A2Off; % A2Off
pkin(2) = A3Off; % A3Off
pkin(3) = A4Off; % A4Off
pkin(4) = L1; % L1
pkin(5) = L2; % L2
pkin_struct.pkin_palh2m2 = pkin;
pkin_struct.theta_offset = theta_offset;
