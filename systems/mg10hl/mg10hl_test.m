% Teste Klassendefinition für Schwerlastroboter Kawasaki MG10HL
% 
% TODO: Dieser Roboter funktioniert noch nicht. Die kinematischen
% Zwangsbedingungen werden nicht geschlossen. Ursache ist wahrscheinlich
% eine geänderte Definition der vorab gespeicherten Werte für die
% Viergelenkkette (fourbar1) oder die Viergelenkkette mit Schubgelenk
% (fourbarpris). In der ursprünglichen Version von Abderahman Bejaoui hat
% es funktioniert.

% * Quelle des Modells: [Bejaoui2018_S749]

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 202-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('mg10hl', 'TE');
RS_DE1 = hybroblib_create_robot_class('mg10hl', 'DE1');
RS_DE2 = hybroblib_create_robot_class('mg10hl', 'DE1');

% Kinematikparameter und Modell initialisieren (ohne models.csv)
% Aus mg10hlTE/kinematic_parameter_values.m
AC=.4; AE=.3; DC=.45; ED=.2; TA=.15; TE=.15; % fuer VGK  A-C-D-E
GP=.3; GK=.5; HP=.2; % fuer Slcr G-P-K 
PM=.2; ML=.2; LW=.2; OT=.3;
CG=.3;
phi3=pi/3; phi23=pi/10; phi34=pi/3;

% Siehe: mg10hlTE_fkine_fixb_rotmat_mdh_sym_varpar
pkin = [AC; AE; CG; DC; ED; GK; GP; HP; LW; ML; OT; PM; TA; TE; phi23; phi3; phi34;];

% Gelenkgrenzen durch Ausprobieren mit Debug-Plot von unten
qlim =  [0.7 2; ...
         -1 0.5; ...
         -1 1; ...
         -3 3; ...
         -2 3; ...
         .05 .25];
% Modelle aktualisieren
RS_DE1.update_mdh(pkin);
RS_DE2.update_mdh(pkin);
RS_TE.update_mdh(pkin);
RS_TE.qlim = qlim;
% Kinematikparameter der offenen Baumstruktur setzen
% Siehe: mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar
pkin_OL = [AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
RS_OL = hybroblib_create_robot_class('mg10hl', 'OL');
RS_OL.update_mdh(pkin_OL);
%% Vergleich der Implementierungen
TSS = RS_TE.gen_testsettings();
n_test = 0;
I_iO = false(TSS.n, 1);
for i = 1:TSS.n
  q=TSS.Q(i,:)';
  try
    T_DE1 = RS_DE1.fkine(q);
    T_DE2 = RS_DE2.fkine(q);
    T_TE = RS_TE.fkine(q);
  catch
    warning('Fehler beim Berechnen der Kinematik. Zwangsbedingungen nicht erfüllt, Winkel außerhalb des Definitionsbereichs.');
    continue; % Gelenkkoordinate gehört nicht zum Definitionsgebiet
  end
  if any(imag(T_TE(:)) > 1e-10)
    warning('Die Transformationsmatrix wird imaginär. Zwangsbedingungen nicht erfüllt, Winkel außerhalb des Definitionsbereichs.');
    continue
  end
  I_iO(i) = true;

  n_test = n_test + 1;
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
  if any(abs([test1(:); test2(:)]) > 1e-10)
    error('Methoden DE1, DE2 und TE stimmen nicht überein');
  end
end
fprintf('Kinematik vom Typ Yaskawa MG10HL für %d Kombinationen in unterschiedlichen Implementierungen getestet. Restliche %d nicht lösbar.\n', ...
  n_test, TSS.n-n_test);
% Funktionierende Gelenkwinkel anzeigen
figure(100);clf;
sgtitle('Testen: Wertebereich für Gelenkwinkel');
for j = 1:length(q)
  subplot(length(q), 1, j); hold on;
  if sum(I_iO) > 0
    plot(TSS.Q(I_iO,j), 1, 'gv');
  end
  if sum(~I_iO) > 0
    plot(TSS.Q(~I_iO,j), 1, 'rx');
  end
  ylabel(sprintf('Achse %d', j));
  xlabel(sprintf('Werte für q%d', j));
  grid on;
end

%% Funktionen direkt testen
% Ausprobieren sinnvoller Gelenkwinkel zum Zeichnen
q = [0 ; ...% Karussel / erste Achse
  -25; ... % zweite Achse
  5; ... % drittletzte Achse (Arm)
  -60; ... % vorletzte Achse (Hand)
  75; ... % letzte Achse
  NaN]*pi/180;
q(6) = 0.2; % Schubachse
% Testweise berechnen der Kinematik
Tc_TE = mg10hlTE_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
T_TE = mg10hlTE_joint_trafo_rotmat_mdh_sym_varpar(q, pkin);
% Roboter in Grundstellung zeichnen
s_plot = struct( 'ks', [1:RS_TE.NJ, RS_TE.NJ+2], 'straight', 0);
Tc_test = RS_TE.fkine(q);

figure(2);clf;
hold on;
grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3);
title('Nullstellung TE');
RS_TE.plot( q, s_plot );
view([0 0])

%% Debugge Zwangsbedingungen mit der offenen Struktur
figure(3);clf;
hold on;
grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3);
title('Nullstellung OL Korrekturen');
q_jv = RS_TE.jointvar(q);
q_OL = q_jv(RS_TE.MDH.sigma~=2);
% TODO: Neuberechnung der Winkel mit fourbar1 und fourbarpris
% Dafür: 
% * Klasse der beiden Mechanismen erstellen
% * Gelenkwinkel für Längen dieses Roboters berechnen
% * Winkel eintragen in die Variable
% * Zeichnung prüfen
% * Maple korrigieren
% Korrekturen
q_OL(3) = q_OL(3) + 45*pi/180;
q_OL(11) = q_OL(11) + 125*pi/180;
q_OL(10) = q_OL(10) - 15*pi/180;
Tc_OL = RS_OL.fkine(q_OL);
RS_OL.plot( q_OL, s_plot );
view([0 0])

% Prüfe die Viergelenkkette als ein Teil des Roboters
% TODO: Das ist noch nicht fertig:
RS_fb1 = hybroblib_create_robot_class('fourbar1', 'TE');
s_plot_fb1 = struct( 'ks', [1:RS_fb1.NJ, RS_fb1.NJ+2], 'straight', 0);
RS_fb1.pkin_names
pkin_fb1 = [TA+TE;AC; DC ; ED];
RS_fb1.update_mdh(pkin_fb1);
RS_fb1.update_base([0;0;OT], [-pi/2;0;pi]);
figure(101);clf;
hold on; grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3); title('Viergelenkkette');
RS_fb1.plot(115*pi/180, s_plot_fb1);
view([0 0])

% Prüfe die Schubkurbel als ein Teil des Roboters
% TODO: Das ist noch nicht fertig:
RS_sc = hybroblib_create_robot_class('fourbarpris', 'TE');
s_plot_sc = struct( 'ks', [1:RS_fb1.NJ, RS_fb1.NJ+2], 'straight', 0);
RS_sc.pkin_names
pkin_sc = [GK; GP; HP];
RS_sc.update_mdh(pkin_sc);
RS_sc.update_base(Tc_OL(1:3,4,10), r2eulxyz(Tc_OL(1:3,1:3,10)));
figure(102);clf;
hold on; grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3); title('Schubkurbel');
RS_sc.plot(q(6)); % 
view([0 0])
