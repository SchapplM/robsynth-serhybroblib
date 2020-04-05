% Teste die Kinematik verschiedener (Ersatz)-Modelle von Palettier-Robotern
% Modelle:
% * (1) Serielles Modell, 5FG
% * (2) Serielles Modell, 4FG (Elimination des abhängigen Winkels)
% * (3) Hybrides Modell, 4FG (drei geschlossene Ketten)
% * (4) Serielles Modell, 4FG (Elimination, zusätzlicher Offset-Parameter)
% * (5) Hybrides Modell, 4FG (zwei geschlossene Ketten)
% 
% Ergebnis:
% Alle Modelle haben die gleiche direkte Kinematik
% 
% TODO:
% * Gelenkwinkelgrenzen für hybrides Modell noch sehr gering

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Pfade initialisieren
if isempty(which('serroblib_path_init.m'))
  warning('Repo mit Robotermodellen ist nicht im Pfad. Beispiel nicht ausführbar.');
  return
end

%% Modelle initialisieren
% Typ des seriellen Roboters auswählen (5FG-Palettierroboter)
SName='S5RRRRR1';
RS(1) = serroblib_create_robot_class(SName);
RS(1).descr = 'S5RRRRR1';
serroblib_addtopath({SName});

% Modell des 5FG-Roboters mit einem abhängigen FG als Ersatzmodell
RS(2) = hybroblib_create_robot_class('palh2m1', 'DE');
RS(2).descr = 'palh2m1';

% Vollständiger hybriden Palettierroboter (drei Parallelogramme)
RS(3) = hybroblib_create_robot_class('palh1m1', 'TE', 'palh1m1Bsp1');
RS(3).descr = 'palh1m1';
% EE-Drehen
RS(3).update_EE([], [pi;0;0]);

% Modell des 5FG-Roboters mit einem abhängigen FG als Ersatzmodell, mit
% Winkel-Offset direkt als Parameter
RS(4) = hybroblib_create_robot_class('palh2m2', 'DE');
RS(4).descr = 'palh2m2';

% Vollständigen hybriden Palettierroboter (zwei Parallelogramme)
RS(5) = hybroblib_create_robot_class('palh3m1', 'TE', 'palh3m1Bsp1');
RS(5).descr = 'palh3m1';
% EE-Drehen
RS(5).update_EE([], [pi;0;0]);
%% Parameter setzen

% Beispiel-Parameter laden
pkin_struct = palletizer_exampleparam();

% Parameter der Modelle aktualisieren
RS(1).update_mdh(pkin_struct.pkin_S5RRRRR1);
RS(2).update_mdh(pkin_struct.pkin_palh2m1);
RS(3).update_mdh(pkin_struct.pkin_palh1m1);
RS(4).update_mdh(pkin_struct.pkin_palh2m2);
RS(5).update_mdh(pkin_struct.pkin_palh3m1);
theta_offset = pkin_struct.theta_offset;

%% Roboter in verschiedenen Posen plotten
for ii = 1:4
  % Roboter in Startpose (Null-Stellung)
  if ii == 1
    q1 = [0; -pi/2 + theta_offset; pi/2 - theta_offset; 0; 0];
    q2 = [0; -pi/2 + theta_offset; pi/2 - theta_offset; 0];
    q3 = [0;     0 - theta_offset;    0 + theta_offset; 0];
    q4 = [0; pi/2;                0;                0];
    q5 = [0; pi/2 - theta_offset;  pi/2 + theta_offset; 0];
  end
  % Variation q2
  if ii == 2
    q1 = [0; -15-90 + theta_offset; 90 - theta_offset; 15; 0]*pi/180;
    q2 = [0; -15-90 + theta_offset; 90 - theta_offset;     0]*pi/180;
    q3 = [0;  15    - theta_offset;  0 + theta_offset;     0]*pi/180;
    q4 = [0;  15+90 - theta_offset; 15 + theta_offset;     0]*pi/180;
    q5 = [0;  15+90 - theta_offset; 90 + theta_offset;     0]*pi/180;
  end
  % Variation q3
  if ii == 3
    q1 = [0; -90 + theta_offset; -30+90 - theta_offset; 30; 0]*pi/180;
    q2 = [0; -90 + theta_offset; -30+90 - theta_offset;     0]*pi/180;
    q3 = [0;   0 - theta_offset;  30    + theta_offset;     0]*pi/180;
    q4 = [0;  90 - theta_offset;  30    + theta_offset;     0]*pi/180;
    q5 = [0;  90 - theta_offset;  30+90 + theta_offset;     0]*pi/180;
  end
  % Variation q1-q4
  if ii == 4
    q1 = [45; -90 + theta_offset; 90 - theta_offset;  0; -45]*pi/180;
    q2 = [45; -90 + theta_offset; 90 - theta_offset;     -45]*pi/180;
    q3 = [45;   0 - theta_offset;  0 + theta_offset;      45]*pi/180;
    q4 = [45;  90 - theta_offset;  0 + theta_offset;     -45]*pi/180;
    q5 = [45;  90 - theta_offset; 90 + theta_offset;      45]*pi/180;
  end
  % Kinematik berechnen-
  T1 = RS(1).fkineEE(q1);
  T2 = RS(2).fkineEE(q2);
  T3 = RS(3).fkineEE(q3);
  T4 = RS(4).fkineEE(q4);
  T5 = RS(5).fkineEE(q5);
  % Plotten
  figure(ii);clf;
  for i = 1:5
    s_plot = struct( 'ks', [1:RS(i).NJ, RS(i).NJ+2], 'straight', 1);
    subplot(2,3,i);
    hold on;grid on;
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');view(3);
    title(RS(i).descr);
    if i == 1, q = q1; end
    if i == 2, q = q2; end
    if i == 3, q = q3; end
    if i == 4, q = q4; end
    if i == 5, q = q5; end
    RS(i).plot( q, s_plot );
    view(0,0);
  end
  sgtitle(sprintf('Modellvergleich in Pose %d', ii));
  if any(abs(T1(:)-T2(:))>1e-10)
    error('Modell 1 stimmt nicht gegen Modell 2');
  end
  if any(abs(T1(:)-T3(:))>1e-10)
    error('Modell 1 stimmt nicht gegen Modell 3');
  end
  % Modell 4 stimmt nicht für diese
  if any(abs(T1(:)-T4(:))>1e-10)
    error('Modell 1 stimmt nicht gegen Modell 4');
  end
  if any(abs(T1(:)-T5(:))>1e-10)
    error('Modell 1 stimmt nicht gegen Modell 5');
  end
end

%% Trajektorie Beispiel
% Beispiel-Trajektorie im Gelenkraum: Fahren zwischen Gelenkwinkelgrenzen
% Trajektorie für vier FG Hybrid

% Anpassung der Grenzen: Hybrides Modell ist anscheinend nur in geringem
% Bereich definiert (sonst klappen die Viergelenkketten um)
RS(3).qlim(2,2) = 25*pi/180;
RS(3).qlim(3,2) = 30*pi/180;

% Für jedes Gelenk Verfahrbewegung an Gesamttrajektorie anhängen
k=1; QE_3 = RS(3).qref';
for i = 1:RS(3).NQJ
  % Fahrt zu Minimalwert
  qi = RS(3).qref;
  qi(i) = RS(3).qlim(i,1);
  k=k+1; QE_3(k,:) = qi;
  % Fahrt zu Maximalwert
  qi(i) = RS(3).qlim(i,2);
  k=k+1; QE_3(k,:) = qi;
  % Fahrt zur Referenzpose
  qi = RS(3).qref;
  k=k+1; QE_3(k,:) = qi;
end
% Trajektorie berechnen
[Q_3,QD_3,QDD_3,T] = traj_trapez2_multipoint(QE_3, 1, 1e-1, 1e-2, 1e-3, 0.25);

% Umrechnen auf 5G seriell
Q_1 = NaN(size(Q_3,1),5);
Q_1(:,1) =  Q_3(:,1);
Q_1(:,2) = -Q_3(:,2) - pi/2;
Q_1(:,3) = -Q_3(:,3) + pi/2;
Q_1(:,4) =  Q_3(:,2)+Q_3(:,3);% Achse 4 so, dass Roboter immer parallel zum Boden
Q_1(:,5) = -Q_3(:,4);

% Umrechnen auf 4FG seriell (mit Elim)
Q_2 = NaN(size(Q_3,1),4);
Q_2(:,1) =  Q_3(:,1);
Q_2(:,2) = -Q_3(:,2) - pi/2;
Q_2(:,3) = -Q_3(:,3) + pi/2;
Q_2(:,4) = -Q_3(:,4);

% Modell 4 ist das gleiche wie 2. Hier müssen nur keine Offsets addiert
% werden
Q_4(:,1) = Q_2(:,1);
Q_4(:,2) = -Q_2(:,2);
Q_4(:,3) = -Q_2(:,3) -Q_2(:,2);
Q_4(:,4) = Q_2(:,4);

% Berücksichtigung des Winkel-Offsets resultierend aus A3off
Q_1 = [Q_1(:,1), Q_1(:,2)+theta_offset, Q_1(:,3)-theta_offset, Q_1(:,4), Q_1(:,5)];
Q_2 = [Q_2(:,1), Q_2(:,2)+theta_offset, Q_2(:,3)-theta_offset, Q_2(:,4)];
Q_3 = [Q_3(:,1), Q_3(:,2)-theta_offset, Q_3(:,3)+theta_offset, Q_3(:,4)];

%% Direkte Kinematik berechnen
% Direkte Kinematik für Trajektorie berechnen
X_1 = NaN(size(Q_1, 1), 6);
X_2 = X_1; X_3 = X_1;
for ii = 1:size(Q_1, 1)
  T_0_Ei = RS(1).fkineEE(Q_1(ii,:)');
  X_1(ii,:) = RS(1).t2x(T_0_Ei);
  T_0_Ei = RS(2).fkineEE(Q_2(ii,:)');
  X_2(ii,:) = RS(2).t2x(T_0_Ei);
  T_0_Ei = RS(3).fkineEE(Q_3(ii,:)');
  X_3(ii,:) = RS(3).t2x(T_0_Ei);
  T_0_Ei = RS(4).fkineEE(Q_4(ii,:)');
  X_4(ii,:) = RS(4).t2x(T_0_Ei);
  
  % Teste direkte Kinematik direkt bei Trajektorienerstellung
  test12 = X_2(ii,:)-X_1(ii,:);
  if any(abs(test12(:)) > 1e-6)
    error('Kinematik 1 gegen 2 stimmt nicht');
  end
  test13 = X_3(ii,:)-X_1(ii,:);
  if any(abs(test13(:)) > 1e-6)
    error('Kinematik 1 gegen 3 stimmt nicht');
  end
  test23 = X_3(ii,:)-X_2(ii,:);
  if any(abs(test23(:)) > 1e-6)
    error('Kinematik 2 gegen 3 stimmt nicht');
  end
  test24 = X_4(ii,:)-X_2(ii,:);
  if any(abs(test24(:)) > 1e-6)
    % TODO: Fehlerursache bei modell palh2m2 beseitigen und Test wieder
    % korrekt durchführen
    error('Kinematik 2 gegen 4 stimmt nicht');
    if ii == 1, warning('Kinematik 2 gegen 4 stimmt nicht'); end
  end
  
  % Zeichnen der Konfiguration (im Fehlerfall)
  continue
  figure(20);clf;
  for i = 1:3
    s_plot = struct( 'ks', [1:RS(i).NJ, RS(i).NJ+2], 'straight', 1);
    subplot(2,2,i);
    hold on;grid on;
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');view(3);
    title(RS(i).descr);
    if i == 1, q0 = Q_1(ii,:)'; end
    if i == 2, q0 = Q_2(ii,:)'; end
    if i == 3, q0 = Q_3(ii,:)'; end
    RS(i).plot( q0, s_plot );
  end
end

%% Daten vergleichen
% Direkte Kinematik für Trajektorie plotten
figure(11);clf;
for k = 1:6
  if k < 4, subplot(3,2,sprc2no(3,2,k,1));hold on;
  else,     subplot(3,2,sprc2no(3,2,k-3,2));hold on; end
  plot(T, X_1(:,k), '-');
  plot(T, X_2(:,k), '-.');
  plot(T, X_3(:,k), '--');
  plot(T, X_4(:,k), '--');
  grid on;
  if k < 4
    ylabel(sprintf('%s_E [m]', char(119+k)));
  else
    ylabel(sprintf('\\phi_{E,%d} [rad]', k));
  end
  if k == 5
    legend({'5FG, seriell', '5FG, elim.', '4FG hybrid', '4FG hybrid M2'});
  end
end
linkxaxes

% Gelenktrajektorie
figure(12);clf;
for k = 1:RS(3).NQJ
  subplot(3,RS(3).NQJ,sprc2no(3,RS(3).NQJ,1,k));hold on;
  plot(T, Q_3(:,k)/RS(3).qunitmult_eng_sci(k));
  plot([0;T(end)], RS(3).qlim(k,1)*[1;1]/RS(3).qunitmult_eng_sci(k), 'r--');
  plot([0;T(end)], RS(3).qlim(k,2)*[1;1]/RS(3).qunitmult_eng_sci(k), 'r--');
  xlabel('t [s]');
  ylabel(sprintf('q_%d / %s', k, RS(3).qunit_eng{k}));
  grid on;
  title(sprintf('Zeitverlauf Gelenkgrößen Achse %d',k));
  subplot(3,RS(3).NQJ,sprc2no(3,RS(3).NQJ,2,k));hold on;
  plot(T, QD_3(:,k)/RS(3).qunitmult_eng_sci(k));
  plot([0;T(end)], RS(3).qDlim(k,1)*[1;1]/RS(3).qunitmult_eng_sci(k), 'r--');
  plot([0;T(end)], RS(3).qDlim(k,2)*[1;1]/RS(3).qunitmult_eng_sci(k), 'r--');
  xlabel('t [s]');
  ylabel(sprintf('qD_%d / %s/s', k, RS(3).qunit_eng{k}));
  grid on;
  subplot(3,RS(3).NQJ,sprc2no(3,RS(3).NQJ,3,k));hold on;
  plot(T, QDD_3(:,k)/RS(3).qunitmult_eng_sci(k));
  xlabel('t [s]');
  ylabel(sprintf('qDD_%d / %s/s^2', k, RS(3).qunit_eng{k}));
  grid on;
end
linkxaxes
