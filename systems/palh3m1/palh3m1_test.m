% Teste Klassendefinition für Palettierroboter mit zwei Parallelogrammen
% Modell: Palh3m1 (Herleitung über allgemeine Viergelenkkette)
% Beispiel: KUKA KR-700PA (Parameter noch nicht validiert)
% Ergebnis:
% * Unterschiedliche Implementierungen haben gleiches Ergebnis
% * Modell sieht plausibel aus

% Quelle:
% * Quelle des Modells: [Shan2019_S828]
% * Vergleich der Implementierungen für diesen Roboter: TE, DE1, DE2
%   Powerpoint-Präsentation der Publikation
%   "Kinematics and Dynamics Model via Explicit Direct and Trigonometric
%   Elimination of Kinematic Constraints" (Schappler et al. 2019)
%   https://www.researchgate.net/publication/337759816_Presentation_Slides_IFToMM_World_Congress

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('palh3m1', 'TE', 'palh3m1Bsp1');
RS_DE1 = hybroblib_create_robot_class('palh3m1', 'DE1', 'palh3m1Bsp1');
RS_DE2 = hybroblib_create_robot_class('palh3m1', 'DE2', 'palh3m1Bsp1');

TSS = RS_TE.gen_testsettings();

%% Vergleich der Implementierungen
for i = 1:TSS.n
  q=TSS.Q(i,:)';
  T_DE1 = RS_DE1.fkine(q);
  T_DE2 = RS_DE2.fkine(q);
  T_TE = RS_TE.fkine(q);
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
 if any(abs([test1(:); test2(:)]) > 1e-10)
   error('Methoden DE1, DE2 und TE stimmen nicht überein');
 end
end

fprintf('Kinematik vom Typ KUKA KR-700PA für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

%% Trajektorie zwischen Gelenkwinkelgrenzen erzeugen
QE = RS_TE.qlim';

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS_TE.fkineEE(Q(ii,:)');
  X(ii,:) = RS_TE.t2x(T_0_Ei);
end

%% CAD-Modell plotten (falls vorhanden)
% s_plot = struct( 'ks', [1:RS_TE.NJ, RS_TE.NJ+2], 'mode', 2);
% q = pi/180*[0; -15; 30; 0];
% figure(5);clf;
% hold on;grid on;
% xlabel('x [m]');ylabel('y [m]');
% zlabel('z [m]');view(3);
% cadhdl=RS_TE.plot( q, s_plot );
% title(sprintf('CAD-Modell (%s)', RS_TE.descr));
%% Gelenkmomentenverlauf berechnen
nt = size(Q,1);
TAU = NaN(nt, RS_TE.NQJ);
for i = 1:nt
  q_i = Q(i,:)';
  qD_i = QD(i,:)';
  qDD_i = QDD(i,:)';
  
  tau_i = RS_TE.invdyn(q_i, qD_i, qDD_i);
  TAU(i,:) = tau_i;
end
%% Plotten
figure(1);clf;
subplot(4,1,1);
plot(t, 180/pi*Q);
xlabel('t [s]');
ylabel('q [deg]');
grid on;
subplot(4,1,2);
plot(t, 180/pi*QD);
xlabel('t [s]');
ylabel('qD [deg/s]');
grid on;
subplot(4,1,3);
plot(t, 180/pi*QDD);
xlabel('t [s]');
ylabel('qDD [deg/s^2]');
grid on;
subplot(4,1,4);
plot(t, TAU);
xlabel('t [s]');
ylabel('tau [Nm]');
grid on;

figure(20);clf;
subplot(2,1,1);
plot(t, X(:,1:3));
ylabel('EE-Position');
subplot(2,1,2);
plot(t, X(:,4:6));
ylabel('EE-Orientierung');
legend({'phi_x', 'phi_y', 'phi_z'});

s_plot = struct( 'ks', [1:RS_TE.NJ, RS_TE.NJ+2], 'straight', 0);
q = [0;60;60;0]*pi/180;
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS_TE.descr);
RS_TE.plot( q, s_plot );

resdir = fileparts(which('palh3m1_test.m'));
s_anim = struct( 'gif_name', fullfile(resdir, 'palh3m1.gif'));
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS_TE.descr);
RS_TE.anim( Q(1:50:end,:), [], s_anim, s_plot);
