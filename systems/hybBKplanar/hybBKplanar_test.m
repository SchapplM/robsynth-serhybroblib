% Teste Klassendefinition für Fünfgelenkkette
%
% Verfügbare Beispiele: Siehe Tabelle models.csv

% Quelle:
% [Brünger2019_M832]
% [Bejaoui2018_S749] S. 23

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Definition der Roboterklasse
RS = hybroblib_create_robot_class('hybBKplanar', '', 'hybBKplanarBsp1');
% RS.update_mdh([1;2;3;2;4])

%% End-Effektor definieren
% Nehme an, dass einer der indirekt bewegten Stäbe den Endeffektor trägt
l = RS.pkin(5);
r_N_E = [0;0;0];
RS.update_EE(r_N_E);

%% Gelenk-Trajektorie im Definitionsgebiet
% Beispiel-Trajektorie, die das Definitionsgebiet bestmöglich abfahren soll
QE = [RS.qlim(1,1), RS.qlim(2,1),0; ...
      RS.qlim(1,2), RS.qlim(2,1),0; ...
      RS.qlim(1,2), RS.qlim(2,2),0; ...
      RS.qlim(1,1), RS.qlim(2,2),0; ...
      RS.qlim(1,1), RS.qlim(2,1),0; ...
      RS.qlim(1,2), RS.qlim(2,2),0];

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS.fkineEE(Q(ii,:)');
  X(ii,:) = RS.t2x(T_0_Ei);
end

%% Definitionsgebiet
% Bestimmung des kompletten Gelenkarbeitsraums durch ausprobieren
n1=100;n2=100;n3=100;
Q1_range = linspace(-pi*1.5, pi*1.5, n1);
Q2_range = linspace(-pi*1.5, pi*1.5, n2);
% Q3_range = linspace(-pi*1.5, pi*1.5, n3);
Q_data = ones(n1*n2,4);
for i = 1:n1
  for j = 1:n2
         q = [Q1_range(i); Q2_range(j) ; 0];
         Tc_ges = RS.fkine(q);
         dg = 1; % Marker für "Teil des Definitionsgebiets"
            if any(~isreal(squeeze(Tc_ges(:))))
              % Wenn es keine Lösung gibt, wird das Ergebnis komplex
              dg = 0;
            end
         Q_data(n2*(i-1)+j, :) = [q', dg];
      end
end
% Gebiet zeichnen
figure(8);clf;hold on;
I_iO = Q_data(:,4) == 1;
I_niO = Q_data(:,4) == 0;
plot(Q_data(I_iO,1)*180/pi, Q_data(I_iO,2)*180/pi, 'gv');
plot(Q_data(I_niO,1)*180/pi, Q_data(I_niO,2)*180/pi, 'rx');
plot(Q(:,1)*180/pi, Q(:,2)*180/pi, 'k-', 'LineWidth', 3);
xlabel('q1 in [deg]');ylabel('q2 in [deg]');
legend({'Lösung', 'keine Lösung', 'Traj'});
grid on;

%% Abstände zwischen den einzelnen Koordinatensystemen
% Hierdurch kann validiert werden, ob die Kinematikparameter zu den
% Koordinatentransformationen passen
% Stimmt überein mit [Bejaoui2018_S749] S. 23
q_test = Q_data(find(I_iO==1,1),1:3)';
Tges = RS.jtraf(q_test);
for i = 1:RS.NL
  iv = RS.MDH.v(i);
  r = norm(Tges(1:3,4,i));
  fprintf('%d -> %d: Abstand %1.1f\n', iv,i,r);
end

%% Gelenkmomentenverlauf berechnen

nt = size(Q,1);
TAU = NaN(nt, RS.NQJ);
for i = 1:nt
  q_i = Q(i,:)';
  qD_i = QD(i,:)';
  qDD_i = QDD(i,:)';
end

%% Plotten
figure(1);clf;
subplot(3,1,1);
plot(t, 180/pi*Q);
xlabel('t [s]');
ylabel('q [deg]');
grid on;
subplot(3,1,2);
plot(t, 180/pi*QD);
xlabel('t [s]');
ylabel('qD [deg/s]');
grid on;
subplot(3,1,3);
plot(t, 180/pi*QDD);
xlabel('t [s]');
ylabel('qDD [deg/s^2]');
grid on;

s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0);
q = [225;90;0]*pi/180;
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view([10, 90]);
title(RS.descr);
RS.plot( q, s_plot );

s_anim = struct( 'gif_name', '');
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS.descr);
RS.anim( Q(1:50:end,:), [], s_anim, s_plot);