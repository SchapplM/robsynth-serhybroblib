% Teste Klassendefinition für 2FG-Palettierer
%
% Verfügbare Beispiele: Siehe Tabelle models.csv
% TODO: Test picker2Dm1 gegen picker2Dm2 fehlt noch

% Quelle:
% http://en.triowin.com/tdr-series2-axisparallelrobot-15268685491267769.html

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS = hybroblib_create_robot_class('picker2Dm1', 'TE', 'TSR2040');
RS.fill_fcn_handles(true); % Nutze Mex-Funktionen, damit es schneller geht
RS.mex_dep();
TSS = RS.gen_testsettings();

% EE-Transformation
RS.update_EE([-0.1;-0.2;0]);
% EE so drehen, dass er immer nach unten zeigt (parallel zum Basis-KS)
T_test = RS.fkineEE(zeros(2,1));
phi_neu = r2eulxyz(T_test(1:3,1:3)');
RS.update_EE([0.1;0.2;0], phi_neu);

%% Beispiel-Trajektorie Gelenkraum
QE = pi/180*[  0,  0; ...
              10, 10; ...
               0, 45; ...
               0,  0; ...
              45,  0; ...
               0, 80; ...
               0,  0; ...
              15,-45; ...
               0,-30];              

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 100e-3, 10e-3, 5e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS.fkineEE(Q(ii,:)');
  X(ii,:) = RS.t2x(T_0_Ei);
end

%% Definitionsgebiet
% Bestimmung des kompletten Gelenkarbeitsraums durch ausprobieren
Q1_range = (-120:5:120)*pi/180;n1=length(Q1_range);
Q2_range = (-120:5:120)*pi/180;n2=length(Q2_range);
Q_data = ones(n1*n2,3);
for i = 1:n1
  for j = 1:n2
    q = [Q1_range(i); Q2_range(j)];
    try % Nutze try catch, da bei Definitionsverletzung bei mex-Funktionen ein Fehler kommt
      Tc_ges = RS.fkine(q);
      dg = 1; % Marker für "Teil des Definitionsgebiets"
    catch
      % Wenn es keine Lösung gibt, wird das Ergebnis komplex
      dg = 0;
    end
    Q_data(n2*(i-1)+j, :) = [q', dg];
  end
end
% Gebiet zeichnen
figure(8);clf;hold on;
I_iO = Q_data(:,3) == 1;
I_niO = Q_data(:,3) == 0;
plot(Q_data(I_iO,1)*180/pi, Q_data(I_iO,2)*180/pi, 'gv');
plot(Q_data(I_niO,1)*180/pi, Q_data(I_niO,2)*180/pi, 'rx');
plot(Q(:,1)*180/pi, Q(:,2)*180/pi, 'k-', 'LineWidth', 3);
xlabel('q1 in [deg]');ylabel('q2 in [deg]');
legend({'Lösung', 'keine Lösung', 'Traj'});
grid on;
title('Definitionsgebiet');

%% Strich-Modell plotten
s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'mode', 1);
q = pi/180*[0; -45];
figure(5);clf;
hold on;grid on;
xlabel('x [m]');ylabel('y [m]');
zlabel('z [m]');view(3);
cadhdl=RS.plot( q, s_plot );
title(sprintf('CAD-Modell (%s)', RS.descr));

%% Animation der Trajektorie
resdir = fileparts(which('picker2Dm1_test.m'));
s_anim = struct( 'gif_name', fullfile(resdir, 'picker2Dm1_trajvis.gif'));
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS.descr);
RS.anim( Q(1:50:end,:), s_anim, s_plot);
