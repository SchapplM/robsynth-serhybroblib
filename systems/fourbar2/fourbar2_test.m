% Teste Klassendefinition für Parallelogramm in unterschiedlichen
% Implementierungen
% 
% Ergebnis: Berechnung der allgemeinen Viergelenkkette aus
% Kreisschnittpunkt stimmt für Geometrie des Parallelogramms überein.
% 
% TOOD: Bei Winkeln außerhalb von (0, pi] klappt die Viergelenkkette um
% (Anti-Parallelogramm)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('fourbar2', 'TE', 'fourbar2Bsp1');
RS_DE1 = hybroblib_create_robot_class('fourbar2', 'DE1', 'fourbar2Bsp1');
RS_DE2 = hybroblib_create_robot_class('fourbar2', 'DE2', 'fourbar2Bsp1');

%% Vergleich der Implementierungen des Parallelogramms
for q = linspace(RS_DE1.qlim(1), RS_DE1.qlim(2), 1000)
  T_DE1 = RS_DE1.fkine(q);
  T_DE2 = RS_DE2.fkine(q);
  T_TE = RS_TE.fkine(q);
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
  if any(abs([test1(:); test2(:)]) > 1e-10)
    error('Methoden DE1, DE2 und TE stimmen nicht überein');
  end
end
fprintf('Kinematik der Viergelenkkette für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

% Ab jetzt nur noch mit einem Modell weiterrechnen und dieses Parametrieren
RS = copy(RS_TE);
%% End-Effektor definieren
% Längen des Parallelogramms
l1 = RS.pkin(1);
l2 = RS.pkin(2);
% Bewegtes Segment ist gegenüber der Basis
r_N_E = [l1/2;0;0];
RS.update_EE(r_N_E);

%% Vergleiche mit allgemeiner Viergelenkkette
RS_m1 = hybroblib_create_robot_class('fourbar1', 'TE');
RS_m1.descr = 'Viergelenkkette';
RS_m1.update_mdh([l1;l2;l1;l2]);
RS_m1.I_EElink = uint8(2);
RS_m1.update_EE(r_N_E);
for q = linspace(0, pi, 1000) % funktioniert für kleiner 0 oder größer pi noch nicht.
  q_vgk = q;
  T_m2 = RS.fkine(q_vgk);
  T_m1 = RS_m1.fkine(q);
  test12 = T_m2-T_m1;
  if any(abs(test12(:)) > 1e-3)
    figure(20);clf;
    for k = 1:2
      if k == 1, RS_plot = RS;q_plot = q; % Parallelogramm
      else,      RS_plot = RS_m1;q_plot = q_vgk; end % Viergelenkkette
      subplot(1,2,k);
      s_plot = struct( 'ks', [1:RS_plot.NJ, RS_plot.NJ+2], 'straight', 0);
      hold on; grid on;
      xlabel('x in m'); ylabel('y in m');zlabel('z in m');
      view([0 90]);
      title(RS_plot.descr);
      RS_plot.plot( q_plot, s_plot );
    end
    error('Die allgemeine Viergelenkkette und das Parallelogramm stimmen nicht überein');
  end
end
fprintf('Kinematik von Viergelenkkette und Parallelogramm für 1000 Kombinationen verglichen\n');

%% Gelenk-Trajektorie mit einem Umlauf 
QE = RS.qlim';

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS.fkineEE(Q(ii,:)');
  X(ii,:) = RS.t2x(T_0_Ei);
end

%% Zeichnen
s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0);
q = pi/3;
figure(2);clf;
hold on;
grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3);
title(RS.descr);
RS.plot( q, s_plot );

s_anim = struct( 'gif_name', '');
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x in m');
ylabel('y in m');
zlabel('z in m');
view(3);
title(RS.descr);
RS.anim( Q(1:50:end,:), [], s_anim, s_plot);