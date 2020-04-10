% Teste Klassendefinition für drehbar aufgestellte Viergelenkkette in
% unterschiedlichen Implementierungen.
% Das System dient zum Testen eines Teils eines Paletierroboters
%
% Verfügbare Beispiele: Siehe Tabelle models.csv
% * Grashoff-Bedingung: Hunt1978, Gl. (3.11) (S. 82)

% Quelle:
% [Hunt1978] Hunt: Kinematic geometry of mechanisms,1978, Oxford Uni. Press
% https://de.wikipedia.org/wiki/Kurbelschwinge

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('fourbar1turn', 'TE', 'fourbar1turnBsp1');
RS_DE1 = hybroblib_create_robot_class('fourbar1turn', 'DE1', 'fourbar1turnBsp1');
RS_DE2 = hybroblib_create_robot_class('fourbar1turn', 'DE2', 'fourbar1turnBsp1');

%% Vergleich der Implementierungen
TSS = RS_TE.gen_testsettings();
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
fprintf('Kinematik der drehbaren Viergelenkkette für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

% Ab jetzt nur noch mit einem Modell weiterrechnen und dieses Parametrieren
RS = copy(RS_TE);
%% End-Effektor definieren
% Kurbelschwinge: Das Endeffektor-Element ist die Schwinge (der Stab
% mit Gestellbefestigung ohne Aktuierung
l = RS.pkin(4);
r_N_E = [-l/2;0;0];
RS.update_EE(r_N_E);

%% Gelenk-Trajektorie mit einem Umlauf 
QE = RS.qlim';

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS.fkineEE(Q(ii,:)');
  X(ii,:) = RS.t2x(T_0_Ei);
end

%% Plotten
s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0);
% q = 0;
q = pi/180*[-15; 30;];
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
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
RS.anim( Q(1:50:end,:), s_anim, s_plot);