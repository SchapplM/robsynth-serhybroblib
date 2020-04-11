% Teste palh4m1

% Quelle:

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('palh4m1', 'TE', 'palh4m1Bsp1');
RS_DE1 = hybroblib_create_robot_class('palh4m1', 'DE1', 'palh4m1Bsp1');
RS_DE2 = hybroblib_create_robot_class('palh4m1', 'DE2', 'palh4m1Bsp1');
RS_OL = hybroblib_create_robot_class('palh4m1', 'OL', 'palh4m1Bsp1');
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
   warning('Methoden DE1, DE2 und TE stimmen nicht überein');
 end
end
fprintf('Kinematik des Palettierers mit einer geschlossenen Kette für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

%% Roboter plotten
s_plot = struct( 'ks', [1:RS_TE.NJ, RS_TE.NJ+2], 'straight', 0);
q = [180;NaN;120;0;40]*pi/180;
q(2) = .3; % d-Parameter
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title([RS_TE.descr, ' - Nullstellung']);
RS_TE.plot( q, s_plot );
Tc_CL = RS_TE.fkine(q);
view([0 0])

% Prüfe die Viergelenkkette als ein Teil des Roboters
% TODO: Das ist noch nicht fertig:
RS_fb1 = hybroblib_create_robot_class('fourbar1', 'TE');
% Trage die Parameter der Viergelenkkette des Roboters in die reine
% Viergelenkkette ein. TODO: Noch unsicher, ob so richtig.
% In Maple stimmen die Variablen so überein.
pkin_fb1 = [RS_TE.pkin(strcmp(RS_TE.pkin_names, 'AB')); ...
            RS_TE.pkin(strcmp(RS_TE.pkin_names, 'AD')); ...
            q(2)+RS_TE.pkin(strcmp(RS_TE.pkin_names, 'HC')); ...
            RS_TE.pkin(strcmp(RS_TE.pkin_names, 'CB'))];
RS_fb1.update_mdh(pkin_fb1);
% RS_fb1.update_base(Tc_CL(1:3,4,2), r2eulxyz(Tc_CL(1:3,1:3,2)));
s_plot_fb1 = struct( 'ks', [1:RS_fb1.NJ, RS_fb1.NJ+2], 'straight', 0);
figure(102);clf;
hold on; grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3); title('Viergelenkkette mit Werten des Roboters');
RS_fb1.plot(q(5)+pi/2, s_plot_fb1)
view([0 90])
% Winkel der Viergelenkkette:
RS_fb1.jointvar(q(5)+pi/2);
Tc_fb1 = RS_fb1.fkine(q(5)+pi/2);
% Vergleiche Längen von Viergelenkkette und Roboter. TODO: Das ist noch falsch
% l_AB_fb1 = norm(Tc_fb1(1:3,4,2)-Tc_fb1(1:3,4,3));
% l_AB_rob = norm(Tc_CL(1:3,4,8)-Tc_CL(1:3,4,9));

% Neuen Roboter mit offener Kette definieren und neue Gelenkwinkel testen
% TODO: Die korrigierten Gelenkwinkel der Viergelenkkette hier einsetzen.
% Das geht aber erst, wenn dort sinnvolle Werte herauskommen.
figure(103);clf;
hold on; grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3); title('Offene Baumstruktur zum Testen');
q_jv = RS_TE.jointvar(q);
q_OL = q_jv(RS_TE.MDH.sigma~=2);
q_OL(7) = q_OL(7)-pi/2;
Tc_OL = RS_OL.fkine(q_OL);
RS_OL.plot( q_OL, s_plot );
view([0 0])