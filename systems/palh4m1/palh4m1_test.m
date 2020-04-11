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

%% Roboter plotten
s_plot = struct( 'ks', [1:RS_DE1.NJ, RS_DE1.NJ+2], 'straight', 1);
q = [0;NaN;120;0;40]*pi/180;
q(2) = 1; % D-Parameter
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title([RS_DE1.descr, ' - Nullstellung']);
RS_TE.plot( q, s_plot );

