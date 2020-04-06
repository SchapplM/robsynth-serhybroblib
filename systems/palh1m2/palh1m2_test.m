% Teste Klassendefinition für Palettierroboter mit drei Parallelogrammen
% Modell: Palh1m2 (Herleitung über Parallelogramm)
% Beispiel: MPL800-Yaskawa (modifizierte Parameter)
% Ergebnis:
% * Herleitung über Allgemeine Viergelenkkette (palh1m1) und Parallelogramm
% (palh1m2) stimmt überein
% * Unterschiedliche Implementierungen haben gleiches Ergebnis

% Quelle:
% * https://www.motoman.com/industrial-robots/mpl800-ii
% * Vergleich der Implementierungen für diesen Roboter:
%   Powerpoint-Präsentation der Publikation
%   "Kinematics and Dynamics Model via Explicit Direct and Trigonometric
%   Elimination of Kinematic Constraints" (Schappler et al. 2019)
%   https://www.researchgate.net/publication/337759816_Presentation_Slides_IFToMM_World_Congress

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('palh1m2', 'TE', 'palh1m2Bsp1');
RS_DE1 = hybroblib_create_robot_class('palh1m2', 'DE1', 'palh1m2Bsp1');
RS_DE2 = hybroblib_create_robot_class('palh1m2', 'DE2', 'palh1m2Bsp1');

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
fprintf('Kinematik vom Typ Yaskawa MPL800 für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

%% Vergleich von palh1m1 und palh1m2
% Das Modell palh1m2 ist einfacher zu berechnen. Wird aber gegen palh1m1
% getestet, weil dieses zuerst erstellt wurde
RS_palh1m1 = hybroblib_create_robot_class('palh1m1', 'TE', 'palh1m1Bsp1');
for i = 1:TSS.n
  q=TSS.Q(i,:)';
  T_m1 = RS_DE1.fkine(q);
  T_m2 = RS_palh1m1.fkine(q);
  test_m1m2 = T_m1-T_m2;
  if any(abs(test_m1m2(:)) > 1e-10)
    error('Palettierer palh1m1 stimmt nicht gegen palh1m2');
  end
end