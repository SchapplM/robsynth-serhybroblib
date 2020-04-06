% Teste Klassendefinition für Palettierroboter mit zwei Parallelogrammen
% Modell: Palh3m1 (Herleitung über allgemeine Viergelenkkette)
% Beispiel: KUKA KR-700PA (Parameter noch nicht validiert)
% Ergebnis:
% * Herleitung über Allgemeine Viergelenkkette (palh3m1) und Parallelogramm
% (palh3m2) stimmt überein
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
RS_TE = hybroblib_create_robot_class('palh3m2', 'TE', 'palh3m2KR1');
RS_DE1 = hybroblib_create_robot_class('palh3m2', 'DE1', 'palh3m2KR1');
RS_DE2 = hybroblib_create_robot_class('palh3m2', 'DE2', 'palh3m2KR1');

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

%% Vergleich von palh3m1 und palh3m2
% Das Modell palh3m2 ist einfacher zu berechnen. Wird aber gegen palh3m1
% getestet, weil dieses zuerst erstellt wurde
RS_palh3m1 = hybroblib_create_robot_class('palh3m1', 'TE', 'palh3m1Bsp1');
for i = 1:TSS.n
  q=TSS.Q(i,:)';
  T_m1 = RS_DE1.fkine(q);
  T_m2 = RS_palh3m1.fkine(q);
  test_m1m2 = T_m1-T_m2;
  if any(abs(test_m1m2(:)) > 1e-4) % grobe Toleranz, da wenige Nachkommastellen in models.csv
    error('Palettierer palh3m1 stimmt nicht gegen palh3m2');
  end
end