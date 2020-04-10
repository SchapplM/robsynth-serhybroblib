% Teste Klassendefinition für 2D-Palettierroboter
% Modell: picker2Dm2 (Herleitung über Parallelogramm)
% Ergebnis:
% * Herleitung über Allgemeine Viergelenkkette (picker2Dm1) und Parallelogramm
% (picker2Dm2) stimmt überein
% * Unterschiedliche Implementierungen haben gleiches Ergebnis

% Quelle:
%   Powerpoint-Präsentation der Publikation
%   "Kinematics and Dynamics Model via Explicit Direct and Trigonometric
%   Elimination of Kinematic Constraints" (Schappler et al. 2019)
%   https://www.researchgate.net/publication/337759816_Presentation_Slides_IFToMM_World_Congress

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('picker2Dm2', 'TE', 'TSR2040');
RS_DE1 = hybroblib_create_robot_class('picker2Dm2', 'DE1', 'TSR2040');
RS_DE2 = hybroblib_create_robot_class('picker2Dm2', 'DE2', 'TSR2040');

TSS = RS_TE.gen_testsettings();
%% Vergleich der Implementierungen
n_test = 0;
I_io = false(TSS.n, 1);
for i = 1:TSS.n
  q=TSS.Q(i,:)';
  try
    T_DE1 = RS_DE1.fkine(q);
    T_DE2 = RS_DE2.fkine(q);
    T_TE = RS_TE.fkine(q);
  catch
    continue; % Gelenkkoordinate gehört nicht zum Definitionsgebiet
  end
  I_io(i) = true;
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
  n_test = n_test + 1;
  if any(abs([test1(:); test2(:)]) > 1e-10)
    error('Methoden DE1, DE2 und TE stimmen nicht überein');
  end
end
fprintf('Kinematik vom Typ RoboWin TSR2040 für %d Kombinationen in unterschiedlichen Implementierungen getestet. Restliche %d nicht lösbar.\n', ...
  n_test, TSS.n-n_test);

%% Vergleich von picker2Dm1 und picker2Dm2
% Das Modell picker2Dm2 ist einfacher zu berechnen. Wird aber gegen picker2Dm1
% getestet, weil dieses zuerst erstellt wurde
RS_picker2Dm1 = hybroblib_create_robot_class('picker2Dm1', 'TE', 'TSR2040');
for i = find(I_io)'
  q=TSS.Q(i,:)';
  T_m1 = RS_DE1.fkine(q);
  T_m2 = RS_picker2Dm1.fkine(q);
  test_m1m2 = T_m1-T_m2;
  if any(abs(test_m1m2(:)) > 1e-4) % grobe Toleranz, da wenige Nachkommastellen in models.csv
    error('Palettierer picker2Dm1 stimmt nicht gegen picker2Dm2');
  end
end
fprintf('Modelle picker2Dm1 und picker2Dm2 stimmen überein\n');