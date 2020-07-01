% Starte die Modultests für alle generierten Robotermodelle
% 
% Ergebnis:
% * Die Modelle werden wenn möglich geneinander getestet und mit Bildern plausibilisiert
% 
% Die Pfadvariable muss immer neu definiert werden, da in den Testskripten
% Der `clear`-Befehl benutzt wird.
% 
% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

%% Viergelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar1', 'fourbar1_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar2', 'fourbar2_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbarpris', 'fourbarpris_test.m')); close all;
%% Fünfgelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fivebar1', 'fivebar1_test.m')); close all;

%% Minimalbeispiele und Testsysteme
% Als Teil der größeren Roboter zum Testen und zur Entwicklung
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar1turn', 'fourbar1turn_test.m')); close all;

%% Palettierroboter mit drei geschlossenen Ketten (Typ Yaskawa MPL800)
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m1', 'palh1m1_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m2', 'palh1m2_test.m')); close all;

%% Palettierroboter mit zwei geschlossenen Ketten (Typ KUKA KR-700PA)
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh3m1', 'palh3m1_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh3m2', 'palh3m2_test.m')); close all;

%% Modellvergleich für Palettierroboter
test_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'test');
addpath(test_path);
run(fullfile(test_path, 'palletizer_mdlcomp')); close all;

%% 2D-Palettierroboter
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'picker2Dm1', 'picker2Dm1_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'picker2Dm2', 'picker2Dm2_test.m')); close all;

%% Beinketten für hybride PKM
% Aus MA Brünger
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'hybBKplanar', 'hybBKplanar_test.m')); close all;
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'hybBKspatial', 'hybBKspatial_test.m')); close all;

%% Alle Tests aus HybridDyn-Repo
% Starte die Tests aus dem HybridDyn-Repo mit den abgespeicherten
% Funktionen aus diesem Repo
% In diesem Repo liegen die Testskripte:
hybriddyn_path = fileparts(which('hybrdyn_path_init.m'));
testscript_names = {'varpar_fixbase_kinematics_test', ...
  'varpar_fixbase_invdyn_test', 'varpar_fixbase_paramlin_test'};
testscript_names_IC = {'test_constr_kinematics', 'test_constr_dynamics'};
roblist_kinematics_only = {'hybBKplanar', 'hybBKspatial'};
% Alle Systeme und Modellierungen durchgehen
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
roblist = dir(systems_path);
testscript_list = {}; % für zu sammelnde existierende Testskripte
% Debug: Nummer des manuell zu untersuchenden Systems herausfinden:
% find(strcmp({roblist.name}, 'palh2m1')), {roblist.name}
for i = 1:length(roblist)
  if roblist(i).isdir == 0, continue; end
  if roblist(i).name(1) == '.', continue; end
  mfcndirlist = dir(fullfile(systems_path, roblist(i).name, 'matlabfcn_*'));
  % Alle Matlab-Funktions-Verzeichnisse hinzufügen
  for j = 1:length(mfcndirlist)
    mdlext = mfcndirlist(j).name(11+length(roblist(i).name):end);
    if strcmp(mdlext, 'IC') % untere Vorgehensweise geht nicht.
      addpath(fullfile(mfcndirlist(j).folder, mfcndirlist(j).name));
    else % Benutze Klasseninitialisierung zum Erstellen der Template-Dateien
      hybroblib_create_robot_class(roblist(i).name, mdlext);
    end
  end
  % Alle Modelle durchgehen und Testskripte sammeln (Start hier nicht
  % möglich wegen `clear`-Befehl.
  for j = 1:length(mfcndirlist)
    mdlext = mfcndirlist(j).name(11+length(roblist(i).name):end);
    % Starte die Testskripte zu diesem System
    if strcmp(mdlext, 'IC') % Verschiedene Testskripte für beide Fälle
      testscript_names_j = testscript_names_IC;
    else
      % Prüfe, ob Funktion nur mit Kinematik generiert wird
      if any(strcmp(roblist_kinematics_only, roblist(i).name))
        testscript_names_j = testscript_names(1);
      else
        testscript_names_j = testscript_names;
      end
    end
    for k = 1:length(testscript_names_j)
      rn = [roblist(i).name,mdlext];
      ts = fullfile(hybriddyn_path, 'codeexport', rn, 'testfcn', ...
        [rn, '_', testscript_names_j{k}, '.m']);
      if exist(ts, 'file')
        testscript_list = [testscript_list, ts]; %#ok<AGROW>
      else
        warning('%s: Testskript %s existiert aktuell nicht. Muss kein Fehler sein.', rn, testscript_names{k});
      end
    end
  end
end
fprintf('%d Testskripte gesammelt. Führe alle aus.\n', length(testscript_list));
% Alle Testskripte abarbeiten. Muss am Ende passieren, da jedes Testskript
% am Anfang `clear` ausführt. Nur die Schleifenvariable `ts` bleibt.
for ts = testscript_list
  run(ts{1});
end
