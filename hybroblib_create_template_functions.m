% Erstelle die Vorlagen-Funktionen für Robotermodelle aus der Bibliothek
% 
% Eingabe:
% Names [nx2 cell)
%   Cell array mit Namen der `n` zu erstellenden Robotermodelle
%   Optional: Wenn nicht angegeben, werden alle erstellt.
%   Zeilenweise Angabe von Robotermodellen. 
%   Spalte 1: Name des Systems (palh1, fourbar1, ...)
%   Spalte 2: Modellierungsmethode (DE, TE, ...)
% skip_existing
%   true: Bereits existierende tpl-Ordner überspringen (Standard)
%   false: Alle Dateien immer neu erstellen
% mex_results
%   true: Die aus Vorlagen generierten Funktionen werden zum testen direkt
%   kompiliert.
% 
% Beispiele:
% hybroblib_create_template_functions({}, false, true)
% hybroblib_create_template_functions({'mg10hl','TE'}, false, true);
% hybroblib_create_template_functions([{'palh3m2','TE'}; {'palh3m2','DE1'}], false, true);
% 
% Siehe auch: serroblib_create_template_functions.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function hybroblib_create_template_functions(Names, skip_existing, mex_results)

%% Initialisierung
old_dir = pwd();
% Prüfe Eingabeargument
repopath=fileparts(which('hybroblib_path_init.m'));
if nargin < 1 || isempty(Names)
  % Stelle Liste aller Roboter zusammen
  Names = {};
  roblist = dir(fullfile(repopath, 'systems'));
  for i = 1:length(roblist)
    if roblist(i).isdir == 0, continue; end
    if roblist(i).name(1) == '.', continue; end
    mfcndirlist = dir(fullfile(repopath, 'systems', roblist(i).name, 'matlabfcn_*'));
    for j = 1:length(mfcndirlist)
      mdlext = mfcndirlist(j).name(11+length(roblist(i).name):end);
      if strcmp(mdlext, 'IC'), continue; end % Hierfür keine Template-Funktionen definiert
      Names = [Names; {roblist(i).name, mdlext}]; %#ok<AGROW>
    end
  end
end
if nargin < 2
  skip_existing = true;
end
if nargin < 3
  mex_results = false;
end
% Kopier-/Ersetz-Befehl vorbereiten (zwei Alternativen
if ispc() % Windows
  % Der sed-Befehl unter Windows ist möglich (wenn Git for Windows
  % installiert ist), ist aber sehr langsam und wird nicht benutzt.
  sedcmd = '"C:\\Program Files\\Git\\usr\\bin\\sed.exe" -i ''s/%s/%s/g'' %s';
  copycommand = 'linewise';
else % Linux
  sedcmd = 'sed -i ''s/%s/%s/g'' %s';
  copycommand = 'sed';
end
%% Alle Vorlagen-Funktionen aus HybrDyn-Repo kopieren
% Alle Funktionen dieser Liste werden roboterspezifisch aus der Liste
% erstellt. Die Anpassung sind nur geringfügig und ermöglichen Kompilierung
function_list_copy_hybrdyn = {...
  'constr2.m', ...
  'constr2grad.m'};
function_list_copy_robotics = {...
  {'kinematics', 'invkin_traj.m'}, ...
  {'kinematics', 'invkin_eulangresidual.m'}};

% Pfad zur Maple-Dynamik-Toolbox, in der einige Vorlagen liegen
mrp = fileparts(which('hybrdyn_path_init.m'));
if isempty(mrp)
  warning('Die HybridDyn-Toolbox muss im Pfad sein (siehe README.MD)');
  return
end
% Alle Vorlagen-Funktionen aus Robotik- und Maple-Repo in Vorlagen-Ordner kopieren
for tmp = function_list_copy_hybrdyn
  tplf = tmp{1};
  copyfile(fullfile(mrp, 'robot_codegen_scripts', 'templates_num', ['robot_',tplf,'.template']), ...
           fullfile(repopath, 'template_functions') );
end
rtp = fileparts(which('robotics_toolbox_path_init.m'));
if isempty(rtp)
  warning('Die Robotik-Toolbox muss im Pfad sein (siehe README.MD)');
  return
end
for tmp = function_list_copy_robotics
  tplf = tmp{1};
  copyfile(fullfile(rtp, tplf{1}, ['robot_',tplf{2},'.template']), ...
           fullfile(repopath, 'template_functions') );
end

% Generiere die Liste der Vorlagen-Funktionen aus den tatsächlich
% existierenden Dateien. Dadurch können durch Benutzer hinzugefügte
% Funktionen auch generiert werden
fl_tmp = dir(fullfile(repopath, 'template_functions', '*.template'));
function_list = cell(1,length(fl_tmp));
for i = 1:length(fl_tmp)
  % Präfix "robot_" und Suffix ".template" entfernen
  function_list{i} = strrep(fl_tmp(i).name(7:end), '.template', '');
end
%% Gehe alle Modellnamen durch und erstelle alle Vorlagen-Funktionen
for i = 1:size(Names,1)
  Name_i = Names{i,1};
  mdlsuffix_i = Names{i,2};
  robopath=fullfile(repopath, 'systems', Name_i);
  % Bestimme Ziel-Ordner des Codes für die Vorlagen-Funktionen
  fcn_dir = fullfile(robopath, sprintf('tplfcn_%s%s', Name_i, mdlsuffix_i));  
  % Definitionsdatei lesen
  def_file = fullfile(robopath, sprintf('matlabfcn_%s%s', Name_i, mdlsuffix_i), ...
    'robot_env.sh');
  if ~exist(def_file, 'file')
    warning('%s existiert nicht. Überspringe.', def_file);
    continue;
  end
  filetext = fileread(def_file);
  % Platzhalter-Ausdrücke für diesen Roboter erhalten (aus robot_env.sh)
  subsexp_array = {'robot_name', 'RN', {}; ...
                   'robot_NQJ', 'NQJ', {}; ...
                   'robot_NJ',  'NJ',  {}; ...
                   'robot_NL',  'NL',  {}; ...
                   'robot_NMPVFIXB',  'NMPVFIXB',  {}; ...
                   'robot_NMPVFLOATB',  'NMPVFLOATB',  {}; ...
                   'robot_NKP',  'NKP',  {}; ...
                   'robot_KP',  'KPDEF',  {}; ...
                   'robot_NTAUJFIXBREGNN',  'NTAUJFIXBREGNN',  {}};
  for ii = 1:size(subsexp_array,1)
    expr = [subsexp_array{ii,1}, '=(.*)'];
    tokens = regexp(filetext,expr,'tokens','dotexceptnewline');
    subsexp_array(ii,3) = strrep(tokens{1},'"','');
    subsexp_array(ii,3) = strrep(subsexp_array(ii,3),sprintf('\r'),''); % CR-Zeichen entfernen (Windows/Linux-Problem)
    if strcmp(subsexp_array{ii,2}, 'KPDEF')
      subsexp_array{ii,3} = sprintf('pkin: %s', subsexp_array{ii,3});
    end
  end
  subsexp_array{10,2} = 'VERSIONINFO';
  subsexp_array{10,3} = 'Generated in SerHybRobLib from template file';
  
  % Kopiere alle Vorlagen-Funktionen an die Ziel-Orte und Ersetze die
  % Platzhalter-Ausdrücke
  if exist(fcn_dir, 'file')
    if skip_existing
      continue % Ordner existiert schon. Überspringe.
    end
  else
    mkdir(fcn_dir);  % Ordner existiert noch nicht. Neu erstellen.
  end
  cd(fcn_dir); % In Ordner wechseln für kürzeren sed-Befehl (und zum Finden der Dateien)
  for tmp = function_list
    tplf = tmp{1};
    file1=fullfile(repopath, 'template_functions', ['robot_',tplf,'.template']);
    file2=fullfile(fcn_dir, [Name_i,mdlsuffix_i,'_',tplf]);
    % Ersetzungsausdruck für Dateinamen vorbereiten
    subsexp_array{11,2} = 'FN';
    [~,subsexp_array{11,3},~] = fileparts(file2);
    % Kopieren und Ausdrücke ersetzen
    if strcmp(copycommand, 'sed')
      % Variante für Linux: Vorlagen-Datei kopieren und anschließend
      % Text-Ersetzung mit sed-Befehl durchführen
      copyfile(file1, file2);
      for ii = 1:size(subsexp_array,1)
        sedcmd_ii = sprintf(sedcmd, ['%',subsexp_array{ii,2},'%'], subsexp_array{ii,3}, file2);
        system(sedcmd_ii);
      end
    elseif strcmp(copycommand, 'linewise')
      % Variante für Windows: Text ersetzen mit sed ist sehr langsam.
      % Nutze daher die Zeilenweise kopier-Funktion von Matlab.
      fid1 = fopen(file1);
      fid2 = fopen(file2, 'w');
      tline = fgetl(fid1);
      while ischar(tline)
        % Zeile weiterverarbeiten: Platzhalter-Ausdrücke ersetzen
        for ii = 1:size(subsexp_array,1)
          tline = strrep(tline, ['%',subsexp_array{ii,2},'%'], subsexp_array{ii,3});
        end
        fwrite(fid2, [tline, newline()]); % Zeile in Zieldatei schreiben
        tline = fgetl(fid1); % nächste Zeile
      end
      fclose(fid1);fclose(fid2);
    else
      error('Befehl nicht definiert');
    end
    % fprintf('%d/%d: Vorlagen-Funktion %s erstellt.\n', i, length(Names), tplf);
    
    % Prüfe, ob die Erstellung erfolgreich war (falls die Code-Generierung
    % nicht vollständig durchgeführt wurde, gehen nicht alle Funktionen)
    fid2 = fopen(file2, 'r');
    tline = fgetl(fid2);
    function_invalid = false;
    while ischar(tline)
      if contains(tline, 'NOTDEFINED') % Wird in HybrDyn eingesetzt, falls Ausdruck nicht gefunden.
        function_invalid = true;
        break;
      end
      tline = fgetl(fid2);
    end
    fclose(fid2);
    if function_invalid
      fprintf('%d/%d: Datei %s konnte nicht erzeugt werden (Einige Variablen nicht definiert)\n', ...
        i, size(Names,1), tplf);
      delete(file2);
    end
  end
  fprintf('%d/%d: Vorlagen-Funktionen für %s erstellt.\n', i, size(Names,1), Name_i);
  
  % Testen: Kompilieren aller Funktionen im Zielordner
  if mex_results
    addpath(fileparts(def_file)); % Symbolisch generierter Code muss im Pfad liegen
    mex_all_matlabfcn_in_dir(fcn_dir);
  end
end
% Zurückwechseln in vorheriges Verzeichnis
cd(old_dir);
