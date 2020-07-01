% Instanz der Roboterklasse SerRob für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des Robotermodells, wie es im Ordner `systems` gespeichert ist
% mdlsuffix
%   Suffix für die Implementierung des Robotermodells, entsprechend den
%   verwendeten Namenskonventionen aus den Unterordnern von `systems`
% RobName
%   Name der Roboterparameter entsprechend der ersten Spalte der Tabelle
%   models.csv des jeweiligen Robotermodells
%   Eine Initialisierung ohne die mit `RobName` verbundenen Parameter ist
%   möglich, aber nicht sinnvoll, da die hybriden Roboter für beliebige
%   Bemaßung und Gelenkstellung nicht berechenbar sind (Zwangsbedingungen)
% 
% Ausgabe:
% RS [SerRob]
%   Instanz der SerRob-Klasse mit den Eigenschaften und Funktionen des zu
%   ladenden Roboters
% 
% Siehe auch: serroblib_create_robot_class.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
% (C) Institut für mechatronische Systeme, Universität Hannover

function RS = hybroblib_create_robot_class(Name, mdlsuffix, RobName)

if nargin < 2 || isempty(mdlsuffix)
  mdlsuffix = '';
end
if nargin < 3
  RobName = '';
end

%% Initialisierung
repopath=fileparts(which('hybroblib_path_init.m'));
robopath=fullfile(repopath, 'systems', Name);
if exist(robopath, 'file')
  addpath(robopath);
end
% Pfad für gewünschte Modellimplementierung hinzufügen
% Viergelenkkette
apd = fullfile(robopath, sprintf('matlabfcn_%s%s',Name, mdlsuffix));
if exist(apd, 'file')
  addpath(apd);
end
tpd = fullfile(robopath, sprintf('tplfcn_%s%s',Name, mdlsuffix));
if ~exist(tpd, 'file')
  % Template-Funktionen existieren nicht. Erzeuge neu.
  hybroblib_create_template_functions({Name, mdlsuffix});
end
addpath(tpd);
%% Klassen-Instanz initialisieren
% Strukturinformationen über das Modell holen
eval(sprintf('[v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = %s%s_structural_kinematic_parameters();', Name, mdlsuffix));
NJ = length(v_mdh);

% Dummy-Parameterstruktur erstellen
PS = struct('beta',  NaN(NJ,1), 'b', NaN(NJ,1), ...
            'alpha', NaN(NJ,1), 'a', NaN(NJ,1), ...
            'theta', NaN(NJ,1), 'd', NaN(NJ,1), ...
            'sigma', sigma_mdh, 'offset', NaN(NJ,1), ...
            'pkin', NaN(NKP,1), 'v', v_mdh, ...
            'mu', mu_mdh, ...
            'm', NaN(NL,1), 'mrSges', NaN(NL,3), 'Ifges', NaN(NL,6), ...
            'NJ', NJ, 'NL', NL, 'NQJ', NQJ, ...
            'qmin', NaN(NQJ,1), 'qmax', NaN(NQJ,1), 'vmax', NaN(NQJ,1), 'qref', NaN(NQJ,1));

%% Kinematik-Parameter aus Tabelle der Modellparameter laden
% In diesem Abschnitt belegte Variablen mit Standardwerten initialisieren:
descr = '';
value_EElink = NL-1; % Standardmäßig letztes Segment
pkin_names_csv = {}; % leer initialisieren
% Parameter aus Tabelle holen
if ~isempty(RobName) % Falls Name des parametrierten Modells gegeben
  found = false; % Marker, ob Parametersatz gefunden wurde
  pardat = fullfile(robopath, 'models.csv');
  if ~exist(pardat, 'file')
    error('Parameterdatei zu %s existiert nicht: %s', Name, pardat);
  end
  % Suche nach gefordertem Roboter in Parameterdatei
  fid = fopen(pardat);
  tline = fgetl(fid); % csv-Datei zeilenweise einlesen
  lc = 0;
  while ischar(tline)
    lc = lc + 1;
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    if lc == 2
      csv_head2 = csvline; % zweite Kopfzeile. Hier steht "pkin1, pkin2, ..."
    elseif lc == 3
      csv_head3 = csvline; % dritte Kopfzeile. Hier stehen die Namen der Par.
      pkin_names_csv = csv_head3(contains(csv_head2, 'pkin'));
      pkin_table = NaN(length(pkin_names_csv),1);
    end
    tline = fgetl(fid); % nächste Zeile
    if isempty(csvline) || strcmp(csvline{1}, '')
      continue
    end
    if strcmp(csvline{1}, RobName)
      found = true; % Erste Spalte der Zeile entspricht dem Roboternamen
      break;
    end
  end
  fclose(fid);
  if ~found
    warning('Parameter für %s nicht gefunden', RobName);
  else
    % Daten aus csv-Zeile extrahieren
    c = 3; % erste Drei Zeilen nicht betrachten (sind allgemeine Felder des Roboters)
    descr = sprintf('%s %s', csvline{2}, csvline{3});
    for kk = 1:length(pkin_table) % über alle Kinematik-Parameter in pkin
      c=c+1; pkin_table(kk) = str2double(csvline{c});
    end
    for kk = 1:NQJ % über alle Minimalkoordinaten
      c=c+1; value_qmin   = str2double(csvline{c});
      c=c+1; value_qmax   = str2double(csvline{c});
      c=c+1; value_vmax   = str2double(csvline{c});
      c=c+1; % Spalte für Achs-Vorzeichen nach Hersteller-Norm (noch nicht implementiert)
      c=c+1; value_qref = str2double(csvline{c});
      
      % Werte belegen, wenn sie in der Tabelle nicht gegeben sind. Ist im
      % Allgemeinen nicht sinnvoll, da aufgrund der geschlossenen
      % kinematischen Ketten dann ungültige Lösungen im Gelenkarbeitsraum
      % liegen
      if isnan(value_qmin)
        value_qmin = -pi;
      end
      if isnan(value_qmax)
        value_qmax = pi;
      end
      if isnan(value_vmax)
        value_vmax = 1;
      end
      PS.qmin(kk) = value_qmin;
      PS.qmax(kk) = value_qmax;
      PS.vmax(kk) = value_vmax;
      PS.qref(kk) = value_qref;
    end
    % Allgemeine Roboterdaten
    c=c+1; value_mload   = str2double(csvline{c}); %#ok<NASGU>
    c=c+1; value_mass   = str2double(csvline{c}); %#ok<NASGU>
    c=c+1; value_EElink   = str2double(csvline{c});
  end
end
%% Kinematik-Parameter prüfen
if ~isempty(RobName) % Falls Name des parametrierten Modells gegeben
  % Gleiche die Kinematik-Parameter aus der Tabelle mit denen des Roboters ab
  if length(pkin_names_csv) ~= length(pkin_names) %#ok<USENS>
    warning('In Matlab-Klasse %d Kinematikparameter, in Tabelle aber %d.', ...
      length(pkin_names), length(pkin_names_csv));
  else
    for i = 1:length(pkin_names_csv)
      if ~strcmp(pkin_names_csv{i}, pkin_names{i})
        warning('pkin %d in Tabelle ist %s. In Matlab-Klasse aber %s', ...
          i, pkin_names_csv{i}, pkin_names{i});
      end
    end
  end
  % Ordne die Parameter automatisch anhand der Namen zu und trage sie in die
  % Klasse ein
  for i = 1:length(pkin_names)
    I_table = strcmp(pkin_names{i}, pkin_names_csv);
    if sum(I_table) == 0
      warning('Kinematik-Parameter %s nicht in Tabelle gefunden', pkin_names{i})
      continue
    elseif sum(I_table) > 1
      warning('Kinematik-Parameter %s mehrfach Tabelle gefunden', pkin_names{i})
      continue
    else
      PS.pkin(i) = pkin_table(I_table);
    end
  end
end
%% Klasse initialisieren

% Klassen-Instanz erstellen
RS = SerRob(PS, [Name,mdlsuffix]);

% Klassen-Instanz vorbereiten
RS.fill_fcn_handles(false);

% Kinematik-Parameter initialisieren
RS.update_mdh(PS.pkin)

RS.qlim = [PS.qmin(:), PS.qmax(:)];
RS.qDlim = [-PS.vmax(:), PS.vmax(:)];
RS.qref = PS.qref(:);
RS.descr = descr;

% Nummer des EE-Segmentes setzen. Bei hybriden Robotern ist das nicht
% unbedingt das letzte in der MDH-Tabelle.
RS.I_EElink = value_EElink;

%% CAD-Modelle initialisieren, falls vorhanden
% TODO: Anpassung an neues Repo mit CAD-Modelldateien
