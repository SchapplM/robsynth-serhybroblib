% CAD-Modell für Roboter Kuka LBR 4+ vorbereiten
%
% Eingabe:
% RS
%   Instanz der Roboterklasse, initialisiert für Roboter
% Name
%   Name des Robotermodells (entsprechend der Unterordner von `systems`)
% RobName
%   Name der Roboterparameter. Das CAD-Modell bezieht sich auf die
%   Roboterparameter mit dem gleichen Namen in der Datei `models.csv`
% 
% Ausgabe:
% RS
%   Roboterklasse mit initialisierten Pfaden zum CAD-Modell
%
% Siehe auch: hybroblib_create_robot_class.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RS = palh1m1YR1_init_CAD(RS, Name, RobName)

% Pfad zu diesem Ordner bestimmen (wo die CAD-Dateien drin sind)
repopath = fileparts(which('hybroblib_path_init.m'));
CAD_basepath = fullfile(repopath, 'systems', Name, sprintf('CAD_%s',RobName));
% Liste der CAD-Dateinamen
CAD_filenames = {'basis.STL', 'AB.STL'};

% Transformation zwischen Ursprung des MDH-Körper-KS und dem Referenz-KS
% des STL-Modells
T_mdh_visual = NaN(4,4,2);
T_mdh_visual(:,:,1) = eye(4);
T_mdh_visual(:,:,2) = eye(4);

RS.CAD_add(fullfile(CAD_basepath, CAD_filenames{1}), 0, T_mdh_visual(:,:,2), 1e-3);
