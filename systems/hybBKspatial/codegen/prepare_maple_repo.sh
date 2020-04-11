#!/bin/bash -ex
# Kopiere alle benötigten Dateien für die Code-Generierung ins HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
#
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
# (C) Institut für Mechatronische Systeme, Universität Hannover

maplerepopath=$1
this_path=$(pwd)

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;
defpath=$maplerepopath/robot_codegen_definitions
constrpath=$maplerepopath/robot_codegen_constraints

# Definitionen  kopieren
cp $this_path/robot_env $defpath/robot_env

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/hybBKspatial_kinematic_constraints.mpl $constrpath/hybBKspatial_kinematic_constraints.mpl

# Werte für Kinematikparameter (für Modultests)
cp $this_path/hybBKspatial_kinematic_parameter_values.m $constrpath/hybBKspatial_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
cp $this_path/hybBKspatial_kinematic_constraints_matlab.m $constrpath/hybBKspatial_kinematic_constraints_matlab.m

