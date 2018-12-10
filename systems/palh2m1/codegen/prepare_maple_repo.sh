#!/bin/bash -ex
# Kopiere alle benötigten Dateien für die Code-Generierung ins HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
#
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
# (C) Institut für mechatronische Systeme, Universität Hannover

maplerepopath=$1
this_path=$(pwd)

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;
defpath=$maplerepopath/robot_codegen_definitions
constrpath=$maplerepopath/robot_codegen_constraints

## Definitionen kopieren
cp $this_path/robot_env_palh2m1 $defpath/examples/robot_env_palh2m1

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/palh2m1_kinematic_constraints.mpl $constrpath/palh2m1_kinematic_constraints.mpl
cp $this_path/palh2m1_kinematic_constraints.mw  $constrpath/palh2m1_kinematic_constraints.mw 
