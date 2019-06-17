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
cp $this_path/robot_env_palh2m2 $defpath/robot_env_palh2m2DE
sed -i "s/palh2m2/palh2m2DE/g" $defpath/robot_env_palh2m2DE

 cp $this_path/robot_env_palh2m2OL $defpath/robot_env_palh2m2OL
 
 cp $this_path/robot_env_palh2m2IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/palh2m2_kinematic_constraints.mpl $constrpath/palh2m2DE_kinematic_constraints.mpl
cp $this_path/palh2m2_kinematic_constraints.mw  $constrpath/palh2m2DE_kinematic_constraints.mw 
cp $this_path/palh2m2IC_kinematic_constraints_implicit.mpl $constrpath/palh2m2IC_kinematic_constraints_implicit.mpl
cp $this_path/palh2m2IC_kinematic_constraints_implicit.mw  $constrpath/palh2m2IC_kinematic_constraints_implicit.mw 