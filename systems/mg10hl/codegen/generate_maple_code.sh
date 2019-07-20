#!/bin/bash -e
# Starte die Code-Generierung für den Roboter

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-07
# (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

## Definitionen für Unterschiedliche Modellierungen kopieren

maplerepopath=$1

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;

deflist="
robot_env_mg10hlTE
robot_env_mg10hlDE1
robot_env_mg10hlDE2
"
##robot_env_mg10hlDE2
##robot_env_mg10hlIC
## robot_env_mg10hlDE1
## robot_env_mg10hlTE

cd $maplerepopath
for df in $deflist; do
  echo "Starte Generierung für $df"
  cp robot_codegen_definitions/$df robot_codegen_definitions/robot_env
  ./robot_codegen_start.sh --fixb_only -p


cp robot_codegen_definitions/robot_env_mg10hlOL robot_codegen_definitions/robot_env
./robot_codegen_start.sh --fixb_only --ic