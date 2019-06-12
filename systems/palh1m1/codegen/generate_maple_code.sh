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
robot_env_palh1m1TE
robot_env_palh1m1DE1
robot_env_palh1m1DE2
"
##robot_env_palh1m1DE2
##robot_env_palh1m1IC
## robot_env_palh1m1DE1
## robot_env_palh1m1TE

cd $maplerepopath
for df in $deflist; do
  echo "Starte Generierung für $df"
  cp robot_codegen_definitions/$df robot_codegen_definitions/robot_env
  ./robot_codegen_start.sh --fixb_only --minimal --notest
done


cp robot_codegen_definitions/robot_env_palh1m1OL robot_codegen_definitions/robot_env
./robot_codegen_start.sh --fixb_only --ic --minimal