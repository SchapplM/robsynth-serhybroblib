#!/bin/bash -e
# Starte die Code-Generierung für den Roboter
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-07
# (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

## Definitionen für Unterschiedliche Modellierungen kopieren

maplerepopath=$1

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;

deflist="
robot_env_fourbar1turnTE
robot_env_fourbar1turnDE1
robot_env_fourbar1turnDE2
"

cd $maplerepopath
for df in $deflist; do
  echo "Starte Generierung für $df"
  cp robot_codegen_definitions/$df robot_codegen_definitions/robot_env
  ./robot_codegen_start.sh -p --fixb_only
done

cp robot_codegen_definitions/robot_env_fourbar1turnOL robot_codegen_definitions/robot_env
./robot_codegen_start.sh --fixb_only --ic -p
