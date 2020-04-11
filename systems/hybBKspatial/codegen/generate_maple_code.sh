#!/bin/bash -e
# Starte die Code-Generierung für den Roboter
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-07
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

## Definitionen für Unterschiedliche Modellierungen kopieren

maplerepopath=$1

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;

cd $maplerepopath
# Starte direkt. Die einzige Definitionsdatei hat schon den richtigen Namen
./robot_codegen_start.sh -p --fixb_only --kinematics_only

