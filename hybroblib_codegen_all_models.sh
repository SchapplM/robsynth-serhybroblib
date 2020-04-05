#!/bin/bash -ex
# Starte die Code-Generierung für alle validierten Roboter dieses Repos
#
# Das Skript muss in dem Verzeichnis ausgeführt werden, in dem es liegt
# Es muss eine Datei hybrdyn_repo_path im selben Ordner angelegt werden, die den Pfad zum HybrDyn-Repo enthält (und nichts anderes)

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

hybroblib_repo_path=$(pwd)
hybrdyn_repo_path=`cat hybrdyn_repo_path`

systems="
fourbar1
fourbar2
fivebar1
palh1m1
palh1m2
palh2m1
palh2m2
palh3m1
palh3m2
picker2Dm1
picker2Dm2
"

for mdl in $systems; do
  echo "Starte Code-Generierung für alle Teilmodelle von $mdl"
  cd $hybroblib_repo_path/systems/$mdl/codegen
  ./prepare_maple_repo.sh $hybrdyn_repo_path
  
  cd $hybroblib_repo_path/systems/$mdl/codegen
  ./generate_maple_code.sh $hybrdyn_repo_path
  
  cd $hybroblib_repo_path/systems/$mdl/codegen
  ./copy_generated_code.sh $hybrdyn_repo_path
done


