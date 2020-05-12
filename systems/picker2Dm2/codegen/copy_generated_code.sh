#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mdlext="DE1 DE2 TE OL IC"
for ext in $mdlext; do
  echo "Kopiere Code für picker2Dm2$ext"
  mkdir -p ../matlabfcn_picker2Dm2${ext}
  cp -u $maplerepopath/codeexport/picker2Dm2${ext}/matlabfcn/*.* ../matlabfcn_picker2Dm2${ext}
  cp -u $maplerepopath/codeexport/picker2Dm2${ext}/tmp/robot*.sh ../matlabfcn_picker2Dm2${ext}
done
