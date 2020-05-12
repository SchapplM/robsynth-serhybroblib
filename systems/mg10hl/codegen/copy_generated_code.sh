#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mdlext="DE1 DE2 TE OL IC"
for ext in $mdlext; do
  echo "Kopiere Code für mg10hl$ext"
  mkdir -p ../matlabfcn_mg10hl${ext}
  cp -u $maplerepopath/codeexport/mg10hl${ext}/matlabfcn/*.* ../matlabfcn_mg10hl${ext}
  cp -u $maplerepopath/codeexport/mg10hl${ext}/tmp/robot*.sh ../matlabfcn_mg10hl${ext}
done
