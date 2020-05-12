#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mdlext="TE DE1 DE2 IC OL"
for ext in $mdlext; do
  echo "Kopiere Code für palh3m1$ext"
  mkdir -p ../matlabfcn_palh3m1${ext}
  cp -u $maplerepopath/codeexport/palh3m1${ext}/matlabfcn/*.* ../matlabfcn_palh3m1${ext}
  cp -u $maplerepopath/codeexport/palh3m1${ext}/tmp/robot*.sh ../matlabfcn_palh3m1${ext}
done
