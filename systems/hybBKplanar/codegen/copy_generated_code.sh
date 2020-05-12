#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_hybBKplanar
echo "Kopiere Code für hybBKplanar"
cp -u $maplerepopath/codeexport/hybBKplanar/matlabfcn/*.* ../matlabfcn_hybBKplanar
cp $maplerepopath/codeexport/hybBKplanar/tmp/robot_env.sh ../matlabfcn_hybBKplanar
