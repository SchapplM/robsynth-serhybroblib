#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_palh1m2TE ../matlabfcn_palh1m2DE1 ../matlabfcn_palh1m2DE2 ../matlabfcn_palh1m2OL ../matlabfcn_palh1m2IC
cp -u $maplerepopath/codeexport/palh1m2TE/matlabfcn/*.* ../matlabfcn_palh1m2TE
cp -u $maplerepopath/codeexport/palh1m2DE1/matlabfcn/*.* ../matlabfcn_palh1m2DE1
cp -u $maplerepopath/codeexport/palh1m2DE2/matlabfcn/*.* ../matlabfcn_palh1m2DE2
cp -u $maplerepopath/codeexport/palh1m2OL/matlabfcn/*.* ../matlabfcn_palh1m2OL
cp -u $maplerepopath/codeexport/palh1m2IC/matlabfcn/*.* ../matlabfcn_palh1m2IC

