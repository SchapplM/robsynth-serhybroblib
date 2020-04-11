#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_palh4m1TE ../matlabfcn_palh4m1DE1 ../matlabfcn_palh4m1DE2 ../matlabfcn_palh4m1OL ../matlabfcn_palh4m1IC
cp -u $maplerepopath/codeexport/palh4m1TE/matlabfcn/*.* ../matlabfcn_palh4m1TE
cp -u $maplerepopath/codeexport/palh4m1DE1/matlabfcn/*.* ../matlabfcn_palh4m1DE1
cp -u $maplerepopath/codeexport/palh4m1DE2/matlabfcn/*.* ../matlabfcn_palh4m1DE2
cp -u $maplerepopath/codeexport/palh4m1OL/matlabfcn/*.* ../matlabfcn_palh4m1OL
cp -u $maplerepopath/codeexport/palh4m1IC/matlabfcn/*.* ../matlabfcn_palh4m1IC

