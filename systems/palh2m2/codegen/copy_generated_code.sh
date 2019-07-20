#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_palh2m2TE ../matlabfcn_palh2m2DE1 ../matlabfcn_palh2m2DE2 ../matlabfcn_palh2m2OL ../matlabfcn_palh2m2IC
cp $maplerepopath/codeexport/palh2m2TE/matlabfcn/* ../matlabfcn_palh2m2TE
cp $maplerepopath/codeexport/palh2m2DE1/matlabfcn/* ../matlabfcn_palh2m2DE1
cp $maplerepopath/codeexport/palh2m2DE2/matlabfcn/* ../matlabfcn_palh2m2DE2
cp $maplerepopath/codeexport/palh2m2OL/matlabfcn/* ../matlabfcn_palh2m2OL
cp $maplerepopath/codeexport/palh2m2IC/matlabfcn/* ../matlabfcn_palh2m2IC

