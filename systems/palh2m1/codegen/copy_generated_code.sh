#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_palh2m1DE ../matlabfcn_palh2m1OL ../matlabfcn_palh2m1IC
#cp $maplerepopath/codeexport/palh2m1TE/matlabfcn/* ../matlabfcn_palh2m1TE
#cp $maplerepopath/codeexport/palh2m1DE1/matlabfcn/* ../matlabfcn_palh2m1DE1
cp $maplerepopath/codeexport/palh2m1DE2/matlabfcn/* ../matlabfcn_palh2m1DE
cp $maplerepopath/codeexport/palh2m1OL/matlabfcn/* ../matlabfcn_palh2m1OL
cp $maplerepopath/codeexport/palh2m1IC/matlabfcn/* ../matlabfcn_palh2m1IC

