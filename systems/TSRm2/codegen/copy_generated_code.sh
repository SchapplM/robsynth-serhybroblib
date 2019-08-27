#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_TSRm2TE ../matlabfcn_TSRm2DE1 ../matlabfcn_TSRm2DE2 ../matlabfcn_TSRm2OL ../matlabfcn_TSRm2IC
cp $maplerepopath/codeexport/TSRm2TE/matlabfcn/* ../matlabfcn_TSRm2TE
cp $maplerepopath/codeexport/TSRm2DE1/matlabfcn/* ../matlabfcn_TSRm2DE1
cp $maplerepopath/codeexport/TSRm2DE2/matlabfcn/* ../matlabfcn_TSRm2DE2
cp $maplerepopath/codeexport/TSRm2OL/matlabfcn/* ../matlabfcn_TSRm2OL
cp $maplerepopath/codeexport/TSRm2IC/matlabfcn/* ../matlabfcn_TSRm2IC

