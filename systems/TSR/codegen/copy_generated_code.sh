#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_TSRTE ../matlabfcn_TSRDE1 ../matlabfcn_TSRDE2 ../matlabfcn_TSROL ../matlabfcn_TSRIC
cp $maplerepopath/codeexport/TSRTE/matlabfcn/* ../matlabfcn_TSRTE
cp $maplerepopath/codeexport/TSRDE1/matlabfcn/* ../matlabfcn_TSRDE1
cp $maplerepopath/codeexport/TSRDE2/matlabfcn/* ../matlabfcn_TSRDE2
cp $maplerepopath/codeexport/TSROL/matlabfcn/* ../matlabfcn_TSROL
cp $maplerepopath/codeexport/TSRIC/matlabfcn/* ../matlabfcn_TSRIC

