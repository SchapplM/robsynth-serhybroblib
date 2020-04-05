#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_picker2Dm2TE ../matlabfcn_picker2Dm2DE1 ../matlabfcn_picker2Dm2DE2 ../matlabfcn_picker2Dm2OL ../matlabfcn_picker2Dm2IC
cp -u $maplerepopath/codeexport/picker2Dm2TE/matlabfcn/*.* ../matlabfcn_picker2Dm2TE
cp -u $maplerepopath/codeexport/picker2Dm2DE1/matlabfcn/*.* ../matlabfcn_picker2Dm2DE1
cp -u $maplerepopath/codeexport/picker2Dm2DE2/matlabfcn/*.* ../matlabfcn_picker2Dm2DE2
cp -u $maplerepopath/codeexport/picker2Dm2OL/matlabfcn/*.* ../matlabfcn_picker2Dm2OL
cp -u $maplerepopath/codeexport/picker2Dm2IC/matlabfcn/*.* ../matlabfcn_picker2Dm2IC

