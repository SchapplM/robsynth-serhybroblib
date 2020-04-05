#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_picker2Dm1TE ../matlabfcn_picker2Dm1DE1 ../matlabfcn_picker2Dm1DE2 ../matlabfcn_picker2Dm1OL ../matlabfcn_picker2Dm1IC
cp -u $maplerepopath/codeexport/picker2Dm1TE/matlabfcn/*.* ../matlabfcn_picker2Dm1TE
cp -u $maplerepopath/codeexport/picker2Dm1DE1/matlabfcn/*.* ../matlabfcn_picker2Dm1DE1
cp -u $maplerepopath/codeexport/picker2Dm1DE2/matlabfcn/*.* ../matlabfcn_picker2Dm1DE2
cp -u $maplerepopath/codeexport/picker2Dm1OL/matlabfcn/*.* ../matlabfcn_picker2Dm1OL
cp -u $maplerepopath/codeexport/picker2Dm1IC/matlabfcn/*.* ../matlabfcn_picker2Dm1IC

