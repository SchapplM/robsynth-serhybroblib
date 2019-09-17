#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbar2TE ../matlabfcn_fourbar2DE1 ../matlabfcn_fourbar2DE2 ../matlabfcn_fourbar2OL ../matlabfcn_fourbar2IC
cp -u $maplerepopath/codeexport/fourbar2TE/matlabfcn/*.* ../matlabfcn_fourbar2TE
cp -u $maplerepopath/codeexport/fourbar2DE1/matlabfcn/*.* ../matlabfcn_fourbar2DE1
cp -u $maplerepopath/codeexport/fourbar2DE2/matlabfcn/*.* ../matlabfcn_fourbar2DE2
cp -u $maplerepopath/codeexport/fourbar2OL/matlabfcn/*.* ../matlabfcn_fourbar2OL
cp -u $maplerepopath/codeexport/fourbar2IC/matlabfcn/*.* ../matlabfcn_fourbar2IC

