#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbarprisTE ../matlabfcn_fourbarprisDE1 ../matlabfcn_fourbarprisDE2 ../matlabfcn_fourbarprisOL ../matlabfcn_fourbarprisIC
cp -u $maplerepopath/codeexport/fourbarprisTE/matlabfcn/*.* ../matlabfcn_fourbarprisTE
cp -u $maplerepopath/codeexport/fourbarprisDE1/matlabfcn/*.* ../matlabfcn_fourbarprisDE1
cp -u $maplerepopath/codeexport/fourbarprisDE2/matlabfcn/*.* ../matlabfcn_fourbarprisDE2
cp -u $maplerepopath/codeexport/fourbarprisOL/matlabfcn/*.* ../matlabfcn_fourbarprisOL
cp -u $maplerepopath/codeexport/fourbarprisIC/matlabfcn/*.* ../matlabfcn_fourbarprisIC

