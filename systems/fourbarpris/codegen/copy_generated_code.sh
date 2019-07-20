#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbarprisTE ../matlabfcn_fourbarprisDE1 ../matlabfcn_fourbarprisDE2 ../matlabfcn_fourbarprisOL ../matlabfcn_fourbarprisIC
cp $maplerepopath/codeexport/fourbarprisTE/matlabfcn/* ../matlabfcn_fourbarprisTE
cp $maplerepopath/codeexport/fourbarprisDE1/matlabfcn/* ../matlabfcn_fourbarprisDE1
cp $maplerepopath/codeexport/fourbarprisDE2/matlabfcn/* ../matlabfcn_fourbarprisDE2
cp $maplerepopath/codeexport/fourbarprisOL/matlabfcn/* ../matlabfcn_fourbarprisOL
cp $maplerepopath/codeexport/fourbarprisIC/matlabfcn/* ../matlabfcn_fourbarprisIC

