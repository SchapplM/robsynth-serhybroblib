#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbar1TE ../matlabfcn_fourbar1DE1 ../matlabfcn_fourbar1DE2 ../matlabfcn_fourbar1OL ../matlabfcn_fourbar1IC
cp $maplerepopath/codeexport/fourbar1TE/matlabfcn/* ../matlabfcn_fourbar1TE
cp $maplerepopath/codeexport/fourbar1DE1/matlabfcn/* ../matlabfcn_fourbar1DE1
cp $maplerepopath/codeexport/fourbar1DE2/matlabfcn/* ../matlabfcn_fourbar1DE2
cp $maplerepopath/codeexport/fourbar1OL/matlabfcn/* ../matlabfcn_fourbar1OL
cp $maplerepopath/codeexport/fourbar1IC/matlabfcn/* ../matlabfcn_fourbar1IC

