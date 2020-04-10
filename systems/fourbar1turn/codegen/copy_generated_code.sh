#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbar1turnTE ../matlabfcn_fourbar1turnDE1 ../matlabfcn_fourbar1turnDE2 ../matlabfcn_fourbar1turnOL ../matlabfcn_fourbar1turnIC
cp -u $maplerepopath/codeexport/fourbar1turnTE/matlabfcn/*.* ../matlabfcn_fourbar1turnTE
cp -u $maplerepopath/codeexport/fourbar1turnDE1/matlabfcn/*.* ../matlabfcn_fourbar1turnDE1
cp -u $maplerepopath/codeexport/fourbar1turnDE2/matlabfcn/*.* ../matlabfcn_fourbar1turnDE2
cp -u $maplerepopath/codeexport/fourbar1turnOL/matlabfcn/*.* ../matlabfcn_fourbar1turnOL
cp -u $maplerepopath/codeexport/fourbar1turnIC/matlabfcn/*.* ../matlabfcn_fourbar1turnIC

