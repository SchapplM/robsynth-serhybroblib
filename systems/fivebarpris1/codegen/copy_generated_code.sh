#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fivebarpris1TE ../matlabfcn_fivebarpris1DE1 ../matlabfcn_fivebarpris1DE2 ../matlabfcn_fivebarpris1OL ../matlabfcn_fivebarpris1IC
cp -u $maplerepopath/codeexport/fivebarpris1TE/matlabfcn/*.* ../matlabfcn_fivebarpris1TE
cp -u $maplerepopath/codeexport/fivebarpris1DE1/matlabfcn/*.* ../matlabfcn_fivebarpris1DE1
cp -u $maplerepopath/codeexport/fivebarpris1DE2/matlabfcn/*.* ../matlabfcn_fivebarpris1DE2
cp -u $maplerepopath/codeexport/fivebarpris1OL/matlabfcn/*.* ../matlabfcn_fivebarpris1OL
cp -u $maplerepopath/codeexport/fivebarpris1IC/matlabfcn/*.* ../matlabfcn_fivebarpris1IC

