#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/palh2m1DE/matlabfcn/ ../matlabfcn_palh2m1DE
cp -a $maplerepopath/codeexport/palh2m1IC/matlabfcn/ ../matlabfcn_palh2m1IC
cp -a $maplerepopath/codeexport/palh2m1OL/matlabfcn/ ../matlabfcn_palh2m1OL