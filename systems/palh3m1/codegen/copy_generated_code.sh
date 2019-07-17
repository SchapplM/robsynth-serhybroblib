#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/palh3m1TE/matlabfcn ../matlabfcn_palh3m1TE
cp -a $maplerepopath/codeexport/palh3m1OL/matlabfcn/ ../matlabfcn_palh3m1OL
cp -a $maplerepopath/codeexport/palh3m1IC/matlabfcn/ ../matlabfcn_palh3m1IC
cp -a $maplerepopath/codeexport/palh3m1DE1/matlabfcn/ ../matlabfcn_palh3m1DE1
cp -a $maplerepopath/codeexport/palh3m1DE2/matlabfcn/ ../matlabfcn_palh3m1DE2

## rsync -rv $maplerepopath/codeexport/palh3m1TE/matlabfcn/ ../matlabfcn_palh3m1TE
## rsync -rv $maplerepopath/codeexport/palh3m1DE1/matlabfcn/ ../matlabfcn_palh3m1DE1
## rsync -rv $maplerepopath/codeexport/palh3m1DE2/matlabfcn/ ../matlabfcn_palh3m1DE2
## rsync -rv $maplerepopath/codeexport/palh3m1IC/matlabfcn/ ../matlabfcn_palh3m1IC

