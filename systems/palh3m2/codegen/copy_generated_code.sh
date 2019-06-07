#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

 cp -a $maplerepopath/codeexport/palh3m2TE/matlabfcn ../matlabfcn_palh3m2TE
 cp -a $maplerepopath/codeexport/palh3m2IC/matlabfcn/ ../matlabfcn_palh3m2IC
 cp -a $maplerepopath/codeexport/palh3m2DE1/matlabfcn/ ../matlabfcn_palh3m2DE1
 cp -a $maplerepopath/codeexport/palh3m2DE2/matlabfcn/ ../matlabfcn_palh3m2DE2
 
 
## rsync -rv $maplerepopath/codeexport/palh3m2TE/matlabfcn/ ../matlabfcn_palh3m2TE
## rsync -rv $maplerepopath/codeexport/palh3m2DE1/matlabfcn/ ../matlabfcn_palh3m2DE1
## rsync -rv $maplerepopath/codeexport/palh3m2DE2/matlabfcn/ ../matlabfcn_palh3m2DE2
## rsync -rv $maplerepopath/codeexport/palh3m2IC/matlabfcn/ ../matlabfcn_palh3m2IC

