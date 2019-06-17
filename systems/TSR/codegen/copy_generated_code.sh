#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

 cp -a $maplerepopath/codeexport/TSRTE/matlabfcn ../matlabfcn_TSRTE
 cp -a $maplerepopath/codeexport/TSROL/matlabfcn/ ../matlabfcn_TSROL
 cp -a $maplerepopath/codeexport/TSRIC/matlabfcn/ ../matlabfcn_TSRIC
 cp -a $maplerepopath/codeexport/TSRDE1/matlabfcn/ ../matlabfcn_TSRDE1
 cp -a $maplerepopath/codeexport/TSRDE2/matlabfcn/ ../matlabfcn_TSRDE2
 
 
## rsync -rv $maplerepopath/codeexport/TSRTE/matlabfcn/ ../matlabfcn_TSRTE
## rsync -rv $maplerepopath/codeexport/TSRDE1/matlabfcn/ ../matlabfcn_TSRDE1
## rsync -rv $maplerepopath/codeexport/TSRDE2/matlabfcn/ ../matlabfcn_TSRDE2
## rsync -rv $maplerepopath/codeexport/TSRIC/matlabfcn/ ../matlabfcn_TSRIC

