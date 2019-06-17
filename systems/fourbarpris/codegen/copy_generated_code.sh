#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/fourbarprisTE/matlabfcn ../matlabfcn_fourbarprisTE
cp -a $maplerepopath/codeexport/fourbarprisOL/matlabfcn/ ../matlabfcn_fourbarprisOL
cp -a $maplerepopath/codeexport/fourbarprisIC/matlabfcn/ ../matlabfcn_fourbarprisIC
cp -a $maplerepopath/codeexport/fourbarprisDE1/matlabfcn/ ../matlabfcn_fourbarprisDE1 
cp -a $maplerepopath/codeexport/fourbarprisDE2/matlabfcn/ ../matlabfcn_fourbarprisDE2

##rsync -rv $maplerepopath/codeexport/fourbarprisTE/matlabfcn/ ../matlabfcn_fourbarprisTE
##rsync -rv $maplerepopath/codeexport/fourbarprisDE1/matlabfcn/ ../matlabfcn_fourbarprisDE1
##rsync -rv $maplerepopath/codeexport/fourbarprisDE2/matlabfcn/ ../matlabfcn_fourbarprisDE2
##rsync -rv $maplerepopath/codeexport/fourbarprisIC/matlabfcn/ ../matlabfcn_fourbarprisIC

