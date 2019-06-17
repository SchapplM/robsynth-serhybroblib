#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/mg10hlTE/matlabfcn ../matlabfcn_mg10hlTE
cp -a $maplerepopath/codeexport/mg10hlOL/matlabfcn/ ../matlabfcn_mg10hlOL
cp -a $maplerepopath/codeexport/mg10hlIC/matlabfcn/ ../matlabfcn_mg10hlIC
cp -a $maplerepopath/codeexport/mg10hlDE1/matlabfcn/ ../matlabfcn_mg10hlDE1 
cp -a $maplerepopath/codeexport/mg10hlDE2/matlabfcn/ ../matlabfcn_mg10hlDE2

##rsync -rv $maplerepopath/codeexport/mg10hlTE/matlabfcn/ ../matlabfcn_mg10hlTE
##rsync -rv $maplerepopath/codeexport/mg10hlDE1/matlabfcn/ ../matlabfcn_mg10hlDE1
##rsync -rv $maplerepopath/codeexport/mg10hlDE2/matlabfcn/ ../matlabfcn_mg10hlDE2
##rsync -rv $maplerepopath/codeexport/mg10hlIC/matlabfcn/ ../matlabfcn_mg10hlIC

