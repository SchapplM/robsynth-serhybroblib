#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/fivebarpris1TE/matlabfcn ../matlabfcn_fivebarpris1TE
cp -a $maplerepopath/codeexport/fivebarpris1OL/matlabfcn/ ../matlabfcn_fivebarpris1OL
cp -a $maplerepopath/codeexport/fivebarpris1IC/matlabfcn/ ../matlabfcn_fivebarpris1IC
cp -a $maplerepopath/codeexport/fivebarpris1DE1/matlabfcn/ ../matlabfcn_fivebarpris1DE1 
cp -a $maplerepopath/codeexport/fivebarpris1DE2/matlabfcn/ ../matlabfcn_fivebarpris1DE2

##rsync -rv $maplerepopath/codeexport/fivebarpris1TE/matlabfcn/ ../matlabfcn_fivebarpris1TE
##rsync -rv $maplerepopath/codeexport/fivebarpris1DE1/matlabfcn/ ../matlabfcn_fivebarpris1DE1
##rsync -rv $maplerepopath/codeexport/fivebarpris1DE2/matlabfcn/ ../matlabfcn_fivebarpris1DE2
##rsync -rv $maplerepopath/codeexport/fivebarpris1IC/matlabfcn/ ../matlabfcn_fivebarpris1IC

