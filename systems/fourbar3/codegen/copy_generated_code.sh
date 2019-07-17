#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

cp -a $maplerepopath/codeexport/fourbar3TE/matlabfcn/ ../matlabfcn_fourbar3TE
cp -a $maplerepopath/codeexport/fourbar3DE1/matlabfcn/ ../matlabfcn_fourbar3DE1
cp -a $maplerepopath/codeexport/fourbar3DE2/matlabfcn/ ../matlabfcn_fourbar3DE2
cp -a $maplerepopath/codeexport/fourbar3DE2/matlabfcn/ ../matlabfcn_fourbar3IC

# rsync -rv $maplerepopath/codeexport/fourbar3TE/matlabfcn/ ../matlabfcn_fourbar3TE
# rsync -rv $maplerepopath/codeexport/fourbar3DE1/matlabfcn/ ../matlabfcn_fourbar3DE1
# rsync -rv $maplerepopath/codeexport/fourbar3DE2/matlabfcn/ ../matlabfcn_fourbar3DE2
# rsync -rv $maplerepopath/codeexport/fourbar3DE2/matlabfcn/ ../matlabfcn_fourbar3IC

