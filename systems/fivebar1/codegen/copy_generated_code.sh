#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, schappler@irt.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1

rsync -rv $maplerepopath/codeexport/fivebar1TE/matlabfcn/ ../matlabfcn_fivebar1TE
#rsync -rv $maplerepopath/codeexport/fivebar1DE1/matlabfcn/ ../matlabfcn_fivebar1DE1
#rsync -rv $maplerepopath/codeexport/fivebar1DE2/matlabfcn/ ../matlabfcn_fivebar1DE2
#rsync -rv $maplerepopath/codeexport/fivebar1DE2/matlabfcn/ ../matlabfcn_fivebar1IC

