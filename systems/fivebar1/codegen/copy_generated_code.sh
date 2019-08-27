#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fivebar1TE #../matlabfcn_fivebar1DE1 ../matlabfcn_fivebar1DE2 ../matlabfcn_fivebar1OL ../matlabfcn_fivebar1IC
cp $maplerepopath/codeexport/fivebar1TE/matlabfcn/* ../matlabfcn_fivebar1TE
#cp $maplerepopath/codeexport/fivebar1DE1/matlabfcn/* ../matlabfcn_fivebar1DE1
#cp $maplerepopath/codeexport/fivebar1DE2/matlabfcn/* ../matlabfcn_fivebar1DE2
#cp $maplerepopath/codeexport/fivebar1OL/matlabfcn/* ../matlabfcn_fivebar1OL
#cp $maplerepopath/codeexport/fivebar1IC/matlabfcn/* ../matlabfcn_fivebar1IC

