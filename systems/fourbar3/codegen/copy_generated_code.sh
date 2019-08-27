#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fourbar3TE ../matlabfcn_fourbar3DE1 ../matlabfcn_fourbar3DE2 #../matlabfcn_fourbar3OL ../matlabfcn_fourbar3IC
cp $maplerepopath/codeexport/fourbar3TE/matlabfcn/* ../matlabfcn_fourbar3TE
cp $maplerepopath/codeexport/fourbar3DE1/matlabfcn/* ../matlabfcn_fourbar3DE1
cp $maplerepopath/codeexport/fourbar3DE2/matlabfcn/* ../matlabfcn_fourbar3DE2
#cp $maplerepopath/codeexport/fourbar3OL/matlabfcn/* ../matlabfcn_fourbar3OL
#cp $maplerepopath/codeexport/fourbar3IC/matlabfcn/* ../matlabfcn_fourbar3IC

