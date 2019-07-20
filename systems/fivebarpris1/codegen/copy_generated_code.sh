#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mkdir -p ../matlabfcn_fiveparpris1TE ../matlabfcn_fiveparpris1DE1 ../matlabfcn_fiveparpris1DE2 ../matlabfcn_fiveparpris1OL ../matlabfcn_fiveparpris1IC
cp $maplerepopath/codeexport/fiveparpris1TE/matlabfcn/* ../matlabfcn_fiveparpris1TE
cp $maplerepopath/codeexport/fiveparpris1DE1/matlabfcn/* ../matlabfcn_fiveparpris1DE1
cp $maplerepopath/codeexport/fiveparpris1DE2/matlabfcn/* ../matlabfcn_fiveparpris1DE2
cp $maplerepopath/codeexport/fiveparpris1OL/matlabfcn/* ../matlabfcn_fiveparpris1OL
cp $maplerepopath/codeexport/fiveparpris1IC/matlabfcn/* ../matlabfcn_fiveparpris1IC

