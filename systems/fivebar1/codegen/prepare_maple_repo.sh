#!/bin/bash -ex
# Kopiere alle benötigten Dateien für die Code-Generierung ins HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
#
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
# (C) Institut für mechatronische Systeme, Universität Hannover

maplerepopath=$1
this_path=$(pwd)

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;
defpath=$maplerepopath/robot_codegen_definitions
constrpath=$maplerepopath/robot_codegen_constraints

## Definitionen für Unterschiedliche Modellierungen kopieren
cp $this_path/robot_env_fivebar1_CL $defpath/robot_env_fivebar1TE
sed -i "s/fivebar1/fivebar1TE/g" $defpath/robot_env_fivebar1TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fivebar1TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fivebar1TE

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/fivebar1TE_kinematic_constraints.mpl $constrpath/fivebar1TE_kinematic_constraints.mpl
cp $this_path/fivebar1TE_kinematic_constraints.mw $constrpath/fivebar1TE_kinematic_constraints.mw

# Parameter für Testszenario
cp $this_path/fivebar1_kinematic_constraints_matlab.m $constrpath/fivebar1TE_kinematic_constraints_matlab.m
cp $this_path/fivebar1_kinematic_parameter_values.m $constrpath/fivebar1TE_kinematic_parameter_values.m

