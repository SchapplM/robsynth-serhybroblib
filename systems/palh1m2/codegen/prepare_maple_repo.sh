#!/bin/bash -ex
# Kopiere alle benötigten Dateien für die Code-Generierung ins HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
#
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
# (C) Institut für Mechatronische Systeme, Universität Hannover

maplerepopath=$1
this_path=$(pwd)

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;

defpath=$maplerepopath/robot_codegen_definitions
constrpath=$maplerepopath/robot_codegen_constraints

## Definitionen für Unterschiedliche Modellierungen kopieren
cp $this_path/robot_env_palh1m2_CL $defpath/robot_env_palh1m2TE
sed -i "s/palh1m2/palh1m2TE/g" $defpath/robot_env_palh1m2TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_palh1m2TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_palh1m2TE
echo "# Bei TE sind die Terme für den simplify-Befehl geeignet." >> $defpath/robot_env_palh1m2TE
echo "simplify_options := Vector(10,3):" >> $defpath/robot_env_palh1m2TE
echo "simplify_options(8):=1: # Lagrange kaum vereinfachen" >> $defpath/robot_env_palh1m2TE

cp $this_path/robot_env_palh1m2_CL $defpath/robot_env_palh1m2DE1
sed -i "s/palh1m2/palh1m2DE1/g" $defpath/robot_env_palh1m2DE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_palh1m2DE1
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_palh1m2DE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_palh1m2DE1

cp $this_path/robot_env_palh1m2_CL $defpath/robot_env_palh1m2DE2
sed -i "s/palh1m2/palh1m2DE2/g" $defpath/robot_env_palh1m2DE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_palh1m2DE2
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_palh1m2DE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_palh1m2DE2
echo "simplify_options := Vector(10,-1): # Term-Vereinfachungen mit Standard-Werten" >> $defpath/robot_env_palh1m2DE2
echo "simplify_options(9) := 0: # Hängt sich bei Massenmatrix auf" >> $defpath/robot_env_palh1m2DE2

cp $this_path/robot_env_palh1m2OL $defpath/robot_env_palh1m2OL

cp $this_path/robot_env_palh1m2IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/palh1m2TE_kinematic_constraints.mpl $constrpath/palh1m2TE_kinematic_constraints.mpl
cp $this_path/palh1m2TE_kinematic_constraints.mw $constrpath/palh1m2TE_kinematic_constraints.mw
cp $this_path/palh1m2DE_kinematic_constraints.mpl $constrpath/palh1m2DE1_kinematic_constraints.mpl
cp $this_path/palh1m2DE_kinematic_constraints.mw $constrpath/palh1m2DE1_kinematic_constraints.mw
cp $this_path/palh1m2DE_kinematic_constraints.mpl $constrpath/palh1m2DE2_kinematic_constraints.mpl
cp $this_path/palh1m2DE_kinematic_constraints.mw $constrpath/palh1m2DE2_kinematic_constraints.mw
cp $this_path/palh1m2IC_kinematic_constraints_implicit.mpl $constrpath/palh1m2IC_kinematic_constraints_implicit.mpl
cp $this_path/palh1m2IC_kinematic_constraints_implicit.mw $constrpath/palh1m2IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/palh1m2_kinematic_parameter_values.m $constrpath/palh1m2OL_kinematic_parameter_values.m
cp $this_path/palh1m2_kinematic_parameter_values.m $constrpath/palh1m2IC_kinematic_parameter_values.m
cp $this_path/palh1m2_kinematic_parameter_values.m $constrpath/palh1m2TE_kinematic_parameter_values.m
cp $this_path/palh1m2_kinematic_parameter_values.m $constrpath/palh1m2DE1_kinematic_parameter_values.m
cp $this_path/palh1m2_kinematic_parameter_values.m $constrpath/palh1m2DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests (nur für Modelle mit Elimination,
# Nicht für OL- oder IC-Modell. Dort reichen Zufallswerte)
cp $this_path/palh1m2_kinematic_constraints_matlab.m $constrpath/palh1m2TE_kinematic_constraints_matlab.m
cp $this_path/palh1m2_kinematic_constraints_matlab.m $constrpath/palh1m2DE1_kinematic_constraints_matlab.m
cp $this_path/palh1m2_kinematic_constraints_matlab.m $constrpath/palh1m2DE2_kinematic_constraints_matlab.m

