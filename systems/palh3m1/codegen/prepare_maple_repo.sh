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
cp $this_path/robot_env_palh3m1_CL $defpath/robot_env_palh3m1TE
sed -i "s/palh3m1/palh3m1TE/g" $defpath/robot_env_palh3m1TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_palh3m1TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_palh3m1TE

cp $this_path/robot_env_palh3m1_CL $defpath/robot_env_palh3m1DE1
sed -i "s/palh3m1/palh3m1DE1/g" $defpath/robot_env_palh3m1DE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_palh3m1DE1
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_palh3m1DE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_palh3m1DE1

cp $this_path/robot_env_palh3m1_CL $defpath/robot_env_palh3m1DE2
sed -i "s/palh3m1/palh3m1DE2/g" $defpath/robot_env_palh3m1DE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_palh3m1DE2
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_palh3m1DE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_palh3m1DE2

cp $this_path/robot_env_palh3m1OL $defpath/robot_env_palh3m1OL

cp $this_path/robot_env_palh3m1IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/palh3m1TE_kinematic_constraints.mpl $constrpath/palh3m1TE_kinematic_constraints.mpl
cp $this_path/palh3m1TE_kinematic_constraints.mw $constrpath/palh3m1TE_kinematic_constraints.mw
cp $this_path/palh3m1DE_kinematic_constraints.mpl $constrpath/palh3m1DE1_kinematic_constraints.mpl
cp $this_path/palh3m1DE_kinematic_constraints.mw $constrpath/palh3m1DE1_kinematic_constraints.mw
cp $this_path/palh3m1DE_kinematic_constraints.mpl $constrpath/palh3m1DE2_kinematic_constraints.mpl
cp $this_path/palh3m1DE_kinematic_constraints.mw $constrpath/palh3m1DE2_kinematic_constraints.mw
cp $this_path/palh3m1IC_kinematic_constraints_implicit.mpl $constrpath/palh3m1IC_kinematic_constraints_implicit.mpl
cp $this_path/palh3m1IC_kinematic_constraints_implicit.mw $constrpath/palh3m1IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/palh3m1_kinematic_parameter_values.m $constrpath/palh3m1OL_kinematic_parameter_values.m
cp $this_path/palh3m1_kinematic_parameter_values.m $constrpath/palh3m1IC_kinematic_parameter_values.m
cp $this_path/palh3m1_kinematic_parameter_values.m $constrpath/palh3m1TE_kinematic_parameter_values.m
cp $this_path/palh3m1_kinematic_parameter_values.m $constrpath/palh3m1DE1_kinematic_parameter_values.m
cp $this_path/palh3m1_kinematic_parameter_values.m $constrpath/palh3m1DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests (nur für Modelle mit Elimination,
# Nicht für OL- oder IC-Modell. Dort reichen Zufallswerte)
cp $this_path/palh3m1_kinematic_constraints_matlab.m $constrpath/palh3m1TE_kinematic_constraints_matlab.m
cp $this_path/palh3m1_kinematic_constraints_matlab.m $constrpath/palh3m1DE1_kinematic_constraints_matlab.m
cp $this_path/palh3m1_kinematic_constraints_matlab.m $constrpath/palh3m1DE2_kinematic_constraints_matlab.m

