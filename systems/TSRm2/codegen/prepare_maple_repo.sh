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
cp $this_path/robot_env_TSRm2_CL $defpath/robot_env_TSRm2TE
sed -i "s/TSRm2/TSRm2TE/g" $defpath/robot_env_TSRm2TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_TSRm2TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_TSRm2TE

 cp $this_path/robot_env_TSRm2_CL $defpath/robot_env_TSRm2DE1
 sed -i "s/TSRm2/TSRm2DE1/g" $defpath/robot_env_TSRm2DE1
 echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_TSRm2DE1
 echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_TSRm2DE1
 echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_TSRm2DE1

 cp $this_path/robot_env_TSRm2_CL $defpath/robot_env_TSRm2DE2
 sed -i "s/TSRm2/TSRm2DE2/g" $defpath/robot_env_TSRm2DE2
 echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_TSRm2DE2
 echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_TSRm2DE2
 echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_TSRm2DE2

 cp $this_path/robot_env_TSRm2OL $defpath/robot_env_TSRm2OL
 
 cp $this_path/robot_env_TSRm2IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

 cp $this_path/TSRm2TE_kinematic_constraints.mpl $constrpath/TSRm2TE_kinematic_constraints.mpl
 cp $this_path/TSRm2TE_kinematic_constraints.mw $constrpath/TSRm2TE_kinematic_constraints.mw
 cp $this_path/TSRm2DE_kinematic_constraints.mpl $constrpath/TSRm2DE1_kinematic_constraints.mpl
 cp $this_path/TSRm2DE_kinematic_constraints.mw $constrpath/TSRm2DE1_kinematic_constraints.mw
 cp $this_path/TSRm2DE_kinematic_constraints.mpl $constrpath/TSRm2DE2_kinematic_constraints.mpl
 cp $this_path/TSRm2DE_kinematic_constraints.mw $constrpath/TSRm2DE2_kinematic_constraints.mw
 cp $this_path/TSRm2IC_kinematic_constraints_implicit.mpl $constrpath/TSRm2IC_kinematic_constraints_implicit.mpl
 cp $this_path/TSRm2IC_kinematic_constraints_implicit.mw $constrpath/TSRm2IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)

 cp $this_path/TSRm2_kinematic_parameter_values.m $constrpath/TSRm2OL_kinematic_parameter_values.m
 cp $this_path/TSRm2_kinematic_parameter_values.m $constrpath/TSRm2IC_kinematic_parameter_values.m
 cp $this_path/TSRm2_kinematic_parameter_values.m $constrpath/TSRm2TE_kinematic_parameter_values.m
 cp $this_path/TSRm2_kinematic_parameter_values.m $constrpath/TSRm2DE1_kinematic_parameter_values.m
 cp $this_path/TSRm2_kinematic_parameter_values.m $constrpath/TSRm2DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests

 cp $this_path/TSRm2_kinematic_constraints_matlab.m $constrpath/TSRm2IC_kinematic_constraints_matlab.m
 cp $this_path/TSRm2_kinematic_constraints_matlab.m $constrpath/TSRm2TE_kinematic_constraints_matlab.m
 cp $this_path/TSRm2_kinematic_constraints_matlab.m $constrpath/TSRm2DE1_kinematic_constraints_matlab.m
 cp $this_path/TSRm2_kinematic_constraints_matlab.m $constrpath/TSRm2DE2_kinematic_constraints_matlab.m
