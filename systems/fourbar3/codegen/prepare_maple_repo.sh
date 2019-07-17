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
cp $this_path/robot_env_fourbar3_CL $defpath/robot_env_fourbar3TE
sed -i "s/fourbar3/fourbar3TE/g" $defpath/robot_env_fourbar3TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fourbar3TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar3TE

cp $this_path/robot_env_fourbar3_CL $defpath/robot_env_fourbar3DE1
sed -i "s/fourbar3/fourbar3DE1/g" $defpath/robot_env_fourbar3DE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_fourbar3TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar3DE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_fourbar3DE1

cp $this_path/robot_env_fourbar3_CL $defpath/robot_env_fourbar3DE2
sed -i "s/fourbar3/fourbar3DE2/g" $defpath/robot_env_fourbar3DE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_fourbar3TE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fourbar3DE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_fourbar3DE2

cp $this_path/robot_env_fourbar3OL $defpath/robot_env_fourbar3OL

cp $this_path/robot_env_fourbar3IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/fourbar3TE_kinematic_constraints.mpl $constrpath/fourbar3TE_kinematic_constraints.mpl
cp $this_path/fourbar3TE_kinematic_constraints.mw $constrpath/fourbar3TE_kinematic_constraints.mw
cp $this_path/fourbar3DE_kinematic_constraints.mpl $constrpath/fourbar3DE1_kinematic_constraints.mpl
cp $this_path/fourbar3DE_kinematic_constraints.mw $constrpath/fourbar3DE1_kinematic_constraints.mw
cp $this_path/fourbar3DE_kinematic_constraints.mpl $constrpath/fourbar3DE2_kinematic_constraints.mpl
cp $this_path/fourbar3DE_kinematic_constraints.mw $constrpath/fourbar3DE2_kinematic_constraints.mw

cp $this_path/fourbar3IC_kinematic_constraints_implicit.mpl $constrpath/fourbar3IC_kinematic_constraints_implicit.mpl
cp $this_path/fourbar3IC_kinematic_constraints_implicit.mw $constrpath/fourbar3IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
#cp $this_path/fourbar3_kinematic_parameter_values.m $constrpath/fourbar3IC_kinematic_parameter_values.m
cp $this_path/fourbar3_kinematic_parameter_values.m $constrpath/fourbar3TE_kinematic_parameter_values.m
cp $this_path/fourbar3_kinematic_parameter_values.m $constrpath/fourbar3DE1_kinematic_parameter_values.m
cp $this_path/fourbar3_kinematic_parameter_values.m $constrpath/fourbar3DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
#cp $this_path/fourbar3_kinematic_constraints_matlab.m $constrpath/fourbar3IC_kinematic_constraints_matlab.m
cp $this_path/fourbar3_kinematic_constraints_matlab.m $constrpath/fourbar3TE_kinematic_constraints_matlab.m
cp $this_path/fourbar3_kinematic_constraints_matlab.m $constrpath/fourbar3DE1_kinematic_constraints_matlab.m
cp $this_path/fourbar3_kinematic_constraints_matlab.m $constrpath/fourbar3DE2_kinematic_constraints_matlab.m

