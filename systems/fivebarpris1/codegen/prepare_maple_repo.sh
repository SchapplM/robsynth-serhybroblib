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
cp $this_path/robot_env_fivebarpris1_CL $defpath/robot_env_fivebarpris1TE
sed -i "s/fivebarpris1/fivebarpris1TE/g" $defpath/robot_env_fivebarpris1TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fivebarpris1TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fivebarpris1TE

cp $this_path/robot_env_fivebarpris1_CL $defpath/robot_env_fivebarpris1DE1
sed -i "s/fivebarpris1/fivebarpris1DE1/g" $defpath/robot_env_fivebarpris1DE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_fivebarpris1TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fivebarpris1DE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_fivebarpris1DE1

cp $this_path/robot_env_fivebarpris1_CL $defpath/robot_env_fivebarpris1DE2
sed -i "s/fivebarpris1/fivebarpris1DE2/g" $defpath/robot_env_fivebarpris1DE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_fivebarpris1TE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fivebarpris1DE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_fivebarpris1DE2

cp $this_path/robot_env_fivebarpris1OL $defpath/robot_env_fivebarpris1OL

cp $this_path/robot_env_fivebarpris1IC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/fivebarpris1TE_kinematic_constraints.mpl $constrpath/fivebarpris1TE_kinematic_constraints.mpl
cp $this_path/fivebarpris1TE_kinematic_constraints.mw $constrpath/fivebarpris1TE_kinematic_constraints.mw
cp $this_path/fivebarpris1DE_kinematic_constraints.mpl $constrpath/fivebarpris1DE1_kinematic_constraints.mpl
cp $this_path/fivebarpris1DE_kinematic_constraints.mw $constrpath/fivebarpris1DE1_kinematic_constraints.mw
cp $this_path/fivebarpris1DE_kinematic_constraints.mpl $constrpath/fivebarpris1DE2_kinematic_constraints.mpl
cp $this_path/fivebarpris1DE_kinematic_constraints.mw $constrpath/fivebarpris1DE2_kinematic_constraints.mw

cp $this_path/fivebarpris1IC_kinematic_constraints_implicit.mpl $constrpath/fivebarpris1IC_kinematic_constraints_implicit.mpl
cp $this_path/fivebarpris1IC_kinematic_constraints_implicit.mw $constrpath/fivebarpris1IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/fivebarpris1_kinematic_parameter_values.m $constrpath/fivebarpris1IC_kinematic_parameter_values.m
cp $this_path/fivebarpris1_kinematic_parameter_values.m $constrpath/fivebarpris1TE_kinematic_parameter_values.m
cp $this_path/fivebarpris1_kinematic_parameter_values.m $constrpath/fivebarpris1DE1_kinematic_parameter_values.m
cp $this_path/fivebarpris1_kinematic_parameter_values.m $constrpath/fivebarpris1DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
cp $this_path/fivebarpris1_kinematic_constraints_matlab.m $constrpath/fivebarpris1IC_kinematic_constraints_matlab.m
cp $this_path/fivebarpris1_kinematic_constraints_matlab.m $constrpath/fivebarpris1TE_kinematic_constraints_matlab.m
cp $this_path/fivebarpris1_kinematic_constraints_matlab.m $constrpath/fivebarpris1DE1_kinematic_constraints_matlab.m
cp $this_path/fivebarpris1_kinematic_constraints_matlab.m $constrpath/fivebarpris1DE2_kinematic_constraints_matlab.m

