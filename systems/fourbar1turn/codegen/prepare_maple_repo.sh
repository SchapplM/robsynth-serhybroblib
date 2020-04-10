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
cp $this_path/robot_env_fourbar1turn_CL $defpath/robot_env_fourbar1turnTE
sed -i "s/fourbar1turn/fourbar1turnTE/g" $defpath/robot_env_fourbar1turnTE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fourbar1turnTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar1turnTE

cp $this_path/robot_env_fourbar1turn_CL $defpath/robot_env_fourbar1turnDE1
sed -i "s/fourbar1turn/fourbar1turnDE1/g" $defpath/robot_env_fourbar1turnDE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_fourbar1turnTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar1turnDE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_fourbar1turnDE1

cp $this_path/robot_env_fourbar1turn_CL $defpath/robot_env_fourbar1turnDE2
sed -i "s/fourbar1turn/fourbar1turnDE2/g" $defpath/robot_env_fourbar1turnDE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_fourbar1turnTE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fourbar1turnDE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_fourbar1turnDE2

cp $this_path/robot_env_fourbar1turnOL $defpath/robot_env_fourbar1turnOL

cp $this_path/robot_env_fourbar1turnIC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/fourbar1turnTE_kinematic_constraints.mpl $constrpath/fourbar1turnTE_kinematic_constraints.mpl
cp $this_path/fourbar1turnTE_kinematic_constraints.mw $constrpath/fourbar1turnTE_kinematic_constraints.mw
cp $this_path/fourbar1turnDE_kinematic_constraints.mpl $constrpath/fourbar1turnDE1_kinematic_constraints.mpl
cp $this_path/fourbar1turnDE_kinematic_constraints.mw $constrpath/fourbar1turnDE1_kinematic_constraints.mw
cp $this_path/fourbar1turnDE_kinematic_constraints.mpl $constrpath/fourbar1turnDE2_kinematic_constraints.mpl
cp $this_path/fourbar1turnDE_kinematic_constraints.mw $constrpath/fourbar1turnDE2_kinematic_constraints.mw

cp $this_path/fourbar1turnIC_kinematic_constraints_implicit.mpl $constrpath/fourbar1turnIC_kinematic_constraints_implicit.mpl
cp $this_path/fourbar1turnIC_kinematic_constraints_implicit.mw $constrpath/fourbar1turnIC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/fourbar1turn_kinematic_parameter_values.m $constrpath/fourbar1turnOL_kinematic_parameter_values.m
cp $this_path/fourbar1turn_kinematic_parameter_values.m $constrpath/fourbar1turnIC_kinematic_parameter_values.m
cp $this_path/fourbar1turn_kinematic_parameter_values.m $constrpath/fourbar1turnTE_kinematic_parameter_values.m
cp $this_path/fourbar1turn_kinematic_parameter_values.m $constrpath/fourbar1turnDE1_kinematic_parameter_values.m
cp $this_path/fourbar1turn_kinematic_parameter_values.m $constrpath/fourbar1turnDE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
cp $this_path/fourbar1turn_kinematic_constraints_matlab.m $constrpath/fourbar1turnIC_kinematic_constraints_matlab.m
cp $this_path/fourbar1turn_kinematic_constraints_matlab.m $constrpath/fourbar1turnTE_kinematic_constraints_matlab.m
cp $this_path/fourbar1turn_kinematic_constraints_matlab.m $constrpath/fourbar1turnDE1_kinematic_constraints_matlab.m
cp $this_path/fourbar1turn_kinematic_constraints_matlab.m $constrpath/fourbar1turnDE2_kinematic_constraints_matlab.m

