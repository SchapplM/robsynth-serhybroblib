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
cp $this_path/robot_env_fourbar2_CL $defpath/robot_env_fourbar2TE
sed -i "s/fourbar2/fourbar2TE/g" $defpath/robot_env_fourbar2TE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fourbar2TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar2TE

cp $this_path/robot_env_fourbar2_CL $defpath/robot_env_fourbar2DE1
sed -i "s/fourbar2/fourbar2DE1/g" $defpath/robot_env_fourbar2DE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_fourbar2TE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbar2DE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_fourbar2DE1

cp $this_path/robot_env_fourbar2_CL $defpath/robot_env_fourbar2DE2
sed -i "s/fourbar2/fourbar2DE2/g" $defpath/robot_env_fourbar2DE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_fourbar2TE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fourbar2DE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_fourbar2DE2

cp $this_path/robot_env_fourbar2IC $defpath/robot_env_fourbar2IC
##sed -i "s/fourbar2/fourbar2IC/g" $defpath/robot_env_fourbar2IC
##echo "# Anwendung der Additionstheoreme für parallele Drehachsen" >> $defpath/robot_env_fourbar2IC
##echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fourbar2IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/fourbar2TE_kinematic_constraints.mpl $constrpath/fourbar2TE_kinematic_constraints.mpl
cp $this_path/fourbar2TE_kinematic_constraints.mw $constrpath/fourbar2TE_kinematic_constraints.mw
cp $this_path/fourbar2DE_kinematic_constraints.mpl $constrpath/fourbar2DE1_kinematic_constraints.mpl
cp $this_path/fourbar2DE_kinematic_constraints.mw $constrpath/fourbar2DE1_kinematic_constraints.mw
cp $this_path/fourbar2DE_kinematic_constraints.mpl $constrpath/fourbar2DE2_kinematic_constraints.mpl
cp $this_path/fourbar2DE_kinematic_constraints.mw $constrpath/fourbar2DE2_kinematic_constraints.mw

cp $this_path/fourbar2IC_kinematic_constraints_implicit.mpl $constrpath/fourbar2IC_kinematic_constraints_implicit.mpl
cp $this_path/fourbar2IC_kinematic_constraints_implicit.mw $constrpath/fourbar2IC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
#cp $this_path/fourbar2_kinematic_parameter_values.m $constrpath/fourbar2IC_kinematic_parameter_values.m
cp $this_path/fourbar2_kinematic_parameter_values.m $constrpath/fourbar2TE_kinematic_parameter_values.m
cp $this_path/fourbar2_kinematic_parameter_values.m $constrpath/fourbar2DE1_kinematic_parameter_values.m
cp $this_path/fourbar2_kinematic_parameter_values.m $constrpath/fourbar2DE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
#cp $this_path/fourbar2_kinematic_constraints_matlab.m $constrpath/fourbar2IC_kinematic_constraints_matlab.m
cp $this_path/fourbar2_kinematic_constraints_matlab.m $constrpath/fourbar2TE_kinematic_constraints_matlab.m
cp $this_path/fourbar2_kinematic_constraints_matlab.m $constrpath/fourbar2DE1_kinematic_constraints_matlab.m
cp $this_path/fourbar2_kinematic_constraints_matlab.m $constrpath/fourbar2DE2_kinematic_constraints_matlab.m

