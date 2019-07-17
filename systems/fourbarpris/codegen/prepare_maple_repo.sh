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
cp $this_path/robot_env_fourbarpris_CL $defpath/robot_env_fourbarprisTE
sed -i "s/fourbarpris/fourbarprisTE/g" $defpath/robot_env_fourbarprisTE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_fourbarprisTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbarprisTE

cp $this_path/robot_env_fourbarpris_CL $defpath/robot_env_fourbarprisDE1
sed -i "s/fourbarpris/fourbarprisDE1/g" $defpath/robot_env_fourbarprisDE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_fourbarprisTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_fourbarprisDE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_fourbarprisDE1

cp $this_path/robot_env_fourbarpris_CL $defpath/robot_env_fourbarprisDE2
sed -i "s/fourbarpris/fourbarprisDE2/g" $defpath/robot_env_fourbarprisDE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_fourbarprisTE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_fourbarprisDE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_fourbarprisDE2

cp $this_path/robot_env_fourbarprisOL $defpath/robot_env_fourbarprisOL

cp $this_path/robot_env_fourbarprisIC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/fourbarprisTE_kinematic_constraints.mpl $constrpath/fourbarprisTE_kinematic_constraints.mpl
cp $this_path/fourbarprisTE_kinematic_constraints.mw $constrpath/fourbarprisTE_kinematic_constraints.mw
cp $this_path/fourbarprisDE_kinematic_constraints.mpl $constrpath/fourbarprisDE1_kinematic_constraints.mpl
cp $this_path/fourbarprisDE_kinematic_constraints.mw $constrpath/fourbarprisDE1_kinematic_constraints.mw
cp $this_path/fourbarprisDE_kinematic_constraints.mpl $constrpath/fourbarprisDE2_kinematic_constraints.mpl
cp $this_path/fourbarprisDE_kinematic_constraints.mw $constrpath/fourbarprisDE2_kinematic_constraints.mw

cp $this_path/fourbarprisIC_kinematic_constraints_implicit.mpl $constrpath/fourbarprisIC_kinematic_constraints_implicit.mpl
cp $this_path/fourbarprisIC_kinematic_constraints_implicit.mw $constrpath/fourbarprisIC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/fourbarpris_kinematic_parameter_values.m $constrpath/fourbarprisIC_kinematic_parameter_values.m
cp $this_path/fourbarpris_kinematic_parameter_values.m $constrpath/fourbarprisTE_kinematic_parameter_values.m
cp $this_path/fourbarpris_kinematic_parameter_values.m $constrpath/fourbarprisDE1_kinematic_parameter_values.m
cp $this_path/fourbarpris_kinematic_parameter_values.m $constrpath/fourbarprisDE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
cp $this_path/fourbarpris_kinematic_constraints_matlab.m $constrpath/fourbarprisIC_kinematic_constraints_matlab.m
cp $this_path/fourbarpris_kinematic_constraints_matlab.m $constrpath/fourbarprisTE_kinematic_constraints_matlab.m
cp $this_path/fourbarpris_kinematic_constraints_matlab.m $constrpath/fourbarprisDE1_kinematic_constraints_matlab.m
cp $this_path/fourbarpris_kinematic_constraints_matlab.m $constrpath/fourbarprisDE2_kinematic_constraints_matlab.m

