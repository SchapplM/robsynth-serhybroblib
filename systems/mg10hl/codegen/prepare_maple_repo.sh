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
cp $this_path/robot_env_mg10hl_CL $defpath/robot_env_mg10hlTE
sed -i "s/mg10hl/mg10hlTE/g" $defpath/robot_env_mg10hlTE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_mg10hlTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_mg10hlTE

cp $this_path/robot_env_mg10hl_CL $defpath/robot_env_mg10hlDE1
sed -i "s/mg10hl/mg10hlDE1/g" $defpath/robot_env_mg10hlDE1
echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_mg10hlTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_mg10hlDE1
echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_mg10hlDE1

cp $this_path/robot_env_mg10hl_CL $defpath/robot_env_mg10hlDE2
sed -i "s/mg10hl/mg10hlDE2/g" $defpath/robot_env_mg10hlDE2
echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_mg10hlTE
echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_mg10hlDE2
echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_mg10hlDE2

cp $this_path/robot_env_mg10hlOL $defpath/robot_env_mg10hlOL

cp $this_path/robot_env_mg10hlIC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

cp $this_path/mg10hlTE_kinematic_constraints.mpl $constrpath/mg10hlTE_kinematic_constraints.mpl
cp $this_path/mg10hlTE_kinematic_constraints.mw $constrpath/mg10hlTE_kinematic_constraints.mw
cp $this_path/mg10hlDE_kinematic_constraints.mpl $constrpath/mg10hlDE1_kinematic_constraints.mpl
cp $this_path/mg10hlDE_kinematic_constraints.mw $constrpath/mg10hlDE1_kinematic_constraints.mw
cp $this_path/mg10hlDE_kinematic_constraints.mpl $constrpath/mg10hlDE2_kinematic_constraints.mpl
cp $this_path/mg10hlDE_kinematic_constraints.mw $constrpath/mg10hlDE2_kinematic_constraints.mw

cp $this_path/mg10hlIC_kinematic_constraints_implicit.mpl $constrpath/mg10hlIC_kinematic_constraints_implicit.mpl
cp $this_path/mg10hlIC_kinematic_constraints_implicit.mw $constrpath/mg10hlIC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/mg10hl_kinematic_parameter_values.m $constrpath/mg10hlIC_kinematic_parameter_values.m
cp $this_path/mg10hl_kinematic_parameter_values.m $constrpath/mg10hlTE_kinematic_parameter_values.m
cp $this_path/mg10hl_kinematic_parameter_values.m $constrpath/mg10hlDE1_kinematic_parameter_values.m
cp $this_path/mg10hl_kinematic_parameter_values.m $constrpath/mg10hlDE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests
cp $this_path/mg10hl_kinematic_constraints_matlab.m $constrpath/mg10hlIC_kinematic_constraints_matlab.m
cp $this_path/mg10hl_kinematic_constraints_matlab.m $constrpath/mg10hlTE_kinematic_constraints_matlab.m
cp $this_path/mg10hl_kinematic_constraints_matlab.m $constrpath/mg10hlDE1_kinematic_constraints_matlab.m
cp $this_path/mg10hl_kinematic_constraints_matlab.m $constrpath/mg10hlDE2_kinematic_constraints_matlab.m

