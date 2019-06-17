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
cp $this_path/robot_env_TSR_CL $defpath/robot_env_TSRTE
sed -i "s/TSR/TSRTE/g" $defpath/robot_env_TSRTE
echo "# Bei Trigonometrischer Elimination können die Additionstheoreme nicht angewendet werden" >> $defpath/robot_env_TSRTE
echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_TSRTE

 cp $this_path/robot_env_TSR_CL $defpath/robot_env_TSRDE1
 sed -i "s/TSR/TSRDE1/g" $defpath/robot_env_TSRDE1
 echo "# Ersetzung der Abhängigen Winkel direkt in den Einzel-Transformationen" >> $defpath/robot_env_TSRDE1
 echo "codegen_kinematics_opt := false:" >> $defpath/robot_env_TSRDE1
 echo "codegen_kinematics_subsorder:=1:" >> $defpath/robot_env_TSRDE1

 cp $this_path/robot_env_TSR_CL $defpath/robot_env_TSRDE2
 sed -i "s/TSR/TSRDE2/g" $defpath/robot_env_TSRDE2
 echo "# Anwendung der Additionstheoreme für parallele Drehachsen und dann Ersetzung der Abhängigen Winkel" >> $defpath/robot_env_TSRDE2
 echo "codegen_kinematics_opt := true:" >> $defpath/robot_env_TSRDE2
 echo "codegen_kinematics_subsorder:=2:" >> $defpath/robot_env_TSRDE2

 cp $this_path/robot_env_TSROL $defpath/robot_env_TSROL
 
 cp $this_path/robot_env_TSRIC $defpath/robot_env_IC

# Maple-Skripte (Kinematische Zwangsbedingungen)

 cp $this_path/TSRTE_kinematic_constraints.mpl $constrpath/TSRTE_kinematic_constraints.mpl
 cp $this_path/TSRTE_kinematic_constraints.mw $constrpath/TSRTE_kinematic_constraints.mw
 cp $this_path/TSRDE_kinematic_constraints.mpl $constrpath/TSRDE1_kinematic_constraints.mpl
 cp $this_path/TSRDE_kinematic_constraints.mw $constrpath/TSRDE1_kinematic_constraints.mw
 cp $this_path/TSRDE_kinematic_constraints.mpl $constrpath/TSRDE2_kinematic_constraints.mpl
 cp $this_path/TSRDE_kinematic_constraints.mw $constrpath/TSRDE2_kinematic_constraints.mw
 cp $this_path/TSRIC_kinematic_constraints_implicit.mpl $constrpath/TSRIC_kinematic_constraints_implicit.mpl
 cp $this_path/TSRIC_kinematic_constraints_implicit.mw $constrpath/TSRIC_kinematic_constraints_implicit.mw

# Werte für Kinematikparameter (für Modultests)

 cp $this_path/TSR_kinematic_parameter_values.m $constrpath/TSRIC_kinematic_parameter_values.m
 cp $this_path/TSR_kinematic_parameter_values.m $constrpath/TSRTE_kinematic_parameter_values.m
 cp $this_path/TSR_kinematic_parameter_values.m $constrpath/TSRDE1_kinematic_parameter_values.m
 cp $this_path/TSR_kinematic_parameter_values.m $constrpath/TSRDE2_kinematic_parameter_values.m

# Winkelgrenzen für die Modultests

 cp $this_path/TSR_kinematic_constraints_matlab.m $constrpath/TSRIC_kinematic_constraints_matlab.m
 cp $this_path/TSR_kinematic_constraints_matlab.m $constrpath/TSRTE_kinematic_constraints_matlab.m
 cp $this_path/TSR_kinematic_constraints_matlab.m $constrpath/TSRDE1_kinematic_constraints_matlab.m
 cp $this_path/TSR_kinematic_constraints_matlab.m $constrpath/TSRDE2_kinematic_constraints_matlab.m
