
# Berechne kinematische Zwangsbedingungen für Viergelenkkette
# Einleitung

# Die impliziten kinematischen Zwangsbedingungen werden aufgestellt: Annahme : Alle MDH-Winkel seien unabhängig und die Kette sei offen
# 
# fourbar1_OL -> Viergelenkkette, modellierung impliziten der Zwangsbedingungen mit Ausdrücken einer offenen Kette  und basierend auf der MDH-Tabelle. 
# Bei Matlab: Abhängige Winkel durch die Ausdrücke von TE aus der Datei "fourbar1TE_kinconstr_expl_mdh_sym_varpar" ersetzen.
# Quelle
# SA Bejaoui: Bejaoui2018_S749; "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen"
# Autor
# Abderahman Bejaoui
# Studienarbeit bei: Moritz Schappler, schappler@irt.uni-hannover.de, 2018-08
# (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover
# 
# Initialisierung
restart:
kin_constraints_exist := true: # Für Speicherung
;
with(StringTools): # Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 1: # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
;
read "../helper/proc_MatlabExport":
read "../transformation/proc_rotx":
read "../transformation/proc_roty":
read "../transformation/proc_rotz":
read "../transformation/proc_trotz":
read "../transformation/proc_transl":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
read "../robot_codegen_constraints/proc_subs_kintmp_exp":
read "../helper/proc_intersect_circle":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Ergebnisse der Kinematik laden
read sprintf("../codeexport/%s/tmp/kinematics_floatb_%s_rotmat_maple.m", robot_name, base_method_name);
Trf := Trf:
Trf_c := Trf_c:
Trf
;
# Fourbar1 0-1-2-4-5-3-0
# Schleife 0-1-2-4
T_0_4 := combine(Matrix(Trf(1..4,1..4, 1)) . Matrix(Trf(1..4,1..4,2)). Matrix(Trf(1..4,1..4, 4)));
# Schleife 0-3-5
T_0_5:= combine(Matrix(Trf(1..4,1..4, 3)) . Matrix(Trf(1..4,1..4,5)));
h1t := T_0_5(1..3,4) - T_0_4(1..3,4);
tmp := Transpose( Matrix(T_0_5(1..3,1..3)) ) . Matrix(T_0_4(1..3,1..3));  # nur anzeigen lassen für h1r
;
 combine(tmp); # nur anzeigen lassen für h1r
;
h1r:=Pi-(-qJ3(t)+qJ1(t)+qJ2(t)+qJ4(t));
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := <h1t([1, 2]);h1r>; 
implconstr_s := convert_t_s(implconstr_t):
# Exportiere Code für folgende Skripte
kin_constraints_exist:=true:
save implconstr_t, implconstr_s, kin_constraints_exist, sprintf("../codeexport/%s/tmp/kinematic_constraints_implicit_maple.m", robot_name):
# Exportieren des vollständigen Ausdruckes
if codegen_act then
  MatlabExport(implconstr_s, sprintf("../codeexport/%s/tmp/kinconstr_impl_matlab.m", robot_name), 2):
end if:
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( implconstr_s )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_matlab.m", robot_name), 2):


