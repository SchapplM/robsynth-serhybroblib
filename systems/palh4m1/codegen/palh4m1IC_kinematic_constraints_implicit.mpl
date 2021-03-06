
# fgkpris implizit
# Initialisierung
restart:
kin_constraints_exist := true: # F�r Speicherung
;
with(StringTools): # F�r Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 1: # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: # Zur Code-Generierung auch f�r Nicht-Inert-Ausdr�cke
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
with(RealDomain): # Schr�nkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env_IC":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name_OL):
# Ergebnisse der Kinematik laden
read sprintf("../codeexport/%s/tmp/kinematics_floatb_%s_rotmat_maple.m", robot_name_OL, base_method_name);
read "../robot_codegen_definitions/robot_env_IC":
Trf := Trf:
Trf_c := Trf_c:
Trf
;
# fgkpris1-2-3-4-8-9-7-1
# Schleife 1-2-3-4-8
T_1_8 := combine(Matrix(Trf(1..4,1..4, 2)) . Matrix(Trf(1..4,1..4,3)). Matrix(Trf(1..4,1..4, 4)).Matrix(Trf(1..4,1..4, 8)));
# Schleife -9-7-1
T_1_9:= combine(Matrix(Trf(1..4,1..4, 7)) . Matrix(Trf(1..4,1..4,9)));
h1t := T_1_8(1..3,4) - T_1_9(1..3,4);
tmp := Transpose( Matrix(T_1_8(1..3,1..3)) ) . Matrix(T_1_9(1..3,1..3));  # nur anzeigen lassen f�r h1r
;
 combine(tmp); # nur anzeigen lassen f�r h1r
;
h1r:=5*Pi/2+(qJ2(t)+qJ4(t)+qJ8(t)-qJ7(t));
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := <h1t([1, 3]);h1r>; 
implconstr_s := convert_t_s(implconstr_t):
# Exportiere Code f�r folgende Skripte
kin_constraints_exist:=true:
save implconstr_t, implconstr_s, kin_constraints_exist, sprintf("../codeexport/%s/tmp/kinematic_constraints_implicit_maple.m", robot_name):
# Exportieren des vollst�ndigen Ausdruckes
if codegen_act then
  MatlabExport(implconstr_s, sprintf("../codeexport/%s/tmp/kinconstr_impl_matlab.m", robot_name), 2):
end if:
# Liste mit abh�ngigen konstanten Kinematikparametern erstellen (wichtig f�r Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( implconstr_s )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_matlab.m", robot_name), 2):


