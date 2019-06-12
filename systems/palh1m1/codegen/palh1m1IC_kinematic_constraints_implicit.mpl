
# Vierachsroboter implizit
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
read "../robot_codegen_definitions/robot_env_IC":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name_OL):
# Ergebnisse der Kinematik laden
read sprintf("../codeexport/%s/tmp/kinematics_floatb_%s_rotmat_maple.m", robot_name_OL, base_method_name);
Trf := Trf:
Trf_c := Trf_c:
Trf
;
# VGK Gelb 1-6-11-7-2-1
# Schleife (1-6)-(6-14)
T_1_14 := combine( Matrix(Trf(1..4,1..4, 6)) . Matrix(Trf(1..4,1..4,14))):
# Schleife (1-2)-(2-7)-(7-11)
T_1_11:= combine( Matrix(Trf(1..4,1..4, 2)) . Matrix(Trf(1..4,1..4,7)) . Matrix(Trf(1..4,1..4,11))):
h1t := T_1_14(1..3,4) - T_1_11(1..3,4);

#tmp := Transpose( Matrix(T_1_14(1..3,1..3)) ) . Matrix(T_1_11(1..3,1..3));  nur anzeigen lassen für h1r
;
#combine(tmp);# nur anzeigen lassen für h1r
;
h1r := -(-qJ6(t)+qJ2(t)+qJ7(t)+phi711+qJ11(t))+Pi/2:
# VGK PINK 2-3-12-15-9-8-2
# Schleife (2-3)-(3-12)
T_2_12 := combine(Matrix(Trf(1..4,1..4, 3)) . Matrix(Trf(1..4,1..4,12))):
# Schleife (2-8)-(8-9)-(9-15)
T_2_15:= combine(Matrix(Trf(1..4,1..4, 8)) . Matrix(Trf(1..4,1..4,9)).Matrix(Trf(1..4,1..4, 15))):
h2t := T_2_12(1..3,4) - T_2_15(1..3,4):
tmp := Transpose( Matrix(T_2_12(1..3,1..3)) ) . Matrix(T_2_15(1..3,1..3));  # nur anzeigen lassen für h1r
;
 combine(tmp); # nur anzeigen lassen für h1r
;
h2r:=(qJ3(t)+phi312+qJ12(t)-qJ8(t)-qJ9(t))+3*Pi/2:
# VGK GRÜN 2-3-4-13-16-10-7-2
# Schleife (2-3)-(3-4)-(4-13)
T_2_13 := combine( Matrix(Trf(1..4,1..4, 3)) . Matrix(Trf(1..4,1..4,4)) . Matrix(Trf(1..4,1..4,13))):
# Schleife (2-7)-(7-10)-(10-16)
T_2_16:= combine( Matrix(Trf(1..4,1..4, 7)) . Matrix(Trf(1..4,1..4,10)).Matrix(Trf(1..4,1..4,16)) ):
h3t := T_2_13(1..3,4) - T_2_16(1..3,4):
#tmp := Transpose( Matrix(T_2_13(1..3,1..3)) ) . Matrix(T_2_16(1..3,1..3));  nur anzeigen lassen für h2r
;
#combine(tmp);  nur anzeigen lassen für h3r
;
h3r:=-(qJ3(t)+qJ4(t)+phi413+qJ13(t)-qJ7(t)+phi710-qJ10(t))+5*Pi/2:
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := <h1t([1, 3],1);h2t([1, 2],1);h3t([1, 2],1); h1r; h2r; h3r>; # TODO: In h1r, h2r, h3r muss das richtige drinstehen.
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

