
# Vierachsroboter KUKA KR700PA implizit
# Initialisierung
# Import (Library)
interface(warnlevel=0): 		# Unterdr�cke die folgende Warnung.
restart: 				# Gibt eine Warnung, wenn �ber Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true:		# F�r Speicherung
;
with(StringTools): 			# F�r Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 1: 			# Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0:		        # Zur Code-Generierung auch f�r Nicht-Inert-Ausdr�cke
;
# Import (Hybriddyn)
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
read "../robot_codegen_definitions/robot_env_IC":			#aktuelle Roboter, MDH-Tabelle
;
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name_OL):	#von Fourbar
;
# Ergebnisse der Kinematik laden
read sprintf("../codeexport/%s/tmp/kinematics_floatb_%s_rotmat_maple.m", robot_name_OL, base_method_name);
read "../robot_codegen_definitions/robot_env_IC":
Trf := Trf:
Trf_c := Trf_c:
Trf:
# VGK Gelb 1-2-7-11-9-6-1
# Schleife (1-6)-(6-11)
T_1_11 := combine( Matrix(Trf(1..4,1..4, 2)) . Matrix(Trf(1..4,1..4, 7)) . Matrix(Trf(1..4,1..4,11))):
# Schleife (1-2)(2-7)(7-9)
T_1_9 := combine(Matrix(Trf(1..4,1..4,6)) . Matrix(Trf(1..4, 1..4, 9))):
h1t := T_1_11(1..3,4) - T_1_9(1..3,4);
tmp := Transpose( Matrix(T_1_11(1..3,1..3)) ) . Matrix(T_1_9(1..3,1..3)): #nur anzeigen lassen f�r h1r
;
combine(tmp);# nur anzeigen lassen f�r h1r
;
h1r := -(-qJ6(t)+qJ2(t)+qJ7(t)+phi79+qJ9(t)) + Pi;
# VGK GR�N 2-7-8-10-12-4-3-2
# Schleife (2-3)-(3-4)-(4-10)
T_2_10 := combine( Matrix(Trf(1..4,1..4, 7)) . Matrix(Trf(1..4,1..4,8)) . Matrix(Trf(1..4,1..4,10))):
# Schleife (2-7)-(7-8)-(8-12)
T_2_12:= combine( Matrix(Trf(1..4,1..4, 3)) . Matrix(Trf(1..4,1..4,4)).Matrix(Trf(1..4,1..4,12)) ):
h2t := T_2_10(1..3,4) - T_2_12(1..3,4);
tmp := Transpose( Matrix(T_2_10(1..3,1..3)) )  .  Matrix(T_2_12(1..3,1..3)):  #nur anzeigen lassen f�r h2r
;
combine(tmp);  #nur anzeigen lassen f�r h2r
;
h2r := (qJ3(t)+qJ4(t)+phi410+qJ10(t)-qJ7(t)+phi78-qJ8(t))+Pi;
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := <h1t([1, 3],1);h2t([1, 2],1); h1r; h2r>; # TODO: In h1r, h2r muss das richtige drinstehen.
implconstr_s := convert_t_s(implconstr_t);
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
