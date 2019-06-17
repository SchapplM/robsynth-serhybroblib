
# Vierachsroboter KUKA KR700PA implizit
# Initialisierung
# Import (Library)
interface(warnlevel=0): 		# Unterdrücke die folgende Warnung.
restart: 				# Gibt eine Warnung, wenn über Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true:		# Für Speicherung
;
with(StringTools): 			# Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 1: 			# Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0:		        # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
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
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
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
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := Matrix(<theta(3) +theta(2);theta(4)+theta(5)>): # TODO: In h1r, h2r muss das richtige drinstehen.
implconstr_s := convert_t_s(implconstr_t);
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
