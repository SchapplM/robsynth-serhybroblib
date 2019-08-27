
# TSR2040 
# Voraussetzung:
# Viergelenkkette muss mit der Toolbox berechnet worden sein (Arbeitsblatt "vgk_kinematic_constraints.mw")
# Autor
# Weipu Shan
# Studienarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
# (C) Institut fuer mechatronische Systeme, Leibniz Universitaet Hannover
# Quellen
# SA Bejaoui: Bejaoui2018_S749; "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen"
# SA Shan: Shan2019_S828, "Reduktion der Modellkomplexität von seriell-hybriden Palettierrobotern"
# 
# 
# Initialisierung
# Import (Library)
interface(warnlevel=0): 		# Unterdrücke die folgende Warnung.
restart: 				# Gibt eine Warnung, wenn über Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true: 		# Für Speicherung
;
with(StringTools): 			# Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
#with(ListTools):
codegen_act := true:
codegen_opt := 1: # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
;

# Import (hybriddyn)
read "../helper/proc_MatlabExport": 	# Exportiert einen Term als Matlab-Code in optimierter Form
read "../transformation/proc_rotx":  	# rotation um X
read "../transformation/proc_roty":  	# rotation um Y
read "../transformation/proc_rotz":  	# rotation um Z
;
read "../helper/proc_convert_s_t":   	# zeitabhängige Variable für Ableitung
read "../helper/proc_convert_t_s":   	# konstante Variable 
read "../robot_codegen_constraints/proc_subs_kintmp_exp": #substitute für Spalte 1 mit Spalte 2
;
with(RealDomain): 			# Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":			#aktuelle Roboter, MDH-Tabelle
;
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions",robot_name):	#von Fourbar
;

# Variable Initialisierung
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
#qJ_t:= <qJ1(t),qJ2(t),qJ3(t),qJ4(t)>:			#zeitabhängige Variable
;
#qJ_s:= <qJ1s,qJ2s,qJ3s,qJ4s>:				#konstante variable
;
# Ergebnisse von Trigonometrischer Elimination lesen
read sprintf("../codeexport/TSRm2TE/tmp/kinematic_constraints_maple_inert.m"):	# von aktuellem Roboter
;
kin_constraints_exist := kin_constraints_exist:
kintmp_qs := kintmp_qs:
kintmp_qt := kintmp_qt:
kintmp_subsexp := kintmp_subsexp:
# Variable entfernen
kintmp_subsexp:= Matrix(2*RowDimension(kintmp_s),2):
# Export
kintmp_qt := convert_s_t(kintmp_qs):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
#printf("Ausdrücke für kintmp_subsexp gespeichert (Maple). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
for i from 1 to RowDimension(kintmp_s) do
  tmp := kintmp_qs(i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d",robot_name, i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d.m", robot_name, i):
end do:
save kin_constraints_exist, kintmp_qs, kintmp_qt,kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert" ,robot_name):
save kin_constraints_exist, kintmp_qs, kintmp_qt, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
#printf("Ausdrücke mit Inert-Arctan exportiert (Matlab). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_qs ));
#kc_symbols :=Transpose(kc_symbols);
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);

#printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
kintmp_qs(1,1):
kc_symbols(1,1):
kintmp_qs(2,1):
