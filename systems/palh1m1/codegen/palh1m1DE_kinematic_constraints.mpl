
# Palettierroboter mit drei geschlossenen kinematischen Ketten
# Voraussetzung:
# Viergelenkkette muss mit der Toolbox berechnet worden sein (Arbeitsblatt "fourbar1DE_kinematic_constraints.mw")
# 
# Erklärung des Dateinamens:
# palh1: Hybrider Palettierroboter Variante 1 (drei geschlossene kinematische Ketten)
# m1: Modellierung der geschlossenen Ketten als allgemeine Viergelenkkette (ohne Parallelogramm-Annahme)
# DE: Berechnung mit direkter Elimination
# kinematic_contraints: Berechnung der kinematischen Zwangsbedingungen (explizite Form)
# Autor
# Abderahman Bejaoui
# Studienarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
# (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover
# Quellen
# TODO
# 
# Initialisierung
# Import 
restart:
kin_constraints_exist := true: # Für Speicherung
;
with(StringTools): # Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
#with(ListTools):
codegen_act := true:
codegen_opt := 1: # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
;
read "../helper/proc_MatlabExport":
read "../transformation/proc_rotx":
read "../transformation/proc_roty":
read "../transformation/proc_rotz":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
read "../robot_codegen_constraints/proc_subs_kintmp_exp":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions",robot_name):
# Ergebnisse von Trigonometrischer Elimination lesen
read sprintf("../codeexport/palh1m1TE/tmp/kinematic_constraints_maple_inert.m"):
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
kintmp_qs(1,1);
kc_symbols(1,1);
kintmp_qs(2,1);

