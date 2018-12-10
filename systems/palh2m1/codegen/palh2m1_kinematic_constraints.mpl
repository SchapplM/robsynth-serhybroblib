
# Kinematik-Berechnung Modell PalH2m1 (Hybrider Palettierer, Variante 2, Modellierung 1)
# Beschreibung
# Ein Palettierroboter mit fünf Achsen als Näherung für einen Palettierroboter mit Viergelenkketten
# Die vierte Achse des Ersatzmodells wird als abhängiges Gelenk in Abhängigkeit der Achsen 2 und 3 beschrieben.
# Dadurch wird der Effekt der Parallelogramme kinematisch nachgebildet.
# Die Gelenkmomente sollten genauso wie durch die Parallelogramme übertragen werden (als ob deren Dynamikparameter Null sind).
# Quellen
# Autor
# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
# (C) Institut für Mechatronische Systeme, Universität Hannover
# Initialisierung
interface(warnlevel=0): # Unterdrücke die folgende Warnung.
restart: # Gibt eine Warnung, wenn über Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true: # Für Speicherung
;
with(StringTools): # Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 2: # Hoher Optimierungsgrad.
;
read "../helper/proc_MatlabExport":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Variable mit Winkeln der Zwangsbedingungen nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
kintmp_qt := Matrix(RowDimension(kintmp_s),1):
# Ersetzungsausdrücke definieren.
# Variable zum speichern des Sinus und Cosinus der Winkel. Für dieses System ist das eigentlich nicht notwendig. Erstelle Variable, da sie von den anderen Skripten erwartet wird
kintmp_subsexp := Matrix(2*RowDimension(kintmp_s),2):
# Zwangsbedingung für abhängigen Gelenkwinkel
# Der vierte Winkel ist der passive Winkel nahe des Endeffektors. Er wird durch die Parallelstruktur (aus zwei Parallelogrammen) bei hybriden Palettierrobotern bestimmt.
# Durch die Parallelogramme ergibt sich ein konstantes Verhältnis zu den Gelenkwinkeln der Achsen 2 und 3, die auch in derselben Ebene liegen.
rho4_qt := -theta(2) -theta(3):
kintmp_qt(1,1) := rho4_qt:
kintmp_qs := convert_t_s(kintmp_qt):
# Exportiere Code für folgende Skripte
# Speichere Maple-Ausdruck (Eingabe-Format und internes Format)
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
printf("Ausdrücke für kintmp_subsexp gespeichert (Maple)\n"):
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert", robot_name):
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_qs )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m", robot_name), 2):
printf("Fertig\n"):

