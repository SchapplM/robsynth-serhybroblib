
# Kinematik-Berechnung Modell PalH2m1 (Hybrider Palettierer, Variante 2, Modellierung 1)
# Beschreibung
# Ein Palettierroboter mit f�nf Achsen als N�herung f�r einen Palettierroboter mit Viergelenkketten
# Anpassung der Parameter an Konvention von Lenze: Parallele Verschiebung von Achse 3 zum Boden
# Die vierte Achse des Ersatzmodells wird als abh�ngiges Gelenk in Abh�ngigkeit der Achsen 2 und 3 beschrieben.
# Dadurch wird der Effekt der Parallelogramme kinematisch nachgebildet.
# Die Gelenkmomente sollten genauso wie durch die Parallelogramme �bertragen werden (als ob deren Dynamikparameter Null sind).
# Quellen
# Autor
# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
# (C) Institut f�r Mechatronische Systeme, Universit�t Hannover
# Initialisierung
interface(warnlevel=0): # Unterdr�cke die folgende Warnung.
restart: # Gibt eine Warnung, wenn �ber Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true: # F�r Speicherung
;
with(StringTools): # F�r Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 2: # Hoher Optimierungsgrad.
;
read "../helper/proc_MatlabExport":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
with(RealDomain): # Schr�nkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Variable mit Winkeln der Zwangsbedingungen nur in Abh�ngigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
kintmp_qt := Matrix(RowDimension(kintmp_s),1):
# Ersetzungsausdr�cke definieren.
# Variable zum speichern des Sinus und Cosinus der Winkel. F�r dieses System ist das eigentlich nicht notwendig. Erstelle Variable, da sie von den anderen Skripten erwartet wird
kintmp_subsexp := Matrix(2*RowDimension(kintmp_s),2):
# Zwangsbedingung f�r abh�ngigen Gelenkwinkel
# Der vierte Winkel ist der passive Winkel nahe des Endeffektors. Er wird durch die Parallelstruktur (aus zwei Parallelogrammen) bei hybriden Palettierrobotern bestimmt.
# Durch die Parallelogramme ergibt sich ein konstantes Verh�ltnis zu den Gelenkwinkeln der Achsen 2 und 3, die auch in derselben Ebene liegen.
rho4_qt := -theta(2):
rho5_qt := -theta(4):
kintmp_qt(1,1) := rho4_qt:
kintmp_qt(2,1) := rho5_qt:
# Nachverarbeitung
# Umrechnung in Substitutionsvariablen
kintmp_qs := convert_t_s(kintmp_qt):
# Exportiere Code f�r folgende Skripte
# Speichere Maple-Ausdruck (Eingabe-Format und internes Format)
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
printf("Ausdr�cke f�r kintmp_subsexp gespeichert (Maple)\n"):
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert", robot_name):
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
# Liste mit abh�ngigen konstanten Kinematikparametern erstellen (wichtig f�r Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_qs )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m", robot_name), 2):
printf("Fertig\n"):


