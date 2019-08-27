
# Berechne kinematische Zwangsbedingungen für Viergelenkkette mit Parallelogram
# Einleitung

# Die kinematischen Zwangsbedingungen werden als Ersetzungsausdruck für die abhängigen Winkel aufgestellt.
# 
# fourbar2TE -> Viergelenkkette, modellierung der Zwangsbedingungen mit Ausdrücken für trigonometrische Elimination
# kinematic_constraint -> Kinematische Zwangsbedingungen
# Quelle
# SA Bejaoui: Bejaoui2018_S749; "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen"
# SA Shan: Shan2019_S828, "Reduktion der Modellkomplexität von seriell-hybriden Palettierrobotern"
# Autor
# todo
# Abderahman Bejaoui
# Weipu Shan
# Studienarbeit bei: Moritz Schappler, schappler@irt.uni-hannover.de, 2018-08
# (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover
# 
# Initialisierung
interface(warnlevel=0): 	# Unterdrücke die folgende Warnung.
restart: 			# Gibt eine Warnung, wenn über Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true: 	# Für Speicherung
;
with(StringTools): 		# Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
#with(ListTools):
codegen_act := true:
codegen_opt := 1: 		# Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: 		  # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
;
read "../helper/proc_MatlabExport":
read "../transformation/proc_rotx":
read "../transformation/proc_roty":
read "../transformation/proc_rotz":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
read "../robot_codegen_constraints/proc_subs_kintmp_exp":
read "../helper/proc_intersect_circle":
with(RealDomain): 		# Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
kintmpD_t:=diff~(kintmp_t,t):
kintmpD_s:=Matrix(< etaD_s, xiD_s, rhoD_s>):
qJ_t := Matrix(NQJ, 1, qJ1(t)):
qJ_s := Matrix(NQJ, 1, qJ1s):
# Konstante Winkel bereits hineinschreiben ( Keine für Viergelenkkette)
for i from 1 to RowDimension(kintmp_s) do
  if diff(kintmp_s(i,1), t) = 0 then
    kintmp_qs(i,1) := kintmp_s(i,1):
  end if:
end do:
# Variablen definieren für die Hilfswinkel
# 
# Ersetzungsausdrücke definieren.
# Speichere Sinus und Cosinus der Winkel direkt ab, da diese in den Rotationsmatrizen direkt auftreten.
# Spalte 1: Zu suchender Ausdruck (sin oder cos eines Winkels)
# Spalte 2: Einzusetzender Ausdruck.
# Dadurch werden arctan-Ausdrücke in der direkten Kinematik reduziert.
# Ähnliches Vorgehen wie in [1].
kintmp_subsexp := Matrix(2*RowDimension(kintmp_s),2):
for i from 1 to RowDimension(kintmp_s) do
  kintmp_subsexp(2*i-1, 1) := sin(kintmp_s(i,1)):
  kintmp_subsexp(2*i,   1) := cos(kintmp_s(i,1)):
  # Initialisierung der rechten Spalte mit gleichen Werten. Später nur Ersetzung, wenn Vorteilhaft.
  kintmp_subsexp(2*i-1, 2) := kintmp_subsexp(2*i-1, 1):
  kintmp_subsexp(2*i,   2) := kintmp_subsexp(2*i,   1):
end do:
printf("Beginn der Berechnungen. %s\n", FormatTime("%Y-%m-%d %H:%M:%S")):
st := time():
# Reihenfolge, in der die Winkel der Parallelstruktur berechnet werden. Hilft beim exportieren von Code, wenn beim Debuggen in der gleichen Reihenfolge ausgegeben wird.
Reihenfolge_kintmp := <1;3;2>:
# Bestimme Winkel
sin(rho_s) := sin(Pi-qJ1s):
cos(rho_s) := cos(Pi-qJ1s):
sin(eta_s) := sin(qJ1s):
cos(eta_s) := cos(qJ1s):
sin(xi_s) := sin(qJ1s):
cos(xi_s) := cos(qJ1s):
# kintmp_subsexp
kintmp_subsexp(1,2):=sin(rho_s):
kintmp_subsexp(2,2):=cos(rho_s):
kintmp_subsexp(3,2):=sin(eta_s):
kintmp_subsexp(4,2):=cos(eta_s):
kintmp_subsexp(5,2):=sin(xi_s):
kintmp_subsexp(6,2):=cos(xi_s):
# Kintmp_qs
kintmp_qs(1,1):=arctan(sin(rho_s),cos(rho_s)):
kintmp_qs(2,1):=arctan(sin(eta_s),cos(eta_s)):
kintmp_qs(3,1):=arctan(sin(xi_s),cos(xi_s)):
# Exportiere Code für folgende Skripte
# Die Annahmen sind im Ausdruck bereits in den Variablen gespeichert. Ersetze das "~"-Zeichen.
# Quelle: http://www.mapleprimes.com/questions/207601-Remove-Assumptions
winkel := <<sin_xi|sin(xi_s)>; <cos_xi|cos(xi_s)>; <sin_rho | sin(rho_s)>; <cos_rho | cos(rho_s)>;<sin_eta | sin(eta_s)>; <cos_eta | cos(eta_s)>>:
#mkdir(sprintf("../codeexport/foubar1TE/"));
#mkdir(sprintf("../codeexport/foubar1TE/tmp/"));
save winkel, sprintf("../codeexport/%s/tmp/%s_angles_trig_elim", robot_name, robot_name):
#save winkel, sprintf("../codeexport/fourbar1TE/tmp/kinematic_constraints_maple_inert"):
kintmp_qt := convert_s_t(kintmp_qs):
# Speichere Maple-Ausdruck (Eingabe-Format und internes Format)
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
printf("Ausdrücke für kintmp_subsexp gespeichert (Maple). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
for i from 1 to RowDimension(kintmp_s) do
  tmp := kintmp_qs(i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d", robot_name, i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d.m", robot_name, i):
end do:
save kin_constraints_exist, kintmp_qs, kintmp_qt,kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert", robot_name):
save kin_constraints_exist, kintmp_qs, kintmp_qt, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
printf("Ausdrücke mit Inert-Arctan exportiert (Matlab). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_subsexp ));
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);
printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):


