
# Berechne kinematische Zwangsbedingungen für Fünfgelenkkette

# Einleitung
# Die kinematischen Zwangsbedingungen werden als Ersetzungsausdruck für die abhängigen Winkel aufgestellt.
# 
# hybBKplanar -> Beinkette für hybride PKM; Fünfgelenkkette, Modellierung der Zwangsbedingungen mit trigonometrischer Elimination der Winkel
# kinematic_constraint -> Berechnungen zu kinematischen Zwangsbedingungen
# 
# TODO: Dieses System unterscheidet sich nur minimal von der reinen Fünfgelenkkette (Drehgelenk am Ende). Lade die Ergebnisse der Fünfgelenkkette statt alles identisch neu zu berechnen.

# Quelle
# MA Brünger: Brünger2019_M832, "Vergleich der Leistungsmerkmale und Modellierung von parallelen Robotern mit hybriden und seriellen Beinketten"
# SA Bejaoui: Bejaoui2018_S749; "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen"
# (Fünfgelenkkette S. 23 in [Bejaoui2018_S749])
# Haruhiko Asada;  Kamal Youcef-Toumi: "Analysis and design of a direct-drive arm with a five-bar-link parallel drive mechanism" (1984)
# Autor
# André Brünger, Moritz Schappler
# Masterarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-01
# (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover
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
#with(ListTools):
codegen_act := true:
codegen_opt := 2:
read "../helper/proc_MatlabExport":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
read "../robot_codegen_constraints/proc_subs_kintmp_exp":
read "../helper/proc_intersect_circle":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions",robot_name):
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
for i from 1 to RowDimension(kintmp_s) do
  if diff(kintmp_s(i,1), t) = 0 then
    kintmp_qs(i,1) := kintmp_s(i,1):
  end if:
end do:
# Ersetzungsausdrücke definieren.
# Speichere Sinus und Cosinus der Winkel direkt ab, da diese in den Rotationsmatrizen direkt auftreten.
# Spalte 1: Zu suchender Ausdruck (sin oder cos eines Winkels)
# Spalte 2: Einzusetzender Ausdruck.
# Dadurch werden arctan-Ausdrücke in der direkten Kinematik reduziert.
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
# Zwangsbedingungen mit Kreisschnittpunkt D zur Bestimmung von rho
r_p1:=<AE*cos(qJ1s), AE*sin(qJ1s)>:
r_p2:=<AB+BC*cos(qJ2s), BC*sin(qJ2s)>:
r1:=ED:
r2:=CD:
A:=intersect_circle(r_p1, r_p2, r1, r2):
# Es gibt zwei Lösungen für den Kreisschnittpunkt. Bei der gewählten Definition der Kreise ist die erste Lösung in der Fünfgelenkkette und die zweite außen (für die gewählten Beispiel-Parameter).
# Die zweite Lösung entspricht der gewählten Lösung in [Bejaoui2018_S749]
r_s1x:=simplify(A(1,2),trig):
r_s1y:=simplify(A(2,2),trig):
# Zwangsbedingungen für eta, zweiter Weg zum Kreisschnittpunkt D über B
radx1:=AB+BC*cos(qJ2s)+CD*(cos(qJ2s)*cos(eta_s)-sin(qJ2s)*sin(eta_s)):
rady1:=BC*sin(qJ2s)+CD*(sin(qJ2s)*cos(eta_s)+cos(qJ2s)*sin(eta_s)):
f11:=r_s1x-radx1=0: #GL3
;
f12:=r_s1y-rady1=0: #GL4
;
# Zwangsbedingungen für xi , zwei Wege Punkt E
racx1:=AE*cos(qJ1s):
racy1:=AE*sin(qJ1s):
racx2:=AB+BC*cos(qJ2s)+CD*(cos(qJ2s)*cos(eta_s)-sin(qJ2s)*sin(eta_s))+ED*(cos(qJ2s)*(cos(eta_s)*cos(xi_s)-sin(eta_s)*sin(xi_s))-sin(qJ2s)*(cos(xi_s)*sin(eta_s)+cos(eta_s)*sin(xi_s))):
racy2:=BC*sin(qJ2s)+CD*(sin(qJ2s)*cos(eta_s)+cos(qJ2s)*sin(eta_s))+ED*(sin(qJ2s)*(cos(eta_s)*cos(xi_s)-sin(eta_s)*sin(xi_s))+cos(qJ2s)*(cos(xi_s)*sin(eta_s)+cos(eta_s)*sin(xi_s))):
f21:=racx1-racx2=0:  # GL5
;
f22:=racy1-racy2=0:  # GL6
;
# Bestimme Winkel
# eta bestimmen
sol1:=solve({f11,f12},{cos(eta_s),sin(eta_s)}): 
assign(sol1):
cos_eta_s:=simplify(rhs(sol1[1]),trig):
sin_eta_s:=simplify(rhs(sol1[2]),trig):
# xi bestimmen
sol2:=solve({f21,f22},{cos(xi_s),sin(xi_s)}):
sol2:=op(sol2):
cos_xi_s:=simplify(rhs(sol2[1]),trig):
sin_xi_s:=simplify(rhs(sol2[2]),trig):
# rho bestimmen
cos_psi_phi:=cos(qJ2s)*cos(qJ1s)+sin(qJ2s)*sin(qJ1s): # cos(psi-phi)
;
cos_eta_xi:=cos_eta_s*cos_xi_s-sin_eta_s*sin_xi_s: # cos(eta+xi)
;
sin_psi_phi:=sin(qJ2s)*cos(qJ1s)-cos(qJ2s)*sin(qJ1s): # sin(psi-phi)
;
sin_eta_xi:=sin_eta_s*cos_xi_s+cos_eta_s*sin_xi_s: #sin(eta+xi)
;
cos_rho_s:=cos_psi_phi*cos_eta_xi-sin_psi_phi*sin_eta_xi:
sin_rho_s:=sin_psi_phi*cos_eta_xi+cos_psi_phi*sin_eta_xi:
# Cos-und Sinus-Terme zuweisen
kintmp_subsexp(1,2):=sin_rho_s:
kintmp_subsexp(2,2):=cos_rho_s:
kintmp_subsexp(3,2):=sin_eta_s:
kintmp_subsexp(4,2):=cos_eta_s:
kintmp_subsexp(5,2):=sin_xi_s:
kintmp_subsexp(6,2):=cos_xi_s:
# Winkel zuweisen:
kintmp_qs(1,1):=%arctan(sin_rho_s,cos_rho_s):
kintmp_qs(2,1):=%arctan(sin_eta_s,cos_eta_s):
kintmp_qs(3,1):=%arctan(sin_xi_s,cos_xi_s):
# Exportiere Code für folgende Skripte
# 
# Die Annahmen sind im Ausdruck bereits in den Variablen gespeichert. Ersetze das "~"-Zeichen.
# Quelle: http://www.mapleprimes.com/questions/207601-Remove-Assumptions
winkel := <<sin_xi|sin(xi_s)>; <cos_xi|cos(xi_s)>; <sin_rho | sin(rho_s)>; <cos_rho | cos(rho_s)>;<sin_eta | sin(eta_s)>; <cos_eta | cos(eta_s)>>:
save winkel, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert", robot_name):
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
printf("Ausdrücke mit Inert-Arctan exportiert (Matlab). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_subsexp ));
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);
printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
kintmp_qs:

# 
