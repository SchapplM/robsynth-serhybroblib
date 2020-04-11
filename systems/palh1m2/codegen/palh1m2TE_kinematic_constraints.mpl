
# Berechne kinematische Zwangsbedingungen für den Palettierroboter MPL800-YASKAWA aus 3 Parallelogrammen
# Einleitung

# Die kinematischen Zwangsbedingungen werden als Ersetzungsausdruck für die abhängigen Winkel aufgestellt.
# 
# palh1m2TE -> MPL800 II-Yaskawa, modellierung der Zwangsbedingungen mit Ausdrücken für trigonometrische Elimination
# kinematic_constraint -> Kinematische Zwangsbedingungen
# Datenblatt des Roboters unter :   https://www.motoman.com/industrial-robots/mpl800-ii
# Voraussetzung
# Parallelogram muss mit der Toolbox berechnet worden sein (Arbeitsblatt "fourbar2TE_kinematic_constraints.mw")
# Quelle
# SA Bejaoui: Bejaoui2018_S749; "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen"
# Autor
# Abderahman Bejaoui
# Studienarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
# (C) Institut fuer mechatronische Systeme, Leibniz Universitaet Hannover
# Initialisierung
# Import 
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
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
# Konstante Winkel bereits hineinschreiben 
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

# Sicherung der aktuellen Daten
backup_kintmp_qs := kintmp_qs:
backup_kintmp_qt := kintmp_qt:
backup_kintmp_subsexp := kintmp_subsexp:

read sprintf("../codeexport/fourbar2TE/tmp/kinematic_constraints_maple_inert.m"): # mit .m datei wird es schneller
;
kintmp_subsexp_fourbar := kintmp_subsexp:
kintmp_qs_fourbar := kintmp_qs:
kintmp_qt_fourbar := kintmp_qt:
kintmp_subsexp_fourbar:
winkel(6,2):
kintmp_subsexp_fourbar(6,2):
kintmp_subsexp_fourbar(3,2):
winkel_neu := Matrix(6,2):
winkel_neu(1,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(3,2)):
winkel_neu(2,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(4,2)):
winkel_neu(3,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(1,2)):
winkel_neu(4,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(2,2)):
winkel_neu(5,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(5,2)):
winkel_neu(6,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(6,2)):
# Vorherige Werte für dieses System wieder herstellen
kintmp_qs := backup_kintmp_qs:
kintmp_qt := backup_kintmp_qt:
kintmp_subsexp := backup_kintmp_subsexp:
winkel_neu(1..6,1):
# Viergelenkkette 1: A-M-L-B
winkel1:=<<sin_rho28|sin(rho28_s)>;<cos_rho28|cos(rho28_s)>;<sin_rho312|sin(rho312_s)>;<cos_rho312|cos(rho312_s)>;<sin_rho89|sin(rho89_s)>;<cos_rho89|cos(rho89_s)>>:
for i from 1 to 6 do
	winkel1(i,2):=kintmp_subsexp_fourbar(i,2):  
     cos_phi_s :=-(sin(qJ3s)*cos(phi312)+cos(qJ3s)*sin(phi312)):
	sin_phi_s := -(cos(qJ3s)*cos(phi312)-sin(qJ3s)*sin(phi312)):
	winkel1(i,2):=subs({sin(qJ1s)=sin_phi_s},winkel1(i,2)): 
	winkel1(i,2):=subs({cos(qJ1s)=cos_phi_s},winkel1(i,2)):
     winkel1(i,2) := subs({l1=BL}, winkel1(i,2)):
     winkel1(i,2) := subs({l2=AB}, winkel1(i,2)):
     winkel1(i,2) := subs({l3=AM}, winkel1(i,2)):
     winkel1(i,2) := subs({l4=ML}, winkel1(i,2)):
end do:
winkel1:
# Umrechnung aus Zwangsbedingugen für Winkel aus MDH-Tabelle
sin_rho28_s:=winkel1(1,2): # schon in MDH
;
cos_rho28_s:=winkel1(2,2):  # schon in MDH
;
sin_rho312_s:=winkel1(3,2):  # schon in MDH
;
cos_rho312_s:=winkel1(4,2):  # schon in MDH  
;
sin_rho89_s:=winkel1(5,2):    # für MDH , am 20.11 der neuen fourbar1 angepasst
;
cos_rho89_s:=winkel1(6,2):  # für MDH , am 20.11 der neuen fourbar1 angepasst
;

#sin_rho312_s:=sin_rho3_s:  # schon in MDH
;
#cos_rho312_s:=-cos_rho3_s:  # schon in MDH
;
# Viergelenkkette 2: A-B-C-D
winkel2:=<<sin_eta1|sin(eta1_s)>;<cos_eta1|cos(eta1_s)>;<sin_eta3|sin(eta3_s)>;<cos_eta3|cos(eta3_s)>;<sin_eta4|sin(eta4_s)>;<cos_eta4|cos(eta4_s)>>:
for i from 1 to 6 do
	winkel2(i,2):=kintmp_subsexp_fourbar(i,2):  #
	cos_phi_s :=sin(qJ2s)*cos(phi2)-cos(qJ2s)*sin(phi2):
	sin_phi_s :=cos(qJ2s)*cos(phi2)+sin(qJ2s)*sin(phi2):
	winkel2(i,2):=subs({sin(qJ1s)=sin_phi_s},winkel2(i,2)): 
	winkel2(i,2):=subs({cos(qJ1s)=cos_phi_s},winkel2(i,2)):
     winkel2(i,2) := subs({l1=AB}, winkel2(i,2)):
     winkel2(i,2) := subs({l2=DA}, winkel2(i,2)):
     winkel2(i,2) := subs({l3=DC}, winkel2(i,2)):
     winkel2(i,2) := subs({l4=BC}, winkel2(i,2)):
 end do:
winkel2(1..6,1):
sin_eta1_s:=winkel2(1,2):
cos_eta1_s:=winkel2(2,2):
sin_eta3_s:=winkel2(3,2):
cos_eta3_s:=winkel2(4,2):
sin_eta4_s:=winkel2(5,2): 
cos_eta4_s:=winkel2(6,2):

cos_eta16_s:=cos(phi1)*cos_eta1_s-sin(phi1)*sin_eta1_s:
sin_eta16_s:=cos(phi1)*sin_eta1_s+cos_eta1_s*sin(phi1):
cos_eta711_s:=-cos_eta4_s:
sin_eta711_s:=sin_eta4_s:
cos_eta27_s:=cos(phi711)*cos_eta3_s+sin(phi711)*sin_eta3_s: # am 20.11 der neuen fourbar1 angepasst
;
sin_eta27_s:=cos(phi711)*sin_eta3_s-cos_eta3_s*sin(phi711): # am 20.11 der neuen fourbar1 angepasst
;
# Viergelenkkette 3: B-E-P-G 
winkel3:=<<sin_xi1016|sin(xi1016_s)>;<cos_xi1016|cos(xi1016_s)>;<sin_xi3|sin(xi3_s)>;<cos_xi3|cos(xi3_s)>;<sin_xi2|sin(xi2_s)>;<cos_xi2|cos(xi2_s)>>:
for i from 1 to 6 do
	winkel3(i,2):=kintmp_subsexp_fourbar(i,2):  
	cos_phi_s :=(cos(phi711+phi710)*(sin(qJ3s)*cos_eta3_s-cos(qJ3s)*sin_eta3_s)+sin(phi711+phi710)*(cos(qJ3s)*cos_eta3_s+sin(qJ3s)*sin_eta3_s)):
	sin_phi_s := (cos(phi711+phi710)*(cos(qJ3s)*cos_eta3_s+sin(qJ3s)*sin_eta3_s)-sin(phi711+phi710)*(sin(qJ3s)*cos_eta3_s-cos(qJ3s)*sin_eta3_s)):
	winkel3(i,2):=subs({sin(qJ1s)=sin_phi_s},winkel3(i,2)): 
	winkel3(i,2):=subs({cos(qJ1s)=cos_phi_s},winkel3(i,2)):
     winkel3(i,2) := subs({l1=BG}, winkel3(i,2)):
     winkel3(i,2) := subs({l2=BE}, winkel3(i,2)):
     winkel3(i,2) := subs({l3=EP}, winkel3(i,2)):
     winkel3(i,2) := subs({l4=GP}, winkel3(i,2)):
 end do:
winkel3(1..6,1):
sin_xi1016_s:=winkel3(1,2):
cos_xi1016_s:=winkel3(2,2):
sin_xi3_s:=winkel3(3,2):  # schon in MDH
;
cos_xi3_s:=winkel3(4,2): # schon in MDH
;
sin_xi2_s:=winkel3(5,2):
cos_xi2_s:=winkel3(6,2):

cos_xi34_s:=cos_xi3_s*cos(phi413)+sin_xi3_s*sin(phi413): # am 20.11 der neuen fourbar1 angepasst
;
sin_xi34_s:=-cos_xi3_s*sin(phi413)+sin_xi3_s*cos(phi413): # am 20.11 der neuen fourbar1 angepasst
;
cos_xi413_s:=-cos_xi2_s:
sin_xi413_s:=sin_xi2_s:
# Kintmp_subsexp
kintmp_subsexp(1,2):=sin_xi34_s:
kintmp_subsexp(2,2):=cos_xi34_s:
kintmp_subsexp(3,2):=sin_eta16_s:
kintmp_subsexp(4,2):=cos_eta16_s:
kintmp_subsexp(5,2):=sin_eta27_s:
kintmp_subsexp(6,2):=cos_eta27_s:
kintmp_subsexp(7,2):=sin_rho28_s:
kintmp_subsexp(8,2):=cos_rho28_s:
kintmp_subsexp(9,2):=sin_rho89_s:
kintmp_subsexp(10,2):=cos_rho89_s:
kintmp_subsexp(11,2):=sin_xi1016_s:
kintmp_subsexp(12,2):=cos_xi1016_s:
kintmp_subsexp(13,2):=sin_eta711_s:
kintmp_subsexp(14,2):=cos_eta711_s:
kintmp_subsexp(15,2):=sin_rho312_s:
kintmp_subsexp(16,2):=cos_rho312_s:
kintmp_subsexp(17,2):=sin_xi413_s:
kintmp_subsexp(18,2):=cos_xi413_s:
kintmp_subsexp(11..18,1):
# Kintmp_qs
kintmp_qs(1,1):=%arctan(sin_xi34_s,cos_xi34_s):
kintmp_qs(2,1):=%arctan(sin_eta16_s,cos_eta16_s):
kintmp_qs(3,1):=%arctan(sin_eta27_s,cos_eta27_s):
kintmp_qs(4,1):=%arctan(sin_rho28_s,cos_rho28_s):
kintmp_qs(5,1):=%arctan(sin_rho89_s,cos_rho89_s):
kintmp_qs(6,1):=%arctan(sin_xi1016_s,cos_xi1016_s):
kintmp_qs(7,1):=%arctan(sin_eta711_s,cos_eta711_s):
kintmp_qs(8,1):=%arctan(sin_rho312_s,cos_rho312_s):
kintmp_qs(9,1):=%arctan(sin_xi413_s,cos_xi413_s):
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
kc_symbols := Matrix(list_constant_expressions( kintmp_subsexp ));
#kc_symbols :=Transpose(kc_symbols);
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);
#printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
