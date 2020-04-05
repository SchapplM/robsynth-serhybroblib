
# Berechne kinematische Zwangsbedingungen für den TSR2040 aus 2 Parallelogrammen und einer Fünfgelenkkette
# Einleitung
# 
# Herleitung der Kinematischen Zwangsbedingungen für einen 2-FG hybriden Palettier-/Verpackungsroboter mit zwei Parallelogrammen und einer Fünfgelenkkette
# Roboter: Triowin (Hersteller), "TDR Series 2-Axis Parallel Robot"
# Link: http://en.triowin.com/tdr-series2-axisparallelrobot-15268685491267769.html
# 
# Die kinematischen Zwangsbedingungen werden als Ersetzungsausdruck für die abhängigen Winkel aufgestellt.
# picker2Dm2TE -> modellierung der Zwangsbedingungen mit Ausdrücken für trigonometrische Elimination

# Voraussetzung
# Vier- und Fünf gelenkkette muss mit der Toolbox berechnet worden sein (Arbeitsblätter "fourbar1TE_kinematic_constraints.mw" und fivebar1TE_kinematic_constraints.mw)
# 
# Autor
# Abderahman Bejaoui
#  Hiwi bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-01
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
# Aufruf der Viergelenkkette (als Parallelogramm)
read sprintf("../codeexport/fourbar2TE/tmp/kinematic_constraints_maple_inert.m"):
kintmp_subsexp_fourbar:= kintmp_subsexp:
kintmp_qs_fourbar := kintmp_qs:
kintmp_qt_fourbar := kintmp_qt:
# Aufrud der Fünfgelenkkette
read sprintf("../codeexport/fivebar1TE/tmp/kinematic_constraints_maple_inert.m"):
kintmp_subsexp_fivebar := kintmp_subsexp:
kintmp_qs_fivebar := kintmp_qs:
kintmp_qt_fivebar := kintmp_qt:
# Vorherige Werte für dieses System wieder herstellen
kintmp_qs := backup_kintmp_qs:
kintmp_qt := backup_kintmp_qt:
kintmp_subsexp := backup_kintmp_subsexp:
# Fünfgelenkkettte A-B-C-E-D
winkel1:=<<sin_eta3|sin(eta3_s)>;<cos_eta3|cos(eta3_s)>;<sin_eta12|sin(eta12_s)>;<cos_eta12|cos(eta12_s)>;<sin_eta24|sin(eta24_s)>;<cos_eta24|cos(eta24_s)>>:

for i from 1 to 6 do
  # Daten der allgemeinen Fünfgelenkkette laden
  winkel1(i,2):=kintmp_subsexp_fivebar(i,2):  
  # Kinematikparameter der FGK an die des Roboters anpassen
  winkel1(i,2) := subs({BC=L1}, winkel1(i,2)):
  winkel1(i,2) := subs({AB=e}, winkel1(i,2)):
  winkel1(i,2) := subs({CD=L3}, winkel1(i,2)):
  winkel1(i,2) := subs({ED=L4}, winkel1(i,2)):
  winkel1(i,2) := subs({AE=L3}, winkel1(i,2)):
  # variable qJ bezogen auf FGK durch allgemeine Variable ersetzen
  winkel1(i,2) := subs({sin(qJ1s)= sin_phi_s},winkel1(i,2)): 
  winkel1(i,2) := subs({cos(qJ1s)= cos_phi_s},winkel1(i,2)):
  winkel1(i,2) := subs({sin(qJ2s)= sin_psi_s},winkel1(i,2)): 
  winkel1(i,2) := subs({cos(qJ2s)= cos_psi_s},winkel1(i,2)):	
  # Allgemeine Variable durch qJ bezogen auf Roboter ersetzen
  winkel1(i,2) := subs({sin_phi_s = cos(qJ2s)},winkel1(i,2)): 
  winkel1(i,2) := subs({cos_phi_s = -sin(qJ2s)},winkel1(i,2)):
  winkel1(i,2) := subs({sin_psi_s = sin(qJ1s)},winkel1(i,2)): 
  winkel1(i,2) := subs({cos_psi_s = cos(qJ1s)},winkel1(i,2)):	
end do:
# Umrechnung aus Zwangsbedingugen für Winkel aus MDH-Tabelle
sin_eta3_s:=winkel1(1,2):
cos_eta3_s:=winkel1(2,2):  
sin_eta12_s:=winkel1(3,2):  
cos_eta12_s:=winkel1(4,2):  
sin_eta24_s:=winkel1(5,2):    
cos_eta24_s:=winkel1(6,2):  
sin_eta410_s:=sin_eta3_s:
cos_eta410_s:=-cos_eta3_s:

# Viergelenkkette B-H-L-C
winkel2:=<<sin_rho18|sin(rho18_s)>;<cos_rho18|cos(rho18_s)>;<sin_rho05|sin(rho05_s)>;<cos_rho05|cos(rho05_s)>;<sin_rho811|sin(rho811_s)>;<cos_rho811|cos(rho811_s)>>:
for i from 1 to 6 do
  # Daten des Parallelogramms laden
  winkel2(i,2) := kintmp_subsexp_fourbar(i,2):  
  # Kinematikparameter des Parallelogramms an die des Roboters anpassen
  winkel2(i,2) := subs({l1=L5}, winkel2(i,2)):
  winkel2(i,2) := subs({l2=L1}, winkel2(i,2)):
  # variable qJ bezogen auf Parallelogramm durch allgemeine Variable ersetzen
  winkel2(i,2) := subs({sin(qJ1s)= sin_phi_s},winkel2(i,2)): 
  winkel2(i,2) := subs({cos(qJ1s)= cos_phi_s},winkel2(i,2)):
  # Allgemeine Variable durch qJ bezogen auf Roboter ersetzen
  winkel2(i,2) := subs({cos_phi_s = -cos(phi05)*cos(qJ1s)-sin(phi05)*sin(qJ1s)},winkel2(i,2)): 
  winkel2(i,2) := subs({sin_phi_s = -cos(phi05)*sin(qJ1s)+sin(phi05)*cos(qJ1s)},winkel2(i,2)):	
end do:
sin_rho18_s:=winkel2(1,2):
cos_rho18_s:=winkel2(2,2):
sin_rho05_s:=winkel2(3,2):
cos_rho05_s:=winkel2(4,2):
sin_rho811_s:=winkel2(5,2): 
cos_rho811_s:=winkel2(6,2): 
# Viergelenkkette C-G-F-P
winkel3:=<<sin_xi23|sin(xi23_s)>;<cos_xi23|cos(xi23_s)>;<sin_xi612|sin(xi612_s)>;<cos_xi612|cos(xi612_s)>;<sin_xi39|sin(xi39_s)>;<cos_xi39|cos(xi39_s)>>:
cos_xi26_s:= - ( -cos(phi1)*(cos_rho18_s*cos_eta12_s+sin_rho18_s*sin_eta12_s)-sin(phi1)*(sin_eta12_s*cos_rho18_s-cos_eta12_s*sin_rho18_s)): 
sin_xi26_s:= -cos(phi1)*(sin_eta12_s*cos_rho18_s-cos_eta12_s*sin_rho18_s)+sin(phi1)*(cos_eta12_s*cos_rho18_s+sin_eta12_s*sin_rho18_s):
for i from 1 to 6 do
  # Daten des Parallelogramms laden
  winkel3(i,2) := kintmp_subsexp_fourbar(i,2):  
  # Kinematikparameter des Parallelogramms an die des Roboters anpassen
  winkel3(i,2) := subs({l1=L6}, winkel3(i,2)):
  winkel3(i,2) := subs({l2=L2}, winkel3(i,2)):
  # variable qJ bezogen auf Parallelogramm durch allgemeine Variable ersetzen
  winkel3(i,2) := subs({sin(qJ1s) = sin_phi_s},winkel3(i,2)): 
  winkel3(i,2) := subs({cos(qJ1s) = cos_phi_s},winkel3(i,2)):
  # Allgemeine Variable durch qJ bezogen auf Roboter ersetzen
  winkel3(i,2) := subs({cos_phi_s = -cos(phi1)*(cos_rho18_s*cos_eta12_s+sin_rho18_s*sin_eta12_s)-sin(phi1)*(sin_eta12_s*cos_rho18_s-cos_eta12_s*sin_rho18_s)},winkel3(i,2)): 
  winkel3(i,2) := subs({sin_phi_s = -cos(phi1)*(sin_eta12_s*cos_rho18_s-cos_eta12_s*sin_rho18_s)+sin(phi1)*(cos_eta12_s*cos_rho18_s+sin_eta12_s*sin_rho18_s)},winkel3(i,2)):
end do:
winkel3(1..6,1):
sin_xi23_s:=winkel3(1,2):
cos_xi23_s:=winkel3(2,2):
sin_xi612_s:=winkel3(3,2): # schon in MDH
;
cos_xi612_s:=winkel3(4,2): # schon in MDH
;
sin_xi39_s:=winkel3(5,2):
cos_xi39_s:=winkel3(6,2):
# Kintmp_subsexp
kintmp_subsexp(1,2) := sin_eta12_s:
kintmp_subsexp(2,2) := cos_eta12_s:
kintmp_subsexp(3,2) := sin_xi23_s:
kintmp_subsexp(4,2) := cos_xi23_s:
kintmp_subsexp(5,2) := sin_eta24_s:
kintmp_subsexp(6,2) := cos_eta24_s:
kintmp_subsexp(7,2) := sin_rho05_s:
kintmp_subsexp(8,2) := cos_rho05_s:
kintmp_subsexp(9,2) := sin_xi26_s:
kintmp_subsexp(10,2) := cos_xi26_s:
kintmp_subsexp(11,2) := sin_rho18_s:
kintmp_subsexp(12,2) := cos_rho18_s:
kintmp_subsexp(13,2) := sin_xi39_s:
kintmp_subsexp(14,2) := cos_xi39_s:
kintmp_subsexp(15,2) := sin_eta410_s:
kintmp_subsexp(16,2) := cos_eta410_s:
kintmp_subsexp(17,2) := sin_rho811_s:
kintmp_subsexp(18,2) := cos_rho811_s:
kintmp_subsexp(19,2) := sin_xi612_s:
kintmp_subsexp(20,2) := cos_xi612_s:
# Kintmp_qs
kintmp_qs(1,1):=%arctan(sin_eta12_s,cos_eta12_s):
kintmp_qs(2,1):=%arctan(sin_xi23_s,cos_xi23_s):
kintmp_qs(3,1):=%arctan(sin_eta24_s,cos_eta24_s):
kintmp_qs(4,1):=%arctan(sin_rho05_s,cos_rho05_s):
kintmp_qs(5,1):=%arctan(sin_xi26_s,cos_xi26_s):
kintmp_qs(6,1):=%arctan(sin_rho18_s,cos_rho18_s):
kintmp_qs(7,1):=%arctan(sin_xi39_s,cos_xi39_s):
kintmp_qs(8,1):=%arctan(sin_eta410_s,cos_eta410_s):
kintmp_qs(9,1):=%arctan(sin_rho811_s,cos_rho811_s):
kintmp_qs(10,1):=%arctan(sin_xi612_s,cos_xi612_s):
# Export
kintmp_qt := convert_s_t(kintmp_qs):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
for i from 1 to RowDimension(kintmp_s) do
  tmp := kintmp_qs(i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d",robot_name, i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d.m", robot_name, i):
end do:
save kin_constraints_exist, kintmp_qs, kintmp_qt,kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert" ,robot_name):
save kin_constraints_exist, kintmp_qs, kintmp_qt, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_subsexp ));
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);

