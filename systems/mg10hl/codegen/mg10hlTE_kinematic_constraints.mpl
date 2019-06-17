
# MG10HL von Kawasaki
# Voraussetzung:
# Slider_crank muss mit der Toolbox berechnet worden sein (Arbeitsblatt "slcr_kinematic_constraints.mw")
# Viergelenkkette muss mit der Toolbox berechnet worden sein (Arbeitsblatt "vgk_kinematic_constraints.mw")
# Autor
# Abderahman Bejaoui
# Studienarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-09
# (C) Institut fuer mechatronische Systeme, Leibniz Universitaet Hannover
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
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
qJ_t:= <qJ1(t),qJ2(t),qJ3(t),qJ4(t),qJ5(t),qJ6(t)>:
qJ_s:= <qJ1s,qJ2s,qJ3s,qJ4s,qJ5,qJ6s>:
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
read sprintf("../codeexport/fourbar1TE/tmp/fourbar1TE_angles_trig_elim"):
read sprintf("../codeexport/fourbarprisTE/tmp/fourbarprisTE_angles_trig_elim"):

# Viergelenkkette A-C-D-E
winkel1:=<<sin_rho910|sin(rho910_s)>;<cos_rho910|cos(rho910_s)>;<sin_rho19|sin(rho19_s)>;<cos_rho19|cos(rho19_s)>;<sin_rho3|sin(rho3_s)>;<cos_rho3|cos(rho3_s)>>:
for i from 1 to 6 do
	winkel1(i,2):=winkel(i,2):  
     cos_phi_s :=sin(qJ2s)*cos(phi23)+cos(qJ2s)*sin(phi23):
	sin_phi_s := cos(qJ2s)*cos(phi23)-sin(qJ2s)*sin(phi23):
	winkel1(i,2):=subs({sin(phi_s)=sin_phi_s},winkel1(i,2)): 
	winkel1(i,2):=subs({cos(phi_s)=cos_phi_s},winkel1(i,2)):
     winkel1(i,2) := subs({l1=AC}, winkel1(i,2)):
     winkel1(i,2) := subs({l2=AE}, winkel1(i,2)):
     winkel1(i,2) := subs({l3=ED}, winkel1(i,2)):
     winkel1(i,2) := subs({l4=DC}, winkel1(i,2)):
end do:
winkel1(1..6,1):
sin_rho910_s:=winkel1(1,2): # für kintmp_subsexp
;
cos_rho910_s:=winkel1(2,2): # für kintmp_subsexp
;
sin_rho19_s:=winkel1(3,2): # für kintmp_subsexp
;
cos_rho19_s:=winkel1(4,2):  # für kintmp_subsexp
;
sin_rho3_s:=winkel1(5,2):
cos_rho3_s:=winkel1(6,2):
cos_rho23_s:=-cos_rho3_s: # für kintmp_subsexp
;
sin_rho23_s:=sin_rho3_s:  # für kintmp_subsexp
;
# Schubkurbel G-P-K
winkel2:=<<sin_eta1112|sin(eta1112_s)>;<cos_eta1112|cos(eta1112_s)>;<sin_eta3|sin(eta3_s)>;<cos_eta3|cos(eta3_s)>;<sin_eta45|sin(eta45_s)>;<cos_eta45|cos(eta45_s)>>:
for i from 1 to 6 do
	winkel2(i,2):=winkel_fourbarpris(i,2):  
	winkel2(i,2):=subs({qJ1s=qJ6s},winkel2(i,2)): 
	winkel2(i,2) := subs({GK=GK}, winkel2(i,2)):
     winkel2(i,2) := subs({GP=GP}, winkel2(i,2)):
     winkel2(i,2) := subs({HP=HP}, winkel2(i,2)):
end do:
winkel2
;
sin_eta1112_s:=winkel2(1,2): # für kintmp_subsexp
;
cos_eta1112_s:=winkel2(2,2): # für kintmp_subsexp
;
sin_eta3_s:=winkel2(3,2): 
cos_eta3_s:=winkel2(4,2):
sin_eta45_s:=winkel2(5,2): # für kintmp_subsexp
;
cos_eta45_s:=winkel2(6,2): # für kintmp_subsexp
;
sin_eta411_s:=sin_eta3_s: # für kintmp_subsexp
;
cos_eta411_s:=-cos_eta3_s: # für kintmp_subsexp
;
cos_eta34_s:=-(cos_eta411_s*cos(phi3)-sin_eta411_s*sin(phi3)):  # für kintmp_subsexp
;
sin_eta34_s:=sin_eta411_s*cos(phi3)+cos_eta411_s*sin(phi3):  # für kintmp_subsexp
;
# Kintmp_subsexp
kintmp_subsexp(1,2):=sin_rho23_s:
kintmp_subsexp(2,2):=cos_rho23_s:
kintmp_subsexp(3,2):=sin_eta34_s:
kintmp_subsexp(4,2):=cos_eta34_s:
kintmp_subsexp(5,2):=sin_eta45_s:
kintmp_subsexp(6,2):=cos_eta45_s:
kintmp_subsexp(7,2):=sin_rho19_s:
kintmp_subsexp(8,2):=cos_rho19_s:
kintmp_subsexp(9,2):=sin_rho910_s:
kintmp_subsexp(10,2):=cos_rho910_s:
kintmp_subsexp(11,2):=sin_eta411_s:
kintmp_subsexp(12,2):=cos_eta411_s:
kintmp_subsexp(13,2):=sin_eta1112_s:
kintmp_subsexp(14,2):=cos_eta1112_s:
kintmp_subsexp(11..14,1);

# Kintmp_qs
kintmp_qs(1,1):=arctan(sin_rho23_s,cos_rho23_s):
kintmp_qs(2,1):=arctan(sin_eta34_s,cos_eta34_s):
kintmp_qs(3,1):=arctan(sin_eta45_s,cos_eta45_s):
kintmp_qs(4,1):=arctan(sin_rho19_s,cos_rho19_s):
kintmp_qs(5,1):=arctan(sin_rho910_s,cos_rho910_s):
kintmp_qs(6,1):=arctan(sin_eta411_s,cos_eta411_s):
kintmp_qs(7,1):=arctan(sin_eta1112_s,cos_eta1112_s):
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

kc_symbols(11)
;

