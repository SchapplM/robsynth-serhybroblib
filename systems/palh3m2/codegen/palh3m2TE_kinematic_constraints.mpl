
# Berechne kinematische Zwangsbedingungen f��r den Palettierroboter KUKA KR700PA aus 2 Parallelogram
# Einleitung
# Die kinematischen Zwangsbedingungen werden als Ersatzungsausdruck f��r die abh�ngigen Winkel aufgestellt.
# 
# palh3m2TE->KUKA KR700PA, modellierung der Zwangsbedingungen mit Ausdr��cken f��r trigonometrische Elimination
# kinematic_constraint->Kinematische Zwangsbedingungen
# Datenblatt des Roboters:
# unter :   https://www.kuka.com/en-de/products/robot-systems/industrial-robots/kr-700-pa
# Voraussetzung(sehr wichtig!!!!!): 
# Viergelenkkette muss mit der Toolbox berechnet worden sein (Arbeitsblatt "fourbar2TE_kinematic_constrains.mv")
# Quelle:
# TODO
# Autor:
# Weipu Shan
# Studienarbeit bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
# (C) Institut fuer mechatronische Systeme, Leibniz Universitaet Hannover
# Initialisierung
# Import (Library)
interface(warnlevel=0): 		# Unterdr��cke die folgende Warnung.
restart: 				# Gibt eine Warnung, wenn ��ber Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true:               # f��r Speicherung
;
with(StringTools):                                  # f��r Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
#with(ListTools):
codegen_act := true:
codegen_opt := 1:      # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0:    # Zur Code-Generierung auch f��r Nicht-Inert Ausdr��cke
;
# Import(hybriddyn)
read "../helper/proc_MatlabExport":    #Matlab Export
read "../transformation/proc_rotx":      #Rotation um X
read "../transformation/proc_roty":      #Rotation um Y
read "../transformation/proc_rotz":      #Rotation um Z
read "../helper/proc_convert_s_t":        #zeitabh�ngige Variable f��r Ableitung
read "../helper/proc_convert_t_s":        #konstante Variable
read "../robot_codegen_constraints/proc_subs_kintmp_exp": #substitute Spalte 1 mit Spalte 2
with (RealDomain):            # Schr�nkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":                                                      #aktuelle Roboter, MDH-Tabelle
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions",robot_name):        
# Variable Initialisierung
# Variable mit Winkeln der Nebenstruktur nur in Abh�ngigkeit der verallgemeinerten Koordinaten  (Nur Aktive Gelenke)
kintmp_qs := Matrix(RowDimension(kintmp_s),1):       	#kintmp_s von aktuellen Roboter, MDH
qJ_t := <qJ1(t),qJ2(t),qJ3(t),qJ4(t)>:                                	#zeitabh�ngig Variable
qJ_s := <qJ1s,qJ2s,qJ3s,qJ4s>:                                          #konstante Variable
;
# 
# Konstante Winkel bereits hineinschreiben (Passive Gelenk)
for i from 1 to RowDimension(kintmp_s) do
	if diff(kintmp_s(i,1),t) = 0 then
		#konstante Variable, passive Winkel von aktuellem Roboter, MDH
		kimtmp_qs(i,1) := kintmp_s(i,1):  
  	end if:
end do
;
# 
# Variable definieren f��r die Hilfswinkel
# Ersetzungsausdr��cke definieren
# Speichere Sinus und Cosinus der Winkel direkt ab, da diese in den Rotationsmatrizen direkt auftreten.
# Spalte 1: Zusuchender Ausdruck (sin oder cos eines Winkels)	#��һ����Ҫ�ҵ������Ǻ������
# Spalte 2: Einzusetzender Ausdruck			#�ڶ����Ǵ���ı��
# Dadurch werden arctan-Ausdr��cke in der direkten Kinematik reduziert.
# �hnliches Vorgehen wie in [1].
# 
# #�½�kintmp_subsexp 12X2, ��һ�б�������sin cos ��ʽ���ڶ��г�ʼ���͵�һ����ͬ����ǰ�����˵�sin cos
kintmp_subsexp := Matrix(2*RowDimension(kintmp_s),2):
for i from 1 to RowDimension(kintmp_s) do 			#kintmp_s von aktuellen Roboter
	kintmp_subsexp(2*i-1,1) := sin(kintmp_s(i,1)):
	kintmp_subsexp(2*i,  1) := cos(kintmp_s(i,1)):
	# Initialisierung der rechten Spalte mit gleichen Werten. Sp�ter nur Ersetzung, wenn Vorteilhaft.
	kintmp_subsexp(2*i-1,2) := kintmp_subsexp(2*i-1, 1):
	kintmp_subsexp(2*i,  2) := kintmp_subsexp(2*i,   1):
end do:
# 
# Sicherung der aktuelle Daten
backup_kintmp_qs := kintmp_qs:		#Backup von aktuellem Roboter, konstante Variable
backup_kintmp_qt := kintmp_qt:		#Backup von aktuellem Roboter, zeitabh�ngige Variable
backup_kintmp_subsexp := kintmp_subsexp:   	#Backup von aktuellem Roboter, sin, cos von 4 Winkeln
 
;
# Daten von Fourbar ��bernehmen
read sprintf("../codeexport/fourbar2TE/tmp/kinematic_constraints_maple_inert.m"): 	
# Winkeln von Fourbar���ȸ���ԭ�ȵı��� kintmp_subsexp, kintmp_qs, kintmp_qt,����һ��ʼ��Ҫ����

;
# Speichern WInkeln von Fourbar, mit Fourbar Variable
kintmp_subsexp_fourbar := kintmp_subsexp:  	#sin,cos von 4 Winkeln, Fourbar
kintmp_qs_fourbar := kintmp_qs:		#konstante Variable, Fourbar
kintmp_qt_fourbar := kintmp_qt:		#zeitabh�ngige Variable, Fourbar
;
# 
# �������������� Was ist der Sinn darunten
#kintmp_subsexp_fourbar:
#winkel(6,2):
#kintmp_subsexp_fourbar(6,2):
#kintmp_subsexp_fourbar(3,3):
#winkel_neu := Matrix(6,2):
winkel_neu(1,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(3,2)):
winkel_neu(2,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(4,2)):
winkel_neu(3,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(1,2)):
winkel_neu(4,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(2,2)):
winkel_neu(5,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(5,2)):
winkel_neu(6,2) := subs({qJ1s=phi_s}, kintmp_subsexp_fourbar(6,2)):
# 
# Vorherige Werte f��r dieses System wieder herstellen
kintmp_qs := backup_kintmp_qs:		#aktuellem Roboter, konstante Variable
kintmp_qt := backup_kintmp_qt:		#aktuellem Roboter, zeitabh�ngige Variable
kintmp_subsexp := backup_kintmp_subsexp:		#aktuellem Roboter, sin, cos von 4 Winkeln

;
#winkel_neu(1..6,1):
# Viergelenkkette 1: A-B-C-D
# Von ABCD: eta1, eta3, eta4  -----> Von Fourbar: rho, eta, xi
winkel1 := <<sin_eta1|sin(eta1_s)>;<cos_eta1|cos(eta1_s)>;<sin_eta3|sin(eta3_s)>;<cos_eta3|cos(eta3_s)>;
<sin_eta4|sin(eta4_s)>;<cos_eta4|cos(eta4_s)>>:
# Ausdr��cke von Fourbar in aktuellem Roboter mitbringen
for i from 1 to 6 do
	winkel1(i,2) := kintmp_subsexp_fourbar(i,2):       		#Subsexpress von Fourbar ��bernehmen
	cos_phi_s := sin(qJ2s)*sin(phi2) - cos(qJ2s)*cos(phi2):	#phi_s ist unabh�ngige Winkel von VGK, hier ist phi_s = eta2_s
	sin_phi_s := sin(qJ2s)*cos(phi2) + cos(qJ2s)*sin(phi2):	#qJ2s ist echte aktive Winkel von Roboter, eta2 ^= phi = Pi - qJ2 - phi2
	winkel1(i,2) := subs({sin(qJ1s)=sin_phi_s},winkel1(i,2)):	#����qJ1s ��fourbar�����ᣬ������ı��ʽ�����棬����������ﵱǰroboter
	winkel1(i,2) := subs({cos(qJ1s)=cos_phi_s},winkel1(i,2)):
	winkel1(i,2) := subs({l1=AB}, winkel1(i,2)):
	winkel1(i,2) := subs({l2=DA}, winkel1(i,2)):
	winkel1(i,2) := subs({l3=DC}, winkel1(i,2)):
	winkel1(i,2) := subs({l4=BC}, winkel1(i,2)):
end do:
winkel1(1..6,1):
# eta1, eta3, eta4 wird wie rho, eta, xi in Fourbar darstellen
sin_eta1_s := winkel1(1,2):	#eta1 = qJ2s = rho
;
cos_eta1_s := winkel1(2,2):	
sin_eta3_s := winkel1(3,2):	#eta3 = qJ3s = eta
;
cos_eta3_s := winkel1(4,2):	
sin_eta4_s := winkel1(5,2):	#eta4 = qJ4s = xi
;
cos_eta4_s := winkel1(6,2):	
# eta16, eta79, eta27, abh�ngige Winkeln von aktuellem Roboter
cos_eta16_s := cos_eta1_s*cos(phi1) + sin_eta1_s*sin(phi1):	#eta16 = eta1 - phi1
;
sin_eta16_s := sin_eta1_s*cos(phi1) - cos_eta1_s*sin(phi1):
cos_eta79_s := -cos_eta4_s:				#eta79 = Pi-eta4
;
sin_eta79_s := sin_eta4_s:
cos_eta27_s := cos_eta3_s*cos(phi79) + sin_eta3_s*sin(phi79):	#eta27 = eta3 - phi79
;
sin_eta27_s := sin_eta3_s*cos(phi79) - cos_eta3_s*sin(phi79):
# Viergelenkkette 2: B-E-F-G
winkel2 := <<sin_xi1|sin(xi1_s)>;<cos_xi1|cos(xi1_s)>;<sin_xi3|sin(xi3_s)>;<cos_xi3|cos(xi3_s)>;
<sin_xi4|sin(xi4_s)>;<cos_xi4|cos(xi4_s)>>:
# Ausdr��cke von Fourbar in aktuellem Roboter mitbringen
for i from 1 to 6 do
	winkel2(i,2) := kintmp_subsexp_fourbar(i,2):       			#Subsexpress von Fourbar ��bernehmen
	#phi_s ist unabh�ngige Winkel von VGK, hier ist phi_s = xi2_s
	#qJ3s ist echte aktive Winkel von Roboter, xi2 ^= phi = Pi - qJ3 - (phi78-eta27)
	cos_phi_s := -cos(phi78+phi79)*(cos_eta3_s*cos(qJ3s)+sin_eta3_s*sin(qJ3s))-sin(phi78+phi79)*(sin_eta3_s*cos(qJ3s)-cos_eta3_s*sin(qJ3s)):
	sin_phi_s := sin(phi78+phi79)*(cos_eta3_s*cos(qJ3s)+sin_eta3_s*sin(qJ3s))-cos(phi78+phi79)*(sin_eta3_s*cos(qJ3s)-cos_eta3_s*sin(qJ3s)):
	#����qJ1s ��fourbar�����ᣬ������ı��ʽ�����棬����������ﵱǰroboter
	winkel2(i,2) := subs({sin(qJ1s)=sin_phi_s},winkel2(i,2)):	
	winkel2(i,2) := subs({cos(qJ1s)=cos_phi_s},winkel2(i,2)):
	winkel2(i,2) := subs({l1=BG}, winkel2(i,2)):
	winkel2(i,2) := subs({l2=BE}, winkel2(i,2)):
	winkel2(i,2) := subs({l3=EP}, winkel2(i,2)):
	winkel2(i,2) := subs({l4=GP}, winkel2(i,2)):
end do:
winkel2(1..6,1):
# xi1, xi3, xi4 wird wie rho, eta, xi in Fourbar darstellen
sin_xi1_s := winkel2(1,2):	#xi1 = qJ2s = rho
;
cos_xi1_s := winkel2(2,2):
sin_xi3_s := winkel2(3,2):	#xi3 = qJ3s = eta
;
cos_xi3_s := winkel2(4,2):
sin_xi4_s := winkel2(5,2):	#xi4 = qJ4s = xi
;
cos_xi4_s := winkel2(6,2):
# xi78, xi34, xi410, abh�ngige Winkeln von aktuellem Roboter
cos_xi78_s := cos_xi1_s:				#xi78 = xi1
;
sin_xi78_s := sin_xi1_s:
cos_xi34_s := cos_xi3_s*cos(phi410) + sin_xi3_s*sin(phi410):	#xi34 = xi3 - phi410
;
sin_xi34_s := sin_xi3_s*cos(phi410) - cos_xi3_s*sin(phi410):
cos_xi410_s := -cos_xi4_s:				#xi410 = Pi - xi4
;
sin_xi410_s := sin_xi4_s:
# Kintmp_subsexp
kintmp_subsexp(1, 2) := sin_xi34_s:
kintmp_subsexp(2, 2) := cos_xi34_s:
kintmp_subsexp(3, 2) := sin_eta16_s:
kintmp_subsexp(4, 2) := cos_eta16_s:
kintmp_subsexp(5, 2) := sin_eta27_s:
kintmp_subsexp(6, 2) := cos_eta27_s:
kintmp_subsexp(7, 2) := sin_xi78_s:
kintmp_subsexp(8, 2) := cos_xi78_s:
kintmp_subsexp(9, 2) := sin_eta79_s:
kintmp_subsexp(10, 2) := cos_eta79_s:
kintmp_subsexp(11, 2) := sin_xi410_s:
kintmp_subsexp(12, 2) := cos_xi410_s:
kintmp_subsexp(1..12,1):

# Kintmp_qs
kintmp_qs(1,1) := %arctan(sin_xi34_s,cos_xi34_s):
kintmp_qs(2,1) := %arctan(sin_eta16_s,cos_eta16_s):
kintmp_qs(3,1) := %arctan(sin_eta27_s,cos_eta27_s):
kintmp_qs(4,1) := %arctan(sin_xi78_s,cos_xi78_s):
kintmp_qs(5,1) := %arctan(sin_eta79_s,cos_eta79_s):
kintmp_qs(6,1) := %arctan(sin_xi410_s,cos_xi410_s):
# Export
kintmp_qt := convert_s_t(kintmp_qs):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple", robot_name):
save kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_subsexp_maple.m", robot_name):
#printf("Ausdr��cke f��r kintmp_subsexp gespeichert (Maple). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
for i from 1 to RowDimension(kintmp_s) do
  tmp := kintmp_qs(i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d",robot_name, i):
  save tmp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert_kintmpq_%d.m", robot_name, i):
end do:
save kin_constraints_exist, kintmp_qs, kintmp_qt,kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert" ,robot_name):
save kin_constraints_exist, kintmp_qs, kintmp_qt, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
save kintmp_qs, sprintf("../codeexport/%s/tmp/kinematic_constraints_kintmp_qs_maple_inert", robot_name):
#printf("Ausdr��cke mit Inert-Arctan exportiert (Matlab). %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):
# Liste mit abh�ngigen konstanten Kinematikparametern erstellen (wichtig f��r Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
interface(rtablesize=20):
kc_symbols := Matrix(list_constant_expressions( kintmp_subsexp ));
#kc_symbols :=Transpose(kc_symbols);
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
kc_symbols;
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);
#printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs. \n", FormatTime("%Y-%m-%d %H:%M:%S"), time() - st):

