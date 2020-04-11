
# Slider Crank / Schubkurbel
# Initialisierung
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
read "../helper/proc_intersect_circle":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Variable mit Winkeln der Nebenstruktur nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
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
Schalter_Opt := <1;1;0;0;0;0;0;1>:
printf("Beginn der Berechnungen. %s\n", FormatTime("%Y-%m-%d %H:%M:%S")):
st := time():
# Zwangsbedingung mit Kreisschnittpunkt G zur Bestimmung von rho
r_p1:= < 0 , 0 >: # Punkt K
;
r_p2:= < qJ1s+HP , 0> : # Punkt P
;
r1:=GK:
r2:=GP:
A:=intersect_circle(r_p1, r_p2, r1, r2):
r_s1x:=simplify(A(1,1),trig):
r_s1y:=simplify(A(2,1),trig):
r_agx:=GK*cos(rho_s): # aus Homogener Trafo
;
r_agy:=-GK*sin(rho_s): # aus Homogener Trafo
;
f1:=r_s1x-r_agx=0:# GL1
;
f2:=r_s1y-r_agy=0: # GL2
;
# Zwangsbedingungen für eta , zweiter Weg zum Kreisschnittpunkt G
r_agx2:=qJ1s+HP+GP*cos(xi_s):
r_agy2:=-GP*sin(xi_s):
f11:=r_s1x-r_agx2=0: # GL3
;
f12:=r_s1y-r_agy2=0: # GL4
;
# Bestimme Winkel
sol:=simplify(solve({f1,f2},[cos(rho_s),sin(rho_s)]),trig):
sol:=op(sol):
cos(rho_s):=simplify(rhs(sol[1]),trig):
sin(rho_s):=simplify(rhs(sol[2]),trig):
# xi bestimmen
sol1:=solve({f11,f12},{cos(xi_s),sin(xi_s)}): 
assign(sol1):
cos(xi_s):=simplify(rhs(sol1[1]),trig): 
sin(xi_s):=simplify(rhs(sol1[2]),trig):
# eta bestimmen
cos(eta_s):=-(cos(rho_s)*cos(xi_s)+sin(rho_s)*sin(xi_s)):
sin(eta_s):=-(sin(rho_s)*cos(xi_s)-cos(rho_s)*sin(xi_s)):
# Kintmp_subsexp
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

# Export und Save ( falls als Arbeitsblatt fuer grosse Roboter)
winkel_fourbarpris := <<sin_rho|sin(rho_s)>; <cos_rho|cos(rho_s)>; <sin_eta | sin(eta_s)>; <cos_eta | cos(eta_s)>;<sin_xi | sin(xi_s)>; <cos_xi | cos(xi_s)>>:
#mkdir(sprintf("../codeexport/slcr/")): falls Ordner noch nicht vorhanden
;
#mkdir(sprintf("../codeexport/slcr/tmp/")): falls Ordner noch nicht vorhanden
;
save winkel_fourbarpris, sprintf("../codeexport/%s/tmp/%s_angles_trig_elim", robot_name, robot_name):
# Exportiere Code für folgende Skripte
kintmp_qt := convert_s_t(kintmp_qs):
#kintmp_subsexp:=convert_s_t(kintmp_subsexp):
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
#kc_symbols :=Transpose(kc_symbols);
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m",robot_name),2);
printf("Fertig. %s. CPU-Zeit bis hier: %1.2fs.\n", FormatTime("%Y-%m-%d %H:%M:%S"), time()-st):


