
# Implizite kinematische Zwangsbedingungen für 2D-Picker

# 
# Erklärung des Dateinamens
# picker2D: 2D-Roboter für Pick&Place
# m1: Modellierung der geschlossenen Ketten als allgemeine Viergelenkkette
# IC: Modellierung des Roboters mit impliziten Zwangsbedingungen
# kinematic_contraints_implizit: Berechnung der kinematischen Zwangsbedingungen (implizite Form)
# 
# Quelle
# Aufzeichnungen Schappler, 6.5.2019
# Ursprüngliche Notizen von Abderahman Bejaoui (Hiwi bei Moritz Schappler)
# 
# Autor
# Abderahman Bejaoui
# Hiwi bei: Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
# (C) Institut fuer mechatronische Systeme, Leibniz Universitaet Hannover
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
codegen_opt := 1: # Geringerer Optimierungsgrad. Sonst zu lange.
codegen_debug := 0: # Zur Code-Generierung auch für Nicht-Inert-Ausdrücke
;
read "../helper/proc_MatlabExport":
read "../transformation/proc_rotx":
read "../transformation/proc_roty":
read "../transformation/proc_rotz":
read "../transformation/proc_trotz":
read "../transformation/proc_transl":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
read "../robot_codegen_constraints/proc_subs_kintmp_exp":
read "../helper/proc_intersect_circle":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env_IC":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name_OL):
# Ergebnisse der Kinematik laden
read sprintf("../codeexport/%s/tmp/kinematics_floatb_%s_rotmat_maple.m", robot_name_OL, base_method_name);
read "../robot_codegen_definitions/robot_env_IC":
Trf := Trf:
Trf_c := Trf_c:
# VGK Gelb 0-7-13-10-4-2-1-0
# Schleife (0-7)-(7-13)
T_1_13 := combine( Matrix(Trf(1..4,1..4, 7)) . Matrix(Trf(1..4,1..4,13))):
# Schleife (0-1)(1-2)-(2-4)-(4-10)
T_1_10:= combine( Matrix(Trf(1..4,1..4, 1)) . Matrix(Trf(1..4,1..4,2)) . Matrix(Trf(1..4,1..4,4)) . Matrix(Trf(1..4,1..4,10))):
h1t := T_1_13(1..3,4) - T_1_10(1..3,4):

#tmp := Transpose( Matrix(T_1_13(1..3,1..3)) ) . Matrix(T_1_10(1..3,1..3));  nur anzeigen lassen für h1r
;
#combine(tmp);# nur anzeigen lassen für h1r
;
h1r := (-qJ7(t)+qJ1(t)+qJ2(t)+qJ4(t)+qJ10(t))+Pi/2:
# VGK PINK 0-5-14-11-8-1-0
# Schleife (0-5)-(5-14)
T_2_14 := combine(Matrix(Trf(1..4,1..4, 5)) . Matrix(Trf(1..4,1..4,14))):
# Schleife (0-1)-(1-8)-(8-11)
T_2_11:= combine(Matrix(Trf(1..4,1..4, 1)) . Matrix(Trf(1..4,1..4,8)).Matrix(Trf(1..4,1..4, 11))):
h2t := T_2_11(1..3,4) - T_2_14(1..3,4):
#tmp := Transpose( Matrix(T_2_11(1..3,1..3)) ) . Matrix(T_2_14(1..3,1..3));  # nur anzeigen lassen für h1r
;
 #combine(tmp); # nur anzeigen lassen für h1r
;
h2r:=(-qJ1(t)-qJ8(t)-qJ11(t)+phi05+qJ5(t)):
# VGK GRÜN 2-6-12-15-9-3
# Schleife (2-6)-(6-12)
T_2_12 := combine( Matrix(Trf(1..4,1..4, 6)) . Matrix(Trf(1..4,1..4,12))):
# Schleife (2-3)-(3-9)-(9-15)
T_2_15:= combine( Matrix(Trf(1..4,1..4, 3)) . Matrix(Trf(1..4,1..4,9)).Matrix(Trf(1..4,1..4,15)) ):
h3t := T_2_12(1..3,4) - T_2_15(1..3,4):
#tmp := Transpose( Matrix(T_2_12(1..3,1..3)) ) . Matrix(T_2_15(1..3,1..3));  nur anzeigen lassen für h2r
;
#combine(tmp);  nur anzeigen lassen für h3r
;
h3r:= qJ6(t) + qJ12(t) - qJ3(t) - qJ9(t):
# Zusätzliche Zwangsbedingung (6 und 8 sind gleiche Körper)
# TODO: In einer neuen Modellierung sollten beide Körper zusammengefasst werden.
h4:=Pi/2 - qJ8(t) - qJ9(t) + qJ2(t):
# Zusammenstellen aller Zwangsbedingungen
implconstr_t := <h1t([1, 2],1);h2t([1, 2],1);h3t([1, 2],1); h1r; h2r; h3r; h4>; # TODO: In h1r, h2r, h3r muss das richtige drinstehen.
implconstr_s := convert_t_s(implconstr_t):
# Exportiere Code für folgende Skripte
kin_constraints_exist:=true:
save implconstr_t, implconstr_s, kin_constraints_exist, sprintf("../codeexport/%s/tmp/kinematic_constraints_implicit_maple.m", robot_name):
# Exportieren des vollständigen Ausdruckes
if codegen_act then
  MatlabExport(implconstr_s, sprintf("../codeexport/%s/tmp/kinconstr_impl_matlab.m", robot_name), 2):
end if:
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( implconstr_s )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_implicit_constraints_symbols_list_matlab.m", robot_name), 2):

