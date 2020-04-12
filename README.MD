# Bibliothek hybrider Roboter

Dieses Repo enthält eine Sammlung von Modellen hybrider Roboter (serielle Roboter mit parallelen Teilstrukturen bzw. sonstigen kinematischen Zwangsbedingungen).
Die Modelle werden analytisch hergeleitet und liegen mit Kinematik und Dynamik in symbolischer Form und als Matlab-Code vor. Daher eignen sie sich für Vorsteuerungen und modellbasierte Regler in Echtzeit-Steuerungen (z.B. SPS oder Simulink).

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-04  
(C) Institut für Mechatronische Systeme, Leibniz Universität Hannover


## Benutzung

Für die Benutzung des Repos muss mit dem Skript `hybroblib_path_init.m` aus dem Hauptordner die Pfadinitialisierung durchgeführt werden.
Anschließend können die einzelnen abgespeicherten Modelle mit Hilfe der Testskripte in den jeweiligen Ordnern geprüft werden.
Der Matlab-Code für Kinematik und Dynamik kann entweder direkt in einer Simulink-basierten Vorsteuerung/Regelung für die Roboter eingesetzt werden oder automatisch in SPS-Code übersetzt werden.

### Abhängigkeiten

Die Matlab-Testskripte der Robotermodelle setzen folgende Toolboxen mit initialisiertem Pfad voraus:

* Trajektorien-Toolbox: [github](https://github.com/SchapplM/robotics-dep-traj), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-material/trajektorien-toolbox)
* Robotik-Toolbox: [github](https://github.com/SchapplM/robotics-toolbox), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-gruppe-robotik/matlab-robotik-toolbox) 
  * Matlab-Toolbox: [github](https://github.com/SchapplM/matlab_toolbox), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-material/matlab-toolbox)
  * Externe Abhängigkeiten: [github](https://github.com/SchapplM/robotics-dep-ext), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-material/matlab-ext)
  * Geometrie: [github](https://github.com/SchapplM/robotics-dep-geometry), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-material/geometrie-toolbox)

Zur Generierung der Modelle mit Maple sind folgende Schritte notwendig
* Einrichtung der Modelgenerierungs-Toolbox: [github](https://github.com/SchapplM/robsynth-modelgen), [Uni Hannover intern](https://gitlab.projekt.uni-hannover.de/imes-gruppe-robotik/hybriddyn)
* Erstellung der Pfad-Definition in der Datei `hybrdyn_repo_path` (ohne Endung) gemäß Hinweisen in `hybrdyn_repo_path.template`

Hinweis: Der symbolisch generierte Code selbst ist ohne Abhängigkeiten ausführbar.

## Aufbau

Im Ordner `systems` befinden sich Unterordner für Systemmodelle verschiedener Roboter:

* systems/fourbar1
* systems/palh1m1
* systems/palh2m1
* ...

Für jeden der Roboter gibt es unterschiedliche Implementierungen der kinematischen Zwangsbedingungen in einem weiteren Unterordner. Die Art der Implementierung ist dem Roboternamen als Kürzel angehängt:

* ...IC (Herleitung der Dynamik durch implizite Zwangsbedingungsgleichungen und Projektion der vollständigen Dynamik auf die Minimalkoordinaten, "Implicit Constraints")
* ...TE (Herleitung der Dynamik in Minimalkoordinaten durch trigonometrische Elimination der abhängigen (passiven) Gelenkwinkel in der Lagrange-Funktion)
* ...DE1 (Herleitung der Dynamik in Minimalkoordinaten durch direkte Elimination der abhängigen Koordinaten)
* ...DE2 (Herleitung wie in DE1, die Minimalkoordinaten werden aber erst substituiert, nachdem trigonometrische Additionstheoreme für parallele Drehachsen angewendet wurden)

Die Ansätze werden in [SchapplerLilHad2019](#quellen) und in der Dokumentation des [HybrDyn-Repos](https://github.com/SchapplM/robsynth-modelgen) beschrieben. Grundsätzlich wird die Dynamik mit der Lagrange'schen Gleichung 2. Art hergeleitet.

### Test-Funktionen

Die jeweiligen Ordner der Roboter enthalten Testskripte mit Beispiel-Aufrufen und Visualisierungen.
Das Skript `hybroblib_test_all_models.m` ruft alle Testfunktionen nacheinander auf und dient zum Gesamt-Test dieses Repos.

### Code-Generierung

Für jedes Robotermodell gibt es den Unterordner `codegen` mit manuell erstellten Definitionen und Herleitungen der kinematischen Zwangsbedingungen. Für jeden Roboter gibt es die Bash-Skripte
* `prepare_maple_repo.sh` zum Kopieren der Definitionsdateien in die Toolbox `HybrDyn`,
* `generate_maple_code.sh` zum Erzeugen des Codes in der HybrDyn-Toolbox für alle Implementierungen,
* `copy_generated_code.sh` zum Zurückkopieren des generierten Codes in diese Datenbank.

Die Code-Generierung muss nur bei Modelländerung oder -erstellung durchgeführt werden, da die fertigen Matlab-Funktionen in dieses Repo eingecheckt werden.

Für jeden Roboter ist eine Datei `...kinematic_parameter_values.m` und `...kinematic_constraints_matlab.m` abgelegt. Diese dienen nur dazu, im HybrDyn-Repo die Testfunktionen automatisch ausführen zu können.

Es gibt für jeden Roboter zwei Definitionsdateien `robot_env_...` mit Endung OL ("Open Loop", ohne kinematische Zwangsbedingungen) und CL ("Closed Loop", mit eingesetzten Minimalkoordinaten für die Eliminierung der Zwangsbedingungen). Diese Dateien werden durch das Kopierskript `prepare_maple_repo.sh` an die konkreten Implementierungen angepasst, damit die Informationen nur an einer Stelle gespeichert sind.
Die Definitionsdatei mit Endung IC ("Implicit Constraints") generiert die Dynamik der geschlossenen Struktur aus der offenen Struktur mit den implizit definierten kinematischen Zwangsbedingungen.

Die Maple-Arbeitsblätter `...TE_kinematic_constraints.mw`, `...DE_kinematic_constraints.mw` und `...IC_kinematic_constraints_implicit.mw` enthalten Definitionen der jeweiligen kinematischen Zwangsbedingungen. Die gleichnamigen Dateien mit Endung `.mpl` enthalten den Kopie des Codes (in Maple exportiert) und werden in das HybrDyn-Repo kopiert. Zum Aufbau der Dateien sei auf die Dokumentation des [HybrDyn-Repos](https://github.com/SchapplM/robsynth-modelgen) verwiesen.

### Funktionen

Benutzung der Datenbank:

* `hybroblib_create_robot_class.m`: Instanz einer Matlab-Klasse für einen Roboter erstellen, damit sind dann weitere Berechnungen mit dem Modell des Roboters möglich

### Übersichtsdateien

Alle Informationen der Roboterbibliothek liegen in csv-Dateien (Trennzeichen: Semikolon) vor. Dadurch ist eine automatische Verarbeitung und eine Versionsverwaltung möglich.

#### Liste der Roboterparameter für ein Robotermodell

Die Datei `systems/Robotername/models.csv` enthält mögliche Roboterparameter für das entsprechende Robotermodell.

Die Tabelle enthält Zahlenwerte für die einzelnen Kinematik-Parameter aus der Variable `pkin` des HybrDyn-Repos sowie weitere Roboterdaten wie Gelenkwinkelgrenzen und ggf. Traglast oder Masse.

## Liste der Robotermodelle

### Übersicht

Folgende Roboter und Mechanismen wurden bisher erfolgreich getestet:

| Roboter | Name des Ordners | Art des Tests | Implementierte Methoden |
| ---     |                  |---            |---                  |
| [Viergelenkkette aus Drehgelenken](#system_fourbar1) | fourbar1 | Plausibilität/Visualisierung, Abgleich Expl/Impl | TE, DE, IC |
| [Parallelogramm](#system_fourbar2) | fourbar2 | Abgleich über palh3m2 | TE, DE, IC |
| [Viergelenkkette mit Schubgelenk](#system_fourbarpris) | fourbarpris | Plausibilität/Visualisierung, Abgleich Expl/Impl | TE, DE, IC |
| [Drehbar aufgerichtete Viergelenkkette](#system_fourbar1turn) | fourbar1turn | Validierung SA Shan | TE, DE, IC |
| [Fünfgelenkkette](#system_fivebar1) | fivebar1 | Plausibilität/Visualisierung/Abgleich über picker2Dm1 | TE, DE, IC |
| [2T0R-Palettierroboter](#system_picker2D) (2 Parallelogramme, 1 Fünfgelenkkette) | picker2Dm1, picker2Dm2 | Plausibilität/Visualisierung | TE, DE, IC |
| [3T1R-Palettierroboter](#system_palh1) (3 Parallelogramme) | palh1m1, palh1m2 | Plausibilität/Visualisierung, Kinematik | TE, DE, IC |
| [3T1R-Palettierroboter](#system_palh2) (serielles Ersatzmodell) | palh2m1, palh2m2 | Kinematik, Abgleich Implementierungen DE | Elim.+IC |
| [3T1R-Palettierroboter](#system_palh3) (2 Parallelogramme) | palh3m1, palh3m2 | Plausibilität | TE, DE, IC |
| [Planare Hybrid-PKM-Beinkette](#system_hybBK) (Parallelogramm-Fünfgelenkkette) | hybBKplanar | MA Brünger | TE |
| [Räumliche Hybrid-PKM-Beinkette](#system_hybBK) (Parallelogramm-Fünfgelenkkette) | hybBKspatial | MA Brünger | TE |

### Viergelenkkette mit Drehgelenken <a name="system_fourbar1"></a> 

Die Viergelenkkette (mit Drehgelenken) dient als allgemeines Beispiel und als Komponente komplexerer Roboter.

| Viergelenkkette |  |
| ---     | |
| Modell-Name | fourbar1 |
| Minimalkoordinaten / Freiheitsgrad | 1 |
| Anzahl geschlossener Ketten | 1 |

### Parallelogramm <a name="system_fourbar2"></a> 

Das Parallelogramm dient als allgemeines Beispiel und als Komponente komplexerer Roboter.

| Parallelogramm |  |
| ---     | |
| Modell-Name | fourbar2 |
| Minimalkoordinaten / Freiheitsgrad | 1 |
| Anzahl geschlossener Ketten | 1 |

### Viergelenkkette mit Schubgelenk <a name="system_fourbarpris"></a> 

Diese Viergelenkkette (mit einem Schubgelenk) dient als Komponente komplexerer Roboter. Es handelt sich bei dem Mechanismus nicht um eine allgemeine Schubkurbel, da das feststehende Glied zwischen den Drehgelenken liegt.

| Viergelenkkette |  |
| ---     | |
| Modell-Name | fourbarpris |
| Minimalkoordinaten / Freiheitsgrad | 1 |
| Anzahl geschlossener Ketten | 1 |

### Fünfgelenkkette mit Drehgelenken <a name="system_fivebar1"></a> 

Die Fünfgelenkkette (mit Drehgelenken) dient als allgemeines Beispiel und als Komponente komplexerer Roboter. Die symbolische Lösung ist noch sehr rechenaufwändig und damit auch die daraus berechneten Roboter.

| Fünfgelenkkette |  |
| ---     | |
| Modell-Name | fivebar1 |
| Minimalkoordinaten / Freiheitsgrad | 2 |
| Anzahl geschlossener Ketten | 1 |


### Drehbar aufgerichtete Viergelenkkette <a name="system_fourbar1turn"></a> 

Die allgemeine Viergelenkkette wird senkrecht mit einem zusätzlichen Drehgelenk aufgestellt. Dieses Modell dient zum Testen kinematischer Zwangsbedingungen des Teils eines räumlichen Palettierroboters.

| Drehbar aufgerichtete Viergelenkkette |  |
| ---     | |
| Modell-Name | fourbar1turn |
| Minimalkoordinaten / Freiheitsgrad | 2 |
| Anzahl geschlossener Ketten | 1 |


### Planarer Palettierroboter <a name="system_picker2D"></a> 

Dieser planare Palettierroboter ("Picker") hat zwei kartesische FG mit fester Orientierung.

| 2T0R-Palettierer |  | Anmerkung |
| ---     | | |
| Modell-Name | picker2Dm1 | Basierend auf allgemeiner Viergelenkkette `fourbar1` |
| Modell-Name | picker2Dm2 | Basierend auf Parallelogramm `fourbar2` |
| Minimalkoordinaten / Freiheitsgrad | 2 |
| Anzahl geschlossener Ketten | 3 | 2 Parallelogramme, 1 Fünfgelenkkette |
| Endeffektor-FG | 2 | nur translatorisch in der Ebene |

Diese Kinematik wird u.a. in den folgenden realen Robotern eingesetzt:

* [Codian D2](https://www.codian-robotics.com/de/roboter/d2-roboter/)
* [TrioWin TDR2040](http://en.triowin.com/tdr-series2-axisparallelrobot-15268685491267769.html)

### Palettierroboter mit drei Parallelogrammen <a name="system_palh1"></a> 

Dieser Palettierroboter hat vier kartesische FG (drei Translationen, eine Rotation um die Hochachse). Die Konstruktion mit drei Parallelogrammen wird überwiegend im Bereich großer Traglasten eingesetzt.

| 3T1R-Palettierer |  | Anmerkung |
| ---     | | |
| Modell-Name | palh1m1 | Basierend auf allgemeiner Viergelenkkette `fourbar1` |
| Modell-Name | palh1m2 | Basierend auf Parallelogramm `fourbar2` |
| Minimalkoordinaten / Freiheitsgrad | 4 |
| Anzahl geschlossener Ketten | 3 | 3 Parallelogramme |
| Antriebe | 4 | Durch das dritte Parallelogramm wird ein Antrieb nach unten gesetzt |
| Endeffektor-FG | 4 | drei translatorisch, 1 rotatorisch |

Diese Kinematik wird u.a. in den folgenden realen Robotern eingesetzt:

* [Yaskawa MPL 800](https://www.motoman.com/industrial-robots/mpl800-ii) und [MPL500](https://www.motoman.com/en-us/products/robots/industrial-robots/packaging-palletizing/mpl-series)

### Palettierroboter mit zwei Parallelogrammen <a name="system_palh3"></a> 

Dieser Palettierroboter hat vier kartesische FG (drei Translationen, eine Rotation um die Hochachse). Die Konstruktion mit zwei Parallelogrammen wird überwiegend im mittleren Traglastbereich eingesetzt.

| 3T1R-Palettierer |  | Anmerkung |
| ---     | | |
| Modell-Name | palh3m1 | Basierend auf allgemeiner Viergelenkkette `fourbar1` |
| Modell-Name | palh3m2 | Basierend auf Parallelogramm `fourbar2` |
| Minimalkoordinaten / Freiheitsgrad | 4 |
| Anzahl geschlossener Ketten | 2 | 2 Parallelogramme |
| Antriebe | 4 |  |
| Endeffektor-FG | 4 | drei translatorisch, 1 rotatorisch |

Diese Kinematik wird u.a. in den folgenden realen Robotern eingesetzt:

* [KUKA KR700 PA](https://www.kuka.com/de-de/produkte-leistungen/robotersysteme/industrieroboter/kr-700-pa)

### Ersatzmodelle für Palettierroboter <a name="system_palh2"></a> 

Für Palettierroboter ohne parallele Zwangsführung durch Parallelogramme existiert eine mathematisch kinematische Zwangsbedingung zwischen den fünf Gelenkwinkeln der Antriebsachsen, damit die Endeffektor-FG auf vier beschränkt bleiben.

Durch diese Annahme lassen sich auch Ersatzmodelle für die komplexeren Palettierroboter mit Parallelogrammen definieren.

| 3T1R-Palettierer |  | Anmerkung |
| ---     | | |
| Modell-Name | palh2m1 | Direkte Einsetzung eines Gelenkwinkels als abhängig |
| Modell-Name | palh2m2 | Mit Zusätzlichem Gelenk |
| Minimalkoordinaten / Freiheitsgrad | 4 |
| Anzahl Gelenke | 5 bzw. 6 | |
| Anzahl geschlossener Ketten | 0 |  |
| Zusätzliche kinematische Zwangsbed. | ja |  |
| Endeffektor-FG | 4 | drei translatorisch, 1 rotatorisch (mathematisch) |

Diese Ersatzmodelle werden in folgendem Testskript gegen die vollständigen Modelle validiert: `test/palletizer_mdlcomp.m`

### Beinkette für Hybrid-PKM <a name="system_hybBK"></a> 

Für PKM mit Beinketten, die ihrerseits aus seriell-hybriden Kinematiken bestehen, gibt es zwei Modelle:

| Hybrid-PKM Beinketten |  | Anmerkung |
| ---     | | |
| Modell-Name | hybBKplanar | Planare Beinkette |
| Modell-Name | hybBKspatial | Räumliche Beinkette |
| Anzahl geschlossener Ketten | 1 |  |
| Minimalkoordinaten / Freiheitsgrad | 3 bzw. 6 |  |
| Endeffektor-FG | 3 bzw. 6 |  |

Die Beinketten eignen sich zum Aufbau einer planaren PKM (2T1R) oder einer räumlichen PKM (3T3R). Die Beinketten basieren jeweils auf einer Fünfgelenkkette mit einem zusätzlichen Plattform-Koppelgelenk. Für Details siehe [MA Brünger](#quellen).

### Palettierroboter mit einer geschlossenen Kette und Schubgelenk <a name="system_palh4"></a> 

Dieser Palettierroboter stammt aus einem Beispiel von Khalil und Bennis ([KhalilBen1995](#quellen), Fig. 4, Beispiel 1). 
Dieses Modell ist noch unfertig und die Kinematik ist noch nicht richtig implementiert.

| Palettierroboter |  | Anmerkung |
| ---     | | |
| Modell-Name | palh4m1 | Modellierung noch falsch|
| Anzahl geschlossener Ketten | 1 |  |
| Minimalkoordinaten / Freiheitsgrad | 5 |  |
| Antriebe | 5 |  |
| Endeffektor-FG | 4 | drei translatorisch, 1 rotatorisch (mathematisch) |


### Schwerlastroboter mit zwei geschlossenen Ketten <a name="system_mg10hl"></a> 

Durch zusätzliche geschlossene Ketten wird die Kraftübertragung der Antriebe für einen Schwerlastroboter verbessert.
Dieses Modell ist noch unfertig und die Kinematik ist noch nicht richtig implementiert.

| Schwerlastroboter |  | Anmerkung |
| ---     | | |
| Modell-Name | mg10hl | Modellierung noch falsch|
| Minimalkoordinaten / Freiheitsgrad | 6 |  |
| Anzahl geschlossener Ketten | 2 | eine Viergelenkkette mit Drehgelenken und eine mit Schubgelenk |
| Antriebe | 6 |  |
| Endeffektor-FG | 6 |  |

Diese Kinematik wird in dem folgenden realen Roboter eingesetzt:

* [Kawasaki MG10HL](https://robotics.kawasaki.com/de1/products/robots/extra-large-payloads/MG10HL/)

## Quellen <a name="quellen"></a> 

* [Bejaoui2018_S749] A. Bejaoui: "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen", 2018 (Studienarbeit bei Moritz Schappler)
* [Brünger2019_M832] A. Brünger: "Vergleich der Leistungsmerkmale und Modellierung von parallelen Robotern mit hybriden und seriellen Beinketten", 2019 (Masterarbeit bei Moritz Schappler)
* [KhalilBen1995] W. Khalil, F. Bennis: "Symbolic Calculation of the Base Inertial Parameters of Closed-Loop Robots", 1995
* [SchapplerLilHad2019] Schappler, M. and Lilge, T. and Haddadin, S.: "Kinematics and Dynamics Model via Explicit Trigonometric Elimination of Kinematic Constraints for a Force Assistance Exoskeleton", 2019
* [Shan2019_S828] W. Shan: "Reduktion der Modellkomplexität von seriell-hybriden Palettierrobotern", 2019 (Studienarbeit bei Moritz Schappler)

