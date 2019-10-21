# Bibliothek hybrider Roboter

Dieses Repo enthält eine Sammlung von Modellen hybrider Roboter (serielle Roboter mit parallelen Teilstrukturen bzw. sonstigen kinematischen Zwangsbedingungen).


Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-04  
(C) Institut für Mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss mit dem Skript `hybroblib_path_init.m` aus dem Hauptordner die Pfadinitialisierung durchgeführt werden.

### Abhängigkeiten

Die Robotermodelle setzen die Toolbox imes-matlab und imes-matlab-robotics mit initialisiertem Pfad voraus.

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

Die Ansätze werden in [SchapplerLilHad2019](#quellen) und in der Dokumentation des HybrDyn-Repos beschrieben.

### Test-Funktionen

Die jeweiligen Ordner der Roboter enthalten Testskripte mit Beispiel-Aufrufen und Visualisierungen.
Das Skript `hybroblib_test_all_models.m` ruft alle Testfunktionen nacheinander auf und dient zum Gesamt-Test dieses Repos.

### Code-Generierung

Für jedes Robotermodell gibt es den Unterordner `codegen` mit manuell erstellten Definitionen und Herleitungen der kinematischen Zwangsbedingungen. Für jeden Roboter gibt es die Bash-Skripte
* `prepare_maple_repo.sh` zum Kopieren der Definitionsdateien in die Toolbox `HybrDyn`,
* `generate_maple_code.sh` zum Erzeugen des Codes in der HybrDyn-Toolbox für alle Implementierungen,
* `copy_generated_code.sh` zum Zurückkopieren des generierten Codes in diese Datenbank.

Die Code-Generierung muss nur bei Modelländerung oder -erstellung durchgeführt werden, da die fertigen Funktionen ins Repo eingecheckt werden.

Für jeden Roboter ist eine Datei `...kinematic_parameter_values.m` und `...kinematic_constraints_matlab.m` abgelegt. Diese dienen nur dazu, im HybrDyn-Repo die Testfunktionen automatisch ausführen zu können.

Es gibt für jeden Roboter zwei Definitionsdateien `robot_env_...` mit Endung OL ("Open Loop", ohne kinematische Zwangsbedingungen) und CL ("Closed Loop", mit eingesetzten Minimalkoordinaten für die Eliminierung der Zwangsbedingungen). Diese Dateien werden durch das Kopierskript `prepare_maple_repo.sh` an die konkreten Implementierungen angepasst, damit die Informationen nur an einer Stelle gespeichert sind.

### Funktionen

Benutzung der Datenbank:

* `hybroblib_create_robot_class.m`: Instanz einer Matlab-Klasse für einen Roboter erstellen, damit sind dann weitere Berechnungen mit dem Modell des Roboters möglich

### Übersichtsdateien

Alle Informationen der Roboterbibliothek liegen in csv-Dateien (Trennzeichen: Semikolon) vor. Dadurch ist eine automatische Verarbeitung und eine Versionsverwaltung möglich.

#### Liste der Roboterparameter für ein Robotermodell

Die Datei `systems/Robotername/models.csv` enthält mögliche Roboterparameter für das entsprechende Robotermodell.

Die Tabelle enthält Zahlenwerte für die einzelnen Kinematik-Parameter aus der Variable `pkin` des HybrDyn-Repos sowie weitere Roboterdaten wie Gelenkwinkelgrenzen und ggf Traglast oder Masse.

## Liste von erfolgreich getesteten Robotern

Folgende Roboter und Mechanismen wurden bisher erfolgreich getestet:

| Roboter | Name des Ordners | Art des Tests | Implementierte Methoden |
| ---     |                  |---            |---                  |
| Viergelenkkette | fourbar1 | Plausibilität/Visualisierung, Abgleich Expl/Impl | TE, DE, IC |
| Parallelogramm | fourbar2 | Abgleich über palh3m2 | TE, DE, IC |
| Drehbar aufgerichtete Viergelenkkette | fourbar3 | Validierung SA Shan | nur Elim. |
| Fünfgelenkkette | fivebar1 | Plausibilität/Visualisierung/Abgleich über TSR | TE, DE, IC |
| 2T0R-Palettierroboter (2 Parallelogramme, 1 Fünfgelenkkette, Herl. über allgemeine Viergelenkkette) | TSR | Plausibilität/Visualisierung | TE, DE, IC |
| 2T0R-Palettierroboter (2 Parallelogramme, 1 Fünfgelenkkette, Herl. über Parallelogramm) | TSRm2 | Plausibilität/Visualisierung | TE, DE, IC |
| 3T1R-Palettierroboter (3 Parallelogramme, Herl. über allgemeine Viergelenkkette) | palh1m1 | Plausibilität/Visualisierung, Kinematik | TE, DE, IC |
| 3T1R-Palettierroboter (3 Parallelogramme, Herl. über Parallelogramm) | palh1m2 | Plausibilität/Visualisierung, Kinematik | TE, DE, IC |
| 3T1R-Palettierroboter (serielles Ersatzmodell) | palh2m1 | Kinematik, Abgleich Implementierungen DE | Elim.+IC |
| 3T1R-Palettierroboter (serielles Ersatzmodell, mit anderen Parametern) | palh2m2 | Validiert gegen Lenze-Daten | Elim.+IC |
| 3T1R-Palettierroboter (2 Parallelogramme, Herleitung über allgemeine Viergelenkkette) | palh3m1 | Plausibilität | TE, DE, IC |
| 3T1R-Palettierroboter (2 Parallelogramme, Herleitung über Parallelogramm) | palh3m2 | Plausibilität | TE, DE, IC |
| Planare Hybrid-PKM-Beinkette (Parallelogramm-Fünfgelenkkette) | hybBKplanar | MA Brünger | TE |
| Räumliche Hybrid-PKM-Beinkette (Parallelogramm-Fünfgelenkkette) | hybBKspatial | MA Brünger | TE |
| 3. Arm KAS5 Modell 5| KAS5m5 | Elim., expl. ZB (noch nicht in neuer Dateistruktur) | Elim. |
| 3. Arm KAS5 Modell 6| KAS5m6 | Elim., expl. ZB (noch nicht in neuer Dateistruktur) | Elim. |
| 3. Arm KAS5 Modell 7| KAS5m7 | Elim., expl. ZB, impl. ZB (noch nicht in neuer Dateistruktur) | Elim. |
| 3. Arm KAS7 Modell 1| KAS7m1 | Elim., expl. ZB (noch nicht in neuer Dateistruktur) | Elim. |

## Quellen <a name="quellen"></a> 

* [Shan2019_S828] W. Shan: "Reduktion der Modellkomplexität von seriell-hybriden Palettierrobotern" (Studienarbeit bei Moritz Schappler)
* [Brünger2019_M832] A. Brünger: "Vergleich der Leistungsmerkmale und Modellierung von parallelen Robotern mit hybriden und seriellen Beinketten" (Masterarbeit bei Moritz Schappler)
* [Bejaoui2018_S749] A. Bejaoui: "Modellierung kinematischer Zwangsbedingungen für hybride serielle Roboter mit planaren Parallelmechanismen" (Studienarbeit bei Moritz Schappler)
* [NuelleSchTapLil2017] K. Nülle, M. Schappler et al.: Projektabschlussbericht „3. Arm“ - Handwerker-Kraftassistenzsystem mit adaptiver Mensch-Technik-Interaktion
* [SchapplerLilHad2019] Schappler, M. and Lilge, T. and Haddadin, S.: Kinematics and Dynamics Model via Explicit Trigonometric Elimination of Kinematic Constraints for a Force Assistance Exoskeleton (2019)
