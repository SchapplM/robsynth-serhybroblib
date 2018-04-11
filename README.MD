# Funktionssammlung für Modelle hybrider Roboter

Dieses Repo enthält Modelle und Modultests für hybride Roboter (serielle Roboter mit parallelen Teilstrukturen oder kinematischen Zwangsbedingungen).


Moritz Schappler, schappler@imes.uni-hannover.de, 2018-04  
(C) Institut für mechatronische Systeme, Universität Hannover

## Aufbau

Im Ordner `tb_code` befinden sich Unterordner für Systemmodelle verschiedener Roboter. Für jeden Roboter gibt es ein zusätzliches Systemmodell mit Endung "u". Bei diesem sind die kinematischen Zwangsbedingungen entfernt, also z.B. geschlossene Ketten geöffnet.

Der Ordner `test` enthält Testfunktionen für die Berechnung der Dynamik basierend auf verschiedenen Darstellungen der kinematischen Zwangsbedingungen.

## Voraussetzungen

Es werden Folgende Repositorien als Abhängigkeit benötigt (Pfadinitialisierung notwendig):

* IMES-Matlab-Toolbox
* IMES-Robotik-Toolbox
* Sammlung externer Matlab Toolboxen

Folgende Repositorien dienen zur Herleitung der Dynamikgleichungen

* HybridDyn (ehem. IRT-Maple-Toolbox)

## Liste von erfolgreich getesteten Robotern

Folgende Roboter und Mechanismen wurden bisher erfolgreich getestet


| Roboter | Name des Ordners | Art des Tests | Implementierung mit |
| ---     |                  |---            |---                  |
| 3. Arm KAS5 Modell 5| KAS5m5 | Elim., expl. ZB | HybridDyn |
| 3. Arm KAS5 Modell 6| KAS5m6 | Elim., expl. ZB | HybridDyn |
| 3. Arm KAS5 Modell 7| KAS5m7 | Elim., expl. ZB, impl. ZB | HybridDyn |
| 3. Arm KAS7 Modell 1| KAS7m1 | Elim., expl. ZB | HybridDyn |

Erklärung zu den Tests:
* Elim.: Abhängige Gelenke aus den kinematischen Zwangsbedingungen werden in Dynamikgleichungen eliminiert
* expl. ZB: Die explizite Form der Zwangsbedingungen dient zur Projektion der Gelenke der offenen Kette auf die geschlossene Kette
* impl. ZB: Die implizite Form der Zwangsbedingungen wird dafür benutzt

## Quellen <a name="quellen"></a> 

* [NuelleSchTapLil2017] K. Nülle et al.: Projektabschlussbericht „3. Arm“ – Handwerker-Kraftassistenzsystem mit adaptiver Mensch-Technik-Interaktion
