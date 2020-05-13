# Vollständige Liste der Ersetzungsausdrücke für palh4m1OL
Format: Kennzeichnung des Ausdrucks als Überschrift "##" gefolgt von dessen Inhalt in neuer Zeile (umschlossen von Code-Hervorhebung ("Backtick"-Zeichen).
Zur Ersetzung im Quelltext werden die Ausdrücke mit "%"-Zeichen umklammert.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07  
(C) Institut für Mechatronische Systeme, Universität Hannover

## RN

```
palh4m1OL
```

## NQJ

```
8
```

## NJ

```
9
```

## NL

```
8
```

## NMPVFIXB

```
NOTDEFINED
```

## NMPVFLOATB

```
NOTDEFINED
```

## NKP

```
7
```

## NKCP

```
palh4m1OL
```

## KPDEF

```
pkin=[AB,CB,CE,EP,OT,TA,TD]';
```

## NTAUJFIXBREGNN

```
NOTDEFINED
```

## VERSIONINFO

```
% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
```

## INPUT_M

```
% m [8x1]
%   mass of all robot links (including the base)
```

## INPUT_R

```
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
```

## INPUT_MR

```
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
```

## INPUT_IC

```
% Icges [8x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
```

## INPUT_IF

```
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
```

## INPUT_MDPFIXB

```
% MDP [NOTDEFINEDx1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh4m1OL_convert_par2_MPV_fixb.m
```

## INPUT_PKIN

```
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
```

## INPUT_QJ

```
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
```

## INPUT_QJD

```
% qJD [8x1]
%   Generalized joint velocities
```

## INPUT_QJDD

```
% qJDD [8x1]
%   Generalized joint accelerations
```

## INPUT_RB

```
% r_base [3x1]
%   Base position in world frame
```

## INPUT_PHIB

```
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
```

## INPUT_PHIBD

```
% phiD_base [3x1]
%   Time Derivative of Base Orientation in world frame.
%   Expressed with XYZ-Euler angles ("rpy")
```

## INPUT_XDB

```
% xD_base [6x1]
%   time derivative of r_base and phi_base
```

## INPUT_XDDB

```
% xDD_base [6x1]
%   second time derivative of r_base and phi_base
```

