# Vollständige Liste der Ersetzungsausdrücke für palh3m1DE1
Format: Kennzeichnung des Ausdrucks als Überschrift "##" gefolgt von dessen Inhalt in neuer Zeile (umschlossen von Code-Hervorhebung ("Backtick"-Zeichen).
Zur Ersetzung im Quelltext werden die Ausdrücke mit "%"-Zeichen umklammert.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07  
(C) Institut für Mechatronische Systeme, Universität Hannover

## RN

```
palh3m1DE1
```

## NQJ

```
4
```

## NJ

```
12
```

## NL

```
9
```

## NMPVFIXB

```
52
```

## NMPVFLOATB

```
NOTDEFINED
```

## NKP

```
19
```

## NKCP

```
palh3m1DE1
```

## KPDEF

```
pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
```

## NTAUJFIXBREGNN

```
130
```

## VERSIONINFO

```
% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
```

## INPUT_M

```
% m [9x1]
%   mass of all robot links (including the base)
```

## INPUT_R

```
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
```

## INPUT_MR

```
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
```

## INPUT_IC

```
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
```

## INPUT_IF

```
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
```

## INPUT_MDPFIXB

```
% MDP [52x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh3m1DE1_convert_par2_MPV_fixb.m
```

## INPUT_PKIN

```
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
```

## INPUT_QJ

```
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
```

## INPUT_QJD

```
% qJD [4x1]
%   Generalized joint velocities
```

## INPUT_QJDD

```
% qJDD [4x1]
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

