# Vollständige Liste der Ersetzungsausdrücke für palh2m2DE
Format: Kennzeichnung des Ausdrucks als Überschrift "##" gefolgt von dessen Inhalt in neuer Zeile (umschlossen von Code-Hervorhebung ("Backtick"-Zeichen).
Zur Ersetzung im Quelltext werden die Ausdrücke mit "%"-Zeichen umklammert.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07  
(C) Institut für Mechatronische Systeme, Universität Hannover

## RN

```
palh2m2DE
```

## NQJ

```
4
```

## NJ

```
6
```

## NL

```
7
```

## NMPVFIXB

```
22
```

## NMPVFLOATB

```
NOTDEFINED
```

## NKP

```
5
```

## NKCP

```
palh2m2DE
```

## KPDEF

```
pkin=[A2Off,A3Off,A4Off,L1,L2]';
```

## NTAUJFIXBREGNN

```
46
```

## VERSIONINFO

```
% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
```

## INPUT_M

```
% m [7x1]
%   mass of all robot links (including the base)
```

## INPUT_R

```
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
```

## INPUT_MR

```
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
```

## INPUT_IC

```
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
```

## INPUT_IF

```
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
```

## INPUT_MDPFIXB

```
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2DE_convert_par2_MPV_fixb.m
```

## INPUT_PKIN

```
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
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

