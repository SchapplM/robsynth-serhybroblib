# Vollständige Liste der Ersetzungsausdrücke für fourbarprisDE1
Format: Kennzeichnung des Ausdrucks als Überschrift "##" gefolgt von dessen Inhalt in neuer Zeile (umschlossen von Code-Hervorhebung ("Backtick"-Zeichen).
Zur Ersetzung im Quelltext werden die Ausdrücke mit "%"-Zeichen umklammert.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07  
(C) Institut für Mechatronische Systeme, Universität Hannover

## RN

```
fourbarprisDE1
```

## NQJ

```
1
```

## NJ

```
5
```

## NL

```
4
```

## NMPVFIXB

```
9
```

## NMPVFLOATB

```
NOTDEFINED
```

## NKP

```
3
```

## NKCP

```
fourbarprisDE1
```

## KPDEF

```
pkin=[GK,GP,HP]';
```

## NTAUJFIXBREGNN

```
9
```

## VERSIONINFO

```
% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
```

## INPUT_M

```
% m [4x1]
%   mass of all robot links (including the base)
```

## INPUT_R

```
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
```

## INPUT_MR

```
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
```

## INPUT_IC

```
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
```

## INPUT_IF

```
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
```

## INPUT_MDPFIXB

```
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisDE1_convert_par2_MPV_fixb.m
```

## INPUT_PKIN

```
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
```

## INPUT_QJ

```
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
```

## INPUT_QJD

```
% qJD [1x1]
%   Generalized joint velocities
```

## INPUT_QJDD

```
% qJDD [1x1]
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

