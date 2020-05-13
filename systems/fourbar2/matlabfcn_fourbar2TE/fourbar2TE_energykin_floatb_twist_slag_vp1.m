% Calculate kinetic energy for
% fourbar2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:17
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar2TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2TE_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:17:23
% EndTime: 2020-04-24 20:17:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (131->81), mult. (194->100), div. (0->0), fcn. (84->2), ass. (0->32)
t58 = Icges(2,4) + Icges(4,4);
t55 = Icges(4,2) / 0.2e1 + Icges(2,2) / 0.2e1;
t54 = Icges(4,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t49 = Icges(2,5) + Icges(4,5);
t48 = Icges(2,6) + Icges(4,6);
t38 = cos(qJ(1));
t47 = t58 * t38 / 0.2e1;
t37 = sin(qJ(1));
t46 = -t58 * t37 / 0.2e1;
t34 = V_base(6) + qJD(1);
t45 = pkin(2) * t34;
t42 = t55 * t38 - t46;
t41 = t55 * t37 - t47;
t40 = t54 * t37 + t47;
t39 = t54 * t38 + t46;
t33 = t38 * rSges(2,1) - t37 * rSges(2,2);
t32 = t38 * rSges(4,1) - t37 * rSges(4,2);
t31 = t37 * rSges(2,1) + t38 * rSges(2,2);
t30 = t37 * rSges(4,1) + t38 * rSges(4,2);
t17 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t16 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t15 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t14 = -V_base(6) * rSges(3,2) + V_base(5) * rSges(3,3) - t37 * t45 + V_base(1);
t13 = V_base(6) * rSges(3,1) - V_base(4) * rSges(3,3) + t38 * t45 + V_base(2);
t12 = V_base(5) * rSges(2,3) - t34 * t31 + V_base(1);
t11 = V_base(5) * rSges(4,3) - t34 * t30 + V_base(1);
t10 = -V_base(4) * rSges(2,3) + t34 * t33 + V_base(2);
t9 = -V_base(5) * rSges(3,1) + V_base(4) * rSges(3,2) + V_base(3) + (t37 * V_base(4) - t38 * V_base(5)) * pkin(2);
t8 = -V_base(4) * rSges(4,3) + V_base(6) * pkin(1) + t34 * t32 + V_base(2);
t7 = V_base(4) * t31 - V_base(5) * t33 + V_base(3);
t6 = V_base(4) * t30 + V_base(3) + (-pkin(1) - t32) * V_base(5);
t1 = m(1) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + m(2) * (t10 ^ 2 + t12 ^ 2 + t7 ^ 2) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t9 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t6 ^ 2 + t8 ^ 2) / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(3,3) / 0.2e1) * V_base(6) ^ 2 + (Icges(2,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t34 ^ 2 + ((Icges(1,2) / 0.2e1 + Icges(3,2) / 0.2e1 + t42 * t38 + t40 * t37) * V_base(5) + (Icges(1,6) + Icges(3,6)) * V_base(6) + (t49 * t37 + t48 * t38) * t34) * V_base(5) + ((Icges(1,1) / 0.2e1 + Icges(3,1) / 0.2e1 + t39 * t38 + t41 * t37) * V_base(4) + (Icges(1,5) + Icges(3,5)) * V_base(6) + (Icges(1,4) + Icges(3,4) + (t40 - t41) * t38 + (t39 - t42) * t37) * V_base(5) + (-t48 * t37 + t49 * t38) * t34) * V_base(4);
T = t1;
