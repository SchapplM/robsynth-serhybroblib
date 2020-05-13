% Calculate kinetic energy for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1OL_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:18
% EndTime: 2020-04-24 20:10:19
% DurationCPUTime: 0.32s
% Computational Cost: add. (203->108), mult. (234->149), div. (0->0), fcn. (124->6), ass. (0->51)
t53 = qJ(1) + qJ(2);
t50 = cos(t53);
t66 = -t50 / 0.2e1;
t56 = cos(qJ(3));
t65 = t56 / 0.2e1;
t57 = cos(qJ(1));
t64 = t57 / 0.2e1;
t48 = V_base(6) + qJD(1);
t63 = pkin(2) * t48;
t55 = sin(qJ(1));
t62 = Icges(2,4) * t55;
t49 = sin(t53);
t61 = Icges(3,4) * t49;
t60 = Icges(3,4) * t50;
t54 = sin(qJ(3));
t59 = Icges(4,4) * t54;
t52 = Icges(2,4) * t57;
t51 = Icges(4,4) * t56;
t47 = V_base(6) + qJD(3);
t46 = qJD(2) + t48;
t45 = rSges(4,1) * t56 - rSges(4,2) * t54;
t44 = rSges(4,1) * t54 + rSges(4,2) * t56;
t43 = t57 * rSges(2,1) - t55 * rSges(2,2);
t42 = t55 * rSges(2,1) + t57 * rSges(2,2);
t41 = Icges(2,1) * t57 - t62;
t40 = Icges(2,1) * t55 + t52;
t39 = Icges(4,1) * t56 - t59;
t38 = Icges(4,1) * t54 + t51;
t37 = -Icges(2,2) * t55 + t52;
t36 = Icges(2,2) * t57 + t62;
t35 = -Icges(4,2) * t54 + t51;
t34 = Icges(4,2) * t56 + t59;
t29 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t28 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t27 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t26 = -t50 * rSges(3,1) + t49 * rSges(3,2);
t25 = -t49 * rSges(3,1) - t50 * rSges(3,2);
t24 = -Icges(3,1) * t50 + t61;
t23 = -Icges(3,1) * t49 - t60;
t22 = Icges(3,2) * t49 - t60;
t21 = -Icges(3,2) * t50 - t61;
t18 = V_base(5) * rSges(2,3) - t48 * t42 + V_base(1);
t17 = V_base(5) * rSges(4,3) - t47 * t44 + V_base(1);
t16 = -V_base(4) * rSges(2,3) + t48 * t43 + V_base(2);
t15 = -V_base(4) * rSges(4,3) + V_base(6) * pkin(1) + t47 * t45 + V_base(2);
t14 = V_base(4) * t42 - V_base(5) * t43 + V_base(3);
t13 = V_base(4) * t44 + V_base(3) + (-pkin(1) - t45) * V_base(5);
t12 = V_base(5) * rSges(3,3) - t46 * t25 - t55 * t63 + V_base(1);
t11 = -V_base(4) * rSges(3,3) + t46 * t26 + t57 * t63 + V_base(2);
t10 = V_base(4) * t25 - V_base(5) * t26 + V_base(3) + (t55 * V_base(4) - t57 * V_base(5)) * pkin(2);
t1 = m(1) * (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t14 ^ 2 + t16 ^ 2 + t18 ^ 2) / 0.2e1 + Icges(2,3) * t48 ^ 2 / 0.2e1 + m(3) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + Icges(3,3) * t46 ^ 2 / 0.2e1 + m(4) * (t13 ^ 2 + t15 ^ 2 + t17 ^ 2) / 0.2e1 + Icges(4,3) * t47 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (-Icges(3,5) * t49 - Icges(3,6) * t50) * t46 + (Icges(4,5) * t54 + Icges(4,6) * t56) * t47 + (Icges(2,5) * t55 + Icges(2,6) * t57) * t48 + (Icges(1,2) / 0.2e1 + t36 * t64 + t55 * t40 / 0.2e1 + t21 * t66 - t49 * t23 / 0.2e1 + t34 * t65 + t54 * t38 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t57 - Icges(2,6) * t55) * t48 + (-Icges(3,5) * t50 + Icges(3,6) * t49) * t46 + (Icges(4,5) * t56 - Icges(4,6) * t54) * t47 + ((t37 + t40) * t57 + (t35 + t38) * t56 + (-t36 + t41) * t55 + (-t34 + t39) * t54 + (-t22 - t23) * t50 + (t21 - t24) * t49) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t55 * t37 / 0.2e1 + t41 * t64 + t49 * t22 / 0.2e1 + t24 * t66 - t54 * t35 / 0.2e1 + t39 * t65) * V_base(4)) * V_base(4);
T = t1;
