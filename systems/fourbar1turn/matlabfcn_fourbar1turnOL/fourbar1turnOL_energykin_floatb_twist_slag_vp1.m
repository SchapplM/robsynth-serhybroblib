% Calculate kinetic energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnOL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:47
% DurationCPUTime: 1.07s
% Computational Cost: add. (480->191), mult. (722->295), div. (0->0), fcn. (574->8), ass. (0->94)
t98 = qJ(2) + qJ(3);
t94 = sin(t98);
t134 = Icges(4,4) * t94;
t95 = cos(t98);
t133 = Icges(4,4) * t95;
t99 = sin(qJ(4));
t132 = Icges(5,4) * t99;
t103 = cos(qJ(2));
t92 = V_base(6) + qJD(1);
t131 = t103 * t92;
t101 = sin(qJ(1));
t130 = Icges(2,4) * t101;
t100 = sin(qJ(2));
t129 = Icges(3,4) * t100;
t128 = Icges(3,4) * t103;
t102 = cos(qJ(4));
t127 = Icges(5,4) * t102;
t126 = V_base(5) * pkin(5) + V_base(1);
t91 = qJD(2) * t101 + V_base(4);
t123 = -rSges(4,1) * t95 + rSges(4,2) * t94;
t122 = rSges(5,1) * t102 - rSges(5,2) * t99;
t121 = rSges(3,1) * t103 - rSges(3,2) * t100;
t120 = -Icges(4,1) * t95 + t134;
t119 = Icges(4,2) * t94 - t133;
t118 = -Icges(4,5) * t95 + Icges(4,6) * t94;
t117 = -V_base(4) * pkin(5) + V_base(2);
t116 = Icges(5,1) * t102 - t132;
t115 = -Icges(5,2) * t99 + t127;
t114 = Icges(5,5) * t102 - Icges(5,6) * t99;
t113 = Icges(3,1) * t103 - t129;
t112 = -Icges(3,2) * t100 + t128;
t111 = Icges(3,5) * t103 - Icges(3,6) * t100;
t104 = cos(qJ(1));
t67 = V_base(5) + (-qJD(2) - qJD(3)) * t104;
t68 = qJD(3) * t101 + t91;
t110 = (-Icges(4,3) * t104 + t101 * t118) * t67 + (Icges(4,3) * t101 + t104 * t118) * t68 + (-Icges(4,5) * t94 - Icges(4,6) * t95) * t92;
t88 = -qJD(4) * t104 + V_base(5);
t90 = qJD(4) * t101 + V_base(4);
t109 = (-Icges(5,3) * t104 + t101 * t114) * t88 + (Icges(5,3) * t101 + t104 * t114) * t90 + (Icges(5,5) * t99 + Icges(5,6) * t102) * t92;
t89 = -qJD(2) * t104 + V_base(5);
t108 = (-Icges(3,3) * t104 + t101 * t111) * t89 + (Icges(3,3) * t101 + t104 * t111) * t91 + (Icges(3,5) * t100 + Icges(3,6) * t103) * t92;
t41 = -Icges(4,6) * t104 + t101 * t119;
t42 = Icges(4,6) * t101 + t104 * t119;
t43 = -Icges(4,5) * t104 + t101 * t120;
t44 = Icges(4,5) * t101 + t104 * t120;
t64 = -Icges(4,2) * t95 - t134;
t65 = -Icges(4,1) * t94 - t133;
t107 = (t42 * t94 - t44 * t95) * t68 + (t41 * t94 - t43 * t95) * t67 + (t64 * t94 - t65 * t95) * t92;
t51 = -Icges(5,6) * t104 + t101 * t115;
t52 = Icges(5,6) * t101 + t104 * t115;
t55 = -Icges(5,5) * t104 + t101 * t116;
t56 = Icges(5,5) * t101 + t104 * t116;
t76 = Icges(5,2) * t102 + t132;
t80 = Icges(5,1) * t99 + t127;
t106 = (t102 * t56 - t52 * t99) * t90 + (t102 * t55 - t51 * t99) * t88 + (t102 * t80 - t76 * t99) * t92;
t53 = -Icges(3,6) * t104 + t101 * t112;
t54 = Icges(3,6) * t101 + t104 * t112;
t57 = -Icges(3,5) * t104 + t101 * t113;
t58 = Icges(3,5) * t101 + t104 * t113;
t77 = Icges(3,2) * t103 + t129;
t81 = Icges(3,1) * t100 + t128;
t105 = (-t100 * t54 + t103 * t58) * t91 + (-t100 * t53 + t103 * t57) * t89 + (-t100 * t77 + t103 * t81) * t92;
t96 = Icges(2,4) * t104;
t87 = rSges(2,1) * t104 - rSges(2,2) * t101;
t86 = rSges(2,1) * t101 + rSges(2,2) * t104;
t85 = rSges(3,1) * t100 + rSges(3,2) * t103;
t84 = rSges(5,1) * t99 + rSges(5,2) * t102;
t83 = Icges(2,1) * t104 - t130;
t82 = Icges(2,1) * t101 + t96;
t79 = -Icges(2,2) * t101 + t96;
t78 = Icges(2,2) * t104 + t130;
t71 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t70 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t69 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t66 = -rSges(4,1) * t94 - rSges(4,2) * t95;
t62 = t101 * rSges(3,3) + t104 * t121;
t61 = t101 * rSges(5,3) + t104 * t122;
t60 = -rSges(3,3) * t104 + t101 * t121;
t59 = -rSges(5,3) * t104 + t101 * t122;
t46 = t101 * rSges(4,3) + t104 * t123;
t45 = -rSges(4,3) * t104 + t101 * t123;
t38 = V_base(5) * rSges(2,3) - t86 * t92 + t126;
t37 = t87 * t92 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t36 = t86 * V_base(4) - t87 * V_base(5) + V_base(3);
t35 = -t60 * t92 + t85 * t89 + t126;
t34 = t62 * t92 - t85 * t91 + t117;
t33 = t60 * t91 - t62 * t89 + V_base(3);
t32 = t84 * t88 + (-pkin(1) * t101 - t59) * t92 + t126;
t31 = -t90 * t84 + (pkin(1) * t104 + t61) * t92 + t117;
t30 = t90 * t59 - t88 * t61 + V_base(3) + (t101 * V_base(4) - V_base(5) * t104) * pkin(1);
t29 = -t45 * t92 + t66 * t67 + (t100 * t89 - t101 * t131) * pkin(2) + t126;
t28 = t92 * t46 - t68 * t66 + (-t100 * t91 + t104 * t131) * pkin(2) + t117;
t27 = t68 * t45 - t67 * t46 + V_base(3) + (t101 * t91 - t104 * t89) * t103 * pkin(2);
t1 = m(1) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(2) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + t91 * (t101 * t108 + t104 * t105) / 0.2e1 + t89 * (t101 * t105 - t104 * t108) / 0.2e1 + m(4) * (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) / 0.2e1 + t68 * (t101 * t110 + t104 * t107) / 0.2e1 + t67 * (t101 * t107 - t104 * t110) / 0.2e1 + m(5) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + t90 * (t101 * t109 + t104 * t106) / 0.2e1 + t88 * (t101 * t106 - t109 * t104) / 0.2e1 + ((-t101 * t78 + t104 * t82 + Icges(1,4)) * V_base(5) + (-t101 * t79 + t104 * t83 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t82 + t104 * t78 + Icges(1,2)) * V_base(5) + (t101 * t83 + t104 * t79 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t58 + t103 * t54) * t91 + (t100 * t57 + t103 * t53) * t89 + (-t42 * t95 - t44 * t94) * t68 + (-t41 * t95 - t43 * t94) * t67 + (t102 * t52 + t56 * t99) * t90 + (t102 * t51 + t55 * t99) * t88 + (t100 * t81 + t102 * t76 + t103 * t77 - t95 * t64 - t94 * t65 + t99 * t80 + Icges(2,3)) * t92) * t92 / 0.2e1 + V_base(4) * t92 * (Icges(2,5) * t104 - Icges(2,6) * t101) + V_base(5) * t92 * (Icges(2,5) * t101 + Icges(2,6) * t104) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
