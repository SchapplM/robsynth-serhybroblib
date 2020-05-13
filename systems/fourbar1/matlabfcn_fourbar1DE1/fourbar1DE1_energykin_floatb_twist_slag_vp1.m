% Calculate kinetic energy for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
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
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1DE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:56:48
% EndTime: 2020-04-24 19:56:49
% DurationCPUTime: 1.24s
% Computational Cost: add. (2099->137), mult. (2886->232), div. (156->4), fcn. (780->4), ass. (0->97)
t122 = cos(qJ(1));
t153 = pkin(2) * t122;
t146 = (-0.2e1 * t153 + pkin(1)) * pkin(1);
t138 = pkin(2) ^ 2 + t146;
t145 = pkin(3) ^ 2 - pkin(4) ^ 2;
t103 = t138 + t145;
t121 = sin(qJ(1));
t147 = t121 * t103;
t117 = -pkin(1) + t153;
t160 = -pkin(3) + pkin(4);
t161 = -pkin(3) - pkin(4);
t129 = sqrt(-((pkin(2) - t161) * (pkin(2) + t161) + t146) * ((pkin(2) - t160) * (pkin(2) + t160) + t146));
t148 = t117 * t129;
t85 = pkin(2) * t147 + t148;
t165 = t85 / 0.2e1;
t104 = t138 - t145;
t152 = t121 * pkin(2);
t86 = -t104 * t152 + t148;
t164 = t86 / 0.2e1;
t114 = 0.1e1 / t138;
t126 = 0.1e1 / pkin(3);
t149 = t114 * t126;
t134 = t149 / 0.2e1;
t173 = t134 * t85;
t124 = 0.1e1 / pkin(4);
t150 = t114 * t124;
t135 = t150 / 0.2e1;
t172 = t135 * t86;
t171 = Icges(3,4) * t165;
t170 = Icges(4,4) * t164;
t119 = V_base(6) + qJD(1);
t139 = qJD(1) * t114 / t129;
t75 = pkin(1) * t139 * t86 + t119;
t169 = t149 * t75;
t76 = V_base(6) - pkin(2) * ((pkin(1) * t122 - pkin(2)) * t129 - pkin(1) * t147) * t139;
t168 = t150 * t76;
t91 = t129 * t152;
t83 = -t117 * t103 + t91;
t167 = t83 / 0.2e1;
t84 = t117 * t104 + t91;
t166 = t84 / 0.2e1;
t163 = -t85 / 0.2e1;
t162 = -t86 / 0.2e1;
t157 = Icges(3,5) / 0.2e1;
t156 = Icges(4,5) / 0.2e1;
t155 = Icges(3,6) / 0.2e1;
t154 = Icges(4,6) / 0.2e1;
t151 = Icges(2,4) * t121;
t142 = V_base(4) / 0.2e1;
t141 = -V_base(5) / 0.2e1;
t137 = t114 * V_base(4);
t136 = t114 * V_base(5);
t133 = t124 * t137;
t132 = t124 * t136;
t131 = t126 * t137;
t130 = t126 * t136;
t120 = Icges(2,4) * t122;
t116 = t122 * rSges(2,1) - t121 * rSges(2,2);
t115 = t121 * rSges(2,1) + t122 * rSges(2,2);
t113 = Icges(2,1) * t122 - t151;
t112 = Icges(2,1) * t121 + t120;
t111 = -Icges(2,2) * t121 + t120;
t110 = Icges(2,2) * t122 + t151;
t109 = Icges(2,5) * t122 - Icges(2,6) * t121;
t107 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t106 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t105 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t102 = -t117 * rSges(3,1) + rSges(3,2) * t152;
t101 = rSges(3,1) * t152 + t117 * rSges(3,2);
t100 = rSges(4,1) * t152 + t117 * rSges(4,2);
t99 = -t117 * rSges(4,1) + rSges(4,2) * t152;
t96 = V_base(5) * rSges(2,3) - t119 * t115 + V_base(1);
t95 = -V_base(4) * rSges(2,3) + t119 * t116 + V_base(2);
t94 = V_base(4) * t115 - V_base(5) * t116 + V_base(3);
t82 = Icges(3,4) * t83 * t134;
t81 = Icges(4,4) * t84 * t135;
t80 = -t103 * t101 + t102 * t129;
t79 = t101 * t129 + t103 * t102;
t78 = t100 * t129 - t104 * t99;
t77 = t104 * t100 + t129 * t99;
t74 = -Icges(3,1) * t173 + t82;
t73 = (Icges(3,1) * t167 + t171) * t149;
t72 = -Icges(4,1) * t172 + t81;
t71 = (Icges(4,1) * t166 + t170) * t150;
t70 = (Icges(3,2) * t167 - t171) * t149;
t69 = Icges(3,2) * t173 + t82;
t68 = (Icges(4,2) * t166 - t170) * t150;
t67 = Icges(4,2) * t172 + t81;
t65 = (t85 * t155 + t83 * t157) * t149;
t63 = (t86 * t154 + t84 * t156) * t150;
t62 = -V_base(5) * pkin(1) + V_base(3) + (t78 * t141 + t77 * t142) * t150;
t61 = V_base(3) + (t79 * t141 + t80 * t142) * t149 + (V_base(4) * t121 - V_base(5) * t122) * pkin(2);
t60 = V_base(1) + V_base(5) * rSges(4,3) - t77 * t168 / 0.2e1;
t59 = t76 * t78 * t135 - V_base(4) * rSges(4,3) + V_base(6) * pkin(1) + V_base(2);
t58 = V_base(1) - t119 * t152 + V_base(5) * rSges(3,3) - t80 * t169 / 0.2e1;
t57 = t75 * t79 * t134 - V_base(4) * rSges(3,3) + t119 * t153 + V_base(2);
t1 = m(1) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + V_base(6) * (Icges(1,5) * V_base(4) + Icges(1,3) * V_base(6)) / 0.2e1 + m(2) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t119 * (Icges(2,3) * t119 + t109 * V_base(4)) / 0.2e1 + m(3) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + t75 * (Icges(3,3) * t75 + t65 * V_base(4)) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t62 ^ 2) / 0.2e1 + t76 * (Icges(4,3) * t76 + t63 * V_base(4)) / 0.2e1 + (Icges(1,5) * V_base(6) + t109 * t119 + (t69 * t165 + t73 * t167) * t131 + (t70 * t165 + t74 * t167) * t130 + t65 * t75 + (t67 * t164 + t71 * t166) * t133 + (t68 * t164 + t72 * t166) * t132 + t63 * t76 + (-t121 * t110 + t122 * t112 + Icges(1,4)) * V_base(5) + (-t121 * t111 + t122 * t113 + Icges(1,1)) * V_base(4)) * t142 + ((t73 * t163 + t69 * t167) * t131 + (t74 * t163 + t70 * t167) * t130 + (t71 * t162 + t67 * t166) * t133 + (t72 * t162 + t68 * t166) * t132 + (t122 * t110 + t121 * t112 + Icges(1,2)) * V_base(5) + (t122 * t111 + t121 * t113 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * (t83 * t155 - t157 * t85) * t169 + V_base(5) * (t84 * t154 - t156 * t86) * t168 + V_base(5) * t119 * (Icges(2,5) * t121 + Icges(2,6) * t122) + V_base(5) * V_base(6) * Icges(1,6);
T = t1;
