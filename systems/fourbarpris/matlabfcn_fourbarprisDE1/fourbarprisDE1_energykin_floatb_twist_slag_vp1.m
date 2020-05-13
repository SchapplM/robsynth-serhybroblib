% Calculate kinetic energy for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:38
% EndTime: 2020-05-07 09:09:44
% DurationCPUTime: 3.04s
% Computational Cost: add. (1454->172), mult. (1372->291), div. (176->7), fcn. (116->2), ass. (0->104)
t186 = Icges(3,4) + Icges(2,6);
t111 = qJ(1) + pkin(3);
t159 = -pkin(2) + t111;
t160 = -pkin(2) - t111;
t142 = (pkin(1) + t160) * (pkin(1) + t159) * (pkin(1) - t159) * (pkin(1) - t160);
t122 = sqrt(-t142);
t106 = 0.1e1 / t111;
t165 = t106 / pkin(1);
t147 = -t165 / 0.2e1;
t138 = t122 * t147;
t118 = (pkin(2) ^ 2);
t120 = (pkin(1) ^ 2);
t117 = (pkin(3) ^ 2);
t114 = (qJ(1) ^ 2);
t171 = (pkin(3) * qJ(1));
t145 = t114 + 2 * t171;
t136 = -t117 - t145;
t90 = -t118 + t120 - t136;
t169 = Icges(3,5) * t90;
t173 = -t122 / 0.2e1;
t141 = t90 * t147;
t83 = Icges(2,4) * t141;
t185 = (t169 / 0.2e1 + Icges(3,3) * t173) * t165 + Icges(2,1) * t138 + t83;
t166 = Icges(2,4) * t122;
t176 = -t90 / 0.2e1;
t80 = Icges(3,5) * t138;
t184 = Icges(3,3) * t141 + (Icges(2,1) * t176 + t166 / 0.2e1) * t165 + t80;
t172 = t122 / 0.2e1;
t175 = t90 / 0.2e1;
t183 = (Icges(2,5) * t176 + Icges(3,6) * t175 + t186 * t172) * t165;
t182 = (Icges(2,5) * t173 + Icges(3,6) * t172 + t186 * t176) * t165;
t119 = t120 ^ 2;
t115 = rSges(4,2) ^ 2;
t116 = rSges(4,1) ^ 2;
t131 = (t116 / 0.4e1 - t115 / 0.4e1) * m(4) + Icges(4,2) / 0.4e1 - Icges(4,1) / 0.4e1;
t112 = pkin(2) + pkin(3);
t113 = pkin(2) - pkin(3);
t164 = t113 ^ 2 * t112 ^ 2;
t92 = (-t115 + t116) * m(4) + Icges(4,2) - Icges(4,1);
t181 = -t92 * t164 / 0.4e1 - t131 * t119;
t107 = V_base(5) ^ 2;
t108 = V_base(4) ^ 2;
t137 = V_base(4) * V_base(5);
t96 = -m(4) * rSges(4,1) * rSges(4,2) + Icges(4,4);
t180 = (-t131 * t108 + t96 * t137 + t92 * t107 / 0.4e1) * t114;
t162 = t118 + t120;
t98 = -t117 + t162;
t179 = t98 - t145;
t178 = pkin(1) * pkin(2);
t174 = m(4) * rSges(4,3);
t170 = -qJD(1) / 0.2e1;
t168 = t92 * t117;
t167 = t98 * t108;
t163 = t120 * t118;
t104 = V_base(6) * pkin(1);
t158 = t104 + V_base(2);
t156 = V_base(4) / 0.2e1;
t155 = -V_base(5) / 0.2e1;
t154 = V_base(5) / 0.2e1;
t93 = (t115 + t116) * m(4) + Icges(4,3);
t153 = t93 * t163;
t151 = 0.2e1 * rSges(4,3) * t178;
t95 = -rSges(4,1) * t174 + Icges(4,5);
t150 = t95 * V_base(5);
t149 = t96 * V_base(5);
t148 = t98 * V_base(5);
t146 = t165 / 0.2e1;
t144 = pkin(2) * t104;
t143 = qJD(1) * t111 * t178;
t140 = V_base(4) * t165;
t139 = V_base(5) * t165;
t94 = rSges(4,2) * t174 - Icges(4,6);
t135 = t94 * t144;
t134 = t95 * t144;
t133 = -V_base(5) * pkin(1) + V_base(3);
t132 = 0.1e1 / t142 / 0.8e1;
t130 = t120 + t136;
t129 = -rSges(4,1) * V_base(1) - rSges(4,2) * V_base(2);
t128 = V_base(5) * t135;
t127 = t179 * V_base(6);
t99 = -3 * t117 + t162;
t91 = -t118 - t130;
t89 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t88 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t87 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t85 = t96 * t148;
t84 = ((t116 / 0.2e1 + t115 / 0.2e1 + rSges(4,3) ^ 2) * m(4) + Icges(4,2) / 0.2e1 + Icges(4,1) / 0.2e1) * t118;
t79 = t92 * t148 - 0.2e1 * t135;
t77 = -rSges(2,1) * t122 - (rSges(2,2) * t90);
t76 = -rSges(3,1) * t122 - (rSges(3,3) * t90);
t75 = -rSges(2,1) * t90 + rSges(2,2) * t122;
t74 = rSges(3,1) * t90 - rSges(3,3) * t122;
t73 = V_base(6) + (-t118 + t130) * qJD(1) * t106 / t122;
t69 = (Icges(3,1) * t173 - t169 / 0.2e1) * t165;
t68 = Icges(3,1) * t90 * t146 + t80;
t67 = (-t166 / 0.2e1 + Icges(2,2) * t176) * t165;
t66 = Icges(2,2) * t122 * t146 + t83;
t59 = t73 * t77 * t147 + V_base(5) * rSges(2,3) + V_base(1);
t58 = t73 * t75 * t146 - V_base(4) * rSges(2,3) + t158;
t57 = (t75 * t155 + t77 * t156) * t165 + t133;
t56 = (t74 * t156 + t76 * t155 + (t90 * t154 + V_base(4) * t173) * qJ(1)) * t165 + t133;
t55 = -V_base(5) * rSges(3,2) + V_base(1) + (t90 * t170 + (qJ(1) * t172 - t74 / 0.2e1) * t73) * t165;
t54 = V_base(4) * rSges(3,2) + (t122 * t170 + (qJ(1) * t176 + t76 / 0.2e1) * t73) * t165 + t158;
t1 = m(1) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + V_base(6) * (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + (Icges(1,3) * V_base(6))) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(3) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + ((-0.8e1 * (-0.2e1 * t93 * t144 + (rSges(4,1) * V_base(2) - rSges(4,2) * V_base(1)) * t91 * m(4) + t179 * (t94 * V_base(5) - t95 * V_base(4))) * t143 - 0.2e1 * (-t96 * t167 - t79 * V_base(4) + t85 * V_base(5) + t145 * (t92 * t137 + (-t107 + t108) * t96) + 0.2e1 * (V_base(6) * t150 + ((V_base(4) * rSges(4,1) + V_base(5) * rSges(4,2)) * V_base(3) + t129 * V_base(6)) * m(4)) * t178) * t142) * t122 - 0.16e2 * t111 ^ 2 * qJD(1) ^ 2 * t153 + (t129 * m(4) + t94 * V_base(4) + t150) / t132 * t143 + 0.4e1 * (((t84 + t168 / 0.2e1) * t120 + t181) * t108 + ((-(2 * t120 * t117) + t119 + t164) * t149 + t98 * t134) * V_base(4) + ((t84 - t168 / 0.2e1) * t120 - t181) * t107 - t98 * t128 + (V_base(6) ^ 2) * t153 + ((-0.2e1 * t99 * t149 - t134) * V_base(4) + t128 + (t108 / 0.2e1 - t107 / 0.2e1) * t99 * t92 + t180) * t114 + (0.4e1 * t180 + t92 * t167 - 0.4e1 * (t85 + t134 / 0.2e1) * V_base(4) - V_base(5) * t79) * t171 + ((-t91 * (-V_base(5) * rSges(4,1) + V_base(4) * rSges(4,2)) * V_base(3) + (-rSges(4,2) * t127 + V_base(5) * t151) * V_base(1) - (-rSges(4,1) * t127 + V_base(4) * t151) * V_base(2)) * t178 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * t163) * m(4)) * t142) / pkin(1) ^ 2 / pkin(2) ^ 2 * t132 + (t182 * V_base(5) + t183 * V_base(4) + (Icges(3,2) + Icges(2,3)) * t73) * t73 / 0.2e1 + (Icges(1,1) * V_base(4) + Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + t183 * t73 + (t66 * t172 + t69 * t173 + t184 * t176) * t140 + (t67 * t172 + t68 * t173 + t185 * t176) * t139) * t156 + (Icges(1,4) * V_base(4) + Icges(1,2) * V_base(5) + Icges(1,6) * V_base(6) + t182 * t73 + (t184 * t173 + t69 * t175 + t66 * t176) * t140 + (t185 * t173 + t68 * t175 + t67 * t176) * t139) * t154;
T = t1;
