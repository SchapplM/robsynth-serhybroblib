% Calculate joint inertia matrix for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:07:24
% EndTime: 2020-05-03 00:07:25
% DurationCPUTime: 0.92s
% Computational Cost: add. (696->214), mult. (988->282), div. (0->0), fcn. (275->42), ass. (0->133)
t107 = qJ(4) + qJ(5);
t91 = qJ(3) + t107;
t161 = -qJ(5) + qJ(3);
t92 = qJ(4) + t161;
t197 = cos(t91) + cos(t92);
t133 = pkin(3) ^ 2;
t117 = cos(qJ(5));
t128 = 0.2e1 * qJ(5);
t180 = (rSges(6,1) * pkin(4));
t152 = 2 * m(6) * t180;
t112 = sin(qJ(5));
t168 = t112 * rSges(6,2);
t186 = pkin(4) * m(6);
t195 = -0.2e1 * t168 * t186 + Icges(5,3) + (rSges(5,1) ^ 2 + rSges(5,2) ^ 2) * m(5);
t130 = rSges(6,2) ^ 2;
t131 = rSges(6,1) ^ 2;
t46 = m(6) * (-t130 + t131) - Icges(6,1) + Icges(6,2);
t179 = rSges(6,2) * m(6);
t73 = rSges(6,1) * t179 - Icges(6,4);
t139 = t46 * cos(t128) / 0.2e1 - t73 * sin(t128) + t117 * t152 + Icges(6,2) / 0.2e1 + Icges(6,1) / 0.2e1 + t195;
t196 = t139 + Icges(4,3) + t133 * m(5) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4);
t106 = -qJ(5) + qJ(2);
t153 = t179 / 0.2e1;
t148 = pkin(2) * t153;
t183 = m(6) * rSges(6,1);
t155 = t183 / 0.2e1;
t149 = pkin(2) * t155;
t194 = sin(t106) * t148 + cos(t106) * t149;
t178 = rSges(6,2) * pkin(4);
t104 = 0.2e1 * t178;
t129 = 2 * Icges(6,6);
t175 = rSges(5,3) * m(5);
t122 = rSges(6,3) + pkin(6);
t170 = rSges(6,2) * t122;
t181 = m(6) * (-t170 + t180);
t182 = m(6) * (t170 + t180);
t77 = qJ(2) + t91;
t55 = sin(t77);
t78 = qJ(2) + t92;
t56 = sin(t78);
t57 = cos(t77);
t58 = cos(t78);
t109 = qJ(3) + qJ(4);
t95 = qJ(2) + t109;
t64 = sin(t95);
t68 = cos(t95);
t76 = t128 + t95;
t79 = -0.2e1 * qJ(5) + t95;
t96 = t122 * rSges(6,1);
t193 = (t129 - 0.2e1 * t182) * t57 / 0.4e1 + ((t104 + 0.2e1 * t96) * m(6) - (2 * Icges(6,5))) * t56 / 0.4e1 + ((t104 - 0.2e1 * t96) * m(6) + (2 * Icges(6,5))) * t55 / 0.4e1 + (rSges(5,2) * t175 - Icges(5,6)) * t68 + (t129 + 0.2e1 * t181) * t58 / 0.4e1 + (cos(t79) / 0.4e1 - cos(t76) / 0.4e1) * t46 + (sin(t76) + sin(t79)) * t73 / 0.2e1 + t64 * (rSges(5,1) * t175 - Icges(5,5));
t118 = cos(qJ(4));
t159 = 0.2e1 * pkin(3) * t118;
t108 = qJ(4) - qJ(5);
t48 = pkin(3) * sin(t108) * t179;
t160 = pkin(3) * t183;
t49 = cos(t107) * t160;
t50 = cos(t108) * t160;
t184 = m(5) * rSges(5,1);
t69 = t184 + t186;
t191 = t69 * t159 + t48 + t49 + t50;
t154 = -t179 / 0.2e1;
t28 = t69 * pkin(2) * cos(t109);
t61 = sin(t91);
t62 = sin(t92);
t190 = pkin(2) * t61 * t154 + t62 * t148 + t197 * t149 + t28;
t40 = t117 * rSges(6,1) - t168;
t173 = pkin(4) + t40;
t189 = t122 * t64 + t173 * t68;
t54 = rSges(5,2) * m(5) - t122 * m(6);
t188 = -0.2e1 * t54;
t187 = m(6) / 0.2e1;
t177 = rSges(3,3) * m(3);
t176 = rSges(4,3) * m(4);
t83 = sin(t109);
t172 = t54 * t83;
t147 = pkin(3) * t153;
t93 = qJ(2) + t161;
t171 = sin(t93) * t147 + pkin(3) * cos(t93) * t155;
t169 = t112 * rSges(6,1);
t167 = t117 * rSges(6,2);
t119 = cos(qJ(3));
t166 = (rSges(4,1) * m(4) + (m(5) + m(6)) * pkin(3)) * t119;
t113 = sin(qJ(4));
t165 = t54 * t113;
t162 = t96 * m(6) - Icges(6,5);
t110 = qJ(2) + qJ(5);
t97 = t130 + t131;
t158 = 0.2e1 * t68;
t81 = sin(t107);
t157 = t81 * t179;
t114 = sin(qJ(3));
t156 = m(4) * rSges(4,2) * t114;
t151 = m(6) * t97 + Icges(6,3);
t150 = t162 * t112;
t105 = t122 ^ 2;
t132 = pkin(4) ^ 2;
t146 = 0.2e1 * t105 + (2 * t132) + t97;
t145 = rSges(5,1) * t68 - rSges(5,2) * t64;
t143 = 0.2e1 * t133 + t146;
t120 = cos(qJ(2));
t111 = qJ(2) + qJ(3);
t85 = sin(t111);
t90 = cos(t111);
t142 = rSges(4,1) * t90 - rSges(4,2) * t85 + t120 * pkin(2) + pkin(1);
t115 = sin(qJ(2));
t141 = rSges(3,1) * t120 - rSges(3,2) * t115 + pkin(1);
t140 = pkin(2) * t114 * t118 + t113 * (t119 * pkin(2) + pkin(3));
t137 = t146 * t187 + t139;
t136 = -0.2e1 * pkin(3) * t165 + t143 * t187 + t196;
t94 = qJ(3) + t110;
t63 = sin(t94);
t67 = cos(t94);
t135 = (rSges(4,2) * t176 - Icges(4,6)) * t90 + (rSges(4,1) * t176 + pkin(3) * t175 - Icges(4,5)) * t85 + t171 + t63 * t147 - t67 * t160 / 0.2e1 + t193;
t134 = pkin(2) ^ 2;
t121 = cos(qJ(1));
t116 = sin(qJ(1));
t89 = cos(t110);
t84 = sin(t110);
t41 = t121 * rSges(2,1) - t116 * rSges(2,2);
t39 = t167 + t169;
t38 = -t116 * rSges(2,1) - t121 * rSges(2,2);
t20 = t173 * m(6) + t184;
t17 = (t119 * pkin(3) + pkin(2)) * t120 - t115 * t114 * pkin(3) + pkin(1);
t13 = t17 * t121;
t8 = -t116 * rSges(3,3) + t141 * t121;
t7 = -t121 * rSges(3,3) - t141 * t116;
t6 = -t116 * rSges(4,3) + t142 * t121;
t5 = -t121 * rSges(4,3) - t142 * t116;
t4 = -t116 * rSges(5,3) + t145 * t121 + t13;
t3 = -t121 * rSges(5,3) + (-t145 - t17) * t116;
t2 = -t116 * t39 + t189 * t121 + t13;
t1 = -t121 * t39 + (-t17 - t189) * t116;
t9 = [t120 ^ 2 * Icges(3,2) + t90 ^ 2 * Icges(4,2) + Icges(2,3) + (Icges(4,1) * t85 + 0.2e1 * Icges(4,4) * t90) * t85 + (Icges(3,1) * t115 + 0.2e1 * Icges(3,4) * t120) * t115 + (Icges(6,3) + Icges(5,2)) * t68 ^ 2 + ((t117 ^ 2 * Icges(6,1) + Icges(5,1) + (-0.2e1 * Icges(6,4) * t117 + Icges(6,2) * t112) * t112) * t64 + (-Icges(6,5) * t117 + Icges(6,6) * t112 + Icges(5,4)) * t158) * t64 + m(3) * (t7 ^ 2 + t8 ^ 2) + m(2) * (t38 ^ 2 + t41 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2); t120 * (rSges(3,2) * t177 - Icges(3,6)) + t135 + ((t175 + t176) * pkin(2) + rSges(3,1) * t177 - Icges(3,5)) * t115 + (-rSges(6,1) * t89 / 0.2e1 + rSges(6,2) * t84 / 0.2e1) * pkin(2) * m(6) + t194; (0.2e1 * t134 + t143) * t187 + 0.2e1 * t28 + (m(4) + m(5)) * t134 + (t113 * t188 - t157) * pkin(3) + (-0.2e1 * t156 + 0.2e1 * t166 + t83 * t188 + ((-t61 + t62) * rSges(6,2) + t197 * rSges(6,1)) * m(6)) * pkin(2) + Icges(3,3) + (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) + t191 + t196; t135; (-t156 + t166 - t172) * pkin(2) - pkin(3) * t157 + t136 + t190 + t191; t20 * t159 + t136; t193; t48 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 - pkin(2) * t172 + (t69 * t118 + t81 * t154 - t165) * pkin(3) + t137 + t190; (t118 * t20 - t165) * pkin(3) + t137; (t105 + t130 + t132) * m(6) + Icges(6,1) + (-0.2e1 * t112 * t73 + t46 * t117 + t152) * t117 + t195; (Icges(6,6) + t181) * t58 / 0.2e1 + ((t96 + t178) * m(6) - Icges(6,5)) * t56 / 0.2e1 + (-Icges(6,6) + t182) * t57 / 0.2e1 + ((t96 - t178) * m(6) - Icges(6,5)) * t55 / 0.2e1 + t151 * t158 / 0.2e1 + m(6) * (pkin(1) * t40 + (-pkin(3) * t63 / 0.2e1 - pkin(2) * t84 / 0.2e1) * rSges(6,2) + (pkin(3) * t67 / 0.2e1 + pkin(2) * t89 / 0.2e1) * rSges(6,1)) + t171 + t194; Icges(6,6) * t117 - t150 + (-t140 * t169 + (-t122 - t140) * t167) * m(6); (Icges(6,6) + (-pkin(3) * t113 - t122) * t179) * t117 - t112 * (t113 * t160 + t162); (-m(6) * t170 + Icges(6,6)) * t117 - t150; t151;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
