% Calculate joint inertia matrix for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnIC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnIC_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnIC_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:33
% EndTime: 2020-05-07 11:33:34
% DurationCPUTime: 1.02s
% Computational Cost: add. (975->141), mult. (1487->240), div. (16->5), fcn. (1319->11), ass. (0->83)
t112 = sin(qJ(1));
t115 = cos(qJ(1));
t163 = t112 * t115;
t111 = sin(qJ(2));
t114 = cos(qJ(2));
t162 = rSges(3,1) * t114 - rSges(3,2) * t111;
t106 = t112 ^ 2;
t107 = t115 ^ 2;
t161 = t112 / 0.2e1;
t160 = -t115 / 0.2e1;
t159 = pkin(2) * t114;
t108 = qJ(2) + qJ(3);
t101 = sin(t108);
t154 = rSges(4,2) * t101;
t102 = cos(t108);
t157 = rSges(4,1) * t102;
t135 = t154 - t157;
t149 = t112 * rSges(4,3) + t115 * t154;
t152 = rSges(4,3) * t115;
t46 = t112 * (t135 * t112 - t152) + t115 * (-t115 * t157 + t149);
t113 = cos(qJ(4));
t156 = rSges(5,1) * t113;
t109 = sin(qJ(4));
t153 = rSges(5,2) * t109;
t151 = t115 * rSges(5,3);
t150 = t112 * rSges(5,3) + t115 * t156;
t148 = Icges(3,4) * t111;
t147 = Icges(3,4) * t114;
t146 = Icges(4,4) * t101;
t145 = Icges(4,4) * t102;
t144 = Icges(5,4) * t109;
t143 = Icges(5,4) * t113;
t142 = -qJ(4) + qJ(2);
t141 = t106 + t107;
t100 = sin(qJ(3) + t142);
t99 = 0.1e1 / t100;
t140 = (-pkin(3) * t100 + pkin(2) * sin(t142)) / pkin(3) * t99;
t83 = -rSges(4,1) * t101 - rSges(4,2) * t102;
t139 = -pkin(2) * t111 - t83;
t118 = -Icges(4,5) * t102 + Icges(4,6) * t101;
t55 = -Icges(4,3) * t115 + t118 * t112;
t56 = Icges(4,3) * t112 + t118 * t115;
t138 = -t115 * (t107 * t55 - t56 * t163) + t112 * (t106 * t56 - t55 * t163);
t121 = Icges(4,2) * t101 - t145;
t124 = -Icges(4,1) * t102 + t146;
t81 = -Icges(4,2) * t102 - t146;
t82 = -Icges(4,1) * t101 - t145;
t131 = t101 * t81 - t102 * t82;
t80 = -Icges(4,5) * t101 - Icges(4,6) * t102;
t137 = (-t101 * (Icges(4,5) * t112 + t124 * t115) - t102 * (Icges(4,6) * t112 + t121 * t115) + t112 * t80 + t131 * t115) * t161 + (-t101 * (-Icges(4,5) * t115 + t124 * t112) - t102 * (-Icges(4,6) * t115 + t121 * t112) + t131 * t112 - t115 * t80) * t160;
t134 = -t153 + t156;
t89 = Icges(3,2) * t114 + t148;
t91 = Icges(3,1) * t111 + t147;
t126 = -t111 * t89 + t114 * t91;
t125 = Icges(3,1) * t114 - t148;
t123 = Icges(5,1) * t113 - t144;
t122 = -Icges(3,2) * t111 + t147;
t120 = -Icges(5,2) * t109 + t143;
t119 = Icges(3,5) * t114 - Icges(3,6) * t111;
t117 = Icges(5,5) * t113 - Icges(5,6) * t109;
t110 = sin(qJ(3));
t95 = rSges(2,1) * t115 - rSges(2,2) * t112;
t94 = -rSges(2,1) * t112 - rSges(2,2) * t115;
t93 = rSges(3,1) * t111 + rSges(3,2) * t114;
t92 = rSges(5,1) * t109 + rSges(5,2) * t113;
t87 = Icges(3,5) * t111 + Icges(3,6) * t114;
t77 = rSges(3,3) * t112 + t162 * t115;
t75 = -t115 * rSges(3,3) + t162 * t112;
t66 = Icges(3,3) * t112 + t119 * t115;
t65 = -Icges(3,3) * t115 + t119 * t112;
t64 = Icges(5,3) * t112 + t117 * t115;
t63 = -Icges(5,3) * t115 + t117 * t112;
t62 = (pkin(1) - t153) * t115 + t150;
t61 = t151 + (-pkin(1) - t134) * t112;
t54 = t139 * t115;
t53 = t139 * t112;
t50 = (-t157 + t159) * t115 + t149;
t49 = t152 + (-t135 - t159) * t112;
t48 = t112 * t75 + t115 * t77;
t47 = t115 * (-t115 * t153 + t150) + (t134 * t112 - t151) * t112;
t41 = t141 * t159 + t46;
t37 = (t111 * (Icges(3,5) * t112 + t125 * t115) + t114 * (Icges(3,6) * t112 + t122 * t115)) * t161 + (t112 * t87 + t115 * t126) * t161 + (t111 * (-Icges(3,5) * t115 + t125 * t112) + t114 * (-Icges(3,6) * t115 + t122 * t112)) * t160 + (t112 * t126 - t115 * t87) * t160 + m(3) * (-t112 * t77 + t115 * t75) * t93 + m(4) * (t49 * t54 + t50 * t53) + (m(4) * (-t112 * t50 - t115 * t49) * t83 + t137) * t140 + t110 * pkin(2) / pkin(4) * t99 * ((t109 * (Icges(5,5) * t112 + t123 * t115) + t113 * (Icges(5,6) * t112 + t120 * t115)) * t161 + (t109 * (-Icges(5,5) * t115 + t123 * t112) + t113 * (-Icges(5,6) * t115 + t120 * t112)) * t160 + m(5) * (-t112 * t62 - t115 * t61) * t92 + (t106 / 0.2e1 + t107 / 0.2e1) * (Icges(5,5) * t109 + Icges(5,6) * t113)) + t137;
t1 = [-t101 * t82 - t102 * t81 + t109 * (Icges(5,1) * t109 + t143) + t111 * t91 + t113 * (Icges(5,2) * t113 + t144) + t114 * t89 + Icges(2,3) + m(2) * (t94 ^ 2 + t95 ^ 2) + m(3) * (t75 ^ 2 + t77 ^ 2) + m(4) * (t49 ^ 2 + t50 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2), t37; t37, m(3) * (t141 * t93 ^ 2 + t48 ^ 2) + t112 * (t106 * t66 - t65 * t163) - t115 * (t107 * t65 - t66 * t163) + m(4) * (t41 ^ 2 + t53 ^ 2 + t54 ^ 2) + t110 ^ 2 * pkin(2) ^ 2 / pkin(4) ^ 2 / t100 ^ 2 * (m(5) * (t141 * t92 ^ 2 + t47 ^ 2) + t112 * (t106 * t64 - t63 * t163) - t115 * (t107 * t63 - t64 * t163)) + t138 + (0.2e1 * t138 + 0.2e1 * m(4) * (t41 * t46 + (-t112 * t53 - t115 * t54) * t83) + (m(4) * (t141 * t83 ^ 2 + t46 ^ 2) + t138) * t140) * t140;];
Mq = t1;
