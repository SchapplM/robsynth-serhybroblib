% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
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
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1DE1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1DE1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:10
% EndTime: 2020-04-24 19:57:15
% DurationCPUTime: 2.16s
% Computational Cost: add. (5517->124), mult. (7779->279), div. (189->9), fcn. (2037->4), ass. (0->128)
t69 = cos(qJ(1));
t158 = pkin(2) * t69;
t78 = pkin(1) ^ 2;
t132 = -0.2e1 * pkin(1) * t158 + t78;
t77 = pkin(2) ^ 2;
t59 = t77 + t132;
t52 = 0.1e1 / t59;
t75 = 0.1e1 / pkin(3);
t135 = t52 * t75;
t118 = t135 / 0.2e1;
t159 = pkin(1) * t77;
t174 = 0.2e1 * t159;
t162 = t52 / 0.2e1;
t68 = sin(qJ(1));
t154 = t68 * pkin(2);
t173 = rSges(4,1) * t68 + rSges(4,2) * t69;
t172 = rSges(3,1) * t68 + rSges(3,2) * t69;
t171 = pkin(4) ^ 2;
t81 = t59 ^ 2;
t53 = 0.1e1 / t81;
t170 = 0.2e1 * t53;
t161 = -pkin(3) - pkin(4);
t35 = (pkin(2) - t161) * (pkin(2) + t161) + t132;
t160 = -pkin(3) + pkin(4);
t36 = (pkin(2) - t160) * (pkin(2) + t160) + t132;
t136 = t35 * t36;
t137 = 0.1e1 / t136;
t79 = sqrt(-t136);
t29 = 0.1e1 / t79;
t123 = pkin(2) * t29 * t53;
t131 = pkin(3) ^ 2 - t171;
t46 = t59 + t131;
t61 = pkin(1) * t69 - pkin(2);
t25 = -pkin(1) * t46 * t68 + t61 * t79;
t116 = t25 * t123;
t105 = qJD(1) * t116;
t100 = -t105 / 0.2e1;
t70 = qJD(1) ^ 2;
t169 = t70 * t154;
t168 = 0.1e1 / t81 ^ 2 / t171;
t117 = t29 * t68 * t159;
t109 = t25 * t117;
t140 = t29 * t52;
t124 = pkin(2) * t140;
t138 = t29 * t61;
t125 = pkin(1) * t154;
t27 = (-t35 - t36) * t125;
t30 = t29 * t137;
t145 = t27 * t30;
t133 = t68 * t79;
t67 = t68 ^ 2;
t86 = -0.2e1 * t78 * t67 * pkin(2) + (-t46 * t69 - t133) * pkin(1);
t4 = (-(t27 * t138 + t86) * t124 + (-pkin(2) * t52 * t145 + t117 * t170) * t25) * qJD(1);
t54 = t52 * t53;
t167 = qJD(1) * t54 * t109 + t4 * t162;
t165 = -t168 / 0.2e1;
t47 = t59 - t131;
t164 = pkin(2) * (t47 * t69 + t133) + t67 * t174;
t72 = 0.1e1 / pkin(4);
t163 = t72 ^ 2;
t60 = pkin(1) - t158;
t139 = t29 * t60;
t26 = qJD(1) * t27;
t146 = t26 * t30;
t24 = t47 * t154 + t60 * t79;
t91 = (-(qJD(1) * t164 + t26 * t139) * t29 - t24 * t146) * t52 * pkin(1);
t98 = 0.2e1 * t24 * t68 * t78 * t123;
t2 = qJD(1) * t91 + t70 * t98;
t157 = t2 * t52;
t50 = t172 * pkin(2);
t43 = -rSges(3,2) * pkin(1) + t50;
t150 = rSges(3,2) * t68;
t44 = rSges(3,1) * t60 + pkin(2) * t150;
t21 = t43 * t79 + t44 * t46;
t113 = t21 * t118;
t129 = qJD(1) * pkin(2);
t120 = t69 * t129;
t16 = -pkin(1) * t24 * t140 + 0.1e1;
t15 = t16 * qJD(1);
t6 = t113 * t15 + t120;
t155 = t6 * t21;
t22 = -t43 * t46 + t44 * t79;
t119 = -t135 / 0.2e1;
t112 = t22 * t119;
t121 = t68 * t129;
t7 = t112 * t15 - t121;
t153 = t7 * t22;
t148 = rSges(4,2) * t68;
t41 = rSges(4,1) * t60 + pkin(2) * t148;
t144 = t29 * t41;
t51 = t173 * pkin(2);
t42 = -rSges(4,2) * pkin(1) + t51;
t143 = t29 * t42;
t142 = t29 * t43;
t141 = t29 * t44;
t127 = qJD(1) * t52;
t126 = qJD(1) * t77;
t115 = pkin(1) * t121;
t114 = t126 * t137;
t111 = -0.2e1 * t115;
t110 = 0.2e1 * t115;
t106 = t53 * t115;
t102 = t70 * t109;
t23 = t25 ^ 2;
t101 = t23 * t26 / t35 ^ 2 / t36 ^ 2 * t126;
t99 = t105 / 0.2e1;
t49 = (rSges(3,1) * t69 - t150) * pkin(2);
t48 = (rSges(4,1) * t69 - t148) * pkin(2);
t97 = t16 * t75 * t106;
t94 = t23 * t137 * t169 * t174;
t93 = t114 * t165;
t13 = qJD(1) * t86 + t26 * t138;
t1 = t102 * t170 + (-t13 * t29 - t25 * t146) * pkin(2) * t127;
t87 = t1 * t162 + t54 * t102;
t40 = t173 * t129;
t39 = t172 * t129;
t38 = qJD(1) * t49;
t37 = qJD(1) * t48;
t20 = -t41 * t47 + t42 * t79;
t19 = t41 * t79 + t42 * t47;
t14 = qJD(1) * t98;
t11 = t110 * t42 + t26 * t144 + t47 * t37 + t40 * t79;
t10 = t111 * t43 + t26 * t141 - t46 * t38 + t39 * t79;
t9 = t110 * t44 + t26 * t142 + t38 * t79 + t46 * t39;
t8 = t111 * t41 + t26 * t143 + t37 * t79 - t47 * t40;
t5 = t14 + (-t24 * t145 + (-t27 * t139 - t164) * t29) * pkin(1) * t127;
t3 = t14 + t91;
t12 = [m(4) * ((t19 * t11 + t20 * t8) * t23 * t93 + ((t11 * t99 - t19 * t87) * t19 - (t100 * t8 + t20 * t87) * t20) * t163 * t116 + (t19 ^ 2 + t20 ^ 2) * (t25 * t13 * t93 + t52 * t94 * t168 + t101 * t165)) / 0.2e1 + m(4) * (-t20 * ((-0.2e1 * t41 * t125 + t27 * t143 - t47 * t51 + t48 * t79) * t100 + t167 * t20) + ((0.2e1 * t42 * t125 + t27 * t144 + t47 * t48 + t51 * t79) * t99 - t167 * t19) * t19) * t163 * t100 + ((-t70 * t158 + (-t22 * t157 / 0.2e1 + (-t10 * t52 / 0.2e1 + t22 * t106) * t15) * t75) * (t112 * t16 - t154) + t7 * (t16 * t10 * t119 + t22 * t97 - t120) + (-t169 + (t21 * t157 / 0.2e1 + (-t21 * t106 + t9 * t162) * t15) * t75) * (t113 * t16 + t158) + t6 * (t16 * t9 * t118 - t21 * t97 - t121) + (-t153 + t155) * t3 * t118 - (-t6 * t68 - t69 * t7) * t129 - ((-t153 / 0.2e1 + t155 / 0.2e1) * t52 * t5 + ((-t7 * (t27 * t141 - t46 * t49 + t50 * t79) / 0.2e1 + t6 * (t27 * t142 + t46 * t50 + t49 * t79) / 0.2e1) * t52 + (t7 * (t22 * t53 + t43 * t52) + t6 * (-t21 * t53 + t44 * t52)) * t125) * t15) * t75) * m(3) + (t2 * t16 + (t3 - t5) * t15) * Icges(3,3) + (-t53 * t101 + t54 * t94 + (-t53 * t13 * t114 + (t4 * qJD(1) - t1) * t124) * t25) * Icges(4,3);];
tauc = t12(:);
