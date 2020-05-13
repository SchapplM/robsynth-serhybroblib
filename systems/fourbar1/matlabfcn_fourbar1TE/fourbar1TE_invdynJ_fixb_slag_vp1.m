% Calculate vector of inverse dynamics joint torques for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1TE_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1TE_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:17
% EndTime: 2020-04-24 19:49:23
% DurationCPUTime: 2.42s
% Computational Cost: add. (6066->146), mult. (8560->310), div. (211->9), fcn. (2259->4), ass. (0->139)
t75 = cos(qJ(1));
t72 = pkin(2) * t75;
t84 = pkin(1) ^ 2;
t144 = -0.2e1 * pkin(1) * t72 + t84;
t83 = pkin(2) ^ 2;
t64 = t83 + t144;
t57 = 0.1e1 / t64;
t81 = 0.1e1 / pkin(3);
t146 = t57 * t81;
t129 = t146 / 0.2e1;
t174 = t57 / 0.2e1;
t74 = sin(qJ(1));
t185 = rSges(4,1) * t74 + rSges(4,2) * t75;
t184 = rSges(3,1) * t74 + rSges(3,2) * t75;
t183 = pkin(4) ^ 2;
t87 = t64 ^ 2;
t58 = 0.1e1 / t87;
t182 = 0.2e1 * t58;
t173 = -pkin(3) - pkin(4);
t39 = (pkin(2) - t173) * (pkin(2) + t173) + t144;
t172 = -pkin(3) + pkin(4);
t40 = (pkin(2) - t172) * (pkin(2) + t172) + t144;
t147 = t39 * t40;
t148 = 0.1e1 / t147;
t170 = pkin(1) * t74;
t143 = pkin(3) ^ 2 - t183;
t51 = t64 + t143;
t66 = pkin(1) * t75 - pkin(2);
t85 = sqrt(-t147);
t29 = -t51 * t170 + t66 * t85;
t33 = 0.1e1 / t85;
t123 = pkin(2) * t29 * t33 * t58;
t111 = qJD(1) * t123;
t106 = -t111 / 0.2e1;
t181 = -t72 + pkin(1);
t180 = 0.1e1 / t87 ^ 2 / t183;
t169 = pkin(1) * t83;
t124 = t33 * t74 * t169;
t115 = t29 * t124;
t151 = t33 * t57;
t135 = pkin(2) * t151;
t149 = t33 * t66;
t164 = t74 * pkin(2);
t136 = pkin(1) * t164;
t31 = (-t39 - t40) * t136;
t34 = t33 * t148;
t156 = t31 * t34;
t168 = pkin(2) * t57;
t145 = t74 * t85;
t73 = t74 ^ 2;
t92 = -0.2e1 * t84 * t73 * pkin(2) + (-t51 * t75 - t145) * pkin(1);
t5 = (-(t31 * t149 + t92) * t135 + (t124 * t182 - t156 * t168) * t29) * qJD(1);
t59 = t57 * t58;
t179 = qJD(1) * t59 * t115 + t5 * t174;
t177 = -t180 / 0.2e1;
t52 = t64 - t143;
t176 = (t52 * t75 + t145) * pkin(2) + 0.2e1 * t73 * t169;
t78 = 0.1e1 / pkin(4);
t175 = t78 ^ 2;
t171 = pkin(1) * t57;
t55 = t184 * pkin(2);
t47 = -rSges(3,2) * pkin(1) + t55;
t161 = rSges(3,2) * t74;
t48 = rSges(3,1) * t181 + pkin(2) * t161;
t25 = t47 * t85 + t48 * t51;
t119 = t25 * t129;
t140 = pkin(2) * qJD(1);
t131 = t75 * t140;
t28 = t52 * t164 + t181 * t85;
t20 = -pkin(1) * t28 * t151 + 0.1e1;
t19 = t20 * qJD(1);
t6 = t19 * t119 + t131;
t166 = t6 * t25;
t26 = -t47 * t51 + t48 * t85;
t130 = -t146 / 0.2e1;
t118 = t26 * t130;
t132 = t74 * t140;
t7 = t19 * t118 - t132;
t165 = t7 * t26;
t159 = rSges(4,2) * t74;
t30 = qJD(1) * t31;
t157 = t30 * t34;
t45 = rSges(4,1) * t181 + pkin(2) * t159;
t155 = t33 * t45;
t56 = t185 * pkin(2);
t46 = -rSges(4,2) * pkin(1) + t56;
t154 = t33 * t46;
t153 = t33 * t47;
t152 = t33 * t48;
t150 = t33 * t181;
t141 = pkin(1) * qJD(1);
t138 = qJD(1) * t83;
t137 = t33 * qJDD(1);
t134 = t58 * t164;
t128 = -0.2e1 * t136;
t127 = 0.2e1 * t136;
t125 = pkin(1) * t134;
t122 = pkin(1) * t132;
t121 = t58 * t81 * t141;
t120 = t138 * t148;
t117 = -0.2e1 * t122;
t116 = 0.2e1 * t122;
t112 = t19 * t121;
t63 = rSges(2,1) * t75 - rSges(2,2) * t74;
t62 = rSges(2,1) * t74 + rSges(2,2) * t75;
t76 = qJD(1) ^ 2;
t109 = t76 * t115;
t27 = t29 ^ 2;
t108 = t27 * t30 / t39 ^ 2 / t40 ^ 2 * t138;
t107 = -t28 * t157 - (t176 * qJD(1) + t30 * t150) * t33;
t105 = t111 / 0.2e1;
t104 = 0.2e1 * t28 * t33 * t84 * t134;
t54 = (rSges(3,1) * t75 - t161) * pkin(2);
t53 = (rSges(4,1) * t75 - t159) * pkin(2);
t103 = t20 * t121 * t164;
t101 = 0.2e1 * t27 * t76 * pkin(2) * t83 * t148 * t170;
t100 = t120 * t177;
t98 = -t25 * t125 + (t48 * t127 + t31 * t153 + t51 * t55 + t54 * t85) * t174;
t97 = -t26 * t125 + (t47 * t128 + t31 * t152 - t51 * t54 + t55 * t85) * t174;
t17 = t92 * qJD(1) + t30 * t149;
t2 = t109 * t182 + (-t29 * t137 + (-t29 * t157 - t17 * t33) * qJD(1)) * t168;
t93 = t59 * t109 + t2 * t174;
t44 = t185 * t140;
t43 = t184 * t140;
t42 = qJD(1) * t54;
t41 = qJD(1) * t53;
t24 = -t45 * t52 + t46 * t85;
t23 = t45 * t85 + t46 * t52;
t18 = qJD(1) * t104;
t15 = t46 * t127 + t31 * t155 + t52 * t53 + t56 * t85;
t12 = t45 * t128 + t31 * t154 - t52 * t56 + t53 * t85;
t11 = t46 * t116 + t30 * t155 + t52 * t41 + t44 * t85;
t10 = t47 * t117 + t30 * t152 - t51 * t42 + t43 * t85;
t9 = t48 * t116 + t30 * t153 + t42 * t85 + t51 * t43;
t8 = t45 * t117 + t30 * t154 + t41 * t85 - t52 * t44;
t4 = t18 + (-t28 * t156 + (-t31 * t150 - t176) * t33) * t57 * t141;
t3 = t107 * t171 + t18;
t1 = t76 * t104 + qJDD(1) + (t107 * qJD(1) - t28 * t137) * t171;
t13 = [-m(4) * (g(1) * (t12 * t174 - t125 * t24) + g(2) * (-t125 * t23 + t15 * t174)) * t78 + Icges(2,3) * qJDD(1) + m(4) * (-t24 * (t12 * t106 + t179 * t24) + (t15 * t105 - t179 * t23) * t23) * t175 * t106 + m(4) * ((t23 * t11 + t24 * t8) * t27 * t100 + ((t105 * t11 - t23 * t93) * t23 - (t106 * t8 + t24 * t93) * t24) * t175 * t123 + (t23 ^ 2 + t24 ^ 2) * (t29 * t17 * t100 + t57 * t101 * t180 + t108 * t177)) / 0.2e1 + ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t74 - (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t75 + (t62 ^ 2 + t63 ^ 2) * qJDD(1)) * m(2) + (g(1) * t164 - g(2) * t72 - (-t6 * t74 - t7 * t75) * t140 + ((-t1 * t26 / 0.2e1 - t19 * t10 / 0.2e1) * t146 + (-t76 * t75 + (t112 * t26 - qJDD(1)) * t74) * pkin(2)) * (t118 * t20 - t164) + t7 * (t10 * t130 * t20 + t103 * t26 - t131) + ((t1 * t25 / 0.2e1 + t19 * t9 / 0.2e1) * t146 + (qJDD(1) * t75 + (-t112 * t25 - t76) * t74) * pkin(2)) * (t119 * t20 + t72) + t6 * (t129 * t20 * t9 - t103 * t25 - t132) + (-t165 + t166) * t3 * t129 + (-g(1) * t98 - g(2) * t97 - (-t165 / 0.2e1 + t166 / 0.2e1) * t57 * t4 - (t6 * t98 - t7 * t97) * t19) * t81) * m(3) + (t1 * t20 + (t3 - t4) * t19) * Icges(3,3) + (t59 * t101 - t58 * t108 + (-t58 * t17 * t120 + (t5 * qJD(1) - t2) * t135) * t29) * Icges(4,3);];
tau = t13;
