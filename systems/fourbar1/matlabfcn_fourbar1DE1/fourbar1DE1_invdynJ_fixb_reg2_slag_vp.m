% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% fourbar1DE1
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
% 
% Output:
% tau_reg [1x(1*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1DE1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:18
% EndTime: 2020-04-24 19:57:22
% DurationCPUTime: 0.78s
% Computational Cost: add. (1852->122), mult. (3015->225), div. (76->10), fcn. (586->4), ass. (0->108)
t60 = cos(qJ(1));
t136 = pkin(2) * t60;
t116 = pkin(1) * t136;
t75 = pkin(1) ^ 2;
t120 = -0.2e1 * t116 + t75;
t72 = pkin(2) ^ 2;
t41 = t72 + t120;
t38 = 0.1e1 / t41 ^ 2;
t128 = t38 * t72;
t140 = pkin(4) - pkin(3);
t61 = pkin(3) + pkin(4);
t124 = (t140 ^ 2 * t61 ^ 2);
t66 = pkin(4) ^ 2;
t68 = pkin(3) ^ 2;
t48 = t66 + t68;
t143 = 0.2e1 * t41 * t48 - (2 * t124);
t59 = sin(qJ(1));
t139 = pkin(1) * t59;
t65 = qJD(1) ^ 2;
t102 = t65 * pkin(2) * t139 * t128;
t37 = 0.1e1 / t41;
t127 = (pkin(2) + t61) * (pkin(2) - t61);
t29 = t120 + t127;
t126 = (pkin(2) - t140) * (pkin(2) + t140);
t30 = t120 + t126;
t129 = t29 * t30;
t76 = sqrt(-t129);
t20 = 0.1e1 / t76;
t142 = -2 * qJD(1);
t141 = t37 / 0.2e1;
t138 = pkin(1) * t60;
t137 = pkin(2) * t59;
t64 = 0.2e1 * t72;
t70 = t72 ^ 2;
t73 = t75 ^ 2;
t23 = t73 + (t64 - 0.2e1 / 0.3e1 * t68 - 0.2e1 / 0.3e1 * t66) * t75 + t70 + (-0.4e1 / 0.3e1 * t68 - 0.4e1 / 0.3e1 * t66) * t72 + t124 / 0.3e1;
t135 = t23 * pkin(2);
t117 = pkin(1) * t137;
t17 = (-t29 - t30) * t117;
t131 = t17 * t20;
t35 = -g(1) * t137 + g(2) * t136;
t31 = g(2) * pkin(1) - t35;
t99 = g(1) * t60 + g(2) * t59;
t36 = t99 * pkin(2);
t134 = -t31 * t131 - t36 * t76;
t32 = -g(1) * pkin(1) + t36;
t133 = -t32 * t131 - t35 * t76;
t13 = qJD(1) * t17;
t132 = t13 * t20;
t45 = -pkin(2) + t138;
t130 = t20 * t45;
t125 = t59 ^ 2 * t75;
t123 = t59 * t76;
t122 = t60 * t75;
t121 = t61 * t140;
t119 = -t66 / 0.2e1 + t72;
t118 = t68 - t66;
t115 = pkin(2) * t125;
t69 = 0.1e1 / pkin(3);
t114 = t69 * t128;
t54 = t60 ^ 2;
t77 = pkin(1) * t75;
t113 = (-t68 / 0.2e1 + t119) * t54 * t77;
t112 = t45 * t121;
t43 = t68 / 0.2e1 + t119;
t111 = t72 + t118;
t110 = t69 * t141;
t108 = 0.12e2 * t113;
t24 = (t72 - t68 / 0.6e1 - t66 / 0.6e1) * t75 + t70 + (-0.5e1 / 0.6e1 * t68 - 0.5e1 / 0.6e1 * t66) * t72 + t124 / 0.6e1;
t107 = -0.12e2 * t24 * t122;
t106 = t65 * t114;
t105 = t59 * t112;
t104 = qJD(1) * t114;
t103 = t38 * t69 * t117;
t101 = t37 * t102;
t100 = t20 * t106;
t98 = t20 * t104;
t97 = qJDD(1) * t20 * t114;
t96 = 0.2e1 * t101;
t95 = 0.1e1 / t30 ^ 2 * t102;
t21 = t20 / t129;
t94 = t17 * t21 * t106;
t93 = t45 * t138 - t125;
t92 = t13 * t21 * t104;
t27 = 0.1e1 / t30;
t91 = 0.1e1 / t29 ^ 2 * t27 * t102;
t90 = t20 * t69 * t101;
t89 = t93 * t76 * t121;
t33 = t111 + t120;
t88 = -0.2e1 * t115 + (-t33 * t60 - t123) * pkin(1);
t34 = t41 - t118;
t87 = 0.2e1 * t115 + (t34 * t60 - t123) * pkin(1);
t42 = t75 + t43;
t86 = 0.4e1 * t48 * t45 * t115 + (0.4e1 * pkin(1) * pkin(2) * t42 - 0.8e1 * t43 * t122) * t123 + t93 * t143;
t67 = 0.1e1 / pkin(4);
t40 = -g(1) * t59 + g(2) * t60;
t25 = 0.1e1 / t29;
t18 = t45 * t76;
t16 = -0.4e1 * t42 * t116 + t73 + t72 * t111 + (0.4e1 * t43 * t54 - t118 + t64) * t75;
t12 = t32 * t76;
t11 = t31 * t76;
t10 = t34 * t139 + t18;
t9 = -t33 * t139 + t18;
t8 = t10 ^ 2;
t7 = t9 ^ 2;
t6 = t17 * t130;
t5 = t13 * t130;
t1 = [0, 0, 0, 0, 0, qJDD(1), -t40, t99, 0, 0, 0, 0, 0, 0, 0, t8 * t91 + (t8 * t95 + (t8 * t96 + (-qJDD(1) * t8 + ((t87 * qJD(1) + t5) * t142 + t65 * (t6 + t87)) * t10) * t128) * t27) * t25, -(t86 * qJD(1) + t16 * t132) * t98 + (t16 * t131 + t86) * t100 / 0.2e1 + (0.2e1 * t32 * t117 + t33 * t35 + t134) * t110 - (t32 * t33 - t11) * t103 + (-t97 - t92 + 0.2e1 * t90 + t94 / 0.2e1) * (t139 * t143 * t45 + t16 * t76), -0.2e1 * (-pkin(1) * t105 * t132 + (-t89 + (t107 + (0.3e1 * pkin(1) * t23 + t108) * pkin(2)) * t59) * qJD(1)) * t98 + (-t89 + (pkin(2) * t108 + t107 + (-t112 * t131 + 0.3e1 * t135) * pkin(1)) * t59) * t100 + (-0.2e1 * t31 * t117 - t33 * t36 + t133) * t110 - (-t31 * t33 - t12) * t103 + (-0.2e1 * t97 - 0.2e1 * t92 + 0.4e1 * t90 + t94) * (-0.4e1 * t113 * t136 + t77 ^ 2 / 0.2e1 + (0.6e1 * t24 * t54 + 0.3e1 / 0.2e1 * t70 - t124 / 0.2e1) * t75 + (0.3e1 / 0.2e1 * t73 - t48 * t75 + t126 * t127 / 0.2e1) * t72 + (-t76 * t105 - 0.3e1 * t60 * t135) * pkin(1)), 0, -pkin(2) * t40 + qJDD(1) * t72, 0, 0, 0, 0, 0, t7 * t91 + (t7 * t95 + (t7 * t96 + (-qJDD(1) * t7 + ((t88 * qJD(1) + t5) * t142 + t65 * (t6 + t88)) * t9) * t128) * t27) * t25, ((-t34 * t35 + t134) * t141 + (-t32 * t37 - (-t32 * t34 - t11) * t38) * t117) * t67, ((t34 * t36 + t133) * t141 + (t31 * t37 - (t31 * t34 - t12) * t38) * t117) * t67, 0, 0;];
tau_reg = t1;
