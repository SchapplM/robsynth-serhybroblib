% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% 
% Output:
% tau_reg [1x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1TE_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar1TE_invdynJ_fixb_regmin_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1TE_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:20:59
% EndTime: 2020-06-26 17:21:03
% DurationCPUTime: 0.79s
% Computational Cost: add. (1850->120), mult. (3010->223), div. (76->10), fcn. (584->4), ass. (0->107)
t59 = cos(qJ(1));
t135 = pkin(2) * t59;
t115 = pkin(1) * t135;
t74 = pkin(1) ^ 2;
t119 = -0.2e1 * t115 + t74;
t71 = pkin(2) ^ 2;
t40 = t71 + t119;
t38 = 0.1e1 / t40 ^ 2;
t127 = t38 * t71;
t139 = pkin(4) - pkin(3);
t60 = pkin(3) + pkin(4);
t123 = (t139 ^ 2 * t60 ^ 2);
t65 = pkin(4) ^ 2;
t67 = pkin(3) ^ 2;
t47 = t65 + t67;
t142 = 0.2e1 * t40 * t47 - (2 * t123);
t58 = sin(qJ(1));
t138 = pkin(1) * t58;
t64 = qJD(1) ^ 2;
t101 = t64 * pkin(2) * t138 * t127;
t37 = 0.1e1 / t40;
t126 = (pkin(2) + t60) * (pkin(2) - t60);
t29 = t119 + t126;
t125 = (pkin(2) - t139) * (pkin(2) + t139);
t30 = t119 + t125;
t128 = t29 * t30;
t75 = sqrt(-t128);
t20 = 0.1e1 / t75;
t141 = -2 * qJD(1);
t140 = t37 / 0.2e1;
t137 = pkin(1) * t59;
t136 = pkin(2) * t58;
t63 = 0.2e1 * t71;
t69 = t71 ^ 2;
t72 = t74 ^ 2;
t23 = t72 + (t63 - 0.2e1 / 0.3e1 * t67 - 0.2e1 / 0.3e1 * t65) * t74 + t69 + (-0.4e1 / 0.3e1 * t67 - 0.4e1 / 0.3e1 * t65) * t71 + t123 / 0.3e1;
t134 = t23 * pkin(2);
t116 = pkin(1) * t136;
t17 = (-t29 - t30) * t116;
t130 = t17 * t20;
t35 = -g(1) * t136 + g(2) * t135;
t31 = g(2) * pkin(1) - t35;
t98 = g(1) * t59 + g(2) * t58;
t36 = t98 * pkin(2);
t133 = -t31 * t130 - t36 * t75;
t32 = -g(1) * pkin(1) + t36;
t132 = -t32 * t130 - t35 * t75;
t13 = qJD(1) * t17;
t131 = t13 * t20;
t44 = -pkin(2) + t137;
t129 = t20 * t44;
t124 = t58 ^ 2 * t74;
t122 = t58 * t75;
t121 = t59 * t74;
t120 = t60 * t139;
t118 = -t65 / 0.2e1 + t71;
t117 = t67 - t65;
t114 = pkin(2) * t124;
t68 = 0.1e1 / pkin(3);
t113 = t68 * t127;
t53 = t59 ^ 2;
t76 = pkin(1) * t74;
t112 = (-t67 / 0.2e1 + t118) * t53 * t76;
t111 = t44 * t120;
t42 = t67 / 0.2e1 + t118;
t110 = t71 + t117;
t109 = t68 * t140;
t24 = (t71 - t67 / 0.6e1 - t65 / 0.6e1) * t74 + t69 + (-0.5e1 / 0.6e1 * t67 - 0.5e1 / 0.6e1 * t65) * t71 + t123 / 0.6e1;
t107 = -0.12e2 * t24 * t121;
t106 = 0.12e2 * t112;
t105 = t64 * t113;
t104 = t58 * t111;
t103 = qJD(1) * t113;
t102 = t38 * t68 * t116;
t100 = t37 * t101;
t99 = t20 * t105;
t97 = t20 * t103;
t96 = qJDD(1) * t20 * t113;
t95 = 0.2e1 * t100;
t94 = 0.1e1 / t30 ^ 2 * t101;
t21 = t20 / t128;
t93 = t17 * t21 * t105;
t92 = t44 * t137 - t124;
t91 = t13 * t21 * t103;
t27 = 0.1e1 / t30;
t90 = 0.1e1 / t29 ^ 2 * t27 * t101;
t89 = t20 * t68 * t100;
t88 = t92 * t75 * t120;
t33 = t110 + t119;
t87 = -0.2e1 * t114 + (-t33 * t59 - t122) * pkin(1);
t34 = t40 - t117;
t86 = 0.2e1 * t114 + (t34 * t59 - t122) * pkin(1);
t41 = t74 + t42;
t85 = 0.4e1 * t47 * t44 * t114 + (0.4e1 * pkin(1) * pkin(2) * t41 - 0.8e1 * t42 * t121) * t122 + t92 * t142;
t66 = 0.1e1 / pkin(4);
t25 = 0.1e1 / t29;
t18 = t44 * t75;
t16 = -0.4e1 * t41 * t115 + t72 + t71 * t110 + (0.4e1 * t42 * t53 - t117 + t63) * t74;
t12 = t32 * t75;
t11 = t31 * t75;
t10 = t34 * t138 + t18;
t9 = -t33 * t138 + t18;
t8 = t10 ^ 2;
t7 = t9 ^ 2;
t6 = t17 * t129;
t5 = t13 * t129;
t1 = [qJDD(1), g(1) * t58 - g(2) * t59, t98, t8 * t90 + (t8 * t94 + (t8 * t95 + (-qJDD(1) * t8 + ((t86 * qJD(1) + t5) * t141 + t64 * (t6 + t86)) * t10) * t127) * t27) * t25, -(t85 * qJD(1) + t16 * t131) * t97 + (t16 * t130 + t85) * t99 / 0.2e1 + (0.2e1 * t32 * t116 + t33 * t35 + t133) * t109 - (t33 * t32 - t11) * t102 + (-t96 - t91 + 0.2e1 * t89 + t93 / 0.2e1) * (t138 * t142 * t44 + t16 * t75), -0.2e1 * (-pkin(1) * t104 * t131 + (-t88 + (t107 + (0.3e1 * pkin(1) * t23 + t106) * pkin(2)) * t58) * qJD(1)) * t97 + (-t88 + (pkin(2) * t106 + t107 + (-t111 * t130 + 0.3e1 * t134) * pkin(1)) * t58) * t99 + (-0.2e1 * t31 * t116 - t33 * t36 + t132) * t109 - (-t31 * t33 - t12) * t102 + (-0.2e1 * t96 - 0.2e1 * t91 + 0.4e1 * t89 + t93) * (-0.4e1 * t112 * t135 + t76 ^ 2 / 0.2e1 + (0.6e1 * t24 * t53 + 0.3e1 / 0.2e1 * t69 - t123 / 0.2e1) * t74 + (0.3e1 / 0.2e1 * t72 - t47 * t74 + t125 * t126 / 0.2e1) * t71 + (-t75 * t104 - 0.3e1 * t59 * t134) * pkin(1)), t7 * t90 + (t7 * t94 + (t7 * t95 + (-t7 * qJDD(1) + ((t87 * qJD(1) + t5) * t141 + t64 * (t6 + t87)) * t9) * t127) * t27) * t25, ((-t34 * t35 + t133) * t140 + (-t32 * t37 - (-t32 * t34 - t11) * t38) * t116) * t66, ((t34 * t36 + t132) * t140 + (t31 * t37 - (t31 * t34 - t12) * t38) * t116) * t66;];
tau_reg = t1;
