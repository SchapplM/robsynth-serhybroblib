% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m2DE_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:45
% EndTime: 2020-05-03 01:06:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (408->110), mult. (867->185), div. (0->0), fcn. (520->8), ass. (0->92)
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t52 = qJD(3) ^ 2;
t67 = qJDD(3) * t45 + t49 * t52;
t47 = sin(qJ(1));
t51 = cos(qJ(1));
t22 = g(1) * t51 + g(2) * t47;
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t53 = qJD(2) ^ 2;
t68 = qJDD(2) * t46 + t50 * t53;
t8 = -pkin(4) * t68 - t22;
t127 = -pkin(5) * t67 + t8;
t125 = g(1) * t47 - g(2) * t51;
t91 = qJD(1) * qJD(2);
t82 = t46 * t91;
t56 = -0.2e1 * pkin(4) * t82 + t125;
t37 = qJD(1) + qJD(4);
t124 = t37 ^ 2;
t118 = pkin(5) * t45;
t119 = pkin(4) * t46;
t20 = qJD(2) * t119 + qJD(3) * t118;
t102 = t20 * qJD(1);
t60 = -0.2e1 * t102 + t125;
t28 = pkin(4) * t50 + pkin(1);
t121 = qJD(2) - qJD(3);
t120 = pkin(4) * pkin(5);
t114 = t49 * g(3);
t113 = t50 * g(3);
t48 = cos(qJ(4));
t112 = t127 * t48;
t36 = qJDD(1) + qJDD(4);
t111 = t36 * t48;
t54 = qJD(1) ^ 2;
t110 = t45 * t54;
t109 = t46 * t54;
t108 = t48 * t20;
t107 = t49 * t54;
t38 = t45 ^ 2;
t40 = t49 ^ 2;
t106 = t38 - t40;
t39 = t46 ^ 2;
t41 = t50 ^ 2;
t105 = t39 - t41;
t23 = pkin(2) + t28;
t19 = pkin(5) * t49 + t23;
t14 = pkin(3) + t19;
t104 = qJD(1) * t14;
t71 = t45 * t50 - t46 * t49;
t103 = qJD(3) * t71;
t101 = -qJD(4) + t37;
t100 = qJDD(1) * t23;
t99 = t14 * qJDD(1);
t98 = t19 * qJDD(1);
t97 = t28 * qJDD(1);
t96 = t45 * qJDD(1);
t95 = t46 * qJDD(1);
t94 = t49 * qJDD(1);
t93 = t50 * qJDD(1);
t4 = t121 * t71;
t92 = t4 * t120;
t90 = qJD(1) * qJD(3);
t70 = t45 * t46 + t49 * t50;
t89 = t70 * t120;
t88 = qJDD(1) * pkin(1);
t87 = t45 * t107;
t86 = t50 * t109;
t85 = pkin(4) ^ 2 * qJDD(2) + t119 * t22;
t84 = pkin(4) * t95;
t44 = sin(qJ(4));
t83 = t44 * t104;
t81 = qJD(2) * t103;
t80 = qJD(4) * (-qJD(1) - t37);
t79 = -0.2e1 * t23 * t90;
t78 = t45 * t49 * t90;
t77 = t50 * t82;
t76 = -0.2e1 * pkin(1) * t91;
t75 = qJD(4) * t108 + (-t102 + t99) * t48 - t127 * t44;
t74 = qJD(3) * t92 + qJDD(3) * t89 + t85;
t73 = pkin(5) ^ 2 * qJDD(3) + qJD(2) * t92 + qJDD(2) * t89 + t118 * t22;
t66 = -t125 + (qJD(1) + t101) * t20;
t64 = pkin(4) * t81 - t114;
t63 = -pkin(5) * t81 - t113;
t62 = t23 * t54 + t22;
t61 = pkin(1) * t54 + t22;
t59 = 0.2e1 * t88 + t125;
t58 = t124 * t48 + t36 * t44;
t57 = -t124 * t44 + t111;
t55 = -t56 - 0.2e1 * t100;
t32 = pkin(4) * t109;
t3 = t121 * t70;
t1 = [0, 0, 0, 0, 0, qJDD(1), t125, t22, 0, 0, qJDD(1) * t39 + 0.2e1 * t77, -0.2e1 * t105 * t91 + 0.2e1 * t46 * t93, t68, qJDD(1) * t41 - 0.2e1 * t77, qJDD(2) * t50 - t46 * t53, 0, t46 * t76 + t50 * t59, -t46 * t59 + t50 * t76, -t22, (t125 + t88) * pkin(1), 0, 0, 0, qJDD(1), 0, 0, 0.2e1 * t97 + t56, 0, t8, (t97 + t56) * t28, qJDD(1) * t38 + 0.2e1 * t78, -0.2e1 * t106 * t90 + 0.2e1 * t45 * t94, t67, qJDD(1) * t40 - 0.2e1 * t78, qJDD(3) * t49 - t45 * t52, 0, t45 * t79 - t49 * t55, t45 * t55 + t49 * t79, t8, (t56 + t100) * t23, 0, 0, 0, qJDD(1), 0, 0, 0.2e1 * t98 + t60, 0, t127, (t60 + t98) * t19, 0, 0, 0, 0, 0, t36, (-t37 * t20 + t125) * t48 + (t44 * t80 + t111) * t14 + t75, t48 * t14 * t80 + ((-qJDD(1) - t36) * t14 + t66) * t44 - t112, 0, (t60 + t99) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t105 * t54, t95, t86, t93, qJDD(2), t46 * t61 - t113, g(3) * t46 + t50 * t61, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t84, (t109 * t28 - t113) * pkin(4) + t85, 0, 0, 0, 0, 0, 0, (t46 * t107 + qJDD(3) * t70 + (-qJD(2) * t71 + t4) * qJD(3)) * pkin(4), (-t45 * t109 - qJDD(3) * t71 + (-qJD(2) * t70 + t3) * qJD(3)) * pkin(4), -t84, (t109 * t23 - t113) * pkin(4) + t85, 0, 0, 0, 0, 0, 0, t32, 0, -t84, (t109 * t19 + t63) * pkin(4) + t74, 0, 0, 0, 0, 0, 0, t58 * t119, t57 * t119, 0, (t109 * t14 + t63) * pkin(4) + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t106 * t54, t96, t87, t94, qJDD(3), -t114 + t62 * t45 + (qJDD(2) * t70 + (t4 + t103) * qJD(2)) * pkin(4), g(3) * t45 + t62 * t49 + (-qJDD(2) * t71 + (qJD(3) * t70 + t3) * qJD(2)) * pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, pkin(5) * t110, 0, -pkin(5) * t96, (t110 * t19 + t64) * pkin(5) + t73, 0, 0, 0, 0, 0, 0, t58 * t118, t57 * t118, 0, (t110 * t14 + t64) * pkin(5) + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -qJD(4) * t83 - t37 * (-t83 + t108) + t48 * t125 + t75, t101 * t48 * t104 + (t66 - t99) * t44 - t112, 0, 0;];
tau_reg = t1;
