% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% palh2m1DE
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m1DE_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:35
% EndTime: 2020-06-30 17:39:37
% DurationCPUTime: 0.57s
% Computational Cost: add. (487->120), mult. (1043->199), div. (0->0), fcn. (803->8), ass. (0->90)
t58 = cos(qJ(2));
t39 = t58 * pkin(2) + pkin(1);
t82 = 2 * qJD(1);
t120 = t39 * t82;
t113 = 2 * qJDD(1);
t55 = sin(qJ(1));
t59 = cos(qJ(1));
t32 = g(1) * t55 - g(2) * t59;
t54 = sin(qJ(2));
t93 = qJD(2) * t54;
t119 = pkin(2) * t82 * t93 - t39 * t113 - t32;
t57 = cos(qJ(3));
t100 = t58 * t57;
t53 = sin(qJ(3));
t105 = t54 * t53;
t27 = t100 - t105;
t65 = t27 * qJD(3);
t66 = t27 * qJD(2);
t13 = t65 + t66;
t101 = t58 * t53;
t38 = t57 * pkin(3) + pkin(2);
t22 = pkin(3) * t101 + t54 * t38;
t33 = g(1) * t59 + g(2) * t55;
t68 = t57 * t54 + t101;
t92 = qJD(2) * t58;
t9 = t38 * t92 + (-t53 * t93 + t65) * pkin(3);
t118 = t9 * qJD(2) + t22 * qJDD(2) + (qJD(3) * t13 + qJDD(3) * t68) * pkin(3) + t33;
t46 = qJD(2) + qJD(3);
t14 = t46 * t68;
t61 = qJD(1) ^ 2;
t109 = t39 * t61;
t69 = t53 * g(3) - t33 * t57;
t77 = t57 * g(3) + t53 * t33;
t117 = t57 * qJDD(2) * pkin(2) + t109 * t68 - t69 * t54 + t77 * t58;
t31 = -0.2e1 * t100 * t105;
t50 = t57 ^ 2;
t51 = t58 ^ 2;
t88 = t51 - 0.1e1 / 0.2e1;
t116 = t31 + (-t53 ^ 2 + t50) * t88;
t41 = t50 - 0.1e1 / 0.2e1;
t49 = t54 ^ 2;
t98 = t49 - t51;
t115 = -t98 * t41 + t31;
t114 = -2 * t61;
t97 = pkin(3) * qJD(3);
t23 = t27 * t97;
t71 = -pkin(3) * t105 + t38 * t58;
t111 = -t71 * qJD(2) - t23 + t9;
t56 = cos(qJ(4));
t110 = t118 * t56;
t47 = qJD(1) + qJD(4);
t108 = t47 * t56;
t45 = qJDD(1) + qJDD(4);
t52 = sin(qJ(4));
t107 = t52 * t45;
t106 = t52 * t47;
t104 = t54 * t58;
t103 = t54 * t61;
t102 = t56 * t45;
t8 = -t38 * t93 + (-t68 * qJD(3) - t53 * t92) * pkin(3);
t96 = t8 * qJD(1);
t95 = qJD(1) * t52;
t94 = qJD(1) * t56;
t91 = qJD(4) * t56;
t90 = -qJD(1) - t47;
t89 = -qJD(4) + t47;
t87 = qJDD(1) * t68;
t67 = pkin(1) + t71;
t19 = pkin(4) + t67;
t86 = t19 * qJDD(1);
t85 = t58 * qJDD(1);
t84 = qJD(1) * qJD(2);
t83 = t61 * t68 * t27;
t79 = t19 * t95;
t76 = -0.2e1 * pkin(1) * t84;
t75 = qJD(2) * (-qJD(3) + t46);
t74 = qJD(3) * (-qJD(2) - t46);
t73 = qJD(4) * t90;
t11 = t22 * qJD(2) + t68 * t97;
t72 = t11 * t91 + t118 * t52 + t56 * t86 + t8 * t94;
t70 = t27 * t109 - t54 * t77 - t58 * t69;
t64 = pkin(1) * t61 + t33;
t63 = pkin(1) * t113 + t32;
t60 = qJD(2) ^ 2;
t44 = qJDD(2) + qJDD(3);
t24 = t68 * pkin(3);
t12 = pkin(3) * t66 + t23;
t5 = qJDD(1) * t27;
t4 = -t87 + (t27 * t46 - t13) * qJD(1);
t1 = [qJDD(1), t32, t33, t49 * qJDD(1) + 0.2e1 * t84 * t104, 0.2e1 * t54 * t85 - 0.2e1 * t98 * t84, -qJDD(2) * t54 - t60 * t58, -qJDD(2) * t58 + t60 * t54, 0, t54 * t76 + t63 * t58, -t63 * t54 + t58 * t76, -(-0.2e1 * qJD(1) * t13 - t87) * t68, 0.4e1 * (t88 * t57 * t53 + t41 * t104) * qJDD(1) + 0.4e1 * t115 * t84 + 0.4e1 * t116 * qJD(1) * qJD(3), -t46 * t13 - t44 * t68, t46 * t14 - t44 * t27, 0, -t119 * t27 - t14 * t120, t119 * t68 - t13 * t120, t67 * t113 + t32 + 0.2e1 * t96, t45, (t47 * t8 + t32) * t56 + (t52 * t73 + t102) * t19 + t72, t56 * t19 * t73 + (-qJD(4) * t11 - t32 + t90 * t8 + (-qJDD(1) - t45) * t19) * t52 + t110; 0, 0, 0, -t58 * t103, t98 * t61, -t54 * qJDD(1), -t85, qJDD(2), g(3) * t58 + t64 * t54, -g(3) * t54 + t64 * t58, -t83, t115 * t114, t4, -t5, t44, (t27 * t103 + t57 * t44 + t53 * t74) * pkin(2) + t117, (-t68 * t103 + (-qJDD(2) - t44) * t53 + t57 * t74) * pkin(2) + t70, t22 * t61, 0, t22 * t107 + (t22 * t108 + t111 * t52) * t47, t22 * t102 + (-t22 * t106 + t111 * t56) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t116 * t114, t4, -t5, t44, t53 * pkin(2) * t75 + t117, (-qJDD(2) * t53 + t57 * t75) * pkin(2) + t70, t24 * t61, 0, -(t52 * t12 - t24 * t94) * t47 + (t13 * t106 - (-t47 * t91 - t107) * t68) * pkin(3), (-t56 * t12 - t24 * t95) * t47 + (t13 * t108 - (qJD(4) * t106 - t102) * t68) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -qJD(4) * t79 - (t56 * t11 - t79) * t47 + t32 * t56 + t72, t89 * t19 * t94 + (t89 * t11 - t32 - t86 - t96) * t52 + t110;];
tau_reg = t1;
