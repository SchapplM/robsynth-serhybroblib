% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x38]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m2OL_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:38
% EndTime: 2020-06-30 18:06:50
% DurationCPUTime: 1.72s
% Computational Cost: add. (3634->178), mult. (8320->315), div. (0->0), fcn. (8939->10), ass. (0->122)
t130 = cos(qJ(4));
t56 = sin(qJ(4));
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t74 = t128 * t57 - t131 * t59;
t75 = t128 * t59 + t131 * t57;
t31 = t130 * t75 - t56 * t74;
t62 = t130 * t74 + t56 * t75;
t58 = cos(qJ(6));
t53 = t58 ^ 2;
t54 = sin(qJ(6));
t123 = t54 ^ 2 - t53;
t94 = qJD(6) * t123;
t136 = pkin(4) * t56;
t129 = cos(qJ(5));
t55 = sin(qJ(5));
t17 = t129 * t31 - t55 * t62;
t119 = qJD(5) * t55;
t32 = (qJD(2) + qJD(3)) * t75;
t70 = t74 * qJD(3);
t63 = -qJD(2) * t74 - t70;
t14 = -qJD(4) * t62 + t130 * t63 - t56 * t32;
t60 = qJD(4) * t31 + t130 * t32 + t56 * t63;
t61 = t129 * t62;
t6 = -qJD(5) * t61 - t119 * t31 + t129 * t14 - t55 * t60;
t135 = t17 * t6;
t120 = qJD(5) * t17;
t7 = t129 * t60 + t55 * t14 + t120;
t134 = t54 * t7;
t133 = t58 * t6;
t132 = t58 * t7;
t112 = t131 * pkin(4);
t48 = t112 + pkin(2);
t105 = t130 * t48;
t89 = qJD(3) * t112;
t98 = t128 * qJD(3);
t28 = -qJD(4) * t105 - t130 * t89 + (qJD(4) * t128 + t98) * t136;
t121 = qJD(4) * t56;
t84 = t130 * t128;
t64 = (-qJD(4) * t84 + (-t131 * t56 - t84) * qJD(3)) * pkin(4);
t29 = -t121 * t48 + t64;
t101 = t129 * t29 + t55 * t28;
t38 = pkin(4) * t84 + t56 * t48;
t104 = t129 * t38;
t36 = -t128 * t136 + pkin(5) + t105;
t20 = t55 * t36 + t104;
t12 = -qJD(5) * t20 + t101;
t127 = t12 * t54;
t99 = t130 * pkin(2);
t47 = t99 + pkin(5);
t100 = qJD(5) * t129;
t103 = t129 * t56;
t71 = (-t130 * t55 - t103) * qJD(4);
t65 = (-t100 * t56 + t71) * pkin(2);
t27 = -t119 * t47 + t65;
t126 = t27 * t54;
t125 = t55 * t38;
t83 = qJD(4) * t99;
t124 = -t47 * t100 - t129 * t83;
t122 = pkin(2) * t56;
t51 = qJD(6) * t54;
t118 = qJD(6) * t58;
t117 = t57 * qJD(2);
t116 = t59 * qJD(2);
t115 = t54 * t133;
t114 = -t36 * t100 + t129 * t28 - t55 * t29;
t50 = pkin(4) * t117;
t113 = pkin(5) * t119;
t111 = t129 * pkin(5);
t16 = t55 * t31 + t61;
t110 = t16 * t51;
t109 = t16 * t118;
t108 = t54 * t118;
t107 = -0.2e1 * pkin(1) * qJD(2);
t106 = -0.2e1 * pkin(3) * qJD(6);
t102 = pkin(2) * t121;
t49 = -t59 * pkin(4) - pkin(1);
t19 = t129 * t36 + pkin(3) - t125;
t35 = -t122 * t55 + t129 * t47 + pkin(3);
t97 = qJD(6) * (-t19 - t35);
t46 = t111 + pkin(3);
t96 = qJD(6) * (-t19 - t46);
t95 = qJD(6) * (-t35 - t46);
t25 = t32 * pkin(2) + t50;
t93 = t58 * t113;
t92 = qJD(6) * (-pkin(3) - t19);
t91 = qJD(6) * (-pkin(3) - t35);
t90 = qJD(6) * (-pkin(3) - t46);
t88 = pkin(4) * t98;
t87 = pkin(5) * t100;
t86 = t58 * t100;
t85 = pkin(2) * t103;
t82 = t16 * t20 + t17 * t19;
t37 = t55 * t47 + t85;
t81 = t16 * t37 + t17 * t35;
t80 = (qJD(4) + qJD(5)) * t122;
t8 = pkin(5) * t60 + t25;
t3 = t7 * pkin(3) + t8;
t34 = pkin(2) * t74 + t49;
t18 = pkin(5) * t62 + t34;
t9 = t16 * pkin(3) + t18;
t79 = t118 * t9 + t3 * t54;
t78 = -t3 * t58 + t51 * t9;
t77 = t118 * t17 + t54 * t6;
t76 = t17 * t51 - t133;
t26 = t55 * t80 + t124;
t11 = t119 * t38 + t114;
t69 = -t11 * t16 + t12 * t17 + t19 * t6 + t20 * t7;
t68 = -t16 * t26 + t17 * t27 + t35 * t6 + t37 * t7;
t45 = 0.2e1 * t108;
t44 = t54 * t113;
t39 = -0.2e1 * t94;
t21 = t27 * t58;
t15 = t17 ^ 2;
t10 = t12 * t58;
t5 = t109 + t134;
t4 = -t110 + t132;
t2 = t17 * t94 - t115;
t1 = 0.4e1 * t108 * t17 + t123 * t6;
t13 = [0, 0, 0, 0.2e1 * t57 * t116, 0.2e1 * (-t57 ^ 2 + t59 ^ 2) * qJD(2), 0, 0, 0, t57 * t107, t59 * t107, 0.2e1 * t75 * t63, -0.2e1 * t32 * t75 - 0.2e1 * t63 * t74, 0, 0, 0, 0.2e1 * t49 * t32 + 0.2e1 * t50 * t74, -0.2e1 * t49 * t70 + 0.2e1 * (pkin(4) * t57 * t75 - t49 * t74) * qJD(2), 0.2e1 * t31 * t14, -0.2e1 * t14 * t62 - 0.2e1 * t31 * t60, 0, 0, 0, 0.2e1 * t25 * t62 + 0.2e1 * t34 * t60, 0.2e1 * t14 * t34 + 0.2e1 * t25 * t31, 0.2e1 * t135, -0.2e1 * t16 * t6 - 0.2e1 * t17 * t7, 0, 0, 0, 0.2e1 * t16 * t8 + 0.2e1 * t18 * t7, 0.2e1 * t17 * t8 + 0.2e1 * t18 * t6, -0.2e1 * t108 * t15 + 0.2e1 * t135 * t53, -0.4e1 * t115 * t17 + 0.2e1 * t15 * t94, -0.2e1 * t132 * t17 + 0.2e1 * t16 * t76, 0.2e1 * t134 * t17 + 0.2e1 * t16 * t77, 0.2e1 * t16 * t7, 0.2e1 * t132 * t9 - 0.2e1 * t16 * t78, -0.2e1 * t134 * t9 - 0.2e1 * t16 * t79; 0, 0, 0, 0, 0, t116, -t117, 0, 0, 0, 0, 0, t63, -t32, 0, 0, 0, 0, 0, t14, -t60, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t2, t1, t5, t4, 0, t118 * t82 + t54 * t69, -t51 * t82 + t58 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, -0.2e1 * t89, 0, 0, 0, 0, 0, 0.2e1 * t29, 0.2e1 * t28, 0, 0, 0, 0, 0, 0.2e1 * t12, 0.2e1 * t11, t45, t39, 0, 0, 0, -0.2e1 * t19 * t51 + 0.2e1 * t10, -0.2e1 * t118 * t19 - 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t32, 0, 0, 0, 0, 0, t14, -t60, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t2, t1, t5, t4, 0, t118 * t81 + t54 * t68, -t51 * t81 + t58 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t89, 0, 0, 0, 0, 0, (-pkin(2) - t48) * t121 + t64, -t83 + t28, 0, 0, 0, 0, 0, pkin(2) * t71 + (-t85 - t104 + (-t36 - t47) * t55) * qJD(5) + t101, (qJD(5) * t38 + t80) * t55 + t114 + t124, t45, t39, 0, 0, 0, t54 * t97 + t10 + t21, (-t12 - t27) * t54 + t58 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t102, -0.2e1 * t83, 0, 0, 0, 0, 0, 0.2e1 * t27, 0.2e1 * t26, t45, t39, 0, 0, 0, -0.2e1 * t35 * t51 + 0.2e1 * t21, -0.2e1 * t118 * t35 - 0.2e1 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t60, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t2, t1, t5, t4, 0, t77 * t46 + (t55 * t109 + (t55 * t7 + (t129 * t16 - t17 * t55) * qJD(5)) * t54) * pkin(5), -t76 * t46 + (t16 * t86 + (-t110 + (t7 - t120) * t58) * t55) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28, 0, 0, 0, 0, 0, (-t104 + (-pkin(5) - t36) * t55) * qJD(5) + t101, (-t111 + t125) * qJD(5) + t114, t45, t39, 0, 0, 0, t54 * t96 + t10 - t93, t58 * t96 - t127 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t83, 0, 0, 0, 0, 0, (-pkin(5) - t47) * t119 + t65, -t87 + t26, t45, t39, 0, 0, 0, t54 * t95 + t21 - t93, t58 * t95 - t126 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113, -0.2e1 * t87, t45, t39, 0, 0, 0, -0.2e1 * t46 * t51 - 0.2e1 * t93, -0.2e1 * t118 * t46 + 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t2, t1, t5, t4, 0, t77 * pkin(3), -t76 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, t45, t39, 0, 0, 0, t54 * t92 + t10, t58 * t92 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, t45, t39, 0, 0, 0, t54 * t91 + t21, t58 * t91 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t87, t45, t39, 0, 0, 0, t54 * t90 - t93, t58 * t90 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t39, 0, 0, 0, t54 * t106, t58 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, -t7, t78, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t51, 0, t11 * t54 - t118 * t20, t11 * t58 + t20 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t51, 0, -t118 * t37 + t26 * t54, t26 * t58 + t37 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t51, 0, (-t100 * t54 - t118 * t55) * pkin(5), (t51 * t55 - t86) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t51, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
