% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m1OL_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:27:53
% EndTime: 2020-05-03 00:28:01
% DurationCPUTime: 1.56s
% Computational Cost: add. (1703->137), mult. (4199->272), div. (0->0), fcn. (4107->8), ass. (0->103)
t55 = sin(qJ(5));
t53 = t55 ^ 2;
t58 = cos(qJ(5));
t54 = t58 ^ 2;
t113 = t53 + t54;
t117 = sin(qJ(3));
t118 = cos(qJ(4));
t56 = sin(qJ(4));
t119 = cos(qJ(3));
t100 = t119 * pkin(2);
t84 = t100 + pkin(3);
t73 = t84 * t118;
t87 = qJD(3) * t100;
t90 = t117 * qJD(3);
t16 = -qJD(4) * t73 - t118 * t87 + (t117 * qJD(4) + t90) * t56 * pkin(2);
t130 = t113 * t16;
t83 = t117 * t118;
t129 = (qJD(3) + qJD(4)) * pkin(2) * (t119 * t56 + t83);
t114 = t53 - t54;
t127 = t114 * qJD(5);
t126 = qJD(2) + qJD(3);
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t94 = t117 * t59;
t35 = t119 * t57 + t94;
t95 = t117 * t57;
t71 = -t119 * t59 + t95;
t68 = t56 * t71;
t23 = -t118 * t35 + t68;
t109 = qJD(4) * t56;
t66 = t71 * qJD(3);
t60 = t71 * qJD(2) + t66;
t61 = t126 * t35;
t64 = t118 * t71;
t7 = -qJD(4) * t64 - t35 * t109 - t118 * t60 - t56 * t61;
t124 = t23 * t7;
t91 = qJD(4) * t118;
t8 = qJD(4) * t68 - t118 * t61 - t35 * t91 + t56 * t60;
t123 = t55 * t8;
t122 = t57 * pkin(2);
t121 = t58 * t7;
t120 = t58 * t8;
t101 = pkin(3) * t109;
t17 = t101 + t129;
t99 = t117 * pkin(2);
t32 = -t56 * t99 + t73;
t28 = -pkin(4) - t32;
t52 = qJD(5) * t58;
t116 = t17 * t55 + t28 * t52;
t31 = pkin(2) * t83 + t84 * t56;
t85 = pkin(3) * t91;
t30 = t113 * t85;
t49 = -t118 * pkin(3) - pkin(4);
t115 = t55 * t101 + t49 * t52;
t112 = pkin(3) * qJD(4);
t111 = qJD(2) * t57;
t110 = qJD(2) * t59;
t107 = qJD(5) * t55;
t22 = -t56 * t35 - t64;
t106 = 0.2e1 * t22 * t8;
t105 = t55 * t121;
t104 = -0.2e1 * t111;
t103 = pkin(4) * t107;
t102 = pkin(4) * t52;
t98 = t55 * t52;
t97 = t57 * t110;
t50 = t119 * pkin(3) + pkin(2);
t26 = -pkin(3) * t95 + t50 * t59 + pkin(1);
t9 = t22 * pkin(4) - t23 * pkin(6) + t26;
t96 = t113 * t9;
t24 = t28 * t107;
t92 = -t17 * t58 + t24;
t89 = -0.2e1 * t101;
t21 = t23 ^ 2;
t88 = t21 * t98;
t86 = pkin(2) * t90;
t82 = t16 * t22 + t17 * t23;
t27 = pkin(6) + t31;
t81 = t22 * t27 - t23 * t28;
t48 = t56 * pkin(3) + pkin(6);
t80 = t22 * t48 - t23 * t49;
t36 = t49 * t107;
t79 = -t58 * t101 + t36;
t67 = t35 * qJD(3);
t15 = -t50 * t111 + (-qJD(2) * t94 - t67) * pkin(3);
t2 = t8 * pkin(4) + t7 * pkin(6) + t15;
t78 = -t55 * t2 - t9 * t52;
t77 = -t9 * t107 + t58 * t2;
t5 = t22 * t52 + t123;
t76 = t22 * t107 - t120;
t75 = -t23 * t52 + t55 * t7;
t74 = t23 * t107 + t121;
t72 = t113 * t118;
t70 = (-t118 * t22 + t23 * t56) * qJD(4);
t65 = -t27 * t8 - t28 * t7 + t82;
t63 = pkin(3) * t70 - t48 * t8 - t49 * t7;
t51 = t59 * pkin(2) + pkin(1);
t46 = -0.2e1 * t98;
t45 = 0.2e1 * t98;
t34 = -0.2e1 * t127;
t3 = t23 * t127 + t105;
t1 = t114 * t7 - 0.4e1 * t23 * t98;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97, 0.2e1 * (-t57 ^ 2 + t59 ^ 2) * qJD(2), 0, -0.2e1 * t97, 0, 0, pkin(1) * t104, -0.2e1 * pkin(1) * t110, 0, 0, -0.2e1 * t35 * t60, 0.2e1 * t126 * (-t35 ^ 2 + t71 ^ 2), 0, 0.2e1 * t71 * t61, 0, 0, -0.2e1 * t51 * t67 + 0.2e1 * (t71 * t122 - t51 * t35) * qJD(2), 0.2e1 * t51 * t66 + 0.2e1 * (t35 * t122 + t51 * t71) * qJD(2), 0, t51 * pkin(2) * t104, -0.2e1 * t124, 0.2e1 * t22 * t7 - 0.2e1 * t23 * t8, 0, t106, 0, 0, 0.2e1 * t15 * t22 + 0.2e1 * t26 * t8, 0.2e1 * t15 * t23 - 0.2e1 * t26 * t7, 0, 0.2e1 * t26 * t15, -0.2e1 * t54 * t124 - 0.2e1 * t88, 0.4e1 * t23 * t105 + 0.2e1 * t21 * t127, 0.2e1 * t23 * t120 - 0.2e1 * t74 * t22, -0.2e1 * t53 * t124 + 0.2e1 * t88, -0.2e1 * t23 * t123 + 0.2e1 * t75 * t22, t106, 0.2e1 * t9 * t120 + 0.2e1 * t77 * t22, -0.2e1 * t9 * t123 + 0.2e1 * t78 * t22, -0.2e1 * t113 * t23 * t2 + 0.2e1 * t7 * t96, 0.2e1 * t2 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, t111, 0, 0, 0, 0, 0, 0, 0, t60, 0, t61, 0, 0, 0, -t35 * t86 + t61 * t99 + (-t126 * t71 + t66) * t100, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, -t31 * t8 + t32 * t7 + t82, 0, -t3, t1, t5, t3, -t76, 0, -t81 * t52 + t65 * t55, t81 * t107 + t65 * t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t87, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, 0.2e1 * t16, 0, -0.2e1 * t16 * t31 - 0.2e1 * t17 * t32, t45, t34, 0, t46, 0, 0, 0.2e1 * t92, 0.2e1 * t116, -0.2e1 * t130, -0.2e1 * t130 * t27 + 0.2e1 * t28 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t61, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, (t118 * t7 - t56 * t8 + t70) * pkin(3), 0, -t3, t1, t5, t3, -t76, 0, -t80 * t52 + t63 * t55, t80 * t107 + t63 * t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t87, 0, 0, 0, 0, 0, 0, 0, 0, t89 - t129, -t85 + t16, 0, (-t118 * t17 - t16 * t56 + (t118 * t31 - t32 * t56) * qJD(4)) * pkin(3), t45, t34, 0, t46, 0, 0, t24 + t36 + (-t17 - t101) * t58, t115 + t116, t30 - t130, t17 * t49 - t48 * t130 + (t72 * t27 + t28 * t56) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -0.2e1 * t85, 0, 0, t45, t34, 0, t46, 0, 0, 0.2e1 * t79, 0.2e1 * t115, 0.2e1 * t30, 0.2e1 * (t72 * t48 + t49 * t56) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, 0, 0, -t3, t1, t5, t3, -t76, 0, t75 * pkin(4) - t5 * pkin(6), t74 * pkin(4) + t76 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t45, t34, 0, t46, 0, 0, t92 - t103, -t102 + t116, -t130, -t17 * pkin(4) - pkin(6) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t85, 0, 0, t45, t34, 0, t46, 0, 0, t79 - t103, -t102 + t115, t30, (-pkin(4) * t56 + t72 * pkin(6)) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t34, 0, t46, 0, 0, -0.2e1 * t103, -0.2e1 * t102, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, t75, t8, t77, t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t107, 0, t55 * t16 - t27 * t52, t27 * t107 + t58 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t107, 0, -t48 * t52 - t55 * t85, t48 * t107 - t58 * t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t107, 0, -pkin(6) * t52, pkin(6) * t107, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
