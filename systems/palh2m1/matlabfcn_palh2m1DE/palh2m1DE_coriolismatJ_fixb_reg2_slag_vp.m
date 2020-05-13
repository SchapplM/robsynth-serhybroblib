% Calculate inertial parameters regressor of coriolis matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = palh2m1DE_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:37
% EndTime: 2020-05-02 23:52:38
% DurationCPUTime: 0.75s
% Computational Cost: add. (354->115), mult. (759->184), div. (0->0), fcn. (485->11), ass. (0->116)
t63 = cos(qJ(2));
t56 = t63 ^ 2;
t130 = 0.2e1 * t56;
t61 = cos(qJ(4));
t101 = qJD(1) * t61;
t59 = sin(qJ(3));
t60 = sin(qJ(2));
t114 = t60 * t59;
t62 = cos(qJ(3));
t27 = t63 * t62 - t114;
t53 = qJD(2) + qJD(3);
t23 = t27 * t53;
t112 = t63 * t59;
t113 = t62 * t60;
t26 = t112 + t113;
t58 = sin(qJ(4));
t129 = pkin(3) * (t26 * t101 - t58 * t23);
t102 = qJD(1) * t58;
t128 = (-t26 * t102 - t61 * t23) * pkin(3);
t64 = pkin(1) + pkin(4);
t127 = t130 - 0.1e1;
t126 = pkin(2) * t60;
t65 = 0.2e1 * qJ(2);
t47 = sin(qJ(3) + t65);
t125 = pkin(3) * t47;
t124 = pkin(3) * t53;
t57 = qJ(2) + qJ(3);
t48 = sin(t57);
t123 = t48 * pkin(1);
t50 = t59 * pkin(2);
t122 = t62 * pkin(2);
t121 = t63 * pkin(2);
t11 = t60 * (t62 ^ 2 - 0.1e1 / 0.2e1) * t63 + (t56 - 0.1e1 / 0.2e1) * t62 * t59;
t120 = t53 * t11;
t119 = t53 * t59;
t118 = t53 * t62;
t54 = qJD(1) + qJD(4);
t117 = t54 * t58;
t116 = t54 * t61;
t55 = sin(t65);
t115 = t55 * pkin(2) ^ 2;
t111 = t64 * t48;
t34 = sin(0.2e1 * t57);
t110 = pkin(3) ^ 2 * t34;
t109 = pkin(2) * qJD(2);
t108 = pkin(3) * qJD(1);
t68 = t55 * pkin(2) / 0.2e1 + t60 * t64;
t88 = pkin(3) * t111;
t107 = (t110 / 0.2e1 + t88 + (t68 + t125) * pkin(2)) * qJD(1);
t93 = pkin(2) * t125;
t106 = (0.2e1 * t93 + t110 + t115 + 0.2e1 * (pkin(3) * t48 + t126) * pkin(1)) * qJD(1);
t105 = qJD(1) * t11;
t83 = t60 * t112;
t18 = (0.4e1 * t83 + (-0.4e1 * t56 + 0.2e1) * t62) * t62 + t127;
t104 = qJD(1) * t18;
t43 = pkin(1) + t121;
t103 = qJD(1) * t43;
t100 = qJD(1) * t63;
t99 = qJD(2) * t60;
t98 = qJD(2) * t63;
t41 = t62 * pkin(3) + pkin(2);
t21 = -pkin(3) * t114 + t41 * t63 + t64;
t97 = qJD(4) * t21;
t30 = pkin(1) * t63 + pkin(2) * t130 - pkin(2);
t42 = pkin(1) + 0.2e1 * t121;
t96 = (t42 * t113 + t30 * t59) * qJD(1);
t86 = pkin(3) * t112;
t24 = t60 * t41 + t86;
t95 = t24 * qJD(1);
t94 = t59 * qJD(3);
t92 = pkin(1) * t118;
t91 = t43 * t126;
t49 = qJD(2) + qJD(3) / 0.2e1;
t90 = t49 * t50;
t89 = pkin(3) * t26 * t54;
t87 = pkin(2) * t119;
t85 = t24 * t117;
t84 = t53 * t114;
t82 = pkin(1) * t99;
t81 = pkin(1) * t60 * qJD(1);
t80 = pkin(1) * t100;
t79 = t26 * t108;
t44 = t59 * t109;
t78 = pkin(3) * t94;
t77 = t26 * t103;
t76 = t27 * t103;
t75 = t60 * t98;
t74 = t108 / 0.2e1;
t73 = -0.2e1 * t49 * t122;
t72 = pkin(3) * t87;
t71 = t58 * t89;
t70 = pkin(2) * t78;
t32 = pkin(3) * t34;
t39 = pkin(2) * t47;
t14 = t32 + t39 + 0.2e1 * t123 + t50;
t69 = qJD(1) * t91;
t29 = pkin(3) * t118 + t109;
t15 = t60 * t29 + t53 * t86;
t46 = pkin(2) * t98;
t45 = t62 * t109;
t38 = pkin(3) * t44;
t35 = t127 * qJD(1);
t33 = t60 * t100;
t31 = cos(t57) * t124;
t25 = -t53 * t110 / 0.2e1;
t22 = t26 * t53;
t20 = t61 * t89;
t19 = t24 * t116;
t17 = -pkin(3) * t84 + t29 * t63;
t10 = t15 * t61;
t9 = t58 * t15;
t8 = -0.2e1 * t105;
t7 = 0.2e1 * t105;
t2 = t17 * t61 + t58 * t95;
t1 = t17 * t58 - t61 * t95;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t127 * qJD(2), 0, -t75, 0, 0, -t82, -pkin(1) * t98, 0, 0, 0.2e1 * t120, -t53 * t18, 0, -0.2e1 * t120, 0, 0, -0.2e1 * t56 * t90 + (-pkin(1) * t119 + t60 * t73) * t63 + t44 - t60 * t92, t56 * t73 + (0.2e1 * t60 * t90 - t92) * t63 + t45 + pkin(1) * t84, 0, -qJD(2) * t91, 0, 0, 0, 0, 0, 0, -t15, 0, 0, t25 - t123 * t124 - qJD(2) * t115 / 0.2e1 + (-t82 + (-t49 * t47 - t94 / 0.2e1) * pkin(3)) * pkin(2), 0, 0, 0, 0, 0, 0, -t58 * t97 - t10, -t61 * t97 + t9, 0, t25 - t49 * t93 - t53 * t88 - pkin(2) * (t78 / 0.2e1 + t68 * qJD(2)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, -t98, -t33, t99, 0, -t81, -t80, 0, 0, t7, -t104, -t23, t8, t22, 0, -t96, -qJD(1) * (-t42 * t114 + t30 * t62), t46, -t69, 0, 0, 0, 0, 0, 0, -t95, 0, t31 + t46, -t106 / 0.2e1, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t104, -t23, t8, t22, 0, -t77, -t76, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, t31, -t14 * t108 / 0.2e1, 0, 0, 0, 0, 0, 0, -t129, -t128, 0, -(t32 / 0.2e1 + t39 / 0.2e1 + t111 + t50 / 0.2e1) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t117, -t21 * t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t35, 0, t33, 0, 0, t81, t80, 0, 0, t8, t104, 0, t7, 0, 0, t96, qJD(1) * (t27 * pkin(1) + (t62 * t130 - t62 - 0.2e1 * t83) * pkin(2)), 0, t69, 0, 0, 0, 0, 0, 0, t95, 0, 0, t106 / 0.2e1, 0, 0, 0, 0, 0, 0, t19, -t85, 0, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(2) * t94, -qJD(3) * t122, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -pkin(2) * t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t104, 0, t7, 0, 0, t77, t76, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, t14 * t74, 0, 0, 0, 0, 0, 0, t20, -t71, 0, (0.2e1 * t48 * pkin(4) + t14) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t102 - t10, t21 * t101 + t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t128, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
