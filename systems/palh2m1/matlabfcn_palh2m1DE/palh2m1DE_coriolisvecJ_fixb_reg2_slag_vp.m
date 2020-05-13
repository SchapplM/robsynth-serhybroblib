% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m1DE_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:30
% EndTime: 2020-05-02 23:52:33
% DurationCPUTime: 0.51s
% Computational Cost: add. (386->94), mult. (1014->197), div. (0->0), fcn. (627->11), ass. (0->102)
t47 = qJD(2) + qJD(3);
t56 = sin(qJ(3));
t60 = cos(qJ(2));
t110 = t60 * t56;
t57 = sin(qJ(2));
t59 = cos(qJ(3));
t72 = t59 * t57 + t110;
t14 = t47 * t72;
t41 = t60 * pkin(2) + pkin(1);
t128 = pkin(2) * t41;
t48 = qJD(1) + qJD(4);
t58 = cos(qJ(4));
t127 = t48 * t58;
t109 = t60 * t59;
t112 = t57 * t56;
t24 = -0.2e1 * t109 * t112;
t52 = t59 ^ 2;
t53 = t60 ^ 2;
t126 = t24 + (-t56 ^ 2 + t52) * (t53 - 0.1e1 / 0.2e1);
t106 = t57 ^ 2 - t53;
t125 = t24 - t106 * (t52 - 0.1e1 / 0.2e1);
t64 = (qJD(1) ^ 2);
t124 = -2 * t64;
t123 = -2 * qJD(1);
t122 = 2 * qJD(1);
t121 = t64 / 0.2e1;
t61 = pkin(1) + pkin(4);
t120 = pkin(3) * t56;
t119 = pkin(3) * t64;
t105 = pkin(3) * qJD(3);
t23 = t109 - t112;
t20 = t23 * t105;
t40 = t59 * pkin(3) + pkin(2);
t75 = -pkin(3) * t112 + t40 * t60;
t101 = qJD(2) * t60;
t102 = qJD(2) * t57;
t69 = t23 * qJD(3);
t9 = t40 * t101 + (-t102 * t56 + t69) * pkin(3);
t118 = -qJD(2) * t75 - t20 + t9;
t17 = t75 + t61;
t117 = t17 * t58;
t54 = qJ(2) + qJ(3);
t116 = sin(0.2e1 * t54) * pkin(3) ^ 2;
t115 = t41 * t64;
t65 = 0.2e1 * qJ(2);
t114 = sin(t65) * pkin(2) ^ 2;
t19 = pkin(3) * t110 + t57 * t40;
t11 = t19 * qJD(2) + t105 * t72;
t55 = sin(qJ(4));
t113 = t55 * t11;
t111 = t57 * t64;
t108 = t60 * t64;
t63 = qJD(2) ^ 2;
t107 = t63 * t60;
t104 = qJD(1) * t55;
t103 = qJD(1) * t58;
t100 = qJD(4) * t58;
t99 = -qJD(1) - t48;
t98 = -qJD(4) + t48;
t70 = t23 * qJD(2);
t13 = t69 + t70;
t3 = t9 * qJD(2) + t105 * t13;
t8 = -t40 * t102 + (-qJD(3) * t72 - t101 * t56) * pkin(3);
t96 = t11 * t100 + t8 * t103 + t55 * t3;
t95 = qJD(1) * qJD(2);
t94 = pkin(2) * t120;
t93 = pkin(2) * t119;
t45 = sin(t54);
t92 = t45 * t119;
t44 = sin(qJ(3) + t65);
t91 = pkin(3) * (qJD(3) + 0.2e1 * qJD(2)) * t44;
t90 = pkin(3) * t47 * t45;
t89 = qJD(3) ^ 2 * t120;
t88 = pkin(1) * t111;
t87 = t57 * t108;
t86 = t64 * t23 * t72;
t85 = -0.2e1 * t102;
t84 = pkin(2) * t102;
t83 = t17 * t104;
t81 = t60 * t95;
t80 = t99 * t55;
t25 = t116 * t121;
t30 = t44 * t93;
t79 = t25 + t30 / 0.2e1 + t56 * t93 / 0.2e1 + t63 * t94;
t35 = -0.2e1 * qJD(3) * t94;
t78 = qJD(2) * t35 + t114 * t121 + t25 + t30;
t77 = qJD(3) * (-qJD(2) - t47);
t76 = qJD(1) * t85;
t74 = t57 * t81;
t73 = pkin(2) * qJD(2) * (-qJD(3) + t47);
t71 = t35 / 0.2e1 - qJD(2) * t114;
t42 = pkin(2) * t107;
t31 = pkin(1) * t92;
t26 = t61 * t92;
t21 = t72 * pkin(3);
t18 = qJD(1) * t47 * t116;
t16 = t23 * t115;
t15 = t72 * t115;
t12 = pkin(3) * t70 + t20;
t4 = (t23 * t47 - t13) * qJD(1);
t2 = t3 * t58;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, -0.2e1 * t106 * t95, -t107, -0.2e1 * t74, t63 * t57, 0, pkin(1) * t76, -0.2e1 * pkin(1) * t81, 0, 0, -t72 * t13 * t123, 0.4e1 * (t125 * qJD(2) + t126 * qJD(3)) * qJD(1), -t47 * t13, t23 * t14 * t123, t47 * t14, 0, (-t14 * t41 - t23 * t84) * t122, (-t13 * t41 + t72 * t84) * t122, t42, t76 * t128, 0, 0, 0, 0, 0, 0, t8 * t122, 0, pkin(3) * t47 ^ 2 * cos(t54) + t42, -t18 + (-pkin(2) * t91 + 0.2e1 * (-t84 - t90) * pkin(1) + t71) * qJD(1), 0, 0, 0, 0, 0, 0, qJD(4) * t17 * t80 + t8 * t127 + t96, t2 + t8 * t80 + (t117 * t99 - t113) * qJD(4), 0, -t18 + (-0.2e1 * t61 * t90 + (t61 * t85 - t91) * pkin(2) + t71) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t106 * t64, 0, t87, 0, 0, t88, pkin(1) * t108, 0, 0, -t86, t125 * t124, t4, t86, 0, 0, t15 + (t111 * t23 + t56 * t77) * pkin(2), t16 + (-t111 * t72 + t59 * t77) * pkin(2), 0, t111 * t128, 0, 0, 0, 0, 0, 0, t19 * t64, 0, 0, t31 + (t88 - t89) * pkin(2) + t78, 0, 0, 0, 0, 0, 0, (t118 * t55 + t19 * t127) * t48, (-t19 * t48 * t55 + t118 * t58) * t48, 0, t26 + (t111 * t61 - t89) * pkin(2) + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t126 * t124, t4, t86, 0, 0, t56 * t73 + t15, t59 * t73 + t16, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t64, 0, 0, t31 + t79, 0, 0, 0, 0, 0, 0, (t21 * t103 - t55 * t12 + (t100 * t72 + t13 * t55) * pkin(3)) * t48, (-t21 * t104 - t12 * t58 + (-qJD(4) * t55 * t72 + t13 * t58) * pkin(3)) * t48, 0, t26 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t83 - (t11 * t58 - t83) * t48 + t96, t2 + t98 * t113 + (t117 * t98 - t55 * t8) * qJD(1), 0, 0;];
tauc_reg = t1;
