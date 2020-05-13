% Calculate inertial parameters regressor of coriolis joint torque vector for
% palh2m2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:20
% EndTime: 2019-06-06 14:45:23
% DurationCPUTime: 0.88s
% Computational Cost: add. (1383->125), mult. (4794->233), div. (0->0), fcn. (2745->6), ass. (0->110)
t60 = sin(qJ(2));
t107 = qJD(2) * t60;
t95 = pkin(4) * t107;
t128 = 0.2e1 * t95;
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t63 = cos(qJ(2));
t72 = t59 * t60 + t62 * t63;
t69 = t72 * qJD(3);
t70 = t72 * qJD(2);
t29 = (t69 - t70) * pkin(4);
t43 = t72 * pkin(4);
t40 = qJD(2) * t43;
t127 = t29 + t40;
t54 = t59 ^ 2;
t56 = t62 ^ 2;
t126 = t54 + t56;
t101 = qJD(1) * qJD(2);
t125 = -0.2e1 * t101;
t65 = qJD(2) ^ 2;
t124 = t72 * t65;
t55 = t60 ^ 2;
t57 = t63 ^ 2;
t53 = t55 + t57;
t48 = qJD(1) * t53;
t68 = t48 * t128;
t116 = t59 * t63;
t94 = qJD(2) * t116;
t51 = pkin(4) * t94;
t93 = t62 * t107;
t85 = pkin(4) * t93;
t41 = -t51 + t85;
t37 = t62 * t41;
t38 = qJD(3) * pkin(5) + t40;
t22 = -t59 * t38 - t37;
t90 = t60 * t101;
t121 = t48 * t59;
t99 = pkin(5) * t121;
t34 = pkin(4) * t90 + qJD(3) * t99;
t58 = sin(qJ(4));
t117 = t59 * t41;
t21 = t62 * t38 - t117;
t106 = qJD(3) * t21;
t25 = (qJD(2) * t69 - t124) * pkin(4);
t114 = t65 * t63;
t100 = pkin(4) * t114;
t115 = t65 * t60;
t26 = qJD(3) * t85 + t59 * t100 + (-qJD(3) * t94 - t62 * t115) * pkin(4);
t75 = t62 * t25 - t59 * t26;
t6 = t75 - t106;
t61 = cos(qJ(4));
t89 = -t63 * pkin(4) - pkin(1);
t49 = t89 * qJD(1);
t122 = pkin(5) * t62;
t88 = -pkin(2) - t122;
t27 = t88 * t48 + t49;
t33 = t126 * t48;
t17 = -t33 * pkin(3) + t27;
t76 = t61 * t17 + t58 * t22;
t1 = -t76 * qJD(4) - t58 * t34 + t61 * t6;
t103 = qJD(4) * t58;
t2 = -t58 * t6 + t17 * t103 + (-qJD(4) * t22 - t34) * t61;
t9 = -t58 * t17 + t61 * t22;
t80 = t58 * t9 - t61 * t76;
t123 = -t80 * qJD(4) + t1 * t61 - t2 * t58;
t120 = t48 * t62;
t119 = t53 * qJD(3) ^ 2;
t73 = -t60 * t62 + t116;
t30 = t51 + (-t73 * qJD(3) - t93) * pkin(4);
t44 = t73 * pkin(4);
t74 = t62 * t43 + t59 * t44;
t10 = -t74 * qJD(3) + t62 * t29 - t59 * t30;
t24 = -t62 * t40 + t117;
t113 = t10 - t24;
t111 = -t38 + t40;
t110 = t54 - t56;
t109 = t55 - t57;
t66 = qJD(1) ^ 2;
t108 = pkin(1) * t66;
t105 = qJD(3) * t53;
t104 = qJD(3) * t62;
t102 = qJD(4) * t61;
t47 = t48 ^ 2;
t98 = t59 * t47 * t62;
t97 = t60 * t66 * t63;
t96 = pkin(4) * qJD(1) * t60;
t92 = t59 * t105;
t7 = t22 * qJD(3) + t59 * t25 + t62 * t26;
t91 = -t21 * t111 * t59 + t7 * t122;
t86 = t48 * t96;
t28 = -t59 * t43 + t62 * t44;
t84 = t7 * t74 + (t28 * qJD(3) + t127 * t59 + t62 * t30 + t37) * t21;
t83 = t63 * t90;
t82 = pkin(1) * t125;
t79 = -t58 * t76 - t61 * t9;
t77 = t92 * t120;
t35 = -t48 * pkin(2) + t49;
t42 = -t53 * pkin(2) + t89;
t71 = qJD(3) * (t35 * t53 + t42 * t48);
t67 = t88 * t53 + t89;
t39 = pkin(5) * t92 + t95;
t36 = t126 * t53;
t32 = qJD(4) + t33;
t20 = -t36 * pkin(3) + t67;
t19 = t61 * t24 - t58 * t96;
t18 = -t58 * t24 - t61 * t96;
t13 = t111 * t62;
t12 = t61 * t13 - t58 * t99;
t11 = -t58 * t13 - t61 * t99;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t83, t109 * t125, t114, -0.2e1 * t83, -t115, 0, t60 * t82, t63 * t82, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, -t53 * t100, t49 * t128, 0.2e1 * t77, -0.2e1 * t110 * t48 * t105, t62 * t119, -0.2e1 * t77, -t59 * t119, 0, t59 * t71 - t62 * t68, t59 * t68 + t62 * t71, (t24 * qJD(3) + t75) * t53, (qJD(1) * t42 + t35) * t95, 0, 0, 0, 0, 0, 0, -t39 * t33 - t34 * t36, 0, t6 * t36, t27 * t39 + t34 * t67, 0, 0, 0, 0, 0, 0, t2 * t36 + (t20 * t103 - t39 * t61) * t32, -t1 * t36 + (t20 * t102 + t39 * t58) * t32, 0, -t80 * t39 + (t79 * qJD(4) - t1 * t58 - t2 * t61) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t109 * t66, 0, t97, 0, 0, t60 * t108, t63 * t108, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, -t49 * t96, 0, 0, 0, 0, 0, 0, t62 * t86 + (t30 + t41) * qJD(3), -qJD(3) * t127 - t59 * t86, (-t117 + (-qJD(3) * t44 - t30) * t59 + (-qJD(3) * t43 + t127) * t62) * t48, t25 * t44 + t26 * t43 - t41 * t29 + t40 * t30 - t35 * t96, 0, 0, 0, 0, 0, 0, t33 * t96, 0, t113 * t33, t113 * t22 - t27 * t96 + t6 * t28 + t84, 0, 0, 0, 0, 0, 0, (-t10 * t58 - t28 * t102 - t18) * t32, (-t10 * t61 + t28 * t103 + t19) * t32, 0, -t79 * t10 + t123 * t28 + t18 * t76 - t9 * t19 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t110 * t47, 0, t98, 0, 0, -t41 * qJD(3) - t35 * t121 + t26, -t35 * t120 + t40 * qJD(3) + (-qJD(3) * t70 + t124) * pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, t33 * t99, 0, (-pkin(5) * t104 - t13) * t33, -t22 * t13 + (-t22 * t104 + (-t27 * t48 - t106 - t6) * t59) * pkin(5) + t91, 0, 0, 0, 0, 0, 0, (-t11 + (t59 * t102 + t58 * t104) * pkin(5)) * t32, (t12 + (-t59 * t103 + t61 * t104) * pkin(5)) * t32, 0, t76 * t11 - t9 * t12 + (t79 * t104 + (-t106 - t123) * t59) * pkin(5) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t32 + t2, -t32 * t76 - t1, 0, 0;];
tauc_reg  = t3;
