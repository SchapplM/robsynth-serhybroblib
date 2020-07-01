% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% MMD_reg [((2+1)*2/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbar1turnTE_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:33
% EndTime: 2020-06-27 16:22:53
% DurationCPUTime: 2.20s
% Computational Cost: add. (23444->139), mult. (32818->351), div. (970->15), fcn. (8852->4), ass. (0->147)
t175 = pkin(4) ^ 2;
t67 = cos(qJ(2));
t156 = pkin(2) * t67;
t75 = pkin(1) ^ 2;
t129 = -0.2e1 * pkin(1) * t156 + t75;
t74 = pkin(2) ^ 2;
t61 = t74 + t129;
t59 = 0.1e1 / t61 ^ 2;
t121 = 0.2e1 * t59;
t174 = 0.6e1 * t67;
t66 = sin(qJ(2));
t163 = -pkin(3) - pkin(4);
t54 = (pkin(2) - t163) * (pkin(2) + t163) + t129;
t162 = pkin(4) - pkin(3);
t55 = (pkin(2) - t162) * (pkin(2) + t162) + t129;
t139 = t54 * t55;
t76 = sqrt(-t139);
t136 = t66 * t76;
t128 = pkin(3) ^ 2 - t175;
t57 = t61 - t128;
t62 = pkin(1) - t156;
t39 = -pkin(2) * t136 + t57 * t62;
t34 = 0.1e1 / t39 ^ 2;
t159 = pkin(2) * t57;
t52 = t66 * t159;
t40 = t62 * t76 + t52;
t145 = t34 * t40;
t100 = pkin(1) * pkin(2) * (-t54 - t55);
t44 = t66 * t100;
t43 = qJD(2) * t44;
t48 = 0.1e1 / t76;
t143 = t48 * t43;
t116 = t66 * t143;
t123 = 0.2e1 * t62 * pkin(1);
t132 = t67 * t76;
t15 = (-t116 + (-t132 + (t57 + t123) * t66) * qJD(2)) * pkin(2);
t33 = 0.1e1 / t39;
t150 = t15 * t33 * t34;
t160 = pkin(1) * t74;
t115 = qJD(2) * t160;
t65 = t66 ^ 2;
t103 = t65 * t115;
t127 = t66 * qJD(2);
t113 = t76 * t127;
t126 = t67 * qJD(2);
t131 = pkin(2) * t113 + t126 * t159;
t141 = t48 * t62;
t16 = t43 * t141 + 0.2e1 * t103 + t131;
t36 = t40 ^ 2;
t28 = t34 * t36 + 0.1e1;
t173 = 0.1e1 / t28 ^ 2 * (t16 * t145 - t36 * t150);
t56 = t61 + t128;
t63 = pkin(1) * t67 - pkin(2);
t38 = -pkin(1) * t136 - t56 * t63;
t31 = 0.1e1 / t38 ^ 2;
t161 = pkin(1) * t56;
t53 = t66 * t161;
t41 = -t63 * t76 + t53;
t146 = t31 * t41;
t169 = -0.2e1 * t63;
t125 = pkin(2) * t169;
t14 = (-t116 + (-t132 + (t56 + t125) * t66) * qJD(2)) * pkin(1);
t30 = 0.1e1 / t38;
t151 = t14 * t30 * t31;
t137 = t65 * t75;
t112 = qJD(2) * t137;
t101 = pkin(2) * t112;
t130 = pkin(1) * t113 + t126 * t161;
t140 = t48 * t63;
t17 = -t43 * t140 + 0.2e1 * t101 + t130;
t37 = t41 ^ 2;
t29 = t31 * t37 + 0.1e1;
t172 = 0.1e1 / t29 ^ 2 * (t17 * t146 - t37 * t151);
t142 = t48 * t44;
t134 = t67 * t41;
t165 = t66 / 0.2e1;
t58 = 0.1e1 / t61;
t171 = (t134 / 0.2e1 + t38 * t165) * t58;
t170 = 0.2e1 * pkin(2);
t168 = -t58 / 0.2e1;
t167 = t58 / 0.2e1;
t166 = -t66 / 0.2e1;
t164 = -t67 / 0.2e1;
t158 = pkin(2) * t58;
t157 = pkin(2) * t59;
t155 = pkin(3) * t61;
t154 = pkin(4) * t61;
t122 = pkin(1) * t157;
t111 = t66 * t122;
t19 = t52 + (-t132 + (t123 - t142) * t66) * pkin(2);
t70 = 0.1e1 / pkin(4);
t10 = (t39 * t111 + t19 * t168) * t70;
t104 = t66 * t115;
t105 = 0.4e1 / t139 * t43 * t142;
t117 = t58 * t143;
t118 = t10 * t145;
t24 = 0.1e1 / t28;
t119 = t24 * t154;
t21 = t44 * t141 + 0.2e1 * t65 * t160 + (t57 * t67 + t136) * pkin(2);
t12 = (-t40 * t111 + t21 * t167) * t70;
t124 = 0.2e1 * t154;
t133 = t67 * t58;
t138 = t59 * t66;
t147 = t16 * t59;
t152 = t12 * t33;
t42 = (t67 * t100 - 0.4e1 * t74 * t137) * qJD(2);
t92 = -t105 / 0.4e1;
t84 = (t66 * t92 + 0.2e1 * (t42 * t166 - t67 * t43) * t48) * t58;
t114 = pkin(2) * t127;
t106 = pkin(1) * t114;
t87 = -0.2e1 * pkin(4) * t24 * t106;
t60 = t58 * t59;
t93 = t60 * t74 * t112;
t1 = 0.2e1 * t87 * t118 + 0.2e1 * (t24 * t150 + t34 * t173) * t10 * t40 * t124 + 0.2e1 * (t124 * t173 + t87) * t152 + 0.2e1 * ((-t10 * t16 + t12 * t15) * t34 + (-((t62 * t105 / 0.4e1 + t42 * t141 + t104 * t174) * t167 + 0.4e1 * t40 * t93 + ((t117 / 0.2e1 - pkin(1) * t147) * t66 + ((t132 + (-t57 + t142) * t66) * t167 + (-t21 * t66 - t40 * t67) * t59 * pkin(1)) * qJD(2)) * pkin(2)) * t33 - ((0.4e1 * t103 + t131) * t168 - 0.4e1 * t39 * t93 + (-t84 / 0.2e1 + (t15 * t138 + (-t62 * t133 + (t19 * t66 + t39 * t67) * t59) * qJD(2)) * pkin(1)) * pkin(2)) * t145) * t70) * t119;
t153 = t1 * t58;
t149 = t15 * t58;
t148 = t16 * t58;
t144 = t40 * t59;
t135 = t67 * t38;
t120 = t30 * t155;
t26 = 0.1e1 / t29;
t110 = t26 * t120;
t109 = t26 * t31 * t155;
t102 = t75 * t114;
t99 = -0.2e1 * t111;
t96 = t41 * t109;
t95 = t59 * t106;
t94 = t60 * t106;
t91 = pkin(3) * t26 * t106;
t90 = t102 * t121;
t18 = t53 + (-t132 + (t125 - t142) * t66) * pkin(1);
t73 = 0.1e1 / pkin(3);
t11 = (t18 * t58 + t38 * t99) * t73;
t20 = -t44 * t140 + t137 * t170 + (t56 * t67 + t136) * pkin(1);
t13 = (t20 * t58 + t41 * t99) * t73;
t4 = -t11 * t96 + t13 * t110 + 0.1e1;
t89 = t4 * t59 * t104;
t88 = 0.8e1 * t93;
t85 = (t41 * t165 - t135 / 0.2e1) * t58;
t71 = 0.1e1 / t175;
t23 = t73 * t85;
t22 = t73 * t171;
t6 = ((t14 * t164 + t17 * t165) * t58 + (t171 + (t66 * t135 - t41 * t65) * t122) * qJD(2)) * t73;
t5 = ((t14 * t166 + t17 * t164) * t58 + (t85 + (t66 * t134 + t38 * t65) * t122) * qJD(2)) * t73;
t3 = 0.2e1 * (-t118 - t152) * t119;
t2 = (((t102 * t174 - t42 * t140 + t63 * t92) * t58 + t41 * t88 + ((-0.2e1 * t17 * t157 + t117) * t66 + ((t132 + (-t56 + t142) * t66) * t58 + (-t20 * t66 - t134) * pkin(2) * t121) * qJD(2)) * pkin(1)) * t110 - ((0.4e1 * t101 + t130) * t58 + t38 * t88 + (t84 + (-0.2e1 * t14 * t138 + (t133 * t169 + (-t18 * t66 - t135) * t121) * qJD(2)) * pkin(2)) * pkin(1)) * t96) * t73 + (-t14 * t109 - 0.2e1 * t120 * t172 + 0.2e1 * t30 * t91) * t13 + (0.2e1 * (t26 * t151 + t31 * t172) * t41 * t155 - t17 * t109 - 0.2e1 * t91 * t146) * t11;
t7 = [0, 0, 0, 0.2e1 * t66 * t126, 0.2e1 * (t67 ^ 2 - t65) * qJD(2), 0, 0, 0, 0, 0, -0.2e1 * t22 * t5, -0.2e1 * t22 * t6 + 0.2e1 * t23 * t5, 0, 0.2e1 * t23 * t6, 0, 0, (-t23 * t127 + t6 * t67) * t170, (-t22 * t127 - t5 * t67) * t170, (t16 * t144 / 0.2e1 - t36 * t94) * t71, (-t15 * t144 / 0.2e1 + (-t147 / 0.2e1 + 0.2e1 * t40 * t94) * t39) * t71, 0, 0, 0, (-pkin(1) * t149 + t39 * t90) * t70, (-pkin(1) * t148 + t40 * t90) * t70; 0, 0, 0, 0, 0, t126, -t127, 0, 0, 0, 0, 0, -t2 * t22 + t4 * t5, 0, t2 * t23 + t4 * t6, 0, 0, 0, 0, 0, (t40 * t153 / 0.2e1 + (t148 / 0.2e1 - t40 * t95) * t3) * t70, (-t39 * t153 / 0.2e1 + (-t149 / 0.2e1 + t39 * t95) * t3) * t70, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t4 * t2, (0.2e1 * t38 * t89 + (-t14 * t4 - t2 * t38) * t158) * t73, (-0.2e1 * t41 * t89 + (t17 * t4 + t2 * t41) * t158) * t73, 0, 0, 0, 0, 0.2e1 * t3 * t1, 0, 0;];
MMD_reg = t7;
