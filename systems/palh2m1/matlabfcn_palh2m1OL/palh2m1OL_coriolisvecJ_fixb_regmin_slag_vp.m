% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m1OL_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:11
% EndTime: 2020-06-30 17:45:20
% DurationCPUTime: 1.95s
% Computational Cost: add. (2872->223), mult. (7498->333), div. (0->0), fcn. (6176->8), ass. (0->139)
t106 = cos(qJ(5));
t157 = qJD(5) * t106;
t103 = sin(qJ(4));
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t161 = qJD(1) * t109;
t146 = t108 * t161;
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t162 = qJD(1) * t105;
t148 = t104 * t162;
t78 = -t146 + t148;
t124 = t104 * t109 + t105 * t108;
t79 = t124 * qJD(1);
t50 = t103 * t79 + t107 * t78;
t196 = t106 * t50;
t199 = t157 - t196;
t86 = qJD(1) * pkin(1) + pkin(2) * t161;
t65 = -pkin(3) * t78 + t86;
t183 = t50 * t65;
t173 = qJD(2) * pkin(2);
t153 = t108 * t173;
t134 = qJD(3) * t153;
t154 = t104 * t173;
t135 = t103 * t154;
t174 = (qJD(3) + qJD(4)) * t135;
t98 = qJD(2) + qJD(3);
t82 = pkin(3) * t98 + t153;
t39 = t107 * (qJD(4) * t82 + t134) - t174;
t198 = -t183 - t39;
t163 = -qJD(5) + t50;
t197 = qJD(5) + t163;
t102 = sin(qJ(5));
t127 = t103 * t78 - t107 * t79;
t96 = qJD(4) + t98;
t36 = t102 * t96 + t106 * t127;
t158 = qJD(5) * t102;
t122 = t124 * qJD(3);
t62 = qJD(2) * t124 + t122;
t114 = t62 * qJD(1);
t159 = qJD(4) * t107;
t160 = qJD(4) * t103;
t156 = qJD(1) * qJD(2);
t144 = t105 * t156;
t55 = qJD(3) * t148 + t104 * t144 - t146 * t98;
t21 = t103 * t114 + t107 * t55 + t159 * t78 + t160 * t79;
t11 = t106 * t21 - t127 * t158 + t157 * t96;
t9 = t11 * t102;
t4 = t199 * t36 + t9;
t22 = qJD(4) * t127 + t103 * t55 - t107 * t114;
t42 = t163 * t157;
t176 = t102 * t22 - t42;
t3 = -t127 * t36 + t163 * t196 + t176;
t195 = t102 * t163;
t20 = t106 * t22;
t34 = t102 * t127 - t106 * t96;
t2 = t127 * t34 - t163 * t195 + t20;
t12 = qJD(5) * t36 + t102 * t21;
t1 = -t102 * t12 + t11 * t106 + t195 * t36 - t199 * t34;
t66 = t107 * t82 - t135;
t63 = -pkin(4) * t96 - t66;
t180 = t63 * t50;
t182 = t127 * t50;
t15 = t127 ^ 2 - t50 ^ 2;
t29 = pkin(4) * t127 - pkin(6) * t50;
t13 = -t50 * t96 + t21;
t184 = t163 * t127;
t23 = -pkin(4) * t50 - pkin(6) * t127 + t65;
t67 = t103 * t82 + t107 * t154;
t64 = pkin(6) * t96 + t67;
t128 = t102 * t64 - t106 * t23;
t145 = t127 * t128 + t158 * t63;
t150 = t82 * t160;
t165 = t104 * t107;
t126 = t103 * t108 + t165;
t189 = qJD(3) * t126 + t104 * t159;
t112 = -t127 * t65 - t173 * t189 - t150;
t14 = t127 * t96 - t22;
t80 = t104 * t105 - t108 * t109;
t59 = -t103 * t124 - t107 * t80;
t61 = t98 * t80;
t25 = -qJD(4) * t59 + t103 * t62 + t107 * t61;
t60 = t103 * t80 - t107 * t124;
t94 = pkin(2) * t109 + pkin(1);
t70 = -pkin(3) * t80 + t94;
t28 = pkin(4) * t59 - pkin(6) * t60 + t70;
t115 = t189 * pkin(2);
t40 = qJD(2) * t115 + t150;
t191 = qJD(5) * t28 * t163 - (qJD(5) * t23 + t39) * t59 + t63 * t25 + t40 * t60;
t17 = t102 * t23 + t106 * t64;
t133 = t102 * t40 + t127 * t17 + t157 * t63;
t188 = pkin(3) * t79;
t186 = t105 * pkin(2);
t185 = t22 * t60;
t181 = t60 * t63;
t179 = t78 * t86;
t178 = t79 * t78;
t177 = t79 * t86;
t168 = -t105 ^ 2 + t109 ^ 2;
t166 = t103 * t104;
t111 = qJD(1) ^ 2;
t164 = t109 * t111;
t152 = t105 * t173;
t151 = pkin(2) * t162;
t149 = t60 * t158;
t27 = -t188 + t29;
t93 = pkin(2) * t108 + pkin(3);
t74 = pkin(2) * t165 + t103 * t93 + pkin(6);
t141 = qJD(5) * t74 - t151 + t27;
t91 = pkin(3) * t103 + pkin(6);
t140 = qJD(5) * t91 + t27;
t137 = -0.2e1 * pkin(1) * t156;
t136 = qJD(3) * (-qJD(2) - t98);
t75 = t126 * t173;
t132 = pkin(3) * t160 - t75;
t131 = (-qJD(3) + t98) * t173;
t26 = qJD(4) * t60 + t103 * t61 - t107 * t62;
t52 = -pkin(3) * t62 - t152;
t130 = t28 * t22 - (pkin(4) * t26 - pkin(6) * t25 + t52) * t163;
t129 = -t163 * t25 + t185;
t125 = t107 * t108 - t166;
t123 = t94 * t124;
t53 = t93 * t159 + (qJD(3) * t125 - t104 * t160) * pkin(2);
t121 = t163 * t53 - t22 * t74 - t180;
t76 = t125 * t173;
t117 = -t91 * t22 - t180 - (-pkin(3) * t159 + t76) * t163;
t41 = (-pkin(3) * t122 + (-pkin(3) * t124 - t186) * qJD(2)) * qJD(1);
t110 = qJD(2) ^ 2;
t92 = -pkin(3) * t107 - pkin(4);
t73 = pkin(2) * t166 - t107 * t93 - pkin(4);
t68 = -t151 - t188;
t54 = t160 * t93 + t115;
t38 = -t78 ^ 2 + t79 ^ 2;
t31 = -t79 * t98 + t114;
t30 = -t78 * t98 + t55;
t6 = t22 * pkin(4) - t21 * pkin(6) + t41;
t5 = t106 * t6;
t7 = [0, 0, 0, 0.2e1 * t109 * t144, 0.2e1 * t168 * t156, -t110 * t109, t110 * t105, 0, t105 * t137, t109 * t137, -t124 * t55 - t61 * t79, -t114 * t124 + t55 * t80 + t61 * t78 - t79 * t62, t61 * t98, t62 * t98, 0, t78 * t152 - t86 * t62 + (-qJD(3) * t123 + (t186 * t80 - t123) * qJD(2)) * qJD(1), 0.2e1 * t152 * t79 + t55 * t94 + t61 * t86, t127 * t25 + t21 * t60, -t127 * t26 - t21 * t59 + t25 * t50 - t185, t25 * t96, -t26 * t96, 0, t22 * t70 + t26 * t65 + t41 * t59 - t50 * t52, t127 * t52 + t21 * t70 + t25 * t65 + t41 * t60, -t36 * t149 + (t11 * t60 + t25 * t36) * t106, (-t102 * t36 - t106 * t34) * t25 + (-t9 - t106 * t12 + (t102 * t34 - t106 * t36) * qJD(5)) * t60, t106 * t129 + t11 * t59 + t149 * t163 + t26 * t36, -t102 * t129 - t12 * t59 - t26 * t34 + t42 * t60, -t163 * t26 + t22 * t59, -t128 * t26 + t5 * t59 + ((-t59 * t64 + t181) * qJD(5) + t130) * t106 + t191 * t102, -t17 * t26 + t191 * t106 + (-(-qJD(5) * t64 + t6) * t59 - qJD(5) * t181 - t130) * t102; 0, 0, 0, -t105 * t164, -t168 * t111, 0, 0, 0, t111 * pkin(1) * t105, pkin(1) * t164, t178, t38, t30, t31, 0, t177 + (t104 * t136 - t162 * t78) * pkin(2), -t179 + (t108 * t136 - t162 * t79) * pkin(2), -t182, t15, t13, t14, 0, t50 * t68 - t54 * t96 + t112, -t127 * t68 - t53 * t96 + t198, t4, t1, t3, t2, t184, t12 * t73 + t34 * t54 + (t141 * t163 - t40) * t106 + t121 * t102 + t145, t106 * t121 + t11 * t73 - t141 * t195 + t36 * t54 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t38, t30, t31, 0, t104 * t131 + t177, t108 * t131 - t179, -t182, t15, t13, t14, 0, t75 * t96 + (-t160 * t96 - t50 * t79) * pkin(3) + t112, t127 * t188 - t183 + t76 * t96 + (-t134 + (-pkin(3) * t96 - t82) * qJD(4)) * t107 + t174, t4, t1, t3, t2, t184, t92 * t12 + t132 * t34 + (t140 * t163 - t40) * t106 + t117 * t102 + t145, t106 * t117 + t92 * t11 + t132 * t36 - t140 * t195 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t15, t13, t14, 0, t67 * t96 + t112, t66 * t96 + t198, t4, t1, t3, t2, t184, -pkin(4) * t12 - t40 * t106 + (-t102 * t66 + t106 * t29) * t163 - t67 * t34 - t102 * t180 - t176 * pkin(6) + t145, -pkin(4) * t11 - (t102 * t29 + t106 * t66) * t163 - t67 * t36 - t63 * t196 + (-t158 * t163 - t20) * pkin(6) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t34, -t34 ^ 2 + t36 ^ 2, -t163 * t34 + t11, -t163 * t36 - t12, t22, -t102 * t39 - t17 * t197 - t36 * t63 + t5, -t102 * t6 - t106 * t39 + t128 * t197 + t34 * t63;];
tauc_reg = t7;
