% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% fourbar1turnDE2
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
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbar1turnDE2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:30
% EndTime: 2020-06-27 16:48:50
% DurationCPUTime: 3.08s
% Computational Cost: add. (34152->156), mult. (48579->397), div. (1709->21), fcn. (13095->8), ass. (0->168)
t205 = pkin(4) ^ 2;
t204 = pkin(3) ^ 2;
t78 = sin(qJ(2));
t79 = cos(qJ(2));
t183 = pkin(2) * t79;
t88 = pkin(1) ^ 2;
t151 = -0.2e1 * pkin(1) * t183 + t88;
t190 = -pkin(3) - pkin(4);
t66 = (pkin(2) - t190) * (pkin(2) + t190) + t151;
t189 = pkin(4) - pkin(3);
t67 = (pkin(2) - t189) * (pkin(2) + t189) + t151;
t162 = t66 * t67;
t91 = sqrt(-t162);
t157 = t78 * t91;
t150 = t204 - t205;
t87 = pkin(2) ^ 2;
t73 = t87 + t151;
t69 = t73 - t150;
t74 = pkin(1) - t183;
t51 = -pkin(2) * t157 + t69 * t74;
t203 = 0.2e1 * t51;
t202 = 0.6e1 * t79;
t46 = 0.1e1 / t51 ^ 2;
t185 = pkin(2) * t69;
t64 = t78 * t185;
t52 = t74 * t91 + t64;
t172 = t46 * t52;
t193 = pkin(1) * pkin(2);
t117 = (-t66 - t67) * t193;
t56 = t78 * t117;
t55 = qJD(2) * t56;
t60 = 0.1e1 / t91;
t166 = t60 * t55;
t135 = t78 * t166;
t145 = 0.2e1 * t74 * pkin(1);
t156 = t79 * t91;
t19 = (-t135 + (-t156 + (t69 + t145) * t78) * qJD(2)) * pkin(2);
t44 = t51 ^ 2;
t45 = 0.1e1 / t51;
t175 = t19 * t45 / t44;
t186 = pkin(1) * t87;
t134 = qJD(2) * t186;
t77 = t78 ^ 2;
t119 = t77 * t134;
t149 = qJD(2) * t78;
t131 = t91 * t149;
t148 = qJD(2) * t79;
t153 = pkin(2) * t131 + t148 * t185;
t164 = t60 * t74;
t20 = t164 * t55 + 0.2e1 * t119 + t153;
t48 = t52 ^ 2;
t38 = t46 * t48 + 0.1e1;
t201 = 0.1e1 / t38 ^ 2 * (t172 * t20 - t175 * t48);
t68 = t73 + t150;
t75 = pkin(1) * t79 - pkin(2);
t50 = -pkin(1) * t157 - t68 * t75;
t42 = 0.1e1 / t50 ^ 2;
t188 = pkin(1) * t68;
t65 = t78 * t188;
t53 = -t75 * t91 + t65;
t173 = t42 * t53;
t194 = -0.2e1 * t75;
t147 = pkin(2) * t194;
t18 = (-t135 + (-t156 + (t68 + t147) * t78) * qJD(2)) * pkin(1);
t40 = t50 ^ 2;
t41 = 0.1e1 / t50;
t176 = t18 * t41 / t40;
t159 = t77 * t88;
t132 = qJD(2) * t159;
t118 = pkin(2) * t132;
t152 = pkin(1) * t131 + t148 * t188;
t163 = t60 * t75;
t21 = -t163 * t55 + 0.2e1 * t118 + t152;
t49 = t53 ^ 2;
t39 = t42 * t49 + 0.1e1;
t200 = (t173 * t21 - t176 * t49) / t39 ^ 2;
t155 = t40 + t49;
t71 = 0.1e1 / t73 ^ 2;
t86 = 0.1e1 / t204;
t33 = t155 * t86 * t71;
t28 = t33 ^ (-0.1e1 / 0.2e1);
t85 = 0.1e1 / pkin(3);
t199 = t28 * t85;
t165 = t60 * t56;
t167 = t53 * t79;
t171 = t50 * t78;
t70 = 0.1e1 / t73;
t198 = (t167 + t171) * t70;
t160 = t71 * t78;
t129 = t160 * t193;
t192 = -t70 / 0.2e1;
t23 = t64 + (-t156 + (t145 - t165) * t78) * pkin(2);
t82 = 0.1e1 / pkin(4);
t14 = (t129 * t51 + t192 * t23) * t82;
t191 = t70 / 0.2e1;
t25 = t56 * t164 + 0.2e1 * t77 * t186 + (t69 * t79 + t157) * pkin(2);
t16 = (-t129 * t52 + t191 * t25) * t82;
t177 = t16 * t45;
t197 = t14 * t172 + t177;
t196 = 0.2e1 * pkin(2);
t154 = t44 + t48;
t83 = 0.1e1 / t205;
t32 = t154 * t83 * t71;
t30 = 0.1e1 / t32;
t26 = t32 ^ (-0.1e1 / 0.2e1);
t195 = -0.2e1 * t26;
t187 = pkin(1) * t70;
t184 = pkin(2) * t70;
t182 = pkin(3) * t73;
t181 = pkin(4) * t73;
t72 = t70 * t71;
t109 = t72 * t87 * t132;
t121 = t78 * t134;
t122 = 0.4e1 / t162 * t55 * t165;
t133 = pkin(2) * t149;
t123 = pkin(1) * t133;
t136 = t70 * t166;
t34 = 0.1e1 / t38;
t140 = t34 * t181;
t161 = t70 * t79;
t174 = t20 * t71;
t54 = (t117 * t79 - 0.4e1 * t159 * t87) * qJD(2);
t108 = -t122 / 0.4e1;
t99 = (-0.2e1 * t79 * t166 + (-t54 * t60 + t108) * t78) * t70;
t1 = -0.4e1 * t197 * pkin(4) * t34 * t123 + 0.4e1 * (t177 * t201 + (t34 * t175 + t46 * t201) * t14 * t52) * t181 + 0.2e1 * ((-t14 * t20 + t16 * t19) * t46 + (-((t74 * t122 / 0.4e1 + t54 * t164 + t121 * t202) * t191 + 0.4e1 * t52 * t109 + ((t136 / 0.2e1 - pkin(1) * t174) * t78 + ((t156 + (-t69 + t165) * t78) * t191 + (-t25 * t78 - t52 * t79) * t71 * pkin(1)) * qJD(2)) * pkin(2)) * t45 - ((0.4e1 * t119 + t153) * t192 - 0.4e1 * t51 * t109 + (-t99 / 0.2e1 + (t19 * t160 + (-t74 * t161 + (t23 * t78 + t51 * t79) * t71) * qJD(2)) * pkin(1)) * pkin(2)) * t172) * t82) * t140;
t180 = t1 * t70;
t110 = t72 * t123;
t103 = 0.4e1 * t110;
t142 = 0.2e1 * t71;
t7 = ((t19 * t51 + t20 * t52) * t142 - t154 * t103) * t83;
t179 = t26 * t30 * t7;
t178 = t28 / t33 * ((t18 * t50 + t21 * t53) * t142 - t155 * t103) * t86;
t170 = t50 * t79;
t169 = t52 * t71;
t168 = t53 * t78;
t158 = t78 * t79;
t36 = 0.1e1 / t39;
t126 = t36 * t42 * t182;
t112 = t53 * t126;
t141 = t41 * t182;
t127 = t36 * t141;
t116 = -0.2e1 * t129;
t22 = t65 + (-t156 + (t147 - t165) * t78) * pkin(1);
t15 = (t116 * t50 + t22 * t70) * t85;
t24 = -t56 * t163 + t159 * t196 + (t68 * t79 + t157) * pkin(1);
t17 = (t116 * t53 + t24 * t70) * t85;
t6 = -t112 * t15 + t127 * t17 + 0.1e1;
t144 = t6 * t178;
t143 = 0.2e1 * t28;
t139 = t70 * t178;
t138 = 0.1e1 / t32 ^ 2 * t7 * t71;
t130 = pkin(2) * t142;
t5 = -0.2e1 * t197 * t140;
t128 = t5 * t70 * t179;
t120 = t88 * t133;
t115 = pkin(1) * t130;
t111 = t71 * t123;
t106 = pkin(3) * t36 * t123;
t105 = 0.8e1 * t109;
t102 = 0.4e1 * t26 * t71 * t120;
t101 = (t168 - t170) * t70;
t100 = t28 * t6 * t71 * t121;
t13 = t101 * t199;
t12 = t198 * t199;
t4 = ((-t168 / 0.2e1 + t170 / 0.2e1) * t139 + ((-t18 * t79 + t21 * t78) * t70 + (t198 + (t158 * t50 - t53 * t77) * t115) * qJD(2)) * t28) * t85;
t3 = ((t171 / 0.2e1 + t167 / 0.2e1) * t139 + ((-t18 * t78 - t21 * t79) * t70 + (t101 + (t158 * t53 + t50 * t77) * t115) * qJD(2)) * t28) * t85;
t2 = (((t75 * t108 + t120 * t202 - t54 * t163) * t70 + t53 * t105 + ((-0.2e1 * t21 * t71 * pkin(2) + t136) * t78 + ((t156 + (-t68 + t165) * t78) * t70 + (-t24 * t78 - t167) * t130) * qJD(2)) * pkin(1)) * t127 - ((0.4e1 * t118 + t152) * t70 + t50 * t105 + (t99 + (-0.2e1 * t18 * t160 + (t161 * t194 + (-t22 * t78 - t170) * t142) * qJD(2)) * pkin(2)) * pkin(1)) * t112) * t85 + (0.2e1 * t41 * t106 - t18 * t126 - 0.2e1 * t141 * t200) * t17 + (-0.2e1 * t106 * t173 - t21 * t126 + 0.2e1 * (t36 * t176 + t42 * t200) * t53 * t182) * t15;
t8 = [0, 0, 0, 0.2e1 * t78 * t148, 0.2e1 * (t79 ^ 2 - t77) * qJD(2), 0, 0, 0, 0, 0, -0.2e1 * t12 * t3, -0.2e1 * t12 * t4 + 0.2e1 * t13 * t3, 0, 0.2e1 * t13 * t4, 0, 0, (-t13 * t149 + t4 * t79) * t196, (-t12 * t149 - t3 * t79) * t196, (-t48 * t138 + (-0.4e1 * t110 * t48 + 0.2e1 * t169 * t20) * t30) * t83, (t52 * t138 * t203 + (-0.2e1 * t19 * t169 + (0.8e1 * t110 * t52 - 0.2e1 * t174) * t51) * t30) * t83, 0, 0, 0, (t51 * t102 + (t179 * t51 + t19 * t195) * t187) * t82, (t52 * t102 + (t179 * t52 + t195 * t20) * t187) * t82; 0, 0, 0, 0, 0, t148, -t149, 0, 0, 0, 0, 0, -t12 * t2 + t3 * t6, 0, t13 * t2 + t4 * t6, 0, 0, 0, 0, 0, (-t52 * t128 / 0.2e1 + (t52 * t180 + (-0.2e1 * t111 * t52 + t20 * t70) * t5) * t26) * t82, (t51 * t128 / 0.2e1 + (-t51 * t180 + (t111 * t203 - t19 * t70) * t5) * t26) * t82, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t6 * t2, (0.4e1 * t50 * t100 + (t50 * t144 + (-t18 * t6 - t2 * t50) * t143) * t184) * t85, (-0.4e1 * t53 * t100 + (-t53 * t144 + (t2 * t53 + t21 * t6) * t143) * t184) * t85, 0, 0, 0, 0, 0.2e1 * t5 * t1, 0, 0;];
MMD_reg = t8;
