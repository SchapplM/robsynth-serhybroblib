% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m1OL_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:27:46
% EndTime: 2020-05-03 00:28:08
% DurationCPUTime: 4.22s
% Computational Cost: add. (5902->327), mult. (14883->495), div. (0->0), fcn. (12228->8), ass. (0->180)
t113 = cos(qJ(5));
t189 = qJD(5) * t113;
t111 = sin(qJ(4));
t232 = cos(qJ(4));
t114 = cos(qJ(2));
t112 = sin(qJ(2));
t231 = sin(qJ(3));
t171 = t231 * t112;
t233 = cos(qJ(3));
t135 = -t233 * t114 + t171;
t90 = t135 * qJD(1);
t170 = t231 * t114;
t134 = -t233 * t112 - t170;
t91 = t134 * qJD(1);
t60 = -t111 * t91 + t232 * t90;
t241 = t113 * t60;
t246 = t189 - t241;
t110 = sin(qJ(5));
t107 = qJD(2) + qJD(3);
t121 = t107 * t134;
t119 = t121 * qJD(1);
t168 = qJD(4) * t232;
t191 = qJD(4) * t111;
t128 = t135 * qJD(2);
t120 = t135 * qJD(3) + t128;
t67 = qJD(1) * t120;
t137 = -t111 * t119 + t90 * t168 - t91 * t191 + t232 * t67;
t139 = -t111 * t90 - t232 * t91;
t187 = -qJD(3) - qJD(4);
t160 = qJD(2) - t187;
t141 = t113 * t160;
t190 = qJD(5) * t110;
t18 = -qJD(5) * t141 - t113 * t137 - t139 * t190;
t47 = t110 * t160 - t113 * t139;
t202 = qJD(5) * t47;
t19 = t110 * t137 + t202;
t239 = qJD(5) - t60;
t245 = t110 * t239;
t45 = -t110 * t139 - t141;
t1 = -t110 * t19 - t18 * t113 - t245 * t47 - t246 * t45;
t118 = t232 * t121;
t162 = qJD(1) * t118 + t111 * t67;
t31 = -t139 * qJD(4) + t162;
t27 = t113 * t31;
t5 = -t139 * t45 - t239 * t245 + t27;
t103 = t233 * pkin(3) + pkin(2);
t82 = -pkin(3) * t171 + t103 * t114 + pkin(1);
t203 = qJD(1) * t82;
t36 = -pkin(4) * t60 + pkin(6) * t139 + t203;
t230 = t111 * pkin(3);
t184 = t233 * pkin(2);
t104 = t184 + pkin(3);
t150 = t232 * t231;
t87 = pkin(2) * t150 + t104 * t111;
t75 = t87 * qJD(2) + qJD(3) * t230;
t73 = t160 * pkin(6) + t75;
t142 = t110 * t73 - t113 * t36;
t244 = t239 * t142;
t24 = t110 * t36 + t113 * t73;
t126 = t134 * qJD(3);
t165 = t231 * qJD(2);
t192 = qJD(2) * t112;
t64 = -t103 * t192 + (-t114 * t165 + t126) * pkin(3);
t10 = pkin(4) * t31 - pkin(6) * t137 + qJD(1) * t64;
t147 = qJD(3) * t168;
t140 = pkin(3) * t147;
t172 = t231 * t111;
t131 = t232 * t233 - t172;
t65 = t104 * t168 + (t131 * qJD(3) - qJD(4) * t172) * pkin(2);
t50 = t65 * qJD(2) + t140;
t4 = -qJD(5) * t24 + t113 * t10 - t110 * t50;
t243 = -t239 * t24 - t4;
t15 = t18 * t110;
t7 = t246 * t47 - t15;
t217 = t110 * t31 + t189 * t239;
t6 = t139 * t47 - t239 * t241 + t217;
t183 = t232 * pkin(3);
t88 = -pkin(2) * t172 + t104 * t232;
t77 = t88 * qJD(2) + qJD(3) * t183;
t74 = -t160 * pkin(4) - t77;
t223 = t60 * t74;
t225 = t239 * t139;
t240 = t60 * t139;
t22 = t139 ^ 2 - t60 ^ 2;
t20 = -t60 * t160 + t137;
t169 = -t139 * t142 + t74 * t190;
t181 = pkin(3) * t191;
t100 = qJD(3) * t181;
t130 = t233 * t111 + t150;
t66 = t104 * t191 + (t130 * qJD(3) + qJD(4) * t150) * pkin(2);
t51 = -t66 * qJD(2) - t100;
t151 = -t51 * t110 - t24 * t139 + t74 * t189;
t21 = -t139 * t107 - t162;
t236 = -pkin(4) * t139 - pkin(6) * t60;
t188 = qJD(1) * qJD(2);
t235 = -0.2e1 * t188;
t93 = t130 * pkin(2);
t197 = t93 * qJD(2);
t138 = t181 - t197;
t234 = t110 * t142 + t113 * t24;
t3 = -t142 * qJD(5) + t10 * t110 + t113 * t50;
t2 = t3 * t113;
t125 = t232 * t135;
t71 = t111 * t134 - t125;
t229 = t31 * t71;
t129 = t111 * t135;
t72 = t134 * t232 + t129;
t228 = t31 * t72;
t227 = t47 * t45;
t226 = t51 * t72;
t220 = t82 * t139;
t219 = t82 * t60;
t218 = t90 * t91;
t213 = pkin(2) * qJD(2);
t211 = t110 * t24;
t209 = t110 * t45;
t116 = qJD(1) ^ 2;
t204 = t116 * t82;
t17 = t19 * t113;
t201 = qJD(5) * t239;
t200 = t112 * t116;
t199 = t114 * t116;
t92 = t131 * pkin(2);
t198 = t92 * qJD(2);
t195 = t112 ^ 2 - t114 ^ 2;
t105 = pkin(2) * t114 + pkin(1);
t194 = qJD(1) * t105;
t193 = qJD(1) * t112;
t186 = -t142 * t241 + t60 * t211 + t2;
t182 = pkin(2) * t192;
t180 = t72 * t190;
t179 = t72 * t189;
t178 = t90 * t194;
t177 = t91 * t194;
t176 = t112 * t199;
t175 = t233 * t107;
t173 = t231 * t107;
t167 = t233 * qJD(2);
t166 = t233 * qJD(3);
t164 = t231 * qJD(3);
t163 = t112 * t188;
t161 = -0.2e1 * qJD(2) + t187;
t86 = -pkin(3) * t170 - t103 * t112;
t37 = qJD(1) * t86 + t236;
t83 = pkin(6) + t87;
t159 = qJD(5) * t83 + t37;
t156 = pkin(1) * t235;
t155 = pkin(2) * t167;
t154 = pkin(2) * t164;
t153 = pkin(3) * t168;
t152 = pkin(2) * t163;
t149 = t114 * t163;
t34 = -qJD(4) * t125 + t111 * t121 - t232 * t120 + t134 * t191;
t146 = -t239 * t34 + t228;
t145 = -t139 * t75 + t77 * t60;
t143 = -t113 * t142 + t211;
t136 = -t239 * t65 - t31 * t83 - t223;
t35 = qJD(4) * t129 + t111 * t120 + t134 * t168 + t118;
t11 = t35 * pkin(4) + t34 * pkin(6) + t64;
t39 = t71 * pkin(4) - t72 * pkin(6) + t82;
t133 = qJD(5) * t72 * t74 + t11 * t239 + t31 * t39;
t132 = -t39 * t201 - t34 * t74 - t226;
t101 = pkin(6) + t230;
t124 = -t101 * t31 - t153 * t239 - t223;
t123 = -t143 * qJD(5) - t4 * t110 + t2;
t117 = pkin(2) * t121 * t231;
t115 = qJD(2) ^ 2;
t102 = -t183 - pkin(4);
t94 = t134 * pkin(3);
t84 = -t88 - pkin(4);
t49 = -t90 ^ 2 + t91 ^ 2;
t42 = t91 * t107 - t119;
t41 = -t90 * t107 + t67;
t38 = qJD(1) * t94 + t236;
t33 = t110 * t236 + t113 * t77;
t32 = -t110 * t77 + t113 * t236;
t29 = t110 * t38 + t113 * t198;
t28 = -t110 * t198 + t113 * t38;
t8 = t245 * t45 - t17;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t149, t195 * t235, -t115 * t114, -0.2e1 * t149, t115 * t112, 0, t112 * t156, t114 * t156, 0, 0, t91 * t120 + t134 * t67, -t134 * t119 + t67 * t135 + t107 * (-t91 * t134 + t135 * t90), t120 * t107, -t135 * t119 - t121 * t90, -t121 * t107, 0, t90 * t182 + 0.2e1 * (t134 * qJD(2) + t126) * t194 + pkin(2) * t128 * t193, 0.2e1 * t105 * t67 - t134 * t152 - t182 * t91, -t128 * t155 + (t134 * t154 - t117) * qJD(2), -0.2e1 * t105 * t152, t137 * t72 + t139 * t34, -t137 * t71 + t139 * t35 - t34 * t60 - t228, -t34 * t160, -t35 * t60 + t229, -t35 * t160, 0, t31 * t82 - t60 * t64 + (t35 * t82 + t64 * t71) * qJD(1), t137 * t82 - t139 * t64 + (-t34 * t82 + t64 * t72) * qJD(1), t34 * t77 - t35 * t75 - t50 * t71 - t226, 0.2e1 * t64 * t203, -t47 * t180 + (-t18 * t72 - t34 * t47) * t113, (t110 * t47 + t113 * t45) * t34 + (t15 - t17 + (-t113 * t47 + t209) * qJD(5)) * t72, t113 * t146 - t18 * t71 - t180 * t239 + t35 * t47, t45 * t179 + (t19 * t72 - t34 * t45) * t110, -t110 * t146 - t179 * t239 - t19 * t71 - t35 * t45, t239 * t35 + t229, t110 * t132 + t113 * t133 - t142 * t35 + t4 * t71, -t110 * t133 + t113 * t132 - t24 * t35 - t3 * t71, (-t11 * t47 + t18 * t39 - t142 * t34 - t4 * t72 + (-t24 * t72 - t39 * t45) * qJD(5)) * t113 + (-t11 * t45 - t19 * t39 + t24 * t34 - t3 * t72 + (-t142 * t72 + t39 * t47) * qJD(5)) * t110, t143 * t11 + (qJD(5) * t234 + t110 * t3 + t113 * t4) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t195 * t116, 0, t176, 0, 0, pkin(1) * t200, pkin(1) * t199, 0, 0, -t218, t49, t41, t218, t42, 0, -t177 + (-t90 * t193 + (-t165 - t173) * qJD(3)) * pkin(2), -t178 + (t91 * t193 + (-t167 - t175) * qJD(3)) * pkin(2), -qJD(1) * t117 + t154 * t91 + t155 * t90 - t184 * t67 + (t165 * t91 + t166 * t90) * pkin(2), t105 * pkin(2) * t200, t240, t22, t20, -t240, t21, 0, -t100 + t161 * t66 + (t60 * t86 + t220) * qJD(1), -t140 + t161 * t65 + (t139 * t86 - t219) * qJD(1), -t137 * t88 - t139 * t66 - t31 * t87 + t60 * t65 + t145, -t204 * t86 + t50 * t87 + t51 * t88 + t65 * t75 - t66 * t77, t7, t1, t6, t8, t5, t225, t19 * t84 + t45 * t66 + (-t159 * t239 + t51) * t113 + t136 * t110 + t169, t113 * t136 + t159 * t245 - t18 * t84 + t47 * t66 + t151, (-t19 * t83 + t37 * t47 - t45 * t65 + (t47 * t83 + t142) * qJD(5)) * t113 + (-t18 * t83 + t37 * t45 + t47 * t65 - t4 + (t45 * t83 - t24) * qJD(5)) * t110 + t186, -t51 * t84 + t66 * t74 + (t142 * t159 + t24 * t65 + t3 * t83) * t113 + (t142 * t65 - t159 * t24 - t4 * t83) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t49, t41, t218, t42, 0, -t177 + (-t164 + t173) * t213, -t178 + (-t166 + t175) * t213, 0, 0, t240, t22, t20, -t240, t21, 0, (t60 * t94 + t220) * qJD(1) + t51 - t138 * t160, (t160 * t92 - t65) * qJD(2) + (t139 * t94 - t219) * qJD(1) + (-t160 * t168 - t147) * pkin(3), (t139 * t93 - t60 * t92) * qJD(2) + (-t232 * t137 - t111 * t31 + (-t111 * t139 + t232 * t60) * qJD(4)) * pkin(3) + t145, -t94 * t204 + (-t75 * t92 + t77 * t93) * qJD(2) + (t232 * t51 + t111 * t50 + (-t111 * t77 + t232 * t75) * qJD(4)) * pkin(3), t7, t1, t6, t8, t5, t225, t102 * t19 - t28 * t239 + t138 * t45 + (-t101 * t201 + t51) * t113 + t124 * t110 + t169, -t102 * t18 + (t101 * t190 + t29) * t239 + t138 * t47 + t124 * t113 + t151, t28 * t47 + t29 * t45 + (-t45 * t153 - t101 * t19 + (t101 * t47 + t142) * qJD(5)) * t113 + (t47 * t153 - t101 * t18 - t4 + (t101 * t45 - t24) * qJD(5)) * t110 + t186, -t74 * t197 - t51 * t102 + t142 * t28 - t24 * t29 + (t111 * t74 + t232 * t234) * qJD(4) * pkin(3) + t123 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, t22, t20, -t240, t21, 0, t139 * t203 + t160 * t75 + t51, t160 * t77 - t203 * t60 - t50, 0, 0, t7, t1, t6, t209 * t239 - t17, t5, t225, -pkin(4) * t19 - pkin(6) * t217 - t110 * t223 + t51 * t113 - t239 * t32 - t45 * t75 + t169, -t74 * t241 + pkin(4) * t18 + t33 * t239 - t47 * t75 + (t190 * t239 - t27) * pkin(6) + t151, t32 * t47 + t33 * t45 + t2 + (t244 + (-t19 + t202) * pkin(6)) * t113 + ((qJD(5) * t45 - t18) * pkin(6) + t243) * t110, pkin(4) * t51 + pkin(6) * t123 + t142 * t32 - t24 * t33 - t74 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, -t45 ^ 2 + t47 ^ 2, t239 * t45 - t18, -t227, t239 * t47 - t19, t31, -t47 * t74 - t243, t45 * t74 - t244 - t3, 0, 0;];
tauc_reg = t9;
