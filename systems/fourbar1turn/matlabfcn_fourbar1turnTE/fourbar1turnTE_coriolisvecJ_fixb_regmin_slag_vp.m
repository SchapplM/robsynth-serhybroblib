% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [2x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1turnTE_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:32
% EndTime: 2020-06-27 16:23:01
% DurationCPUTime: 5.50s
% Computational Cost: add. (67281->252), mult. (94933->627), div. (2796->15), fcn. (25412->4), ass. (0->248)
t107 = pkin(1) ^ 2;
t97 = cos(qJ(2));
t255 = pkin(2) * t97;
t214 = -0.2e1 * pkin(1) * t255 + t107;
t263 = -pkin(3) - pkin(4);
t83 = (pkin(2) - t263) * (pkin(2) + t263) + t214;
t262 = pkin(4) - pkin(3);
t84 = (pkin(2) - t262) * (pkin(2) + t262) + t214;
t231 = t83 * t84;
t108 = sqrt(-t231);
t96 = sin(qJ(2));
t216 = t108 * t96;
t288 = pkin(4) ^ 2;
t209 = pkin(3) ^ 2 - t288;
t106 = pkin(2) ^ 2;
t90 = t106 + t214;
t86 = t90 - t209;
t91 = pkin(1) - t255;
t66 = -pkin(2) * t216 + t86 * t91;
t61 = 0.1e1 / t66 ^ 2;
t256 = pkin(2) * t96;
t81 = t86 * t256;
t67 = t108 * t91 + t81;
t240 = t61 * t67;
t102 = 0.1e1 / pkin(4);
t213 = qJD(2) * t96;
t182 = pkin(2) * t213;
t88 = 0.1e1 / t90 ^ 2;
t160 = t88 * t182;
t141 = pkin(1) * t160;
t159 = pkin(1) * pkin(2) * (-t83 - t84);
t73 = t96 * t159;
t72 = qJD(2) * t73;
t77 = 0.1e1 / t108;
t236 = t72 * t77;
t186 = t96 * t236;
t203 = 0.2e1 * t91 * pkin(1);
t215 = t97 * t108;
t36 = (-t186 + (-t215 + (t86 + t203) * t96) * qJD(2)) * pkin(2);
t87 = 0.1e1 / t90;
t247 = t36 * t87;
t119 = -t247 / 0.2e1 + t66 * t141;
t27 = t119 * t102;
t94 = t96 ^ 2;
t218 = t106 * t94;
t192 = pkin(1) * t218;
t157 = qJD(2) * t192;
t172 = t108 * t213;
t212 = qJD(2) * t97;
t223 = (t212 * t86 + t172) * pkin(2);
t234 = t77 * t91;
t37 = t72 * t234 + 0.2e1 * t157 + t223;
t246 = t37 * t87;
t118 = t246 / 0.2e1 - t67 * t141;
t29 = t118 * t102;
t60 = 0.1e1 / t66;
t131 = t27 * t240 + t29 * t60;
t287 = -0.4e1 * t66;
t286 = 0.6e1 * t97;
t62 = t60 * t61;
t63 = t67 ^ 2;
t239 = t62 * t63;
t21 = -t36 * t239 + t37 * t240;
t51 = t61 * t63 + 0.1e1;
t48 = 0.1e1 / t51 ^ 2;
t285 = t21 * t48;
t284 = t72 * t73;
t253 = pkin(4) * t90;
t47 = 0.1e1 / t51;
t196 = t47 * t253;
t165 = t61 * t196;
t283 = t102 * t67 * t165;
t85 = t90 + t209;
t92 = pkin(1) * t97 - pkin(2);
t65 = -pkin(1) * t216 - t85 * t92;
t58 = 0.1e1 / t65 ^ 2;
t258 = pkin(1) * t96;
t82 = t85 * t258;
t68 = -t108 * t92 + t82;
t242 = t58 * t68;
t105 = 0.1e1 / pkin(3);
t136 = -0.2e1 * t141;
t217 = t107 * t94;
t191 = pkin(2) * t217;
t156 = qJD(2) * t191;
t222 = (t212 * t85 + t172) * pkin(1);
t233 = t77 * t92;
t38 = -t72 * t233 + 0.2e1 * t156 + t222;
t30 = (t68 * t136 + t38 * t87) * t105;
t57 = 0.1e1 / t65;
t248 = t30 * t57;
t276 = -0.2e1 * t92;
t204 = pkin(2) * t276;
t171 = t85 + t204;
t35 = (-t186 + (t171 * t96 - t215) * qJD(2)) * pkin(1);
t28 = (t65 * t136 + t35 * t87) * t105;
t282 = -t28 * t242 + t248;
t277 = 0.2e1 * t90;
t151 = 0.2e1 * t48 * t277;
t237 = t67 * t90;
t168 = 0.4e1 * t62 * t237;
t201 = pkin(1) * t256;
t175 = -0.4e1 * t201;
t199 = t61 * t277;
t232 = 0.4e1 * t77 / t231;
t176 = -t232 / 0.4e1;
t138 = t176 * t284;
t266 = -t96 / 0.2e1;
t179 = t106 * t217;
t117 = t97 * t159 - 0.4e1 * t179;
t70 = t117 * qJD(2);
t116 = (t96 * t138 + 0.2e1 * (t70 * t266 - t97 * t72) * t77) * t87;
t89 = t87 * t88;
t158 = t89 * t179;
t133 = 0.4e1 * t67 * t158;
t137 = qJD(2) * t158;
t195 = t60 * t253;
t139 = t102 * t47 * t195;
t194 = t106 * t258;
t146 = t194 * t286;
t155 = t91 * t232 / 0.4e1;
t230 = t87 * t97;
t185 = t91 * t230;
t187 = t87 * t236;
t229 = t88 * t96;
t235 = t77 * t73;
t245 = t37 * t88;
t259 = pkin(1) * t88;
t267 = t87 / 0.2e1;
t268 = -t87 / 0.2e1;
t40 = t81 + (-t215 + (t203 - t235) * t96) * pkin(2);
t44 = t73 * t234 + 0.2e1 * t192 + (t97 * t86 + t216) * pkin(2);
t260 = -0.2e1 * ((0.4e1 * t157 + t223) * t268 + t137 * t287 + (-t116 / 0.2e1 + (t36 * t229 + (-t185 + (t40 * t96 + t66 * t97) * t88) * qJD(2)) * pkin(1)) * pkin(2)) * t283 - 0.2e1 * ((qJD(2) * t146 + t155 * t284 + t70 * t234) * t267 + qJD(2) * t133 + ((t187 / 0.2e1 - pkin(1) * t245) * t96 + ((t215 + (-t86 + t235) * t96) * t267 + (-t44 * t96 - t67 * t97) * t259) * qJD(2)) * pkin(2)) * t139;
t278 = -0.2e1 * t90;
t1 = (t131 * (-t40 * t239 + t44 * t240) * t151 + ((t60 * t175 + t40 * t199) * t29 + (t40 * t168 + (t67 * t175 + t44 * t278) * t61) * t27) * t47) * pkin(4) + t260;
t148 = 0.2e1 * t196;
t12 = t131 * t148;
t170 = t88 * t201;
t244 = t40 * t87;
t31 = (-t244 / 0.2e1 + t66 * t170) * t102;
t177 = t44 * t267;
t33 = (-t67 * t170 + t177) * t102;
t129 = t31 * t240 + t33 * t60;
t13 = t129 * t148;
t164 = pkin(1) * t182;
t147 = -0.4e1 * t164;
t3 = (t129 * t21 * t151 + ((t60 * t147 + t36 * t199) * t33 + (t36 * t168 + (t67 * t147 + t37 * t278) * t61) * t31) * t47) * pkin(4) + t260;
t281 = (-qJD(2) * t13 + t12) * t170 - (t3 / 0.2e1 - t1 / 0.2e1) * t87;
t280 = 0.2e1 * t68;
t279 = -0.2e1 * t87;
t275 = -t35 / 0.2e1;
t274 = -t36 / 0.2e1;
t273 = t37 / 0.2e1;
t272 = t38 / 0.2e1;
t43 = -t73 * t233 + 0.2e1 * t191 + (t85 * t97 + t216) * pkin(1);
t271 = -t43 / 0.2e1;
t270 = -t65 / 0.2e1;
t269 = t68 / 0.2e1;
t265 = t96 / 0.2e1;
t264 = -t97 / 0.2e1;
t152 = 0.2e1 * t201;
t59 = t57 * t58;
t200 = t59 * t280;
t64 = t68 ^ 2;
t241 = t59 * t64;
t52 = t58 * t64 + 0.1e1;
t50 = 0.1e1 / t52 ^ 2;
t243 = t50 * t90;
t132 = 0.8e1 * t137;
t193 = t107 * t256;
t145 = t193 * t286;
t198 = 0.2e1 * t88;
t257 = pkin(2) * t88;
t205 = -0.2e1 * t257;
t225 = t97 * t68;
t250 = ((qJD(2) * t145 + t92 * t138 - t70 * t233) * t87 + t68 * t132 + ((t38 * t205 + t187) * t96 + ((t215 + (-t85 + t235) * t96) * t87 + (-t43 * t96 - t225) * pkin(2) * t198) * qJD(2)) * pkin(1)) * t105 * t57;
t39 = t82 + (-t215 + (t204 - t235) * t96) * pkin(1);
t49 = 0.1e1 / t52;
t254 = pkin(3) * t90;
t197 = t49 * t254;
t166 = t58 * t197;
t144 = t68 * t166;
t226 = t97 * t65;
t251 = pkin(1) * t105;
t8 = (((0.4e1 * t156 + t222) * t87 + t65 * t132) * t105 + (t116 + (-0.2e1 * t35 * t229 + (t230 * t276 + (-t39 * t96 - t226) * t198) * qJD(2)) * pkin(2)) * t251) * t144;
t2 = -t8 + (-0.2e1 * t282 * (-t39 * t241 + t43 * t242) * t243 + (t282 * t152 + (-t30 * t39 * t58 + t250 + (t39 * t200 - t43 * t58) * t28) * t90) * t49) * pkin(3);
t153 = -0.2e1 * t170;
t32 = (t65 * t153 + t39 * t87) * t105;
t34 = (t68 * t153 + t43 * t87) * t105;
t128 = t32 * t242 - t34 * t57;
t22 = -t35 * t241 + t38 * t242;
t190 = t22 * t243;
t4 = -t8 + (0.2e1 * t128 * t190 + (-t128 * qJD(2) * t152 + (-t34 * t35 * t58 + t250 + (t35 * t200 - t38 * t58) * t32) * t90) * t49) * pkin(3);
t261 = -t2 + t4;
t71 = t72 ^ 2;
t154 = t71 * t176;
t99 = qJD(2) ^ 2;
t69 = t117 * t99;
t120 = t108 * t99 - t77 * t69 + t154;
t180 = qJD(2) * t236;
t162 = 0.2e1 * t180;
t125 = -t86 * t99 + t162;
t181 = t99 * t215;
t224 = t99 * t97;
t5 = -0.2e1 * ((t99 * t146 + t71 * t155 + t69 * t234) * t267 + t99 * t133 + ((t125 * t96 + t181) * t267 + (-0.2e1 * t37 * t213 - t67 * t224) * t259) * pkin(2)) * t139 + 0.4e1 * t29 * t195 * t285 - 0.2e1 * ((t107 * t89 * t287 + pkin(1) * t279) * t99 * t218 + ((t120 * t96 - t125 * t97) * t268 + (-t99 * t185 + (0.2e1 * t36 * t213 + t66 * t224) * t88) * pkin(1)) * pkin(2)) * t283 + 0.2e1 * (-t27 * t37 + t29 * t36) * t165 + 0.2e1 * (0.2e1 * (t62 * t47 * t36 + t61 * t285) * t27 * t237 - 0.2e1 * t131 * t47 * t164) * pkin(4);
t252 = t5 * t87;
t238 = t67 * t88;
t228 = t96 * t65;
t227 = t96 * t97;
t221 = -t97 ^ 2 + t94;
t220 = t105 * t87;
t219 = t105 * t99;
t100 = qJD(1) ^ 2;
t211 = t100 * t102;
t103 = 0.1e1 / t288;
t210 = t100 * t103;
t208 = qJD(1) * t102;
t207 = qJD(1) * t103;
t206 = qJD(1) * t105;
t202 = pkin(1) * t257;
t183 = t4 / 0.2e1 - t2 / 0.2e1;
t178 = 0.2e1 * qJD(1) * qJD(2);
t174 = -t39 / 0.2e1 + t269;
t173 = t270 + t271;
t169 = t89 * t201;
t167 = t57 * t197;
t163 = t88 * t193;
t140 = t89 * t164;
t134 = 0.2e1 * t107 * t160;
t126 = t225 / 0.2e1 + t228 / 0.2e1;
t14 = -t32 * t144 + t34 * t167 + 0.1e1;
t124 = t14 * t88 * t194 * t219;
t123 = (t68 * t265 - t226 / 0.2e1) * t87;
t46 = t105 * t123;
t122 = (t96 * t225 + t65 * t94) * t202;
t121 = (t96 * t226 - t68 * t94) * t202;
t17 = ((t38 * t264 + t35 * t266) * t87 + (t123 + t122) * qJD(2)) * t105;
t18 = ((t35 * t264 + t38 * t265) * t87 + (t126 * t87 + t121) * qJD(2)) * t105;
t45 = t126 * t220;
t42 = qJD(1) * t46;
t41 = (t225 + t228) * t206 * t268;
t20 = (t121 + (-t173 * t96 + t174 * t97) * t87) * t206;
t19 = (t122 + (t173 * t97 + t174 * t96) * t87) * t206;
t16 = qJD(1) * t18;
t15 = qJD(1) * t17;
t11 = -t28 * t144 + t30 * t167 + qJD(2);
t6 = ((t99 * t145 + t92 * t154 - t69 * t233) * t87 + 0.8e1 * t68 * t99 * t158 + ((t181 + (-t85 * t99 + t162) * t96) * t87 + (-0.4e1 * t38 * t213 - 0.2e1 * t68 * t224) * t257) * pkin(1)) * t105 * t167 - ((0.8e1 * t106 * t65 * t89 + 0.4e1 * pkin(2) * t87) * t217 * t219 + ((t180 * t279 + (t171 * t87 + t65 * t205) * t99) * t97 + (-0.4e1 * t35 * qJD(2) * t257 + t120 * t87) * t96) * t251) * t144 + (-t28 * t38 - t30 * t35) * t166 + (t58 * t50 * t22 + t59 * t49 * t35) * t28 * t254 * t280 + 0.2e1 * (t164 * t282 * t49 - t190 * t248) * pkin(3);
t7 = [0, 0, 0, t178 * t227, -t221 * t178, t224, -t99 * t96, 0, 0, 0, -t15 * t45 + t17 * t41, t15 * t46 - t16 * t45 + t17 * t42 + t18 * t41, t11 * t17 - t45 * t6, t16 * t46 + t18 * t42, t11 * t18 + t46 * t6, 0, (-t42 * t213 + t16 * t97 + (t18 * t97 - t46 * t213) * qJD(1)) * pkin(2), (t41 * t213 - t15 * t97 + (-t17 * t97 - t45 * t213) * qJD(1)) * pkin(2), (-t63 * t140 + t238 * t273) * t207, (t238 * t274 + (-t245 / 0.2e1 + 0.2e1 * t67 * t140) * t66) * t207, (t67 * t252 / 0.2e1 - t118 * t12) * t102, (-t66 * t252 / 0.2e1 - t119 * t12) * t102, 0, (-pkin(1) * t247 + t134 * t66) * t208, (-pkin(1) * t246 + t134 * t67) * t208; 0, 0, 0, -t100 * t227, t221 * t100, 0, 0, 0, 0, 0, -t41 * t19, -t19 * t42 - t20 * t41, -t11 * t19 + t14 * t15 + t261 * t41, -t42 * t20, -t11 * t20 + t14 * t16 + t261 * t42, t261 * t11 + t14 * t6, t65 * t124 + ((-t20 * t97 + t42 * t96) * qJD(1) + (t11 * t275 + t6 * t270 + (t14 * t275 + t39 * t11 / 0.2e1 - t183 * t65) * qJD(2)) * t220) * pkin(2), -t68 * t124 + ((t19 * t97 - t41 * t96) * qJD(1) + (t11 * t272 + t6 * t269 + (t11 * t271 + t14 * t272 + t183 * t68) * qJD(2)) * t220) * pkin(2), (-t44 * t238 / 0.4e1 + t63 * t169 / 0.2e1) * t210, (t40 * t238 / 0.4e1 + (t44 * t88 / 0.4e1 - t67 * t169) * t66) * t210, ((-t13 * t273 + t44 * t12 / 0.2e1) * t87 - t281 * t67) * t208, ((-t13 * t274 - t40 * t12 / 0.2e1) * t87 + t281 * t66) * t208, -t13 * t5 - (-t1 + t3) * t12, (pkin(1) * t244 / 0.2e1 - t66 * t163) * t211, (pkin(1) * t177 - t163 * t67) * t211;];
tauc_reg = t7;
