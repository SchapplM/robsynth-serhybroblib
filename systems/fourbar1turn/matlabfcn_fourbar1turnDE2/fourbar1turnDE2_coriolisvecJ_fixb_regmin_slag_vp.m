% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [2x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1turnDE2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:30
% EndTime: 2020-06-27 16:49:06
% DurationCPUTime: 8.31s
% Computational Cost: add. (97843->245), mult. (139742->620), div. (4898->21), fcn. (37534->8), ass. (0->246)
t108 = cos(qJ(2));
t110 = qJD(2) ^ 2;
t242 = t110 * t108;
t325 = pkin(1) * pkin(2);
t119 = pkin(1) ^ 2;
t289 = pkin(2) * t108;
t237 = -0.2e1 * pkin(1) * t289 + t119;
t300 = -pkin(3) - pkin(4);
t94 = (pkin(2) - t300) * (pkin(2) + t300) + t237;
t299 = pkin(4) - pkin(3);
t95 = (pkin(2) - t299) * (pkin(2) + t299) + t237;
t263 = t94 * t95;
t107 = sin(qJ(2));
t183 = (-t94 - t95) * t325;
t84 = t107 * t183;
t83 = qJD(2) * t84;
t122 = sqrt(-t263);
t88 = 0.1e1 / t122;
t324 = 0.4e1 * t88 / t263 * t83 ^ 2;
t102 = pkin(1) - t289;
t244 = t107 * t122;
t118 = pkin(2) ^ 2;
t101 = t118 + t237;
t318 = pkin(3) ^ 2;
t319 = pkin(4) ^ 2;
t235 = t318 - t319;
t97 = t101 - t235;
t77 = -pkin(2) * t244 + t102 * t97;
t72 = 0.1e1 / t77 ^ 2;
t290 = pkin(2) * t107;
t92 = t97 * t290;
t78 = t102 * t122 + t92;
t268 = t72 * t78;
t113 = 0.1e1 / pkin(4);
t221 = pkin(1) * t290;
t184 = qJD(2) * t221;
t99 = 0.1e1 / t101 ^ 2;
t166 = t99 * t184;
t153 = t77 * t166;
t265 = t83 * t88;
t208 = t107 * t265;
t312 = 0.2e1 * pkin(1);
t226 = t102 * t312;
t243 = t108 * t122;
t45 = (-t208 + (-t243 + (t97 + t226) * t107) * qJD(2)) * pkin(2);
t98 = 0.1e1 / t101;
t282 = t45 * t98;
t35 = (-t282 / 0.2e1 + t153) * t113;
t152 = t78 * t166;
t231 = qJD(2) * t107;
t195 = t122 * t231;
t105 = t107 ^ 2;
t247 = t105 * t118;
t206 = pkin(1) * t247;
t230 = qJD(2) * t108;
t257 = t102 * t88;
t306 = 0.2e1 * qJD(2);
t46 = t206 * t306 + t83 * t257 + (t230 * t97 + t195) * pkin(2);
t281 = t46 * t98;
t37 = (t281 / 0.2e1 - t152) * t113;
t74 = t78 ^ 2;
t64 = t72 * t74 + 0.1e1;
t60 = 0.1e1 / t64;
t71 = 0.1e1 / t77;
t140 = t60 * (t35 * t268 + t37 * t71);
t323 = 0.6e1 * t242;
t103 = pkin(1) * t108 - pkin(2);
t96 = t101 + t235;
t76 = -pkin(1) * t244 - t103 * t96;
t66 = t76 ^ 2;
t67 = 0.1e1 / t76;
t69 = t67 / t66;
t292 = pkin(1) * t107;
t93 = t96 * t292;
t79 = -t103 * t122 + t93;
t75 = t79 ^ 2;
t269 = t69 * t75;
t68 = 0.1e1 / t76 ^ 2;
t270 = t68 * t79;
t228 = -0.2e1 * pkin(2) * t103;
t193 = t96 + t228;
t44 = (-t208 + (t193 * t107 - t243) * qJD(2)) * pkin(1);
t246 = t105 * t119;
t205 = pkin(2) * t246;
t256 = t103 * t88;
t47 = t205 * t306 - t83 * t256 + (t230 * t96 + t195) * pkin(1);
t65 = t68 * t75 + 0.1e1;
t63 = 0.1e1 / t65 ^ 2;
t285 = (-t44 * t269 + t47 * t270) * t63;
t322 = -0.2e1 * t285;
t100 = t98 * t99;
t199 = t118 * t246;
t321 = t110 * t100 * t199;
t317 = -0.2e1 * t46;
t287 = pkin(4) * t101;
t227 = 0.2e1 * t287;
t116 = 0.1e1 / pkin(3);
t157 = -0.2e1 * t166;
t36 = (t76 * t157 + t44 * t98) * t116;
t38 = (t79 * t157 + t47 * t98) * t116;
t314 = -t36 * t270 + t38 * t67;
t176 = 0.2e1 * t221;
t18 = t140 * t227;
t192 = t99 * t221;
t302 = -t98 / 0.2e1;
t264 = t88 * t84;
t49 = t92 + (-t243 + (t226 - t264) * t107) * pkin(2);
t39 = (t77 * t192 + t49 * t302) * t113;
t301 = t98 / 0.2e1;
t51 = t84 * t257 + 0.2e1 * t206 + (t108 * t97 + t244) * pkin(2);
t41 = (-t78 * t192 + t51 * t301) * t113;
t139 = t60 * (-t39 * t268 - t41 * t71);
t19 = t139 * t227;
t175 = 0.4e1 * t221;
t61 = 0.1e1 / t64 ^ 2;
t194 = 0.4e1 * t61 * t268;
t70 = t77 ^ 2;
t73 = t71 / t70;
t267 = t73 * t74;
t25 = -t45 * t267 + t46 * t268;
t271 = t61 * t71;
t216 = t25 * t271;
t309 = 0.4e1 * t78;
t223 = t73 * t309;
t225 = 0.2e1 * t60 * t72;
t27 = -t49 * t267 + t51 * t268;
t310 = -0.2e1 * t72;
t298 = ((-qJD(2) * t139 - t140) * t175 + ((t49 * t225 + 0.4e1 * t27 * t271) * t37 + (t27 * t194 + (t49 * t223 + t51 * t310) * t60) * t35 - (t45 * t225 + 0.4e1 * t216) * t41 - (t25 * t194 + (t45 * t223 + t46 * t310) * t60) * t39) * t101) * pkin(4);
t313 = t176 * (qJD(2) * t19 + t18) * t99 + t298 * t98;
t114 = 0.1e1 / t319;
t260 = t70 + t74;
t58 = t260 * t99 * t114;
t56 = 0.1e1 / t58;
t52 = t58 ^ (-0.1e1 / 0.2e1);
t311 = -0.2e1 * t52;
t308 = -0.2e1 * t79;
t307 = 0.2e1 * t79;
t117 = 0.1e1 / t318;
t261 = t66 + t75;
t59 = t261 * t99 * t117;
t54 = t59 ^ (-0.1e1 / 0.2e1);
t158 = t100 * t175;
t137 = t260 * t158;
t222 = 0.2e1 * t99;
t21 = ((t45 * t77 + t46 * t78) * t222 - qJD(2) * t137) * t114;
t305 = t21 / 0.2e1;
t138 = t261 * t158;
t22 = ((t44 * t76 + t47 * t79) * t222 - qJD(2) * t138) * t117;
t304 = t22 / 0.2e1;
t23 = ((t49 * t77 + t51 * t78) * t222 - t137) * t114;
t303 = -t23 / 0.2e1;
t177 = -0.2e1 * t192;
t48 = t93 + (-t243 + (t228 - t264) * t107) * pkin(1);
t279 = t48 * t98;
t40 = (t76 * t177 + t279) * t116;
t50 = -t84 * t256 + 0.2e1 * t205 + (t108 * t96 + t244) * pkin(1);
t278 = t50 * t98;
t42 = (t79 * t177 + t278) * t116;
t149 = t40 * t270 - t42 * t67;
t62 = 0.1e1 / t65;
t297 = ((t149 * qJD(2) + t314) * t62 * t176 + (-0.2e1 * t314 * t63 * (-t48 * t269 + t50 * t270) + t149 * t322 + ((t48 * t36 - t44 * t40) * t69 * t307 + (-t50 * t36 - t38 * t48 + t47 * t40 + t42 * t44) * t68) * t62) * t101) * pkin(3);
t296 = pkin(1) * t98;
t294 = pkin(2) * t99;
t178 = -t324 / 0.4e1;
t240 = t110 * t122;
t80 = (t108 * t183 - 0.4e1 * t199) * t110;
t133 = -t88 * t80 + t178 + t240;
t288 = pkin(3) * t101;
t189 = t62 * t68 * t288;
t170 = t79 * t189;
t203 = qJD(2) * t265;
t186 = 0.2e1 * t203;
t218 = t67 * t288;
t190 = t62 * t218;
t198 = t108 * t240;
t209 = t119 * t290;
t241 = t110 * t116;
t6 = ((t103 * t178 + t209 * t323 - t80 * t256) * t98 + 0.8e1 * t79 * t321 + ((t198 + (-t110 * t96 + t186) * t107) * t98 + (-0.4e1 * t47 * t231 + t242 * t308) * t294) * pkin(1)) * t116 * t190 + t38 * t218 * t322 - ((0.8e1 * t100 * t118 * t76 + 0.4e1 * pkin(2) * t98) * t241 * t246 + ((-0.2e1 * t98 * t203 + (t193 * t98 - 0.2e1 * t294 * t76) * t110) * t108 + (-0.4e1 * t44 * qJD(2) * t294 + t133 * t98) * t107) * pkin(1) * t116) * t170 + (-t36 * t47 - t38 * t44) * t189 + (t69 * t62 * t44 + t68 * t285) * t36 * t288 * t307 + 0.2e1 * t314 * pkin(3) * t62 * t184;
t293 = t54 * t6;
t53 = t52 * t56;
t277 = t53 * t77;
t276 = t53 * t78;
t275 = t53 * t98;
t55 = t54 / t59;
t274 = t55 * t76;
t273 = t55 * t79;
t272 = t55 * t98;
t266 = t78 * t99;
t255 = t105 * t76;
t254 = t107 * t79;
t252 = t108 * t76;
t250 = t116 * t98;
t249 = t76 * t107;
t248 = t79 * t108;
t245 = t107 * t108;
t111 = qJD(1) ^ 2;
t239 = t111 * t113;
t238 = t111 * t114;
t236 = -t108 ^ 2 + t105;
t234 = qJD(1) * t113;
t233 = qJD(1) * t114;
t232 = qJD(1) * t116;
t144 = -t110 * t97 + t186;
t217 = t60 * t287;
t188 = t72 * t217;
t210 = t118 * t292;
t5 = 0.2e1 * (-t35 * t46 + t37 * t45) * t188 - 0.4e1 * pkin(4) * t184 * t140 + 0.2e1 * (t37 * t216 + (t72 * t61 * t25 + t73 * t60 * t45) * t35 * t78) * t227 + 0.2e1 * (-((t102 * t324 / 0.4e1 + t210 * t323 + t80 * t257) * t301 + t309 * t321 + ((t144 * t107 + t198) * t301 + (t231 * t317 - t78 * t242) * t99 * pkin(1)) * pkin(2)) * t71 * t217 - ((-0.4e1 * t100 * t119 * t77 - 0.2e1 * t296) * t110 * t247 + ((t133 * t107 - t144 * t108) * t302 + ((-t102 * t98 + t77 * t99) * t242 + t45 * t107 * t99 * t306) * pkin(1)) * pkin(2)) * t78 * t188) * t113;
t220 = t5 * t52 * t98;
t213 = t77 * t275;
t57 = 0.1e1 / t58 ^ 2;
t212 = t57 * t74 * t99;
t207 = t54 * t250;
t202 = qJD(1) * t306;
t201 = -t276 / 0.2e1;
t200 = t294 * t312;
t191 = t100 * t221;
t187 = t57 * t77 * t266;
t182 = t222 * t325;
t181 = qJD(1) * t207;
t165 = t52 * t99 * t209;
t163 = t100 * t184;
t159 = t248 + t249;
t156 = -0.2e1 * t165;
t146 = -t18 * t303 + t19 * t305;
t17 = -t36 * t170 + t38 * t190 + qJD(2);
t20 = -t40 * t170 + t42 * t190 + 0.1e1;
t24 = ((t48 * t76 + t50 * t79) * t222 - t138) * t117;
t145 = t20 * t304 - t17 * t24 / 0.2e1;
t143 = 0.4e1 * qJD(2) * t165;
t141 = (-t252 + t254) * t98;
t136 = t20 * t54 * t99 * t210 * t241;
t33 = t159 * t207;
t135 = (-t254 / 0.2e1 + t252 / 0.2e1) * t272;
t134 = (t249 / 0.2e1 + t248 / 0.2e1) * t272;
t131 = (-t105 * t79 + t76 * t245) * t182;
t10 = (t22 * t135 + ((t107 * t47 - t108 * t44) * t98 + (t159 * t98 + t131) * qJD(2)) * t54) * t116;
t9 = (t22 * t134 + ((-t107 * t44 - t108 * t47) * t98 + (t141 + (t79 * t245 + t255) * t182) * qJD(2)) * t54) * t116;
t43 = t181 * t252;
t34 = t54 * t116 * t141;
t32 = t181 * t254 - t43;
t31 = qJD(1) * t33;
t12 = (t24 * t135 + (t131 + ((-t48 + t79) * t108 + (t50 + t76) * t107) * t98) * t54) * t232;
t11 = -t43 + (t24 * t134 + (t200 * t255 - t108 * t278 + (-t279 + (t108 * t200 + t98) * t79) * t107) * t54) * t232;
t8 = qJD(1) * t10;
t7 = qJD(1) * t9;
t1 = [0, 0, 0, t202 * t245, -t236 * t202, t242, -t110 * t107, 0, 0, 0, -t31 * t9 - t33 * t7, -t10 * t31 + t32 * t9 - t33 * t8 + t34 * t7, t17 * t9 - t33 * t6, t10 * t32 + t34 * t8, t10 * t17 + t34 * t6, 0, (-t32 * t231 + t108 * t8 + (t10 * t108 - t34 * t231) * qJD(1)) * pkin(2), (-t31 * t231 - t108 * t7 + (-t108 * t9 - t231 * t33) * qJD(1)) * pkin(2), (-t21 * t212 + (-0.4e1 * t163 * t74 + 0.2e1 * t266 * t46) * t56) * t233, (0.2e1 * t21 * t187 + (-0.2e1 * t45 * t266 + (0.8e1 * t163 * t78 + t99 * t317) * t77) * t56) * t233, (t78 * t220 - (t98 * t21 * t201 + (-0.2e1 * t152 + t281) * t52) * t18) * t113, (-t77 * t220 - (t213 * t305 + (0.2e1 * t153 - t282) * t52) * t18) * t113, 0, (t77 * t143 + (t21 * t277 + t45 * t311) * t296) * t234, (t78 * t143 + (t21 * t276 + t46 * t311) * t296) * t234; 0, 0, 0, -t111 * t245, t236 * t111, 0, 0, 0, 0, 0, t31 * t11, -t11 * t32 + t12 * t31, -t11 * t17 + t20 * t7 + t297 * t31, -t32 * t12, -t12 * t17 + t20 * t8 - t297 * t32, -t297 * t17 + t20 * t6, 0.2e1 * t76 * t136 + ((t107 * t32 - t108 * t12) * qJD(1) + (-t76 * t293 + (t274 * t304 - t44 * t54) * t17 + (t145 * t274 + (t48 * t17 - t44 * t20 + t297 * t76) * t54) * qJD(2)) * t250) * pkin(2), t136 * t308 + ((t107 * t31 + t108 * t11) * qJD(1) + (t79 * t293 + (t47 * t54 - t22 * t273 / 0.2e1) * t17 + (-t145 * t273 + (-t50 * t17 + t47 * t20 - t297 * t79) * t54) * qJD(2)) * t250) * pkin(2), (t23 * t212 / 0.2e1 + (0.2e1 * t191 * t74 - t266 * t51) * t56) * t238, (-t23 * t187 + (t49 * t266 + (-0.4e1 * t191 * t78 + t51 * t99) * t77) * t56) * t238, (-t146 * t78 * t275 + ((t18 * t51 + t19 * t46) * t98 - t313 * t78) * t52) * t234, (t146 * t213 + ((-t18 * t49 - t19 * t45) * t98 + t313 * t77) * t52) * t234, t298 * t18 + t19 * t5, (t77 * t156 + (t277 * t303 + t49 * t52) * t296) * t239, (t78 * t156 + (t201 * t23 + t51 * t52) * t296) * t239;];
tauc_reg = t1;
