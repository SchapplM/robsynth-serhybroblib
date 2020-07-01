% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% tau_reg [2x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:30
% EndTime: 2020-06-27 16:49:12
% DurationCPUTime: 9.97s
% Computational Cost: add. (106554->292), mult. (152085->702), div. (5460->21), fcn. (40979->13), ass. (0->270)
t342 = pkin(4) ^ 2;
t120 = 0.1e1 / t342;
t124 = pkin(2) ^ 2;
t125 = pkin(1) ^ 2;
t113 = cos(qJ(2));
t305 = pkin(2) * t113;
t333 = -2 * pkin(1);
t254 = t305 * t333 + t125;
t105 = t124 + t254;
t102 = 0.1e1 / t105;
t103 = 0.1e1 / t105 ^ 2;
t104 = t102 * t103;
t111 = sin(qJ(2));
t306 = pkin(2) * t111;
t235 = pkin(1) * t306;
t192 = 0.4e1 * t235;
t171 = t104 * t192;
t341 = pkin(3) ^ 2;
t252 = t341 - t342;
t101 = t105 - t252;
t106 = pkin(1) - t305;
t324 = -pkin(3) - pkin(4);
t98 = (pkin(2) - t324) * (pkin(2) + t324) + t254;
t323 = pkin(4) - pkin(3);
t99 = (pkin(2) - t323) * (pkin(2) + t323) + t254;
t284 = t98 * t99;
t128 = sqrt(-t284);
t260 = t111 * t128;
t82 = -pkin(2) * t260 + t101 * t106;
t75 = t82 ^ 2;
t96 = t101 * t306;
t83 = t106 * t128 + t96;
t79 = t83 ^ 2;
t280 = t75 + t79;
t146 = t280 * t171;
t233 = 0.2e1 * t103;
t239 = 0.2e1 * t106 * pkin(1);
t259 = t113 * t128;
t352 = pkin(1) * pkin(2);
t197 = (-t98 - t99) * t352;
t88 = t111 * t197;
t92 = 0.1e1 / t128;
t285 = t92 * t88;
t50 = t96 + (-t259 + (t239 - t285) * t111) * pkin(2);
t109 = t111 ^ 2;
t263 = t109 * t124;
t222 = pkin(1) * t263;
t275 = t106 * t92;
t52 = t88 * t275 + 0.2e1 * t222 + (t113 * t101 + t260) * pkin(2);
t23 = ((t50 * t82 + t52 * t83) * t233 - t146) * t120;
t354 = t23 / 0.2e1;
t353 = 0.8e1 * t104;
t77 = 0.1e1 / t82 ^ 2;
t289 = t77 * t83;
t119 = 0.1e1 / pkin(4);
t249 = qJD(2) * t111;
t220 = pkin(2) * t249;
t196 = t103 * t220;
t182 = pkin(1) * t196;
t317 = -t102 / 0.2e1;
t87 = qJD(2) * t88;
t286 = t87 * t92;
t223 = t111 * t286;
t46 = (-t223 + (-t259 + (t101 + t239) * t111) * qJD(2)) * pkin(2);
t35 = (t82 * t182 + t46 * t317) * t119;
t316 = t102 / 0.2e1;
t213 = t128 * t249;
t248 = qJD(2) * t113;
t348 = 0.2e1 * qJD(2);
t47 = t222 * t348 + t87 * t275 + (t101 * t248 + t213) * pkin(2);
t37 = (-t83 * t182 + t47 * t316) * t119;
t68 = t77 * t79 + 0.1e1;
t64 = 0.1e1 / t68;
t76 = 0.1e1 / t82;
t152 = t64 * (t35 * t289 + t37 * t76);
t116 = qJD(2) ^ 2;
t351 = t116 * t124;
t262 = t109 * t125;
t350 = 0.2e1 * pkin(2) * t262;
t100 = t105 + t252;
t107 = pkin(1) * t113 - pkin(2);
t81 = -pkin(1) * t260 - t100 * t107;
t71 = t81 ^ 2;
t72 = 0.1e1 / t81;
t74 = t72 / t71;
t311 = pkin(1) * t111;
t97 = t100 * t311;
t84 = -t107 * t128 + t97;
t80 = t84 ^ 2;
t290 = t74 * t80;
t73 = 0.1e1 / t81 ^ 2;
t291 = t73 * t84;
t241 = -0.2e1 * pkin(2) * t107;
t209 = t100 + t241;
t45 = (-t223 + (t209 * t111 - t259) * qJD(2)) * pkin(1);
t274 = t107 * t92;
t48 = qJD(2) * t350 - t87 * t274 + (t100 * t248 + t213) * pkin(1);
t69 = t73 * t80 + 0.1e1;
t67 = 0.1e1 / t69 ^ 2;
t299 = (-t45 * t290 + t48 * t291) * t67;
t349 = -0.2e1 * t299;
t123 = 0.1e1 / t341;
t281 = t71 + t80;
t147 = t281 * t171;
t22 = ((t45 * t81 + t48 * t84) * t233 - qJD(2) * t147) * t123;
t347 = t22 / 0.2e1;
t304 = pkin(3) * t105;
t66 = 0.1e1 / t69;
t206 = t66 * t73 * t304;
t186 = t84 * t206;
t232 = t72 * t304;
t207 = t66 * t232;
t122 = 0.1e1 / pkin(3);
t166 = -0.2e1 * t182;
t36 = (t102 * t45 + t81 * t166) * t122;
t38 = (t102 * t48 + t84 * t166) * t122;
t17 = -t36 * t186 + t38 * t207 + qJD(2);
t224 = t103 * t306;
t208 = pkin(1) * t224;
t190 = -0.2e1 * t208;
t49 = t97 + (-t259 + (t241 - t285) * t111) * pkin(1);
t40 = (t102 * t49 + t81 * t190) * t122;
t51 = -t88 * t274 + t350 + (t100 * t113 + t260) * pkin(1);
t42 = (t102 * t51 + t84 * t190) * t122;
t20 = -t40 * t186 + t42 * t207 + 0.1e1;
t24 = ((t49 * t81 + t51 * t84) * t233 - t147) * t123;
t63 = t281 * t123 * t103;
t58 = t63 ^ (-0.1e1 / 0.2e1);
t59 = t58 / t63;
t346 = (t17 * t347 + (t20 * t347 - t17 * t24 / 0.2e1) * qJD(2)) * t59;
t345 = t311 * t351;
t344 = -t36 * t291 + t38 * t72;
t303 = pkin(4) * t105;
t240 = 0.2e1 * t303;
t338 = 2 * qJDD(2);
t337 = pkin(1) * t122;
t18 = t152 * t240;
t39 = (t82 * t208 + t50 * t317) * t119;
t41 = (-t83 * t208 + t52 * t316) * t119;
t151 = t64 * (-t39 * t289 - t41 * t76);
t19 = t151 * t240;
t193 = 0.2e1 * t235;
t62 = t280 * t120 * t103;
t56 = t62 ^ (-0.1e1 / 0.2e1);
t336 = (qJD(2) * t19 + t18) * t103 * t56 * t193;
t258 = t116 * t113;
t335 = qJDD(2) * t111 + t258;
t60 = 0.1e1 / t62;
t332 = -0.2e1 * t77;
t331 = 0.4e1 * t83;
t330 = -0.2e1 * t84;
t329 = 0.2e1 * t84;
t327 = -0.4e1 * qJD(2);
t21 = ((t46 * t82 + t47 * t83) * t233 - qJD(2) * t146) * t120;
t325 = t21 / 0.2e1;
t65 = 0.1e1 / t68 ^ 2;
t211 = 0.4e1 * t65 * t289;
t78 = t76 / t75;
t288 = t78 * t79;
t25 = -t46 * t288 + t47 * t289;
t292 = t65 * t76;
t230 = t25 * t292;
t236 = t78 * t331;
t238 = 0.2e1 * t64 * t77;
t27 = -t50 * t288 + t52 * t289;
t322 = ((-qJD(2) * t151 - t152) * t192 + ((t50 * t238 + 0.4e1 * t27 * t292) * t37 + (t27 * t211 + (t50 * t236 + t52 * t332) * t64) * t35 - (t46 * t238 + 0.4e1 * t230) * t41 - (t25 * t211 + (t46 * t236 + t47 * t332) * t64) * t39) * t105) * pkin(4);
t164 = t40 * t291 - t42 * t72;
t321 = ((t164 * qJD(2) + t344) * t66 * t193 + (-0.2e1 * t344 * t67 * (-t49 * t290 + t51 * t291) + t164 * t349 + ((t49 * t36 - t45 * t40) * t74 * t329 + (-t51 * t36 - t38 * t49 + t48 * t40 + t42 * t45) * t73) * t66) * t105) * pkin(3);
t320 = g(3) * t82;
t319 = g(3) * t83;
t155 = qJDD(2) * t128 + t286 * t348;
t141 = -t116 * t101 + t155;
t157 = (t87 ^ 2 / t284 + t197 * t335 - 0.4e1 * t262 * t351) * t92;
t256 = t116 * t128;
t142 = -t157 + t256;
t170 = -t102 * t106 + t103 * t82;
t202 = pkin(1) * t220;
t231 = t64 * t303;
t205 = t77 * t231;
t216 = t116 * t263;
t245 = qJDD(2) * t101;
t250 = qJD(2) * t103;
t257 = t116 * t125;
t313 = pkin(1) * t102;
t5 = 0.2e1 * (-t35 * t47 + t37 * t46) * t205 - 0.4e1 * pkin(4) * t202 * t152 + 0.2e1 * (t37 * t230 + (t77 * t65 * t25 + t78 * t64 * t46) * t35 * t83) * t240 + 0.2e1 * (-((t157 * t106 + 0.6e1 * t113 * t345) * t316 + (t104 * t257 * t331 + qJDD(2) * t313) * t263 + (((t245 + t256) * t113 + t141 * t111) * t316 + (-t83 * t258 + (-0.2e1 * qJD(2) * t47 - qJDD(2) * t83) * t111) * pkin(1) * t103) * pkin(2)) * t76 * t231 - ((-0.4e1 * t104 * t125 * t82 - 0.2e1 * t313) * t216 + ((t170 * t116 * pkin(1) + t141 * t316) * t113 + ((t142 + t245) * t317 + (t170 * qJDD(2) + 0.2e1 * t46 * t250) * pkin(1)) * t111) * pkin(2)) * t83 * t205) * t119;
t318 = t5 * t56;
t114 = cos(qJ(1));
t315 = g(1) * t114;
t112 = sin(qJ(1));
t314 = g(2) * t112;
t117 = (qJD(1) ^ 2);
t310 = pkin(1) * t117;
t308 = pkin(2) * t102;
t307 = pkin(2) * t103;
t57 = t56 * t60;
t300 = t23 * t57;
t296 = t57 * t82;
t295 = t57 * t83;
t294 = t60 * t83;
t61 = 0.1e1 / t62 ^ 2;
t293 = t61 * t79;
t287 = t82 * t83;
t282 = -t49 + t84;
t276 = pkin(1) * qJD(1);
t273 = t111 * t84;
t272 = t113 * t81;
t270 = t81 * t111;
t269 = t84 * t113;
t172 = t269 + t270;
t148 = t102 * t58 * t172;
t33 = t122 * t148;
t31 = qJD(1) * t33;
t268 = t102 * t107;
t266 = t102 * t122;
t261 = t111 * t113;
t255 = t117 * t120;
t253 = -t113 ^ 2 + t109;
t251 = qJD(1) * t122;
t247 = qJDD(1) * t102;
t246 = qJDD(2) * t100;
t243 = t113 * qJDD(1);
t242 = -0.2e1 * t307;
t234 = -2 * t276;
t227 = t61 * t287;
t226 = -0.2e1 * t117 * t125;
t219 = t81 * t266;
t218 = t84 * t266;
t217 = t24 * t59 / 0.2e1;
t212 = qJDD(1) * t103 * t60;
t6 = qJDD(2) + ((-t157 * t268 + (t84 * t216 * t353 + (t109 * t338 + 0.6e1 * t111 * t258) * t308) * t125) * t122 + (((t246 + t256) * t102 + t84 * t116 * t242) * t113 + ((-t100 * t116 + t155) * t102 + (qJDD(2) * t330 + t48 * t327) * t307) * t111) * t337) * t207 + t38 * t232 * t349 - ((t124 * t81 * t353 + 0.4e1 * t308) * t122 * t109 * t257 + ((-t155 * t102 + (t209 * t102 + t81 * t242) * t116) * t113 + ((t142 + t246) * t102 + (-0.4e1 * t45 * t250 + (-t103 * t81 - t268) * t338) * pkin(2)) * t111) * t337) * t186 + (-t36 * t48 - t38 * t45) * t206 + (t74 * t66 * t45 + t73 * t299) * t36 * t304 * t329 + 0.2e1 * t344 * pkin(3) * t66 * t202;
t203 = qJDD(2) * t20 + t6;
t201 = t21 * t57 * t276;
t200 = t56 * t224;
t191 = t19 * t56 * t247;
t187 = t104 * t60 * t235;
t181 = t314 + t315;
t180 = g(1) * t112 - g(2) * t114;
t179 = t56 * t196;
t178 = t58 * t113 * t219;
t177 = t58 * t111 * t218;
t175 = t58 * t233 * t352;
t173 = -t272 + t273;
t167 = t79 * t187;
t161 = t18 * t354 + t19 * t325;
t158 = t187 * t287;
t156 = pkin(1) * t18 * t179;
t154 = t103 * t122 * t58 * t345;
t153 = 0.2e1 * t181;
t150 = t181 + t310;
t149 = qJDD(1) * t333 - t180;
t145 = t59 * (-t273 / 0.2e1 + t272 / 0.2e1);
t144 = 0.4e1 * qJD(1) * t125 * t179;
t143 = -t310 / 0.2e1 - t315 / 0.2e1 - t314 / 0.2e1;
t140 = (t109 * t81 + t84 * t261) * t175;
t139 = (-t109 * t84 + t81 * t261) * t175;
t137 = qJD(2) * t140 + ((t270 / 0.2e1 + t269 / 0.2e1) * t59 * t22 + (t173 * qJD(2) - t111 * t45 - t48 * t113) * t58) * t102;
t136 = qJD(2) * t139 + (t22 * t145 + (t172 * qJD(2) + t48 * t111 - t45 * t113) * t58) * t102;
t55 = qJ(2) + atan2(t218, t219);
t54 = cos(t55);
t53 = sin(t55);
t43 = qJD(1) * t178;
t34 = t177 - t178;
t32 = qJD(1) * t177 - t43;
t12 = (t139 + (t24 * t145 + (t282 * t113 + (t51 + t81) * t111) * t58) * t102) * t251;
t11 = -t43 + (t140 + ((t84 * t217 - t51 * t58) * t113 + (t81 * t217 + t282 * t58) * t111) * t102) * t251;
t10 = t136 * t122;
t9 = t137 * t122;
t8 = (t173 * t58 * t247 + t136 * qJD(1)) * t122;
t7 = (t137 * qJD(1) - qJDD(1) * t148) * t122;
t1 = [qJDD(1), t180, t181, 0.2e1 * qJD(1) * t111 * t248 + qJDD(1) * t109, -0.2e1 * t253 * qJD(2) * qJD(1) + 0.2e1 * t111 * t243, t335, qJDD(2) * t113 - t111 * t116, 0, t180 * t113, -t180 * t111, -t31 * t9 - t33 * t7, -t10 * t31 + t32 * t9 - t33 * t8 + t34 * t7, t17 * t9 - t33 * t6, t10 * t32 + t34 * t8, t10 * t17 + t34 * t6, 0, -t180 * t54 + ((-qJD(1) * t34 - t32) * t249 + (qJD(1) * t10 + qJDD(1) * t34 + t8) * t113) * pkin(2), t180 * t53 + (-0.2e1 * t31 * t249 + (-qJD(1) * t9 + qJDD(1) * t33 - t7) * t113) * pkin(2), (t79 * t212 + (t167 * t327 + (-t21 * t293 + 0.2e1 * t294 * t47) * t103) * qJD(1)) * t120, (-0.2e1 * t212 * t287 + (0.8e1 * qJD(2) * t158 + (t21 * t227 + (-t46 * t83 - t47 * t82) * t60) * t233) * qJD(1)) * t120, (0.2e1 * t83 * t156 + (t83 * t318 - (t47 * t56 - t21 * t295 / 0.2e1) * t18) * t102) * t119, (-0.2e1 * t82 * t156 + (-t82 * t318 - (t296 * t325 - t46 * t56) * t18) * t102) * t119, 0, (t82 * t144 + (t82 * t201 + (t149 * t82 + t234 * t46) * t56) * t102) * t119, (t83 * t144 + (t83 * t201 + (t149 * t83 + t234 * t47) * t56) * t102) * t119; 0, 0, 0, -t117 * t261, t253 * t117, t111 * qJDD(1), t243, qJDD(2), -g(3) * t113 + t111 * t181, g(3) * t111 + t113 * t181, t31 * t11, -t11 * t32 + t12 * t31, -t11 * t17 + t20 * t7 + t31 * t321, -t32 * t12, -t12 * t17 + t20 * t8 - t32 * t321, -t17 * t321 + t20 * t6, (g(3) * t54 + 0.2e1 * t154 * t81 - t181 * t53) * t20 + ((t111 * t32 - t113 * t12) * qJD(1) + (t81 * t346 + (-t45 * t17 - t203 * t81 + (t17 * t49 - t20 * t45 + t321 * t81) * qJD(2)) * t58) * t266) * pkin(2), (-g(3) * t53 + t154 * t330 - t181 * t54) * t20 + ((t113 * t11 + t111 * t31) * qJD(1) + (-t84 * t346 + (t48 * t17 + t203 * t84 + (-t17 * t51 + t20 * t48 - t321 * t84) * qJD(2)) * t58) * t266) * pkin(2), (0.2e1 * t167 + (t293 * t354 - t52 * t294) * t103) * t255, (-0.4e1 * t158 + (-t23 * t227 + (t50 * t83 + t52 * t82) * t60) * t103) * t255, (t83 * t191 + (-t83 * t336 + (-t161 * t295 + (t52 * t18 + t47 * t19 - t322 * t83) * t56) * t102) * qJD(1)) * t119, (-t82 * t191 + (t82 * t336 + (t161 * t296 + (-t50 * t18 - t46 * t19 + t322 * t82) * t56) * t102) * qJD(1)) * t119, t18 * t322 + t19 * t5, ((t82 * t226 + (-t153 * t82 + 0.2e1 * t319) * pkin(1)) * t200 + ((-g(3) * t52 + t150 * t50) * t56 + (t319 / 0.2e1 + t143 * t82) * t300) * t102) * t119, ((t83 * t226 + (-t153 * t83 - 0.2e1 * t320) * pkin(1)) * t200 + ((g(3) * t50 + t150 * t52) * t56 + (-t320 / 0.2e1 + t143 * t83) * t300) * t102) * t119;];
tau_reg = t1;
