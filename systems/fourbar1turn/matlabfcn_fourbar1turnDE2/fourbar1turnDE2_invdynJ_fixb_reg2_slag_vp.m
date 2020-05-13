% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [2x(2*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:34:59
% EndTime: 2020-04-12 19:35:23
% DurationCPUTime: 12.04s
% Computational Cost: add. (119796->334), mult. (171634->793), div. (6421->22), fcn. (46265->13), ass. (0->289)
t381 = pkin(4) ^ 2;
t137 = 0.1e1 / t381;
t142 = pkin(2) ^ 2;
t143 = pkin(1) ^ 2;
t130 = cos(qJ(2));
t341 = pkin(2) * t130;
t372 = -2 * pkin(1);
t286 = t341 * t372 + t143;
t122 = t142 + t286;
t118 = 0.1e1 / t122;
t119 = 0.1e1 / t122 ^ 2;
t120 = t118 * t119;
t128 = sin(qJ(2));
t342 = pkin(2) * t128;
t266 = pkin(1) * t342;
t218 = 0.4e1 * t266;
t380 = pkin(3) ^ 2;
t284 = t380 - t381;
t117 = t122 - t284;
t123 = pkin(1) - t341;
t360 = -pkin(3) - pkin(4);
t114 = (pkin(2) - t360) * (pkin(2) + t360) + t286;
t359 = pkin(4) - pkin(3);
t115 = (pkin(2) - t359) * (pkin(2) + t359) + t286;
t305 = t114 * t115;
t146 = sqrt(-t305);
t296 = t128 * t146;
t98 = -pkin(2) * t296 + t117 * t123;
t91 = t98 ^ 2;
t112 = t117 * t342;
t99 = t123 * t146 + t112;
t95 = t99 ^ 2;
t318 = t91 + t95;
t167 = t318 * t120 * t218;
t263 = 0.2e1 * t119;
t269 = 0.2e1 * t123 * pkin(1);
t295 = t130 * t146;
t386 = pkin(1) * pkin(2);
t220 = (-t114 - t115) * t386;
t104 = t128 * t220;
t108 = 0.1e1 / t146;
t309 = t108 * t104;
t64 = t112 + (-t295 + (t269 - t309) * t128) * pkin(2);
t126 = t128 ^ 2;
t299 = t126 * t142;
t253 = pkin(1) * t299;
t308 = t108 * t123;
t66 = t104 * t308 + 0.2e1 * t253 + (t117 * t130 + t296) * pkin(2);
t29 = ((t64 * t98 + t66 * t99) * t263 - t167) * t137;
t388 = t29 / 0.2e1;
t387 = 0.8e1 * t120;
t93 = 0.1e1 / t98 ^ 2;
t323 = t93 * t99;
t136 = 0.1e1 / pkin(4);
t281 = qJD(2) * t128;
t252 = pkin(2) * t281;
t223 = t119 * t252;
t206 = pkin(1) * t223;
t353 = -t118 / 0.2e1;
t103 = qJD(2) * t104;
t310 = t103 * t108;
t244 = t128 * t310;
t60 = (-t244 + (-t295 + (t117 + t269) * t128) * qJD(2)) * pkin(2);
t41 = (t98 * t206 + t60 * t353) * t136;
t352 = t118 / 0.2e1;
t240 = t146 * t281;
t280 = qJD(2) * t130;
t384 = 0.2e1 * qJD(2);
t61 = t103 * t308 + t253 * t384 + (t117 * t280 + t240) * pkin(2);
t43 = (-t99 * t206 + t61 * t352) * t136;
t84 = t93 * t95 + 0.1e1;
t80 = 0.1e1 / t84;
t92 = 0.1e1 / t98;
t173 = t80 * (t41 * t323 + t43 * t92);
t133 = qJD(2) ^ 2;
t293 = t133 * t142;
t298 = t126 * t143;
t385 = 0.2e1 * t298 * pkin(2);
t347 = pkin(1) * t128;
t383 = t347 * t293;
t116 = t122 + t284;
t113 = t116 * t347;
t124 = pkin(1) * t130 - pkin(2);
t100 = -t124 * t146 + t113;
t97 = -pkin(1) * t296 - t116 * t124;
t89 = 0.1e1 / t97 ^ 2;
t314 = t100 * t89;
t139 = 0.1e1 / pkin(3);
t191 = -0.2e1 * t206;
t271 = -0.2e1 * pkin(2) * t124;
t235 = t116 + t271;
t59 = (-t244 + (t235 * t128 - t295) * qJD(2)) * pkin(1);
t42 = (t118 * t59 + t191 * t97) * t139;
t307 = t108 * t124;
t62 = -t103 * t307 + qJD(2) * t385 + (t116 * t280 + t240) * pkin(1);
t44 = (t100 * t191 + t118 * t62) * t139;
t88 = 0.1e1 / t97;
t189 = t42 * t314 - t44 * t88;
t339 = pkin(4) * t122;
t270 = 0.2e1 * t339;
t376 = 2 * qJDD(2);
t375 = pkin(1) * t139;
t18 = t173 * t270;
t255 = t119 * t342;
t234 = pkin(1) * t255;
t45 = (t98 * t234 + t64 * t353) * t136;
t47 = (-t99 * t234 + t66 * t352) * t136;
t172 = t80 * (-t45 * t323 - t47 * t92);
t19 = t172 * t270;
t219 = 0.2e1 * t266;
t78 = t318 * t137 * t119;
t70 = t78 ^ (-0.1e1 / 0.2e1);
t374 = (qJD(2) * t19 + t18) * t119 * t70 * t219;
t294 = t133 * t130;
t373 = qJDD(2) * t128 + t294;
t74 = 0.1e1 / t78;
t140 = 0.1e1 / t380;
t87 = t97 ^ 2;
t96 = t100 ^ 2;
t319 = t87 + t96;
t79 = t319 * t140 * t119;
t76 = 0.1e1 / t79;
t371 = 0.2e1 * t60;
t370 = -0.2e1 * t93;
t369 = 0.4e1 * t99;
t72 = t79 ^ (-0.1e1 / 0.2e1);
t368 = -0.2e1 * t100;
t367 = 0.2e1 * t100;
t365 = -0.4e1 * qJD(2);
t340 = pkin(3) * t122;
t85 = t89 * t96 + 0.1e1;
t82 = 0.1e1 / t85;
t231 = t82 * t89 * t340;
t208 = t100 * t231;
t262 = t88 * t340;
t232 = t82 * t262;
t17 = -t42 * t208 + t44 * t232 + qJD(2);
t363 = -t17 / 0.2e1;
t27 = ((t60 * t98 + t61 * t99) * t263 - qJD(2) * t167) * t137;
t362 = t27 / 0.2e1;
t233 = t120 * t266;
t187 = t319 * t233;
t168 = -0.4e1 * t187;
t171 = 0.2e1 * t100 * t62 + 0.2e1 * t59 * t97;
t28 = (qJD(2) * t168 + t119 * t171) * t140;
t361 = t28 / 0.2e1;
t81 = 0.1e1 / t84 ^ 2;
t238 = 0.4e1 * t81 * t323;
t94 = t92 / t91;
t322 = t94 * t95;
t31 = -t60 * t322 + t61 * t323;
t325 = t81 * t92;
t260 = t31 * t325;
t267 = t94 * t369;
t268 = 0.2e1 * t80 * t93;
t33 = -t64 * t322 + t66 * t323;
t358 = ((-qJD(2) * t172 - t173) * t218 + ((t64 * t268 + 0.4e1 * t33 * t325) * t43 + (t33 * t238 + (t64 * t267 + t66 * t370) * t80) * t41 - (t60 * t268 + 0.4e1 * t260) * t47 - (t31 * t238 + (t60 * t267 + t61 * t370) * t80) * t45) * t122) * pkin(4);
t215 = -0.2e1 * t234;
t63 = t113 + (-t295 + (t271 - t309) * t128) * pkin(1);
t46 = (t118 * t63 + t97 * t215) * t139;
t65 = -t104 * t307 + t385 + (t130 * t116 + t296) * pkin(1);
t48 = (t100 * t215 + t118 * t65) * t139;
t188 = t46 * t314 - t48 * t88;
t90 = t88 / t87;
t324 = t90 * t96;
t83 = 0.1e1 / t85 ^ 2;
t335 = (t62 * t314 - t59 * t324) * t83;
t357 = ((-t188 * qJD(2) + t189) * t82 * t219 + (0.2e1 * t188 * t335 - 0.2e1 * t189 * t83 * (t65 * t314 - t63 * t324) + ((-t63 * t42 + t59 * t46) * t90 * t367 + (t65 * t42 + t44 * t63 - t62 * t46 - t48 * t59) * t89) * t82) * t122) * pkin(3);
t356 = g(3) * t98;
t355 = g(3) * t99;
t175 = qJDD(2) * t146 + t310 * t384;
t160 = -t117 * t133 + t175;
t176 = (t103 ^ 2 / t305 + t373 * t220 - 0.4e1 * t298 * t293) * t108;
t290 = t146 * t133;
t161 = -t176 + t290;
t195 = -t118 * t123 + t119 * t98;
t227 = pkin(1) * t252;
t261 = t80 * t339;
t230 = t93 * t261;
t245 = t126 * t293;
t275 = qJDD(2) * t117;
t282 = qJD(2) * t119;
t292 = t133 * t143;
t349 = pkin(1) * t118;
t5 = 0.2e1 * (-t41 * t61 + t43 * t60) * t230 - 0.4e1 * pkin(4) * t227 * t173 + 0.2e1 * (t43 * t260 + (t93 * t81 * t31 + t94 * t80 * t60) * t41 * t99) * t270 + 0.2e1 * (-((t123 * t176 + 0.6e1 * t130 * t383) * t352 + (t120 * t292 * t369 + qJDD(2) * t349) * t299 + (((t275 + t290) * t130 + t160 * t128) * t352 + (-t99 * t294 + (-0.2e1 * qJD(2) * t61 - qJDD(2) * t99) * t128) * pkin(1) * t119) * pkin(2)) * t92 * t261 - ((-0.4e1 * t120 * t143 * t98 - 0.2e1 * t349) * t245 + ((t195 * t133 * pkin(1) + t160 * t352) * t130 + ((t161 + t275) * t353 + (qJDD(2) * t195 + t282 * t371) * pkin(1)) * t128) * pkin(2)) * t99 * t230) * t136;
t354 = t5 * t70;
t131 = cos(qJ(1));
t351 = g(1) * t131;
t129 = sin(qJ(1));
t350 = g(2) * t129;
t134 = (qJD(1) ^ 2);
t346 = pkin(1) * t134;
t344 = pkin(2) * t118;
t343 = pkin(2) * t119;
t75 = 0.1e1 / t78 ^ 2;
t337 = t27 * t75;
t71 = t70 * t74;
t336 = t29 * t71;
t196 = t100 * t130 + t128 * t97;
t166 = t118 * t72 * t196;
t39 = t139 * t166;
t334 = t39 * t97;
t331 = t71 * t98;
t330 = t71 * t99;
t73 = t72 * t76;
t329 = t73 * t28;
t328 = t73 * t97;
t327 = t74 * t98;
t326 = t74 * t99;
t321 = t98 * t99;
t317 = t97 + t65;
t316 = pkin(1) * qJD(1);
t302 = t118 * t139;
t247 = t100 * t302;
t198 = t72 * t128 * t247;
t250 = t97 * t302;
t203 = t72 * t130 * t250;
t40 = t198 - t203;
t315 = t100 * t40;
t313 = t97 * t130;
t37 = qJD(1) * t39;
t312 = t100 * t128;
t304 = t118 * t124;
t297 = t128 * t130;
t291 = t134 * t137;
t127 = t130 ^ 2;
t285 = t126 - t127;
t283 = qJD(1) * t139;
t279 = qJD(1) * qJD(2);
t278 = qJDD(1) * t118;
t277 = qJDD(1) * t127;
t276 = qJDD(2) * t116;
t273 = t130 * qJDD(1);
t272 = -0.2e1 * t343;
t265 = -2 * t316;
t258 = t75 * t321;
t257 = -0.2e1 * t134 * t143;
t249 = t75 * t388;
t248 = -t328 / 0.2e1;
t243 = t134 * t297;
t242 = -t100 * t73 / 0.2e1;
t239 = qJDD(1) * t119 * t74;
t236 = t128 * t279;
t20 = -t46 * t208 + t48 * t232 + 0.1e1;
t6 = qJDD(2) + ((-t176 * t304 + (t100 * t245 * t387 + (t126 * t376 + 0.6e1 * t128 * t294) * t344) * t143) * t139 + (((t276 + t290) * t118 + t100 * t133 * t272) * t130 + ((-t133 * t116 + t175) * t118 + (qJDD(2) * t368 + t62 * t365) * t343) * t128) * t375) * t232 - 0.2e1 * t44 * t262 * t335 - ((t142 * t97 * t387 + 0.4e1 * t344) * t139 * t126 * t292 + ((-t175 * t118 + (t235 * t118 + t97 * t272) * t133) * t130 + ((t161 + t276) * t118 + (-0.4e1 * t59 * t282 + (-t119 * t97 - t304) * t376) * pkin(2)) * t128) * t375) * t208 + (-t42 * t62 - t44 * t59) * t231 + (t90 * t82 * t59 + t89 * t335) * t42 * t340 * t367 - 0.2e1 * t189 * pkin(3) * t82 * t227;
t228 = qJDD(2) * t20 + t6;
t226 = t27 * t71 * t316;
t225 = t70 * t255;
t217 = t19 * t70 * t278;
t212 = t130 * t236;
t211 = t74 * t233;
t121 = t350 + t351;
t205 = g(1) * t129 - g(2) * t131;
t204 = t70 * t223;
t202 = t315 + t334;
t201 = t100 * t65 + t63 * t97;
t200 = t72 * t263 * t386;
t197 = 0.2e1 * t211;
t190 = t119 * t72 * t383;
t184 = t18 * t388 + t19 * t362;
t30 = (t201 * t263 + t168) * t140;
t183 = t20 * t361 + t30 * t363;
t56 = qJD(1) * t198;
t180 = t211 * t321;
t179 = t211 * t365;
t178 = pkin(1) * t18 * t204;
t177 = t139 * t190;
t174 = 0.2e1 * t121;
t170 = t121 + t346;
t169 = qJDD(1) * t372 - t205;
t165 = 0.4e1 * qJD(1) * t143 * t204;
t164 = t73 * (-t312 / 0.2e1 + t313 / 0.2e1);
t163 = -t346 / 0.2e1 - t351 / 0.2e1 - t350 / 0.2e1;
t162 = -g(3) * t130 + t121 * t128;
t159 = (-t100 * t126 + t97 * t297) * t200;
t158 = (-t100 * t297 - t126 * t97) * t200;
t156 = qJD(2) * t158 + ((t28 * t248 + t59 * t72) * t128 + (t28 * t242 + (qJD(2) * t97 + t62) * t72) * t130) * t118;
t155 = qJD(2) * t159 + (t28 * t164 + (t196 * qJD(2) + t62 * t128 - t59 * t130) * t72) * t118;
t69 = qJ(2) + atan2(t247, t250);
t68 = cos(t69);
t67 = sin(t69);
t38 = -qJD(1) * t203 + t56;
t12 = (t159 + (t30 * t164 + ((t100 - t63) * t130 + t317 * t128) * t72) * t118) * t283;
t11 = -t56 + (t158 + ((t30 * t248 + t63 * t72) * t128 + (t30 * t242 + t317 * t72) * t130) * t118) * t283;
t10 = t155 * t139;
t9 = -qJD(2) * t198 + t139 * t156;
t8 = ((t312 - t313) * t72 * t278 + t155 * qJD(1)) * t139;
t7 = -qJD(2) * t56 + (qJD(1) * t156 + qJDD(1) * t166) * t139;
t1 = [0, 0, 0, 0, 0, qJDD(1), t205, t121, 0, 0, qJDD(1) * t126 + 0.2e1 * t212, 0.2e1 * t128 * t273 - 0.2e1 * t285 * t279, t373, -0.2e1 * t212 + t277, qJDD(2) * t130 - t133 * t128, 0, t205 * t130, -t205 * t128, -t121, 0, t37 * t9 + t39 * t7, -t10 * t37 - t38 * t9 - t39 * t8 - t40 * t7, -t17 * t9 - t39 * t6, t10 * t38 + t40 * t8, t10 * t17 + t40 * t6, 0, -t205 * t68 + ((-qJD(1) * t40 - t38) * t281 + (qJD(1) * t10 + qJDD(1) * t40 + t8) * t130) * pkin(2), t205 * t67 + (-0.2e1 * t37 * t281 + (qJD(1) * t9 + qJDD(1) * t39 + t7) * t130) * pkin(2), (0.2e1 * t202 * t190 + (-t202 * t72 * qJDD(2) + ((t315 / 0.2e1 + t334 / 0.2e1) * t329 + (-t100 * t10 - t59 * t39 - t62 * t40 - t97 * t9) * t72) * qJD(2)) * t344) * t139 - t121, t142 * t277 + (pkin(2) * t205 - 0.2e1 * t142 * t236) * t130, (t95 * t239 + (t95 * t179 + (0.2e1 * t61 * t326 - t95 * t337) * t119) * qJD(1)) * t137, (-0.2e1 * t239 * t321 + (0.8e1 * qJD(2) * t180 + (t27 * t258 + (-t60 * t99 - t61 * t98) * t74) * t263) * qJD(1)) * t137, (0.2e1 * t99 * t178 + (t99 * t354 - (t61 * t70 - t27 * t330 / 0.2e1) * t18) * t118) * t136, (t91 * t239 + (t91 * t179 + (t327 * t371 - t91 * t337) * t119) * qJD(1)) * t137, (-0.2e1 * t98 * t178 + (-t98 * t354 - (t331 * t362 - t60 * t70) * t18) * t118) * t136, 0, (t98 * t165 + (t98 * t226 + (t169 * t98 + t265 * t60) * t70) * t118) * t136, (t99 * t165 + (t99 * t226 + (t169 * t99 + t265 * t61) * t70) * t118) * t136, -t121, pkin(1) * t205 + qJDD(1) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t285 * t134, t128 * qJDD(1), t243, t273, qJDD(2), t162, g(3) * t128 + t121 * t130, 0, 0, -t37 * t11, t11 * t38 + t12 * t37, t11 * t17 - t20 * t7 - t357 * t37, -t38 * t12, -t12 * t17 + t20 * t8 + t357 * t38, t357 * t17 + t20 * t6, (g(3) * t68 - t121 * t67 + 0.2e1 * t177 * t97) * t20 + ((-t130 * t12 + t128 * t38) * qJD(1) + ((t183 * qJD(2) + t17 * t361) * t328 + (-t59 * t17 - t228 * t97 + (t17 * t63 - t20 * t59 - t357 * t97) * qJD(2)) * t72) * t302) * pkin(2), (-g(3) * t67 - t121 * t68 + t177 * t368) * t20 + ((-t130 * t11 + t128 * t37) * qJD(1) + ((t62 * t17 + (-t17 * t65 + t20 * t62) * qJD(2)) * t72 + (t329 * t363 + t228 * t72 + (-t183 * t73 + t357 * t72) * qJD(2)) * t100) * t302) * pkin(2), ((-t100 * t8 - t59 * t37 - t62 * t38 - t97 * t7 + (t100 * t12 + t11 * t97 + t37 * t63 + t38 * t65) * qJD(2)) * t72 + (-qJD(2) * t30 + t28) * t73 * (t100 * t38 / 0.2e1 + t97 * t37 / 0.2e1)) * pkin(2) * t302, t142 * t243 + t162 * pkin(2) + (-0.2e1 * t76 * t187 * t293 + (((t96 / 0.2e1 + t87 / 0.2e1) * t30 * t133 - t319 * t28 * qJD(2)) / t79 ^ 2 + (qJD(2) * t171 + qJDD(2) * t319 - t133 * t201) * t76) * t142 * t119) * t140, (t95 * t197 + (t249 * t95 - t66 * t326) * t119) * t291, (-0.4e1 * t180 + (-t29 * t258 + (t64 * t99 + t66 * t98) * t74) * t119) * t291, (t99 * t217 + (-t99 * t374 + (-t184 * t330 + (t18 * t66 + t19 * t61 - t358 * t99) * t70) * t118) * qJD(1)) * t136, (t91 * t197 + (t249 * t91 - t64 * t327) * t119) * t291, (-t98 * t217 + (t98 * t374 + (t184 * t331 + (-t18 * t64 - t19 * t60 + t358 * t98) * t70) * t118) * qJD(1)) * t136, t358 * t18 + t19 * t5, ((t98 * t257 + (-t174 * t98 + 0.2e1 * t355) * pkin(1)) * t225 + ((-g(3) * t66 + t170 * t64) * t70 + (t355 / 0.2e1 + t163 * t98) * t336) * t118) * t136, ((t99 * t257 + (-t174 * t99 - 0.2e1 * t356) * pkin(1)) * t225 + ((g(3) * t64 + t170 * t66) * t70 + (-t356 / 0.2e1 + t163 * t99) * t336) * t118) * t136, 0, 0;];
tau_reg = t1;
