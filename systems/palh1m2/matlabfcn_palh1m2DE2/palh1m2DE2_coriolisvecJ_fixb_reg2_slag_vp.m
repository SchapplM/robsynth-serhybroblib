% Calculate inertial parameters regressor of coriolis joint torque vector for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh1m2DE2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:06:26
% EndTime: 2020-05-02 21:07:09
% DurationCPUTime: 7.54s
% Computational Cost: add. (2107->355), mult. (5509->696), div. (0->0), fcn. (4745->22), ass. (0->334)
t203 = cos(pkin(18));
t204 = cos(pkin(17));
t352 = t203 * t204;
t198 = sin(pkin(18));
t199 = sin(pkin(17));
t359 = t198 * t199;
t117 = t352 + t359;
t197 = sin(qJ(2));
t202 = cos(qJ(2));
t435 = -t198 * t204 + t203 * t199;
t245 = -t117 * t197 - t202 * t435;
t57 = -t117 * t202 + t435 * t197;
t458 = t245 * t57;
t328 = qJD(1) * qJD(2);
t457 = -0.2e1 * t328 * t458;
t206 = qJD(2) ^ 2;
t201 = cos(qJ(3));
t411 = pkin(5) * t201;
t325 = pkin(1) * t411;
t456 = t206 * t325;
t193 = cos(pkin(20));
t392 = sin(pkin(20));
t111 = t203 * t193 + t198 * t392;
t196 = sin(qJ(3));
t153 = t203 * t392;
t360 = t198 * t193;
t107 = -t153 + t360;
t383 = t107 * t201;
t443 = -t196 * t111 + t383;
t381 = t111 * t201;
t384 = t107 * t196;
t63 = t381 + t384;
t248 = -t197 * t63 + t202 * t443;
t407 = pkin(5) * qJD(3);
t317 = t248 * t407;
t412 = pkin(5) * t196;
t160 = pkin(1) + t412;
t45 = -pkin(5) * t383 + t111 * t160;
t90 = pkin(5) * t381;
t46 = t107 * t160 + t90;
t23 = t197 * t46 + t202 * t45;
t389 = qJD(2) * t23;
t9 = -t317 + t389;
t194 = cos(pkin(19));
t367 = t196 * t194;
t191 = sin(pkin(19));
t370 = t191 * t201;
t106 = -t367 - t370;
t371 = t191 * t196;
t244 = -t194 * t201 + t371;
t228 = t244 * t197;
t216 = t106 * t202 + t228;
t207 = qJD(1) ^ 2;
t355 = t202 * t207;
t455 = t216 * t355;
t208 = pkin(5) ^ 2;
t364 = t196 * t208;
t143 = 0.2e1 * t201 * t364;
t185 = t202 ^ 2;
t448 = -0.2e1 * pkin(5);
t349 = pkin(1) * t448;
t98 = (t349 - 0.4e1 * t364) * t201 * t185;
t454 = t143 + t98;
t22 = -t197 * t45 + t46 * t202;
t453 = qJD(2) * t22;
t189 = sin(pkin(22));
t192 = cos(pkin(22));
t452 = -t189 * t203 + t192 * t198;
t176 = pkin(22) + pkin(21);
t166 = cos(t176);
t165 = sin(t176);
t385 = t107 * t165;
t39 = t111 * t166 + t385;
t446 = 2 * qJD(1);
t356 = t202 * t196;
t363 = t197 * t201;
t114 = -t356 - t363;
t226 = t114 * qJD(3);
t338 = qJD(2) * t201;
t283 = t197 * t338;
t337 = qJD(2) * t202;
t284 = t160 * t337;
t52 = -t284 + (t226 - t283) * pkin(5);
t315 = pkin(5) * t363;
t92 = t160 * t202 + t315;
t59 = qJD(2) * t92 - t114 * t407;
t451 = qJD(4) * t59 + (t39 * t446 + qJD(4)) * t52;
t178 = t194 ^ 2;
t161 = t178 - 0.1e1 / 0.2e1;
t169 = t178 / 0.2e1;
t181 = t196 ^ 2;
t184 = t201 ^ 2;
t132 = t161 * t184;
t365 = t196 * t201;
t368 = t194 * t191;
t263 = t365 * t368;
t224 = -0.2e1 * t263 + t132;
t344 = t181 - t184;
t362 = t197 * t202;
t296 = t161 * t365;
t430 = t344 * t368 - 0.2e1 * t296;
t450 = (t161 * t181 - t224) * t185 - t344 * (t169 - 0.1e1 / 0.4e1) - t430 * t362 - t263;
t210 = pkin(1) ^ 2;
t373 = t184 * t208;
t112 = t196 * t349 - t208 - t210 + 0.2e1 * t373;
t357 = t201 * t202;
t149 = pkin(5) * t357;
t137 = t149 + pkin(15);
t449 = (-t112 * t197 - t137 * t160) * t202;
t334 = qJD(3) * t202;
t221 = t196 * t334 + t283;
t447 = 0.2e1 * t221;
t445 = pkin(2) * pkin(15);
t24 = t197 * t443 + t63 * t202;
t48 = t443 * qJD(3);
t14 = t24 * qJD(2) + t48 * t197 + t63 * t334;
t327 = qJD(1) * qJD(4);
t128 = -pkin(9) * t198 + pkin(11) * t203;
t129 = pkin(9) * t203 + pkin(11) * t198;
t78 = -t392 * t128 + t129 * t193;
t287 = t78 * t327;
t227 = t63 * qJD(3);
t43 = pkin(5) * t227;
t44 = pkin(5) * t48;
t4 = t197 * t44 + t202 * t43 + t453;
t186 = t203 ^ 2;
t164 = t186 - 0.1e1 / 0.2e1;
t354 = t203 * t198;
t314 = pkin(9) * t354;
t104 = t164 * pkin(11) - t314;
t105 = t164 * pkin(9) + pkin(11) * t354;
t177 = t193 ^ 2;
t274 = t193 * t392;
t434 = 0.1e1 / 0.4e1 - t186 / 0.2e1;
t404 = t165 * (t104 * t177 + t105 * t274 + t314 / 0.2e1 + t434 * pkin(11));
t378 = t160 * t197;
t269 = t149 - t378;
t89 = pkin(15) + t269;
t405 = t111 * t89;
t444 = -qJD(2) * t4 - t14 * t407 + t287 - (0.4e1 * t404 + t405) * t207;
t179 = qJD(2) + qJD(3);
t366 = t196 * t197;
t243 = -t357 + t366;
t50 = t114 * t191 - t194 * t243;
t20 = t179 * t50;
t240 = t24 * t407;
t205 = qJD(4) ^ 2;
t440 = (t165 * t24 - t166 * t248) * t205;
t320 = -0.2e1 * t362;
t120 = t320 * t365;
t439 = t120 - t344 * (t185 - 0.1e1 / 0.2e1);
t163 = t184 - 0.1e1 / 0.2e1;
t182 = t197 ^ 2;
t343 = t182 - t185;
t438 = -t343 * t163 + t120;
t436 = -0.1e1 / 0.4e1 + t263 - t132;
t56 = t245 * qJD(2);
t209 = pkin(2) ^ 2;
t291 = t194 * t356;
t295 = t184 * t368;
t414 = pkin(2) * t197;
t85 = pkin(2) * t161 * t366 - pkin(15) * t194 / 0.2e1;
t433 = -t291 * t445 - 0.4e1 * t209 * (-t178 / 0.2e1 - t436) * t362 + 0.2e1 * (t295 * t414 + t85 * t201 - t191 * (-t196 * pkin(15) + t194 * t414) / 0.2e1) * t414;
t431 = t163 * t368 + t296;
t173 = t197 * pkin(15);
t154 = t173 - pkin(1);
t416 = pkin(1) * t185;
t123 = -t154 - t416;
t174 = t197 * pkin(1);
t155 = t174 - pkin(15);
t254 = -t155 * t197 + t416;
t429 = qJD(2) * t254 - qJD(3) * t123;
t292 = t194 * t363;
t428 = -t292 * t445 + pkin(2) * (-pkin(2) * t367 + t173) * t371 + (0.2e1 * t430 * t185 - (-0.2e1 * t178 + 0.1e1) * t365) * t209;
t13 = qJD(2) * t248 - t197 * t227 + t48 * t202;
t3 = -t197 * t43 + t202 * t44 - t389;
t382 = t107 * t207;
t144 = -t354 / 0.2e1;
t159 = t166 ^ 2;
t235 = -t104 * t274 + t105 * t177;
t393 = (pkin(9) * t434 + pkin(11) * t144 + t235) * t159;
t423 = 0.4e1 * t207;
t427 = t393 * t423 + (-qJD(2) * t3 - t13 * t407 - t89 * t382) * t165;
t152 = t185 * t194;
t293 = t194 * t362;
t294 = t191 * t362;
t426 = -(-t182 * t194 + t152 - 0.2e1 * t294) * t201 - (t343 * t191 - 0.2e1 * t293) * t196;
t425 = -0.2e1 * t202;
t424 = -0.2e1 * t207;
t422 = -2 * qJD(1);
t420 = 4 * qJD(1);
t419 = t59 / 0.2e1;
t418 = -0.2e1 * t185;
t415 = pkin(1) * t202;
t225 = t243 * qJD(3);
t51 = -qJD(2) * t378 + (t201 * t337 - t225) * pkin(5);
t99 = t243 * t407;
t410 = qJD(2) * t269 - t51 - t99;
t409 = 0.2e1 * t92;
t406 = pkin(14) * t207;
t188 = qJ(2) + qJ(3);
t171 = sin(t188);
t109 = t189 * t198 + t192 * t203;
t36 = -pkin(15) + (sin(pkin(21)) * t452 + cos(pkin(21)) * t109) * pkin(4);
t33 = t174 + t36;
t403 = t171 * t33;
t402 = t179 * t20;
t401 = t179 * t33;
t397 = t202 * pkin(15);
t375 = t166 * t107;
t376 = t165 * t111;
t41 = t375 - t376;
t37 = t41 ^ 2;
t396 = t207 * t37;
t395 = t207 * t50;
t394 = t207 * t92;
t391 = qJD(1) * t52;
t55 = qJD(2) * t57;
t102 = t114 * pkin(5);
t387 = t102 * t207;
t379 = t164 * t274;
t172 = cos(t188);
t374 = t172 * t207;
t361 = t197 * t207;
t351 = t206 * t197;
t350 = t206 * t202;
t298 = qJD(3) * pkin(1) * qJD(2);
t146 = t201 * t298;
t335 = qJD(3) * t201;
t306 = pkin(1) * t335;
t348 = t179 * t306 + t146;
t347 = qJD(3) ^ 2 * t325 + 0.2e1 * pkin(5) * t146;
t346 = t171 ^ 2 - t172 ^ 2;
t195 = sin(qJ(4));
t200 = cos(qJ(4));
t345 = t195 ^ 2 - t200 ^ 2;
t342 = qJD(1) * t179;
t341 = qJD(1) * t195;
t340 = qJD(1) * t200;
t339 = qJD(2) * t179;
t332 = qJD(4) * t200;
t330 = -qJD(2) - t179;
t326 = pkin(2) * t425;
t324 = pkin(14) * t446;
t323 = t20 * t422;
t21 = t179 * (t114 * t194 + t191 * t243);
t322 = t21 * t446;
t321 = -0.2e1 * t391;
t318 = pkin(1) * t355;
t316 = pkin(5) * t366;
t290 = t191 * t357;
t49 = t191 * t366 - t290 - t291 - t292;
t313 = t207 * pkin(15) * t49;
t312 = t49 * t395;
t311 = t207 * t458;
t286 = t197 * t337;
t309 = t454 * qJD(3) - 0.2e1 * t112 * t286;
t308 = pkin(1) * t338;
t307 = pkin(1) * t337;
t304 = t165 * (t153 * t360 + t164 * t177 + t434) * t166;
t303 = t33 * t374;
t302 = pkin(15) * t355;
t301 = t39 * t387;
t299 = pkin(15) * t337;
t289 = t207 * t243 * t114;
t288 = qJD(4) * (-t3 - t9);
t285 = t197 * t335;
t282 = -t376 / 0.4e1;
t280 = qJD(4) * (t240 - t4 + t453);
t279 = t197 * t328;
t278 = t202 * t328;
t276 = t200 * t327;
t270 = -0.2e1 * t299;
t268 = pkin(1) * t299;
t266 = -0.2e1 * t285;
t265 = t330 * qJD(3);
t264 = t195 * t200 * t396;
t262 = t171 * t172 * t342;
t261 = t197 * t278;
t256 = 0.2e1 * pkin(1) * t278;
t119 = (-pkin(15) + 0.2e1 * t174) * t202;
t255 = qJD(2) * t119 + t155 * t334;
t53 = -pkin(5) * t384 - t90;
t54 = t443 * pkin(5);
t253 = ((t197 * t53 + t202 * t54) * qJD(2) + t317) * t165 + ((t197 * t54 - t202 * t53) * qJD(2) + t240) * t166;
t252 = -t13 * t165 - t14 * t166;
t250 = t165 * (t128 * t193 + t392 * t129) - t166 * t78;
t249 = t196 * (t185 * t191 - t191 + t293) - t201 * (t152 - t194 - t294);
t247 = t37 * t195 * t276;
t242 = -t179 * t308 + t146;
t241 = (pkin(1) * t339 - t298) * t196;
t231 = t375 / 0.4e1 + t282;
t69 = t177 * t354 + t144 - t379;
t15 = t231 * qJD(4) + (-t304 + (t159 - 0.1e1 / 0.2e1) * t69) * qJD(1);
t239 = (qJD(4) * ((t360 / 0.4e1 - t153 / 0.4e1) * t166 + t282) - t15) * t420;
t238 = 0.4e1 * qJD(4) * (qJD(1) * (t69 * t159 - t304 + t379 / 0.2e1 + (-t177 / 0.2e1 + 0.1e1 / 0.4e1) * t354) + t15);
t32 = -t203 * t129 / 0.2e1 + t235;
t236 = t32 * qJD(1) + t89 * qJD(4) / 0.2e1;
t234 = -t205 * t22 - t382 * t92;
t233 = -t111 * t394 - t205 * t23;
t230 = -t85 * t196 + pkin(15) * t370 / 0.2e1;
t187 = t204 ^ 2;
t223 = (-t164 * t199 * t204 + t187 * t354 + t144) * t320 + t343 * (t164 * t187 + t352 * t359 + t434);
t222 = -t431 * t320 - t343 * (t169 + t436);
t220 = t231 * t423;
t218 = t201 * (t196 * t337 + t285);
t70 = -qJD(2) * t243 - t225;
t26 = qJD(2) * t51 + t70 * t407;
t211 = -t26 + (t250 + t89) * t205 + (-0.8e1 * t393 + (0.8e1 * t404 + 0.2e1 * t405) * t166 + 0.2e1 * t89 * t385 + 0.4e1 * t32) * t327;
t175 = t179 ^ 2;
t162 = pkin(1) * t351;
t148 = t210 * t362;
t145 = t197 * t355;
t142 = t210 * t286;
t127 = -0.2e1 * t261;
t126 = 0.2e1 * t261;
t124 = t343 * t207;
t121 = t171 * t374;
t113 = t316 + t155;
t103 = t346 * t207;
t100 = 0.2e1 * t343 * t328;
t97 = t244 * qJD(3);
t96 = t106 * qJD(3);
t87 = (-pkin(1) * t397 + t148) * t207;
t82 = qJD(2) * pkin(2) * t97;
t81 = (t142 - t268) * t446;
t71 = qJD(2) * t114 + t226;
t64 = (t149 - t316) * qJD(2) - t99;
t42 = pkin(15) + (t106 * t197 - t202 * t244) * pkin(2);
t38 = pkin(15) * t395;
t31 = (-t114 * t179 + t71) * qJD(1);
t30 = (t179 * t243 + t70) * qJD(1);
t18 = (-t179 * t49 + t21) * qJD(1);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t100, -t351, t126, -t350, 0, qJD(1) * t270, 0.2e1 * pkin(15) * t279, 0, 0, t114 * t70 * t422, (t438 * qJD(2) + t439 * qJD(3)) * t420, t179 * t70, -t243 * t71 * t446, t179 * t71, 0, (t255 * t196 - t429 * t201) * t446, (t429 * t196 + t255 * t201) * t446, t162, t81, 0, 0, 0, 0, 0, 0, t39 * t321, -t41 * t321, -t26, (-0.2e1 * t268 + (pkin(5) * t218 * t425 + t113 * t447 - 0.2e1 * t154 * t335 + t196 * t270 + t308 * t418) * pkin(5) + t309) * qJD(1), -0.2e1 * t247, 0.2e1 * t345 * t37 * t327, t195 * t238, 0.2e1 * t247, t200 * t238, 0, t211 * t195 - t451 * t200, t451 * t195 + t211 * t200, -0.8e1 * t231 * t391, (-0.2e1 * t137 * t284 + t309 + (t137 * t266 + t378 * t447 + 0.2e1 * t306) * pkin(5) + pkin(15) * t221 * t448 + 0.2e1 * t250 * t52) * qJD(1), t457, 0.8e1 * t223 * t328, qJD(2) * t56, -t457, qJD(2) * t55, 0, -t55 * t324, t56 * t324, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t256, -t452 * t256, t162, t81, t49 * t323, -0.8e1 * (t222 * qJD(2) + qJD(3) * t450) * qJD(1), t402, t50 * t322, t179 * t21, 0, pkin(15) * t322, pkin(15) * t323, 0, 0, t127, t100, -t351, t126, -t350, 0, (t249 * qJD(3) * pkin(2) + (pkin(2) * t426 - t397) * qJD(2)) * t446, 0.2e1 * t42 * t279 - 0.2e1 * (qJD(2) * t216 + qJD(3) * t228 + t202 * t96) * t202 * pkin(2) * qJD(1), -pkin(2) * t402, (((t161 * t218 + (t163 * t337 + t196 * t266) * t368) * t326 + (pkin(2) * t194 * t335 - t299) * t370) * pkin(2) + t433 * qJD(2) + (t230 * t326 + t428) * qJD(3)) * t446, 0.2e1 * t262, -0.2e1 * t346 * t342, t175 * t172, -0.2e1 * t262, -t175 * t171, 0, (t171 * t401 - t172 * t307) * t446, (t171 * t307 + t172 * t401) * t446, t162, (t307 * t36 + t142) * t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, -t124, 0, -t145, 0, 0, t302, -pkin(15) * t361, 0, 0, -t289, t438 * t424, t30, t289, t31, 0, -(t119 * t196 - t201 * t254) * t207 + t348, -t119 * t201 * t207 + (t155 * t361 + (-t185 * t207 + t265) * pkin(1)) * t196, 0, -t87, 0, 0, 0, 0, 0, 0, -t39 * t394, t41 * t394, t410 * qJD(1), -(t113 * t315 + t449) * t207 + t347, 0, 0, 0, 0, 0, 0, (t195 * t280 + t200 * t233) * t166 + (t195 * t288 + t200 * t234) * t165 + (t195 * t410 - t332 * t409) * qJD(1), (-t195 * t233 + t200 * t280) * t166 + (-t195 * t234 + t200 * t288) * t165 + (qJD(4) * t195 * t409 + t200 * t410) * qJD(1), -t92 * t220, (-t449 - (t160 * t182 - t173) * t411 + t250 * t92) * t207 + t347, t311, -0.4e1 * t207 * t223, 0, -t311, 0, 0, t57 * t406, -t245 * t406, 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t318, t452 * t318, 0, -t87, t312, t222 * t423, 0, -t312, t18, 0, -t313, t38, 0, 0, t145, -t124, 0, -t145, 0, 0, t302 - t82 + (-t179 * t97 - t207 * t426) * pkin(2), -t42 * t361 + (t330 * t96 + t455) * pkin(2), 0, -(t209 * t418 * t431 - t290 * t445 + t433) * t207, -t121, t103, 0, t121, 0, 0, (t172 * t415 - t403) * t207 + t348, -t303 + (-t171 * t355 + t196 * t265) * pkin(1), 0, -(t36 * t415 + t148) * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, t439 * t424, t30, t289, t31, 0, -(t123 * t201 + t155 * t356) * t207 + t242, -(-t123 * t196 + t155 * t357) * t207 + t241, 0, 0, 0, 0, 0, 0, 0, 0, t301, -t41 * t387, (-pkin(5) * t70 + t64) * qJD(1), -t456 + (-(t113 * t412 - t197 * t373) * t202 - t454 / 0.2e1 + t154 * t411) * t207, 0, 0, 0, 0, 0, 0, t64 * t341 + t200 * t301 + (t102 * t340 + t195 * t253) * qJD(4) + (-t70 * t341 - t200 * t440 + (t114 * t340 + t195 * t252) * qJD(4)) * pkin(5), t64 * t340 - t195 * t301 + (-t102 * t341 + t200 * t253) * qJD(4) + (-t70 * t340 + t195 * t440 + (-t114 * t341 + t200 * t252) * qJD(4)) * pkin(5), t102 * t220, -t456 + (-t143 / 0.2e1 - t98 / 0.2e1 - t250 * t102 + ((t137 * t197 - pkin(1)) * t201 + (pkin(15) - t378) * t356) * pkin(5)) * t207, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, t450 * t423, 0, -t312, t18, 0, -t313, t38, 0, 0, 0, 0, 0, 0, 0, 0, -t82 + (-t207 * t249 + t244 * t339) * pkin(2), (t455 + (t106 * t179 - t96) * qJD(2)) * pkin(2), 0, -((t224 * t414 + t230) * t326 + t209 * t295 + t428) * t207, -t121, t103, 0, t121, 0, 0, -t207 * t403 + t242, t241 - t303, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, -t345 * t396, t195 * t239, -t264, t200 * t239, 0, t200 * t419 * t446 - t52 * t340 + ((-t200 * t9 - t341 * t78) * qJD(4) + t332 * t9) * t166 + (t444 * t166 - t236 * t446 + t89 * t327 + t427) * t195, t195 * t419 * t422 + t89 * t276 + t52 * t341 + (t236 * t422 + (-t287 + t444) * t166 + t427) * t200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t1;
