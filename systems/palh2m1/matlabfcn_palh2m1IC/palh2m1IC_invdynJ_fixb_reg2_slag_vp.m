% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% palh2m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tau_reg [4x(5*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m1IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:49
% EndTime: 2020-05-03 01:03:54
% DurationCPUTime: 5.47s
% Computational Cost: add. (6891->423), mult. (16125->598), div. (0->0), fcn. (13593->10), ass. (0->210)
t200 = cos(qJ(3));
t201 = cos(qJ(2));
t175 = t200 * t201;
t196 = sin(qJ(3));
t197 = sin(qJ(2));
t295 = qJD(1) * t197;
t154 = -qJD(1) * t175 + t196 * t295;
t303 = t196 * t201;
t159 = t197 * t200 + t303;
t155 = qJD(1) * t159;
t195 = sin(qJ(4));
t348 = cos(qJ(4));
t109 = t348 * t154 + t155 * t195;
t104 = qJD(5) - t109;
t347 = pkin(2) * t200;
t179 = pkin(3) + t347;
t270 = t348 * t196;
t151 = pkin(2) * t270 + t179 * t195;
t144 = pkin(6) + t151;
t178 = pkin(3) * t200 + pkin(2);
t150 = -pkin(3) * t303 - t178 * t197;
t226 = t195 * t154 - t348 * t155;
t361 = pkin(4) * t226 - pkin(6) * t109;
t47 = qJD(1) * t150 + t361;
t376 = (qJD(5) * t144 + t47) * t104;
t198 = sin(qJ(1));
t343 = g(2) * t198;
t202 = cos(qJ(1));
t344 = g(1) * t202;
t164 = t343 + t344;
t276 = g(3) * t348;
t230 = t195 * t164 + t276;
t186 = t195 * g(3);
t242 = t164 * t348 - t186;
t304 = t196 * t197;
t143 = -pkin(3) * t304 + t178 * t201 + pkin(1);
t297 = qJD(1) * t143;
t46 = -pkin(4) * t109 - pkin(6) * t226 + t297;
t190 = qJDD(2) + qJDD(3);
t250 = qJDD(4) + t190;
t305 = t195 * t196;
t224 = t348 * t200 - t305;
t268 = qJD(4) * t348;
t293 = qJD(4) * t195;
t113 = t179 * t268 + (qJD(3) * t224 - t196 * t293) * pkin(2);
t275 = t348 * pkin(3);
t248 = qJD(3) * t275;
t346 = pkin(3) * t195;
t66 = t113 * qJD(2) + qJD(4) * t248 + t151 * qJDD(2) + qJDD(3) * t346;
t58 = pkin(6) * t250 + t66;
t375 = -qJD(5) * t46 - t197 * (t196 * t242 + t230 * t200) - (t196 * t230 - t242 * t200) * t201 - t58;
t189 = pkin(6) * t195;
t277 = pkin(4) * t348;
t245 = -t189 - t277;
t161 = pkin(3) - t245;
t170 = g(3) * pkin(4) - pkin(6) * t343;
t280 = pkin(4) * t343;
t171 = g(3) * pkin(6) + t280;
t181 = t189 + pkin(3);
t251 = t348 * t343;
t373 = -pkin(6) * t276 + t170 * t195 + t171 * t348 + t181 * t343 + (t161 + t245) * t344 + (-t186 - t251) * pkin(4);
t180 = pkin(2) * t201 + pkin(1);
t296 = qJD(1) * t180;
t307 = t164 * t196;
t258 = g(3) * t200 + t307;
t341 = g(3) * t196;
t351 = t164 * t200 - t341;
t369 = t351 * t197 + t201 * t258;
t152 = -pkin(2) * t305 + t179 * t348;
t132 = t152 * qJD(2) + t248;
t194 = sin(qJ(5));
t199 = cos(qJ(5));
t290 = qJD(5) * t199;
t368 = pkin(6) * t290 - t132 * t194 + t199 * t361;
t291 = qJD(5) * t194;
t367 = pkin(6) * t291 + t132 * t199 + t194 * t361;
t191 = qJD(2) + qJD(3);
t264 = qJD(4) + t191;
t124 = -pkin(4) * t264 - t132;
t366 = t124 * t109;
t309 = t155 * t191;
t364 = -pkin(4) * t276 + pkin(6) * t251 + g(3) * t181 + t170 * t348 + (-t171 + t280) * t195;
t130 = t151 * qJD(2) + qJD(3) * t346;
t123 = pkin(6) * t264 + t130;
t38 = -t123 * t194 + t199 * t46;
t39 = t123 * t199 + t194 * t46;
t239 = t39 * t194 + t38 * t199;
t360 = t373 * t196 + t200 * t364;
t359 = -t364 * t196 + t373 * t200;
t358 = qJDD(2) * t347 + t155 * t296 + t369;
t355 = -t130 * t264 + t226 * t297;
t252 = t154 * t191;
t352 = t175 - t304;
t210 = t191 * t352;
t208 = qJD(1) * t210;
t283 = t159 * qJDD(1);
t80 = -t208 - t283;
t64 = t80 - t252;
t211 = t191 * t159;
t222 = t352 * qJDD(1);
t206 = t211 * qJD(1) - t222;
t85 = t194 * t264 + t199 * t226;
t317 = qJD(5) * t85;
t34 = -t154 * t268 - t155 * t293 - t195 * t206 - t348 * t80;
t17 = -t194 * t34 - t199 * t250 + t317;
t350 = -t194 * t38 + t199 * t39;
t342 = g(2) * t202;
t345 = g(1) * t198;
t243 = -t342 + t345;
t225 = t195 * t200 + t270;
t114 = t179 * t293 + (qJD(3) * t225 + t196 * t268) * pkin(2);
t67 = -pkin(3) * (qJD(3) * t293 - t348 * qJDD(3)) - t114 * qJD(2) + t152 * qJDD(2);
t340 = g(3) * t201;
t319 = qJD(5) * t38;
t285 = qJDD(1) * t143;
t223 = qJD(3) * t159;
t294 = qJD(2) * t197;
t112 = -t178 * t294 + (-qJD(2) * t303 - t223) * pkin(3);
t298 = qJD(1) * t112;
t35 = qJD(4) * t226 + t195 * t80 - t348 * t206;
t6 = pkin(4) * t35 + pkin(6) * t34 + t285 + t298;
t2 = t194 * t6 + t199 * t58 + t319;
t339 = t2 * t199;
t318 = qJD(5) * t39;
t5 = t199 * t6;
t3 = -t194 * t58 - t318 + t5;
t338 = t3 * t194;
t237 = t199 * t264;
t83 = t194 * t226 - t237;
t337 = t85 * t83;
t176 = pkin(6) + t346;
t336 = pkin(6) - t176;
t332 = t109 * t38;
t331 = t109 * t39;
t16 = -qJD(5) * t237 - t194 * t250 + t199 * t34 + t226 * t291;
t330 = t16 * t194;
t33 = qJDD(5) + t35;
t329 = t194 * t33;
t327 = t194 * t83;
t325 = t199 * t85;
t316 = t226 * t143;
t314 = t109 * t143;
t204 = qJD(1) ^ 2;
t311 = t143 * t204;
t310 = t154 * t155;
t308 = t243 * t201;
t192 = t197 ^ 2;
t193 = t201 ^ 2;
t299 = t192 - t193;
t292 = qJD(5) * t176;
t156 = t224 * pkin(2);
t289 = t156 * qJD(2);
t157 = t225 * pkin(2);
t288 = t157 * qJD(2);
t287 = 0.2e1 * pkin(1);
t286 = qJD(1) * qJD(2);
t284 = qJDD(1) * t180;
t282 = t197 * qJDD(1);
t281 = t201 * qJDD(1);
t274 = pkin(2) * t294;
t271 = t197 * t204 * t201;
t269 = t196 ^ 2 + t200 ^ 2;
t267 = t197 * t286;
t266 = t201 * t286;
t265 = pkin(1) * t204 + t164;
t263 = qJD(5) * t83 - t16;
t256 = -qJD(5) * t123 + t243;
t255 = t80 - t283;
t254 = qJD(2) * (-qJD(3) + t191);
t253 = qJD(3) * (-qJD(2) - t191);
t249 = pkin(3) * t268;
t246 = t197 * t266;
t236 = -t113 * t104 - t366;
t125 = t225 * t243;
t126 = t224 * t243;
t235 = t125 * t201 + t126 * t197;
t59 = -pkin(4) * t250 - t67;
t232 = t59 * pkin(4) + pkin(6) * t338 + t124 * t130 + t367 * t39 + t368 * t38;
t229 = pkin(6) * t330 + t239 * t109 - t367 * t83 - t368 * t85;
t228 = pkin(3) * t293 - t288;
t227 = t348 * t352;
t221 = -t154 * t296 - t197 * t258 + t201 * t351;
t219 = -pkin(4) * t16 - t367 * t104 + t130 * t85 + t199 * t366;
t217 = pkin(4) * t17 + pkin(6) * t329 + t368 * t104 + t130 * t83 + t194 * t366;
t216 = t132 * t109 + t130 * t226;
t215 = t109 * t297 - t132 * t264;
t122 = -t348 * t159 - t195 * t352;
t115 = t125 * t197;
t44 = qJD(4) * t227 - t159 * t293 - t195 * t211 + t348 * t210;
t45 = t122 * qJD(4) - t195 * t210 - t348 * t211;
t12 = pkin(4) * t45 + pkin(6) * t44 + t112;
t121 = -t159 * t195 + t227;
t53 = pkin(4) * t121 - pkin(6) * t122 + t143;
t214 = qJD(5) * t122 * t124 + t104 * t12 + t33 * t53 - t115;
t213 = -qJD(5) * t104 * t53 + t122 * t59 - t124 * t44 + t164;
t212 = -0.2e1 * t309;
t203 = qJD(2) ^ 2;
t185 = pkin(1) * t345;
t177 = -t275 - pkin(4);
t158 = t159 * pkin(3);
t145 = -t152 - pkin(4);
t129 = t243 * (-pkin(4) * t195 + pkin(6) * t348);
t120 = t161 * t345 + (-t277 - t181) * t342;
t116 = t224 * t308;
t86 = -t154 ^ 2 + t155 ^ 2;
t65 = t206 - t309;
t48 = -qJD(1) * t158 + t361;
t41 = t194 * t48 + t199 * t289;
t40 = -t194 * t289 + t199 * t48;
t1 = [0, 0, 0, 0, 0, qJDD(1), t243, t164, 0, 0, qJDD(1) * t192 + 0.2e1 * t246, 0.2e1 * t197 * t281 - 0.2e1 * t286 * t299, -qJDD(2) * t197 - t201 * t203, qJDD(1) * t193 - 0.2e1 * t246, -qJDD(2) * t201 + t197 * t203, 0, t308 + (-t267 + t281) * t287, -t243 * t197 + (-t266 - t282) * t287, t164, t185 + (qJDD(1) * pkin(1) - t342) * pkin(1), t155 * t210 - t80 * t159, -t80 * t175 + ((t283 - t252) * t201 + t212 * t197) * t200 + (t212 * t201 + (t252 + t255) * t197) * t196, -t159 * t190 - t191 * t210, t154 * t211 - t206 * t352, -t190 * t352 + t191 * t211, 0, t154 * t274 + 0.2e1 * (-qJD(2) * t159 - t223) * t296 + (-t274 * qJD(1) + t243 + 0.2e1 * t284) * t352, -t243 * t303 + (0.2e1 * qJD(2) * pkin(2) * t155 - t200 * t243) * t197 + (-t208 + t255) * t180, ((t200 * t159 - t196 * t352) * qJDD(2) + ((-t196 * t159 - t200 * t352) * qJD(3) + t191 * t201 * t269) * qJD(2)) * pkin(2) + t164, (-0.2e1 * pkin(2) * t267 + t243 + t284) * t180, -t122 * t34 - t226 * t44, -t109 * t44 + t121 * t34 - t122 * t35 - t226 * t45, t122 * t250 - t264 * t44, -t109 * t45 + t121 * t35, -t121 * t250 - t264 * t45, 0, t116 - t115 + (qJD(1) * t121 - t109) * t112 + (qJD(1) * t45 + qJDD(1) * t121 + t35) * t143, (qJD(1) * t122 + t226) * t112 + (-qJD(1) * t44 + qJDD(1) * t122 - t34) * t143 - t235, -t121 * t66 - t122 * t67 - t130 * t45 + t132 * t44 + t164, (t285 + 0.2e1 * t298 + t243) * t143, -t44 * t325 + (-t16 * t199 - t291 * t85) * t122, (t194 * t85 + t199 * t83) * t44 + (t330 - t17 * t199 + (-t325 + t327) * qJD(5)) * t122, t122 * t199 * t33 - t121 * t16 + t45 * t85 + (-t122 * t291 - t199 * t44) * t104, -t44 * t327 + (t17 * t194 + t290 * t83) * t122, -t122 * t329 - t121 * t17 - t45 * t83 + (-t122 * t290 + t194 * t44) * t104, t104 * t45 + t121 * t33, t3 * t121 + t38 * t45 + (t116 + t214) * t199 + t213 * t194, -t2 * t121 - t39 * t45 + t213 * t199 + (-t126 * t201 - t214) * t194, (-t12 * t85 - t122 * t3 + t16 * t53 + t38 * t44 + (-t122 * t39 - t53 * t83) * qJD(5)) * t199 + (-t12 * t83 - t122 * t2 - t17 * t53 + t39 * t44 + (t122 * t38 + t53 * t85) * qJD(5)) * t194 + t235, (pkin(2) * t243 + t120 * t200 + t129 * t196) * t201 + (-t120 * t196 + t129 * t200) * t197 + t185 - pkin(1) * t342 + t239 * t12 + (t350 * qJD(5) + t2 * t194 + t3 * t199) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, t299 * t204, -t282, t271, -t281, qJDD(2), t197 * t265 + t340, -g(3) * t197 + t201 * t265, 0, 0, t310, t86, t64, -t310, t65, t190, (-t154 * t295 + t200 * t190 + t196 * t253) * pkin(2) + t358, (-t155 * t295 + (-qJDD(2) - t190) * t196 + t200 * t253) * pkin(2) + t221, (-t222 * t196 - t64 * t200) * pkin(2), (t340 + (t180 * t204 + t164) * t197 + t269 * qJDD(2) * pkin(2)) * pkin(2), 0, 0, 0, 0, 0, 0, -t114 * t264 + t152 * t250 + (t109 * t150 - t316) * qJD(1) + t355, -t113 * t264 - t151 * t250 + (-t150 * t226 - t314) * qJD(1) + t215, t109 * t113 + t114 * t226 - t151 * t35 + t152 * t34 + t216, t66 * t151 + t130 * t113 + t67 * t152 - t132 * t114 - t150 * t311 - (pkin(3) * t341 - t164 * t178) * t197 + (pkin(3) * t307 + g(3) * t178) * t201, 0, 0, 0, 0, 0, 0, t114 * t83 + t145 * t17 + (-t144 * t33 + t236) * t194 - t376 * t199 + t217, t114 * t85 - t145 * t16 + ((pkin(6) - t144) * t33 + t236) * t199 + t376 * t194 + t219, (t113 * t85 + t144 * t263 + t47 * t83 - t331) * t194 + (pkin(6) * t17 - t332 - t113 * t83 + t47 * t85 + (-t17 + t317) * t144) * t199 + t229, t124 * t114 + t59 * t145 + (g(3) * pkin(2) + t360) * t201 + (pkin(2) * t164 + t359) * t197 + (-t38 * t113 - t39 * t47 + (-t3 - t318) * t144) * t194 + (-pkin(6) * t2 + t113 * t39 - t38 * t47 + (t2 - t319) * t144) * t199 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t86, t64, -t310, t65, t190, pkin(2) * t196 * t254 + t358, (-qJDD(2) * t196 + t200 * t254) * pkin(2) + t221, 0, 0, 0, 0, 0, 0, 0, 0, t264 * t288 + (-t109 * t158 - t316) * qJD(1) + (t348 * t250 - t264 * t293) * pkin(3) + t355, t264 * t289 + (t158 * t226 - t314) * qJD(1) + (-t195 * t250 - t264 * t268) * pkin(3) + t215, (-t109 * t156 - t157 * t226) * qJD(2) + (t348 * t34 - t195 * t35 + (t109 * t348 + t195 * t226) * qJD(4)) * pkin(3) + t216, t158 * t311 + (-t130 * t156 + t132 * t157) * qJD(2) + (t66 * t195 + t67 * t348 + (t130 * t348 - t132 * t195) * qJD(4) + t369) * pkin(3), 0, 0, 0, 0, 0, 0, t177 * t17 + t228 * t83 + (-t176 * t33 - t366) * t194 + (-t176 * t290 - t194 * t249 - t40) * t104 + t217, -t177 * t16 + t228 * t85 + (t336 * t33 - t366) * t199 + t219 + (t292 * t194 - t249 * t199 + t41) * t104, t40 * t85 + t41 * t83 + (t176 * t263 + t249 * t85 - t331) * t194 + (t336 * t17 - t83 * t249 + t85 * t292 - t332) * t199 + t229, -pkin(6) * t339 - t124 * t288 + t59 * t177 - t38 * t40 - t39 * t41 + t360 * t201 + t359 * t197 + (t124 * t195 + t350 * t348) * qJD(4) * pkin(3) + (-qJD(5) * t239 - t338 + t339) * t176 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, -t83 ^ 2 + t85 ^ 2, t104 * t83 - t16, -t337, t85 * t104 - t17, t33, t39 * t104 - t124 * t85 + t375 * t194 + t256 * t199 + t5, t38 * t104 + t124 * t83 + t375 * t199 + (-t256 - t6) * t194, 0, 0;];
tau_reg = t1;
