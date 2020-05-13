% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tau_reg [4x(6*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m2IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:28
% EndTime: 2020-05-03 06:51:35
% DurationCPUTime: 6.42s
% Computational Cost: add. (12226->463), mult. (28524->684), div. (0->0), fcn. (25230->17), ass. (0->251)
t271 = qJDD(1) * pkin(1);
t170 = sin(qJ(1));
t174 = cos(qJ(1));
t347 = -g(1) * t170 + g(2) * t174;
t365 = -t347 + t271;
t169 = sin(qJ(2));
t276 = qJD(1) * qJD(2);
t244 = t169 * t276;
t173 = cos(qJ(2));
t272 = t173 * qJDD(1);
t364 = (t244 - t272) * pkin(4) - t365;
t337 = cos(qJ(4));
t245 = t337 * pkin(2);
t144 = t245 + pkin(5);
t167 = sin(qJ(5));
t335 = sin(qJ(4));
t254 = pkin(2) * t335;
t336 = cos(qJ(5));
t112 = t144 * t336 - t167 * t254;
t269 = pkin(5) * t336;
t172 = cos(qJ(3));
t145 = pkin(4) * t172 + pkin(2);
t168 = sin(qJ(3));
t258 = t335 * t168;
t117 = -pkin(4) * t258 + t145 * t337;
t110 = pkin(5) + t117;
t260 = t337 * t168;
t116 = pkin(4) * t260 + t145 * t335;
t302 = t116 * t167;
t68 = t110 * t336 - t302;
t57 = qJD(2) * t68 + qJD(3) * t112 + qJD(4) * t269;
t263 = t167 * t335;
t198 = t336 * t337 - t263;
t277 = pkin(2) * qJD(3);
t58 = (t117 * t336 - t302) * qJD(2) + t198 * t277;
t363 = t57 - t58;
t257 = t335 * t172;
t141 = pkin(5) * t257;
t207 = t172 * t337 - t258;
t86 = t141 * t173 + (t169 * t207 + t173 * t260) * pkin(5);
t356 = qJD(1) * t86;
t160 = t173 * pkin(4);
t345 = -pkin(1) - t160;
t131 = t345 * qJD(1);
t162 = qJD(2) + qJD(3);
t361 = t131 * t162;
t135 = g(1) * t174 + g(2) * t170;
t360 = t135 * pkin(4) * t169;
t123 = t207 * pkin(4);
t208 = t260 + t257;
t124 = t208 * pkin(4);
t304 = qJD(2) * (t123 * t336 - t124 * t167);
t251 = qJD(4) * t335;
t252 = qJD(4) * t337;
t82 = t145 * t252 + (qJD(3) * t207 - t168 * t251) * pkin(4);
t83 = -t145 * t251 + (-qJD(3) * t208 - t168 * t252) * pkin(4);
t33 = qJD(5) * t68 + t167 * t83 + t336 * t82;
t250 = qJD(5) * t336;
t78 = t144 * t250 + (qJD(4) * t198 - qJD(5) * t263) * pkin(2);
t322 = t33 - t78;
t231 = t304 + t322;
t354 = t168 * t169;
t219 = -t172 * t173 + t354;
t119 = t219 * qJD(1);
t220 = t168 * t173 + t169 * t172;
t120 = t220 * qJD(1);
t206 = t119 * t335 - t120 * t337;
t359 = qJD(4) * t206;
t76 = -t119 * t337 - t120 * t335;
t358 = -t167 * t76 + t206 * t336;
t357 = pkin(5) * t167;
t234 = t336 * t335;
t199 = -t167 * t337 - t234;
t259 = t336 * t116;
t111 = pkin(2) * t234 + t144 * t167;
t213 = -t110 * t167 - t259;
t268 = qJD(4) * t357;
t55 = -qJD(2) * t213 + qJD(3) * t111 + t268;
t319 = -(-t117 * t167 - t259) * qJD(2) - t199 * t277 - t55;
t311 = t112 - t68;
t309 = -t213 - t111;
t200 = -qJDD(1) * t345 - t364;
t352 = t200 * t172;
t351 = t219 * t337 + t220 * t335;
t350 = -t245 + t117;
t349 = t254 - t116;
t291 = pkin(2) * t172;
t146 = -pkin(4) - t291;
t106 = pkin(2) * t354 + t146 * t173 - pkin(1);
t275 = qJDD(1) * t106;
t285 = qJD(2) * t173;
t286 = qJD(2) * t169;
t80 = -t146 * t286 + (qJD(3) * t220 + t168 * t285) * pkin(2);
t305 = qJD(1) * t80;
t348 = t305 + t275;
t295 = (-t146 - t291) * t169;
t346 = qJD(1) * t295;
t34 = qJD(5) * t213 - t167 * t82 + t336 * t83;
t283 = qJD(5) * t167;
t79 = -t144 * t283 + (qJD(4) * t199 - qJD(5) * t234) * pkin(2);
t230 = t34 - t79 + (-t123 * t167 - t124 * t336) * qJD(2);
t228 = pkin(2) * t251;
t344 = -qJD(2) * t124 + t228 + t83;
t343 = -qJD(2) * t123 + qJD(4) * t245 - t82;
t342 = pkin(4) * (t168 ^ 2 + t172 ^ 2);
t188 = t162 * t219;
t189 = t162 * t220;
t182 = t188 * t335 - t189 * t337;
t88 = -t219 * t335 + t220 * t337;
t180 = qJD(4) * t88 - t182;
t308 = pkin(3) * t358;
t143 = pkin(5) * t337 + pkin(2);
t210 = -pkin(5) * t258 + t143 * t172;
t109 = -pkin(4) - t210;
t115 = t143 * t168 + t141;
t303 = t115 * t173;
t66 = -t109 * t169 + t303;
t223 = -qJD(1) * t66 + qJD(6) * t309 + t308;
t175 = qJD(2) ^ 2;
t341 = qJDD(2) * t169 + t173 * t175;
t51 = pkin(2) * (qJD(3) * t252 + qJDD(3) * t335) + t82 * qJD(2) + t116 * qJDD(2);
t211 = t219 * qJDD(1);
t184 = -qJD(1) * t189 - t211;
t185 = -qJD(1) * t188 + qJDD(1) * t220;
t179 = -t184 * t337 + t185 * t335;
t178 = t179 - t359;
t340 = qJD(4) * t351 + t188 * t337 + t189 * t335;
t329 = g(3) * t173;
t166 = sin(qJ(6));
t171 = cos(qJ(6));
t156 = qJD(4) + t162;
t226 = qJD(5) + t156;
t216 = t171 * t226;
t161 = qJDD(2) + qJDD(3);
t155 = qJDD(4) + t161;
t221 = qJDD(5) + t155;
t282 = qJD(6) * t166;
t25 = t119 * t252 + t120 * t251 - t184 * t335 - t185 * t337;
t9 = t167 * t178 - t206 * t283 + t25 * t336 - t250 * t76;
t4 = qJD(6) * t216 + t166 * t221 + t171 * t9 - t282 * t358;
t328 = t166 * t4;
t242 = -t166 * t9 + t171 * t221;
t42 = -t166 * t226 - t171 * t358;
t5 = qJD(6) * t42 + t242;
t327 = t171 * t5;
t40 = -t166 * t358 + t216;
t326 = t42 * t40;
t325 = t206 * t76;
t324 = t166 * t363 - t171 * t356;
t323 = -t166 * t356 - t171 * t363;
t69 = t169 * t210 + t303;
t318 = -t66 + t69;
t316 = pkin(4) * qJD(2);
t315 = t166 * t40;
t313 = t171 * t42;
t176 = qJD(1) ^ 2;
t64 = t109 * t173 + t115 * t169 - pkin(1);
t312 = t176 * t64;
t253 = qJD(3) * t335;
t284 = qJD(3) * t143;
t32 = -t109 * t286 + (t168 * t284 + (qJD(4) * t208 + t172 * t253) * pkin(5)) * t173 + (t172 * t284 + (qJD(4) * t207 - t168 * t253) * pkin(5)) * t169 + t115 * t285;
t306 = qJD(1) * t32;
t301 = t166 * t170;
t300 = t166 * t174;
t299 = t170 * t171;
t298 = t171 * t174;
t165 = qJ(2) + qJ(3);
t163 = t169 ^ 2;
t164 = t173 ^ 2;
t294 = t163 - t164;
t290 = qJD(1) * t106;
t289 = qJD(1) * t169;
t281 = qJD(6) * t171;
t280 = qJDD(1) * t64;
t159 = qJ(4) + t165;
t153 = qJ(5) + t159;
t140 = cos(t153);
t279 = t140 * pkin(3);
t147 = t269 + pkin(3);
t278 = pkin(3) - t147;
t273 = t169 * qJDD(1);
t270 = t55 + t319;
t193 = t167 * t351;
t54 = t336 * t88 - t193;
t267 = t54 * t282;
t266 = t54 * t281;
t265 = t169 * t176 * t173;
t264 = t166 * t336;
t262 = t171 * t336;
t10 = -qJD(5) * t358 - t167 * t25 + t178 * t336;
t8 = -qJDD(6) + t10;
t261 = t309 * t8;
t158 = cos(t165);
t142 = pkin(2) * t158;
t256 = t142 - t345;
t139 = sin(t153);
t14 = t33 * qJD(2) + t78 * qJD(3) - t213 * qJDD(2) + t111 * qJDD(3) + (qJD(4) * t250 + qJDD(4) * t167) * pkin(5);
t255 = g(3) * t139 - t14;
t3 = pkin(3) * t10 + t280 + t306;
t243 = qJD(6) * t55 + t3;
t46 = -t167 * t206 - t336 * t76;
t241 = qJD(1) * t318;
t238 = t162 ^ 2;
t237 = qJD(1) * t162;
t236 = -0.2e1 * pkin(1) * t276;
t43 = -qJD(6) + t46;
t235 = t43 * t250;
t233 = t173 * t244;
t232 = t347 * t139;
t191 = t336 * t351;
t11 = qJD(5) * t191 + t167 * t180 + t283 * t88 + t336 * t340;
t227 = t11 * t43 - t54 * t8;
t28 = pkin(3) * t46 + qJD(1) * t64;
t16 = -t166 * t55 - t171 * t28;
t17 = -t166 * t28 + t171 * t55;
t224 = -t16 * t171 - t166 * t17;
t150 = sin(t159);
t215 = t135 * t150;
t214 = -t171 * t243 + t28 * t282;
t209 = pkin(1) * t176 + t135;
t205 = -t347 + 0.2e1 * t271;
t89 = pkin(2) * t253 + qJD(2) * t116;
t53 = t167 * t88 + t191;
t36 = pkin(3) * t53 + t64;
t50 = pkin(3) * t226 + t57;
t12 = -qJD(5) * t193 - t167 * t340 + t180 * t336 + t250 * t88;
t6 = pkin(3) * t12 + t32;
t204 = qJD(6) * t50 * t54 + t36 * t8 + t43 * t6;
t15 = qJD(2) * t34 + qJD(3) * t79 - qJD(5) * t268 + qJDD(2) * t68 + qJDD(3) * t112 + qJDD(4) * t269;
t13 = pkin(3) * t221 + t15;
t202 = -qJD(6) * t36 * t43 - t11 * t50 + t13 * t54;
t201 = -t173 * t237 - t273;
t151 = cos(t159);
t197 = -g(3) * t151 + t215;
t52 = qJD(2) * t83 - qJD(3) * t228 + qJDD(2) * t117 + qJDD(3) * t245;
t194 = 0.2e1 * t361;
t192 = -t160 * g(3) + t360;
t186 = -t168 * t200 + t172 * t194;
t138 = pkin(5) * t151;
t103 = t140 * t298 - t301;
t102 = -t140 * t300 - t299;
t101 = -t140 * t299 - t300;
t100 = t140 * t301 - t298;
t96 = t138 + t256;
t91 = qJD(2) * t117 + qJD(3) * t245;
t30 = qJD(1) * t69 - t308;
t21 = -t166 * t30 + t171 * t304;
t20 = -t166 * t304 - t171 * t30;
t2 = -t14 * t166 + t214;
t1 = qJD(6) * t16 + t14 * t171 - t166 * t3;
t7 = [0, 0, 0, 0, 0, qJDD(1), -t347, t135, 0, 0, qJDD(1) * t163 + 0.2e1 * t233, 0.2e1 * t169 * t272 - 0.2e1 * t276 * t294, t341, qJDD(1) * t164 - 0.2e1 * t233, qJDD(2) * t173 - t169 * t175, 0, t169 * t236 + t173 * t205, -t169 * t205 + t173 * t236, -t135, t365 * pkin(1), -t120 * t188 + t185 * t220, (t188 * t219 - t189 * t220) * qJD(1) - 0.2e1 * t220 * t211 + t162 * (t119 * t219 - t120 * t220), t161 * t220 - t162 * t188, t119 * t189 - t184 * t219, -t161 * t219 - t162 * t189, 0, (t168 * t194 + t352) * t173 + (t119 * t316 + t186) * t169, t186 * t173 + (t120 * t316 - t352 + (-t237 * t345 - t361) * t168) * t169, -t341 * t342 - t135, pkin(4) * t131 * t286 + t345 * t364, t206 * t340 - t25 * t88, -t206 * t182 + t25 * t351 - t340 * t76 + (-t179 + 0.2e1 * t359) * t88, t155 * t88 - t156 * t340, t178 * t351 - t180 * t76, -t155 * t351 - t156 * t180, 0, t106 * t178 - t151 * t347 + t180 * t290 + t348 * t351 - t80 * t76, -t106 * t25 + t150 * t347 - t206 * t80 - t290 * t340 + t348 * t88, -t180 * t89 + t340 * t91 - t351 * t51 - t52 * t88 - t135, -t347 * t256 + (t275 + 0.2e1 * t305) * t106, t11 * t358 - t54 * t9, -t10 * t54 + t11 * t46 + t12 * t358 + t53 * t9, -t11 * t226 + t221 * t54, t10 * t53 + t12 * t46, -t12 * t226 - t221 * t53, 0, (qJD(1) * t53 + t46) * t32 - t347 * t140 + (qJD(1) * t12 + qJDD(1) * t53 + t10) * t64, (qJD(1) * t54 - t358) * t32 + (-qJD(1) * t11 + qJDD(1) * t54 - t9) * t64 + t232, t11 * t57 - t12 * t55 - t14 * t53 - t15 * t54 - t135, -t347 * t96 + (t280 + 0.2e1 * t306) * t64, -t42 * t267 + (-t11 * t42 - t4 * t54) * t171, (t166 * t42 + t171 * t40) * t11 + (t328 - t327 + (-t313 + t315) * qJD(6)) * t54, -t12 * t42 + t171 * t227 + t267 * t43 + t4 * t53, t40 * t266 + (-t11 * t40 + t5 * t54) * t166, t12 * t40 - t166 * t227 + t266 * t43 + t5 * t53, t12 * t43 + t53 * t8, -g(1) * t101 - g(2) * t103 - t12 * t16 + t166 * t202 + t171 * t204 - t2 * t53, -g(1) * t100 - g(2) * t102 + t1 * t53 + t12 * t17 - t166 * t204 + t171 * t202, (t11 * t16 - t2 * t54 - t36 * t4 + t42 * t6 + (-t17 * t54 + t36 * t40) * qJD(6)) * t171 + (-t1 * t54 + t11 * t17 + t36 * t5 + t40 * t6 + (t16 * t54 - t36 * t42) * qJD(6)) * t166 + t232, t224 * t6 + (-t1 * t166 - t2 * t171 + (t16 * t166 - t17 * t171) * qJD(6)) * t36 - t347 * (t96 + t279); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t265, t294 * t176, t273, t265, t272, qJDD(2), t169 * t209 - t329, g(3) * t169 + t173 * t209, 0, 0, 0, 0, 0, 0, 0, 0, (-t119 * t289 + t161 * t172 - t168 * t238) * pkin(4), (-t120 * t289 - t161 * t168 - t172 * t238) * pkin(4), ((-t119 * t162 + t172 * t201) * t172 + (t120 * t162 + t168 * t201) * t168) * pkin(4), (-t329 + (-qJD(1) * t131 + t135) * t169 + qJDD(2) * t342) * pkin(4), 0, 0, 0, 0, 0, 0, t155 * t350 + t156 * t344 + t346 * t76, t155 * t349 + t156 * t343 + t206 * t346, t178 * t349 + t206 * t344 + t25 * t350 - t343 * t76, t51 * t116 + t89 * t82 + t52 * t117 + t91 * t83 - g(3) * (t142 + t160) + t360 - t295 * t176 * t106 + (t123 * t89 - t124 * t91) * qJD(2) + (-t52 * t337 - t51 * t335 + g(3) * t158 + (t335 * t91 - t337 * t89) * qJD(4)) * pkin(2), 0, 0, 0, 0, 0, 0, -t221 * t311 + t226 * t230 + t241 * t46, -t221 * t309 - t226 * t231 - t241 * t358, -t10 * t309 + t230 * t358 - t231 * t46 - t311 * t9, t14 * t309 - t15 * t311 + t230 * t57 + t231 * t55 + t312 * t318 + t192, 0, 0, 0, 0, 0, 0, -t311 * t5 + t166 * t261 + t230 * t40 + (t166 * t322 + t171 * t223 - t20) * t43, t311 * t4 + t171 * t261 + t230 * t42 + (-t166 * t223 + t171 * t322 + t21) * t43, -t20 * t42 - t21 * t40 + (t223 * t42 - t309 * t5 - t322 * t40) * t171 + (t223 * t40 - t309 * t4 + t322 * t42) * t166, t16 * t20 + t17 * t21 - t311 * t13 + t230 * t50 + (t1 * t309 - t16 * t223 + t17 * t322) * t171 + (-t16 * t322 - t17 * t223 - t2 * t309) * t166 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, t206 ^ 2 - t76 ^ 2, -t156 * t76 - t25, -t325, -t156 * t206 - t178, t155, t156 * t89 + t206 * t290 + t197 + t52, g(3) * t150 + t135 * t151 + t91 * t156 - t290 * t76 - t51, 0, 0, 0, 0, 0, 0, 0, 0, t221 * t269 - t46 * t356 + (-pkin(5) * t283 + t319) * t226, -t221 * t357 + t358 * t356 + (-pkin(5) * t250 - t363) * t226, -t55 * t358 + t58 * t46 - t57 * t46 + t270 * t358 + (t336 * t9 - t10 * t167 + (-t167 * t358 - t336 * t46) * qJD(5)) * pkin(5), -t86 * t312 - t55 * t58 + t270 * t57 + (t336 * t15 + t14 * t167 + (-t167 * t57 + t336 * t55) * qJD(5) + t197) * pkin(5), 0, 0, 0, 0, 0, 0, -t278 * t5 + t324 * t43 + t319 * t40 + (t166 * t235 + (-qJD(5) * t40 + t166 * t8 + t281 * t43) * t167) * pkin(5), -t323 * t43 + t319 * t42 + t278 * t4 + (t171 * t235 + (-qJD(5) * t42 + t171 * t8 - t282 * t43) * t167) * pkin(5), t324 * t42 + t323 * t40 + ((-t262 * t40 + t264 * t42) * qJD(5) + (-t328 - t327 + (t313 + t315) * qJD(6)) * t167) * pkin(5), t13 * t147 - g(3) * (t138 + t279) + t319 * t50 - t323 * t17 - t324 * t16 + (g(3) * t140 - t13) * pkin(3) + (t215 + (-t16 * t264 + t17 * t262) * qJD(5) + (-qJD(5) * t50 + qJD(6) * t224 + t1 * t171 - t166 * t2) * t167) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t326, -t40 ^ 2 + t42 ^ 2, -t40 * t43 - t4, -t326, -t242 + (-qJD(6) - t43) * t42, -t8, -g(1) * t102 + g(2) * t100 + t166 * t255 - t17 * t43 - t42 * t50 + t214, g(1) * t103 - g(2) * t101 - t16 * t43 + t40 * t50 + t243 * t166 + (qJD(6) * t28 + t255) * t171, 0, 0;];
tau_reg = t7;
