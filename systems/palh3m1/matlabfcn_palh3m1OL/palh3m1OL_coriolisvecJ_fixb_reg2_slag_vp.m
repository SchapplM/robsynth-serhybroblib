% Calculate inertial parameters regressor of coriolis joint torque vector for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% tauc_reg [10x(10*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh3m1OL_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:11:00
% EndTime: 2020-04-20 17:12:17
% DurationCPUTime: 8.06s
% Computational Cost: add. (7470->477), mult. (19451->699), div. (0->0), fcn. (15934->16), ass. (0->268)
t202 = cos(qJ(5));
t199 = sin(qJ(3));
t200 = sin(qJ(2));
t203 = cos(qJ(2));
t354 = cos(qJ(3));
t285 = t354 * t203;
t237 = t199 * t200 - t285;
t228 = t237 * qJD(2);
t213 = t237 * qJD(3) + t228;
t106 = t213 * qJD(1);
t304 = qJD(1) * t200;
t149 = -qJD(1) * t285 + t199 * t304;
t159 = t199 * t203 + t200 * t354;
t151 = qJD(1) * t159;
t198 = sin(qJ(4));
t298 = qJD(2) + qJD(3);
t211 = t298 * t159;
t210 = t211 * qJD(1);
t353 = cos(qJ(4));
t280 = qJD(4) * t353;
t303 = qJD(4) * t198;
t235 = t353 * t106 + t149 * t280 + t151 * t303 + t198 * t210;
t236 = -t198 * t149 + t151 * t353;
t297 = -qJD(3) - qJD(4);
t275 = qJD(2) - t297;
t246 = t202 * t275;
t197 = sin(qJ(5));
t302 = qJD(5) * t197;
t21 = -qJD(5) * t246 - t202 * t235 - t236 * t302;
t17 = t202 * t21;
t72 = t197 * t275 - t202 * t236;
t319 = qJD(5) * t72;
t22 = t197 * t235 + t319;
t301 = qJD(5) * t202;
t70 = -t197 * t236 - t246;
t358 = t22 * t197 + t70 * t301;
t97 = t353 * t149 + t198 * t151;
t366 = t202 * t97;
t363 = qJD(5) - t97;
t371 = t197 * t363;
t1 = t366 * t70 - t371 * t72 - t17 - t358;
t208 = t353 * t211;
t274 = -qJD(1) * t208 + t198 * t106;
t38 = -t236 * qJD(4) + t274;
t34 = t38 * t202;
t6 = -t236 * t70 - t363 * t371 + t34;
t294 = t354 * pkin(1);
t262 = qJD(2) * t294;
t223 = pkin(4) * t298 - t262;
t152 = t198 * t223;
t284 = t353 * t199;
t265 = pkin(1) * t284;
t127 = qJD(2) * t265 - t152;
t123 = pkin(10) * t275 - t127;
t181 = -t203 * pkin(1) - pkin(12);
t163 = qJD(1) * t181;
t125 = -pkin(4) * t149 + t163;
t41 = -pkin(8) * t97 + pkin(10) * t236 + t125;
t244 = t123 * t197 - t202 * t41;
t370 = t363 * t244;
t258 = t353 * t354;
t336 = pkin(1) * qJD(2);
t238 = t258 * t336;
t161 = qJD(3) * t238;
t349 = pkin(1) * t199;
t291 = qJD(3) * t349;
t177 = qJD(2) * t291;
t234 = t177 * t198 - t161;
t332 = t125 * t97;
t369 = -t234 - t332;
t31 = t123 * t202 + t197 * t41;
t300 = qJD(1) * qJD(2);
t277 = t200 * t300;
t179 = pkin(1) * t277;
t84 = -pkin(4) * t210 + t179;
t10 = t38 * pkin(8) - pkin(10) * t235 + t84;
t292 = t199 * t336;
t169 = t198 * t292;
t219 = t353 * t223;
t356 = t169 + t219;
t82 = qJD(4) * t356 + t234;
t4 = -qJD(5) * t31 + t202 * t10 - t197 * t82;
t368 = -t363 * t31 - t4;
t18 = t21 * t197;
t7 = -t18 + (t301 - t366) * t72;
t339 = t38 * t197 + t301 * t363;
t5 = t236 * t72 - t363 * t366 + t339;
t122 = -pkin(8) * t275 - t356;
t333 = t122 * t97;
t343 = t363 * t236;
t365 = t97 * t236;
t352 = cos(qJ(7));
t282 = t352 * t203;
t257 = qJD(1) * t282;
t195 = sin(qJ(7));
t281 = t195 * t304;
t147 = -t257 + t281;
t158 = t195 * t203 + t200 * t352;
t150 = qJD(1) * t158;
t194 = sin(pkin(15));
t337 = cos(pkin(15));
t350 = sin(qJ(8));
t351 = cos(qJ(8));
t156 = t194 * t350 + t337 * t351;
t224 = -t194 * t351 + t337 * t350;
t241 = -t147 * t224 + t150 * t156;
t242 = t147 * t156 + t150 * t224;
t364 = t242 * t241;
t189 = qJD(2) + qJD(7);
t27 = t236 ^ 2 - t97 ^ 2;
t232 = t198 * t354 + t284;
t216 = (qJD(3) * t232 + t199 * t280) * pkin(1);
t83 = qJD(2) * t216 - qJD(4) * t152;
t207 = t125 * t236 + t83;
t25 = -t275 * t97 + t235;
t259 = t122 * t301 - t83 * t197 - t236 * t31;
t276 = t122 * t302 - t236 * t244;
t26 = -t298 * t236 - t274;
t52 = -pkin(8) * t236 - pkin(10) * t97;
t187 = t200 * t336;
t360 = 0.2e1 * t187;
t359 = -0.2e1 * t300;
t278 = t352 * qJD(7);
t279 = t352 * qJD(2);
t355 = -t279 - t278;
t247 = t197 * t244 + t202 * t31;
t145 = t156 * qJD(8);
t146 = t224 * qJD(8);
t293 = t195 * t336;
t348 = pkin(3) * t194;
t153 = t189 * t348 + t293;
t261 = pkin(1) * t279;
t289 = t337 * pkin(3);
t154 = t189 * t289 + t261;
t243 = qJD(7) * t261;
t260 = qJD(7) * t293;
t46 = t145 * t153 + t146 * t154 + t156 * t260 + t224 * t243;
t75 = t163 + (t147 * t337 - t150 * t194) * pkin(3);
t220 = t241 * t75 + t46;
t45 = -t145 * t154 + t146 * t153 - t156 * t243 + t224 * t260;
t217 = -t242 * t75 - t45;
t188 = qJD(8) + t189;
t312 = t195 * t200;
t251 = t189 * t312;
t307 = t189 * t257;
t104 = qJD(1) * t251 - t307;
t121 = t189 * t158;
t105 = t121 * qJD(1);
t23 = t104 * t156 + t105 * t224 + t145 * t147 + t146 * t150;
t13 = -t188 * t242 + t23;
t24 = -t104 * t224 + t156 * t105 + t145 * t150 - t146 * t147;
t14 = -t188 * t241 + t24;
t3 = -qJD(5) * t244 + t10 * t197 + t202 * t82;
t12 = t241 ^ 2 - t242 ^ 2;
t347 = pkin(4) * t151;
t2 = t3 * t202;
t344 = t70 * t72;
t222 = t353 * t237;
t114 = -t159 * t198 - t222;
t335 = t114 * t38;
t229 = t198 * t237;
t115 = -t159 * t353 + t229;
t334 = t115 * t83;
t330 = t197 * t31;
t328 = t197 * t70;
t325 = t202 * t72;
t20 = t22 * t202;
t231 = pkin(1) * (t156 * t195 + t224 * t352);
t322 = qJD(2) * t231 - (t145 * t194 + t146 * t337) * pkin(3);
t230 = pkin(1) * (-t156 * t352 + t195 * t224);
t321 = qJD(2) * t230 - (-t145 * t337 + t146 * t194) * pkin(3);
t318 = qJD(5) * t363;
t317 = t149 * t151;
t316 = t150 * t147;
t315 = t150 * t163;
t314 = t163 * t147;
t311 = t198 * t199;
t201 = cos(qJ(6));
t206 = qJD(1) ^ 2;
t310 = t201 * t206;
t309 = t203 * t206;
t196 = sin(qJ(6));
t306 = t196 ^ 2 - t201 ^ 2;
t305 = t200 ^ 2 - t203 ^ 2;
t299 = qJD(1) * qJD(6);
t296 = -t244 * t366 + t97 * t330 + t2;
t186 = pkin(1) * t304;
t290 = 0.2e1 * t299;
t288 = t196 * t310;
t287 = t200 * t309;
t283 = t352 * t189;
t185 = -t294 + pkin(4);
t142 = t198 * t185 - t265;
t137 = pkin(10) + t142;
t47 = -t347 + t52;
t42 = t186 + t47;
t271 = qJD(5) * t137 + t42;
t270 = qJD(3) * t262 - t163 * t149;
t269 = t163 * t151 + t177;
t268 = pkin(6) * t290;
t267 = pkin(12) * t359;
t266 = t199 * t298;
t264 = t163 * t186;
t263 = pkin(4) * t280;
t256 = t196 * t201 * t299;
t255 = t203 * t277;
t139 = t232 * t336;
t253 = t303 * pkin(4) + t139;
t59 = t153 * t224 - t154 * t156;
t60 = t153 * t156 + t154 * t224;
t252 = t241 * t60 + t242 * t59;
t249 = t127 * t236 + t356 * t97;
t248 = -t202 * t244 + t330;
t245 = t298 * t354;
t239 = t302 * t363 - t34;
t141 = pkin(1) * t311 + t185 * t353;
t102 = t185 * t280 + (t199 * t303 + (-t258 + t311) * qJD(3)) * pkin(1);
t233 = -t102 * t363 - t137 * t38 - t333;
t88 = (t147 * t194 + t150 * t337) * pkin(3);
t43 = -qJD(4) * t222 - t159 * t303 - t198 * t211 - t213 * t353;
t44 = qJD(4) * t229 - t159 * t280 + t198 * t213 - t208;
t99 = -pkin(4) * t211 + t187;
t11 = t44 * pkin(8) + t43 * pkin(10) + t99;
t130 = -pkin(4) * t237 + t181;
t50 = t114 * pkin(8) - t115 * pkin(10) + t130;
t226 = qJD(5) * t115 * t122 + t11 * t363 + t38 * t50;
t225 = -t122 * t43 - t318 * t50 - t334;
t183 = t198 * pkin(4) + pkin(10);
t221 = -t183 * t38 - t263 * t363 - t333;
t218 = -qJD(5) * t248 - t4 * t197 + t2;
t209 = t211 * t349;
t205 = qJD(2) ^ 2;
t204 = qJD(6) ^ 2;
t184 = -pkin(4) * t353 - pkin(8);
t165 = pkin(1) * t352 + t289;
t164 = pkin(1) * t195 + t348;
t157 = -t282 + t312;
t140 = -t238 + t169;
t136 = -pkin(8) - t141;
t128 = t186 - t347;
t124 = t163 * t360;
t120 = t355 * t203 + t251;
t108 = (-t156 * t194 - t224 * t337) * pkin(3);
t107 = (-t156 * t337 + t194 * t224) * pkin(3);
t103 = -t185 * t303 + t216;
t89 = (t157 * t337 - t158 * t194) * pkin(3) + t181;
t81 = t186 + t88;
t80 = -t156 * t164 - t165 * t224;
t79 = -t156 * t165 + t164 * t224;
t78 = -t147 ^ 2 + t150 ^ 2;
t77 = -t149 ^ 2 + t151 ^ 2;
t69 = -t156 * t158 + t157 * t224;
t68 = t156 * t157 + t158 * t224;
t65 = -t151 * t298 + t210;
t64 = -t149 * t298 + t106;
t62 = t307 + (t147 - t281) * t189;
t53 = t187 + (t120 * t194 + t121 * t337) * pkin(3);
t51 = t179 + (t104 * t194 + t337 * t105) * pkin(3);
t49 = qJD(7) * t231 + t145 * t164 + t146 * t165;
t48 = qJD(7) * t230 - t145 * t165 + t146 * t164;
t40 = t197 * t52 + t202 * t356;
t39 = -t197 * t356 + t202 * t52;
t36 = t140 * t202 + t197 * t47;
t35 = -t140 * t197 + t202 * t47;
t29 = -t120 * t224 + t121 * t156 + t145 * t158 - t146 * t157;
t28 = t120 * t156 + t121 * t224 + t145 * t157 + t146 * t158;
t8 = t371 * t70 - t20;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t255, t305 * t359, t205 * t203, -0.2e1 * t255, -t205 * t200, 0, t200 * t267, t203 * t267, 0, 0, -t106 * t159 - t151 * t213, t106 * t237 - t159 * t210 + t298 * (t149 * t237 - t151 * t159), t298 * t213, t149 * t211 + t210 * t237, t211 * t298, 0, -t149 * t187 - 0.2e1 * t163 * t211 - t186 * t228, t181 * t106 - t151 * t187 - t159 * t179 + t163 * t213, -qJD(2) * t209 + t159 * t177 + t228 * t262, t124, t115 * t235 + t236 * t43, -t114 * t235 - t115 * t38 + t236 * t44 - t43 * t97, -t275 * t43, -t44 * t97 + t335, -t275 * t44, 0, t114 * t84 + t125 * t44 + t130 * t38 - t97 * t99, t115 * t84 - t125 * t43 + t130 * t235 - t236 * t99, -t114 * t82 + t127 * t44 + t356 * t43 - t334, t125 * t99 + t130 * t84, -t43 * t325 + (-t302 * t72 - t17) * t115, (t197 * t72 + t202 * t70) * t43 + (t18 - t20 + (-t325 + t328) * qJD(5)) * t115, -t202 * t363 * t43 - t114 * t21 - t115 * t239 + t44 * t72, t358 * t115 - t43 * t328, -t114 * t22 - t115 * t339 + t371 * t43 - t44 * t70, t363 * t44 + t335, t114 * t4 + t197 * t225 + t202 * t226 - t244 * t44, -t114 * t3 - t197 * t226 + t202 * t225 - t31 * t44, (-t11 * t72 - t115 * t4 + t21 * t50 - t244 * t43 + (-t115 * t31 - t50 * t70) * qJD(5)) * t202 + (-t11 * t70 - t115 * t3 - t22 * t50 + t31 * t43 + (-t115 * t244 + t50 * t72) * qJD(5)) * t197, t248 * t11 + (qJD(5) * t247 + t197 * t3 + t202 * t4) * t50, 0.2e1 * t256, -t306 * t290, t204 * t201, -0.2e1 * t256, -t204 * t196, 0, t196 * t268, t201 * t268, 0, 0, -t104 * t158 - t120 * t150, t104 * t157 - t105 * t158 + t120 * t147 - t121 * t150, -t189 * t120, t105 * t157 + t121 * t147, -t189 * t121, 0, t105 * t181 + t121 * t163 + (qJD(1) * t157 + t147) * t187, -t104 * t181 - t120 * t163 + t150 * t360, (t352 * t120 - t121 * t195 + (-t157 * t352 + t158 * t195) * qJD(7)) * t336, t124, t23 * t69 - t241 * t28, t23 * t68 + t24 * t69 - t241 * t29 + t242 * t28, t188 * t28, t24 * t68 + t242 * t29, t188 * t29, 0, -t24 * t89 - t242 * t53 - t29 * t75 - t51 * t68, t23 * t89 - t241 * t53 + t28 * t75 + t51 * t69, -t28 * t59 - t29 * t60 + t45 * t68 - t46 * t69, t51 * t89 + t53 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t305 * t206, 0, t287, 0, 0, t206 * pkin(12) * t200, pkin(12) * t309, 0, 0, t317, t77, t64, -t317, t65, 0, (qJD(3) * t266 + t149 * t304) * pkin(1) + t269, (qJD(3) * t245 + t151 * t304) * pkin(1) + t270, -qJD(1) * t209 - t149 * t262 + (-t149 * qJD(3) + t106) * t294 + (t291 + t292) * t151, -t264, t365, t27, t25, -t365, t26, 0, t103 * t275 + t128 * t97 + t207, -t102 * t298 + t128 * t236 + (-t102 - t356) * qJD(4) + t369, t102 * t97 + t103 * t236 - t141 * t235 - t142 * t38 + t249, -t102 * t127 + t103 * t356 - t125 * t128 + t141 * t83 + t142 * t82, t7, t1, t5, t8, t6, t343, -t103 * t70 + t136 * t22 + (-t271 * t363 + t83) * t202 + t233 * t197 + t276, -t103 * t72 - t136 * t21 + t202 * t233 + t271 * t371 + t259, (-t102 * t70 - t137 * t22 + t42 * t72 + (t137 * t72 + t244) * qJD(5)) * t202 + (t102 * t72 - t137 * t21 + t42 * t70 - t4 + (t137 * t70 - t31) * qJD(5)) * t197 + t296, t102 * t247 - t103 * t122 - t136 * t83 + t137 * t218 - t248 * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t78, t62, -t316, 0, 0, -t315 + (-t147 * t304 + (-qJD(2) - t189) * t195 * qJD(7)) * pkin(1), t314 + (-t150 * t304 + (-t279 - t283) * qJD(7)) * pkin(1), (t352 * t104 + t355 * t147 + (t189 * t150 - t105) * t195) * pkin(1), -t264, t364, t12, t13, -t364, t14, 0, t188 * t49 + t242 * t81 + t220, -t48 * t188 + t241 * t81 + t217, -t23 * t79 + t24 * t80 + t241 * t49 + t242 * t48 + t252, t45 * t80 + t46 * t79 - t48 * t60 + t49 * t59 - t75 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t77, t64, -t317, t65, 0, -t266 * t336 + t269, -t245 * t336 + t270, 0, 0, t365, t27, t25, -t365, t26, 0, -t139 * t275 + (-t151 * t97 - t275 * t303) * pkin(4) + t207, -qJD(4) * t219 + t169 * t297 - t347 * t236 + t161 - t332 + (t140 - t263) * t275, -t139 * t236 - t140 * t97 + (-t353 * t235 - t198 * t38 + (-t198 * t236 + t353 * t97) * qJD(4)) * pkin(4) + t249, -t356 * t139 + t127 * t140 + (t353 * t83 + t125 * t151 + t198 * t82 + (-t127 * t353 - t198 * t356) * qJD(4)) * pkin(4), t7, t1, t5, t8, t6, t343, t184 * t22 - t35 * t363 + t253 * t70 + (-t183 * t318 + t83) * t202 + t221 * t197 + t276, -t184 * t21 + (t183 * t302 + t36) * t363 + t253 * t72 + t221 * t202 + t259, t35 * t72 + t36 * t70 + (-t70 * t263 - t183 * t22 + (t183 * t72 + t244) * qJD(5)) * t202 + (t72 * t263 - t183 * t21 - t4 + (t183 * t70 - t31) * qJD(5)) * t197 + t296, t122 * t139 - t83 * t184 + t244 * t35 - t31 * t36 + (t122 * t198 + t247 * t353) * qJD(4) * pkin(4) + t218 * t183, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365, t27, t25, -t365, t26, 0, -t127 * t275 + t207, t356 * t298 + t369, 0, 0, t7, t1, t5, t328 * t363 - t20, t6, t343, -pkin(8) * t22 - pkin(10) * t339 + t127 * t70 - t197 * t333 + t83 * t202 - t363 * t39 + t276, pkin(8) * t21 + pkin(10) * t239 - t122 * t366 + t127 * t72 + t363 * t40 + t259, t39 * t72 + t40 * t70 + t2 + (t370 + (-t22 + t319) * pkin(10)) * t202 + ((qJD(5) * t70 - t21) * pkin(10) + t368) * t197, pkin(8) * t83 + pkin(10) * t218 + t122 * t127 + t244 * t39 - t31 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t344, -t70 ^ 2 + t72 ^ 2, t363 * t70 - t21, -t344, t363 * t72 - t22, t38, -t122 * t72 - t368, t122 * t70 - t3 - t370, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, t306 * t206, 0, t288, 0, 0, -t206 * pkin(6) * t196, -pkin(6) * t310, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t78, t62, -t316, 0, 0, -t315 + (-qJD(7) + t189) * t293, t314 + (-t278 + t283) * t336, 0, 0, t364, t12, t13, -t364, t14, 0, -t188 * t322 + t242 * t88 + t220, t188 * t321 + t241 * t88 + t217, -t107 * t23 + t108 * t24 - t241 * t322 - t242 * t321 + t252, t107 * t46 + t108 * t45 + t321 * t60 - t322 * t59 - t75 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t364, t12, t13, -t364, t14, 0, -t188 * t60 + t220, t59 * t188 + t217, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t9;
