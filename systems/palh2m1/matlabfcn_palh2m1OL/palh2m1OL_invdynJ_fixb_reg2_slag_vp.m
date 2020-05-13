% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% palh2m1OL
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m1OL_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:27:46
% EndTime: 2020-05-03 00:28:22
% DurationCPUTime: 8.56s
% Computational Cost: add. (7357->464), mult. (17259->640), div. (0->0), fcn. (14349->10), ass. (0->237)
t168 = sin(qJ(5));
t173 = cos(qJ(5));
t174 = cos(qJ(3));
t175 = cos(qJ(2));
t148 = t174 * t175;
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t273 = qJD(1) * t171;
t127 = -qJD(1) * t148 + t170 * t273;
t280 = t170 * t175;
t132 = t171 * t174 + t280;
t128 = qJD(1) * t132;
t169 = sin(qJ(4));
t323 = cos(qJ(4));
t201 = -t169 * t127 + t323 * t128;
t151 = pkin(3) * t174 + pkin(2);
t281 = t170 * t171;
t116 = -pkin(3) * t281 + t151 * t175 + pkin(1);
t275 = qJD(1) * t116;
t80 = t323 * t127 + t128 * t169;
t37 = -pkin(4) * t80 + pkin(6) * t201 + t275;
t322 = pkin(2) * t174;
t152 = pkin(3) + t322;
t247 = t323 * t170;
t124 = pkin(2) * t247 + t152 * t169;
t321 = pkin(3) * t169;
t102 = t124 * qJD(2) + qJD(3) * t321;
t165 = qJD(2) + qJD(3);
t236 = qJD(4) + t165;
t97 = t236 * pkin(6) + t102;
t29 = -t168 * t97 + t173 * t37;
t30 = t168 * t37 + t173 * t97;
t301 = t168 * t30;
t340 = t173 * t80;
t184 = t165 * t132;
t326 = t148 - t281;
t196 = t326 * qJDD(1);
t180 = t184 * qJD(1) - t196;
t242 = qJD(4) * t323;
t269 = qJD(4) * t169;
t183 = t165 * t326;
t182 = qJD(1) * t183;
t262 = t132 * qJDD(1);
t54 = -t182 - t262;
t198 = t127 * t242 + t128 * t269 + t169 * t180 + t323 * t54;
t264 = qJDD(1) * t116;
t238 = t169 * t54 - t323 * t180;
t28 = -qJD(4) * t201 + t238;
t197 = qJD(3) * t132;
t270 = qJD(2) * t171;
t86 = -t151 * t270 + (-qJD(2) * t280 - t197) * pkin(3);
t291 = qJD(1) * t86;
t10 = pkin(4) * t28 - pkin(6) * t198 + t264 + t291;
t164 = qJDD(2) + qJDD(3);
t224 = qJDD(4) + t164;
t254 = t323 * pkin(3);
t222 = qJD(3) * t254;
t246 = t323 * t174;
t282 = t169 * t170;
t199 = t246 - t282;
t87 = t152 * t242 + (qJD(3) * t199 - t170 * t269) * pkin(2);
t47 = t87 * qJD(2) + qJD(4) * t222 + t124 * qJDD(2) + qJDD(3) * t321;
t43 = t224 * pkin(6) + t47;
t3 = qJD(5) * t29 + t168 * t10 + t173 * t43;
t2 = t3 * t173;
t172 = sin(qJ(1));
t318 = g(2) * t172;
t176 = cos(qJ(1));
t319 = g(1) * t176;
t137 = t318 + t319;
t255 = g(3) * t323;
t204 = t169 * t137 + t255;
t160 = t169 * g(3);
t213 = t137 * t323 - t160;
t343 = t170 * t204 - t213 * t174;
t344 = t170 * t213 + t204 * t174;
t350 = t171 * t344 + t343 * t175;
t354 = t2 + t350;
t359 = t29 * t340 + t80 * t301 + t354;
t349 = -t171 * t343 + t344 * t175;
t194 = -t47 - t350;
t149 = pkin(6) + t321;
t337 = qJD(5) - t80;
t289 = qJD(5) * t337;
t358 = t149 * t289 - t349;
t125 = -pkin(2) * t282 + t152 * t323;
t200 = t169 * t174 + t247;
t88 = t152 * t269 + (t200 * qJD(3) + t170 * t242) * pkin(2);
t251 = t88 * qJD(2) - t125 * qJDD(2) - qJDD(3) * t254;
t252 = pkin(3) * t269;
t192 = qJD(3) * t252 + t251;
t44 = -t224 * pkin(4) + t192;
t357 = t349 - t44;
t355 = -qJD(5) * t37 - t350 - t43;
t267 = qJD(5) * t173;
t353 = t267 - t340;
t352 = t349 - t192;
t153 = pkin(2) * t175 + pkin(1);
t274 = qJD(1) * t153;
t208 = t173 * t236;
t268 = qJD(5) * t168;
t16 = -qJD(5) * t208 - t168 * t224 - t173 * t198 - t201 * t268;
t57 = t168 * t236 - t173 * t201;
t290 = qJD(5) * t57;
t17 = t168 * t198 - t173 * t224 + t290;
t346 = t168 * t337;
t55 = -t168 * t201 - t208;
t1 = -t16 * t173 - t168 * t17 - t346 * t57 - t353 * t55;
t26 = qJDD(5) + t28;
t24 = t173 * t26;
t5 = -t201 * t55 - t337 * t346 + t24;
t348 = t29 * t337;
t347 = t30 * t337;
t13 = t16 * t168;
t7 = t353 * t57 - t13;
t284 = t137 * t170;
t230 = g(3) * t174 + t284;
t316 = g(3) * t170;
t325 = t137 * t174 - t316;
t341 = t325 * t171 + t175 * t230;
t6 = t168 * t26 + t201 * t57 + t353 * t337;
t104 = t125 * qJD(2) + t222;
t98 = -t236 * pkin(4) - t104;
t309 = t98 * t80;
t313 = t337 * t201;
t339 = t80 * t201;
t286 = t128 * t165;
t22 = t201 ^ 2 - t80 ^ 2;
t163 = pkin(6) * t169;
t256 = pkin(4) * t323;
t217 = -t163 - t256;
t134 = pkin(3) - t217;
t259 = pkin(4) * t318;
t144 = g(3) * pkin(6) + t259;
t154 = t163 + pkin(3);
t72 = -pkin(4) * t160 + t134 * t319 + t144 * t323 + t154 * t318;
t143 = g(3) * pkin(4) - pkin(6) * t318;
t138 = -pkin(4) * t169 + pkin(6) * t323;
t257 = t138 * t319;
t74 = g(3) * t154 + t143 * t323 + t169 * t259 - t257;
t335 = t170 * t72 + t74 * t174;
t243 = t29 * t201 + t98 * t268;
t333 = -t170 * t74 + t174 * t72;
t18 = -t80 * t236 + t198;
t220 = t44 * t168 - t30 * t201 + t98 * t267;
t19 = -t201 * t165 - t238;
t331 = -pkin(4) * t201 - pkin(6) * t80;
t330 = qJDD(2) * t322 + t128 * t274 + t341;
t117 = pkin(6) + t124;
t123 = -pkin(3) * t280 - t151 * t171;
t38 = qJD(1) * t123 + t331;
t328 = (qJD(5) * t117 + t38) * t337;
t225 = t127 * t165;
t45 = t54 - t225;
t324 = -t168 * t29 + t173 * t30;
t317 = g(2) * t176;
t320 = g(1) * t172;
t215 = -t317 + t320;
t315 = g(3) * t175;
t314 = t57 * t55;
t304 = t116 * t201;
t303 = t116 * t80;
t300 = t168 * t55;
t15 = t17 * t173;
t178 = qJD(1) ^ 2;
t288 = t116 * t178;
t287 = t127 * t128;
t285 = t215 * t175;
t166 = t171 ^ 2;
t167 = t175 ^ 2;
t276 = t166 - t167;
t129 = t199 * pkin(2);
t272 = qJD(2) * t129;
t130 = t200 * pkin(2);
t271 = qJD(2) * t130;
t266 = 0.2e1 * pkin(1);
t265 = qJD(1) * qJD(2);
t263 = qJDD(1) * t153;
t261 = t171 * qJDD(1);
t260 = t175 * qJDD(1);
t253 = pkin(2) * t270;
t96 = -t323 * t132 - t169 * t326;
t250 = t96 * t268;
t249 = t96 * t267;
t248 = t171 * t178 * t175;
t245 = t170 ^ 2 + t174 ^ 2;
t241 = t171 * t265;
t240 = t175 * t265;
t239 = pkin(1) * t178 + t137;
t231 = -qJD(5) * t97 + t215;
t228 = t54 - t262;
t227 = qJD(2) * (-qJD(3) + t165);
t226 = qJD(3) * (-qJD(2) - t165);
t223 = pkin(3) * t242;
t219 = t171 * t240;
t202 = t323 * t326;
t35 = qJD(4) * t202 - t132 * t269 - t169 * t184 + t323 * t183;
t214 = t26 * t96 - t337 * t35;
t210 = -t102 * t201 + t104 * t80;
t209 = t173 * t29 + t301;
t100 = t199 * t215;
t99 = t200 * t215;
t207 = t171 * t100 + t99 * t175;
t203 = t252 - t271;
t195 = -t117 * t26 - t337 * t87 - t309;
t193 = -t127 * t274 - t171 * t230 + t175 * t325;
t36 = t96 * qJD(4) - t169 * t183 - t323 * t184;
t11 = t36 * pkin(4) + t35 * pkin(6) + t86;
t95 = -t132 * t169 + t202;
t40 = pkin(4) * t95 - pkin(6) * t96 + t116;
t89 = t99 * t171;
t190 = qJD(5) * t96 * t98 + t11 * t337 + t26 * t40 - t89;
t189 = -t40 * t289 - t35 * t98 + t44 * t96 + t137;
t188 = -0.2e1 * t286;
t187 = pkin(4) * t255 + t144 * t169 - t257;
t186 = -t149 * t26 - t223 * t337 - t309;
t9 = t173 * t10;
t4 = -qJD(5) * t30 - t168 * t43 + t9;
t185 = -qJD(5) * t209 - t4 * t168 + t2;
t177 = qJD(2) ^ 2;
t159 = pkin(1) * t320;
t150 = -t254 - pkin(4);
t131 = t132 * pkin(3);
t118 = -pkin(4) - t125;
t101 = t215 * t138;
t94 = t134 * t320 + (-t256 - t154) * t317;
t90 = t199 * t285;
t76 = pkin(6) * t255 - t143 * t169 - t217 * t319 + t256 * t318;
t58 = -t127 ^ 2 + t128 ^ 2;
t46 = t180 - t286;
t39 = -qJD(1) * t131 + t331;
t34 = t104 * t173 + t168 * t331;
t33 = -t104 * t168 + t173 * t331;
t32 = t168 * t39 + t173 * t272;
t31 = -t168 * t272 + t173 * t39;
t8 = t346 * t55 - t15;
t12 = [0, 0, 0, 0, 0, qJDD(1), t215, t137, 0, 0, qJDD(1) * t166 + 0.2e1 * t219, 0.2e1 * t171 * t260 - 0.2e1 * t276 * t265, -qJDD(2) * t171 - t175 * t177, qJDD(1) * t167 - 0.2e1 * t219, -qJDD(2) * t175 + t171 * t177, 0, t285 + (-t241 + t260) * t266, -t215 * t171 + (-t240 - t261) * t266, t137, t159 + (pkin(1) * qJDD(1) - t317) * pkin(1), t128 * t183 - t54 * t132, -t54 * t148 + ((t262 - t225) * t175 + t188 * t171) * t174 + (t188 * t175 + (t225 + t228) * t171) * t170, -t132 * t164 - t165 * t183, t127 * t184 - t180 * t326, -t164 * t326 + t165 * t184, 0, t127 * t253 + 0.2e1 * (-qJD(2) * t132 - t197) * t274 + (-t253 * qJD(1) + t215 + 0.2e1 * t263) * t326, -t215 * t280 + (0.2e1 * pkin(2) * qJD(2) * t128 - t174 * t215) * t171 + (-t182 + t228) * t153, ((t174 * t132 - t170 * t326) * qJDD(2) + ((-t170 * t132 - t174 * t326) * qJD(3) + t165 * t175 * t245) * qJD(2)) * pkin(2) + t137, (-0.2e1 * pkin(2) * t241 + t215 + t263) * t153, t198 * t96 + t201 * t35, -t198 * t95 + t201 * t36 - t28 * t96 - t35 * t80, t224 * t96 - t236 * t35, t28 * t95 - t36 * t80, -t224 * t95 - t236 * t36, 0, t90 - t89 + (qJD(1) * t95 - t80) * t86 + (qJD(1) * t36 + qJDD(1) * t95 + t28) * t116, (qJD(1) * t96 - t201) * t86 + (-qJD(1) * t35 + qJDD(1) * t96 + t198) * t116 - t207, -t102 * t36 + t104 * t35 + t192 * t96 - t47 * t95 + t137, (t264 + 0.2e1 * t291 + t215) * t116, -t57 * t250 + (-t16 * t96 - t35 * t57) * t173, (t168 * t57 + t173 * t55) * t35 + (t13 - t15 + (-t173 * t57 + t300) * qJD(5)) * t96, -t16 * t95 + t173 * t214 - t250 * t337 + t57 * t36, t55 * t249 + (t17 * t96 - t35 * t55) * t168, -t168 * t214 - t17 * t95 - t249 * t337 - t55 * t36, t26 * t95 + t337 * t36, t29 * t36 + t4 * t95 + (t190 + t90) * t173 + t189 * t168, -t3 * t95 - t30 * t36 + t189 * t173 + (-t100 * t175 - t190) * t168, (-t11 * t57 + t16 * t40 + t29 * t35 - t4 * t96 + (-t30 * t96 - t40 * t55) * qJD(5)) * t173 + (-t11 * t55 - t17 * t40 - t3 * t96 + t30 * t35 + (t29 * t96 + t40 * t57) * qJD(5)) * t168 + t207, (pkin(2) * t215 + t101 * t170 + t174 * t94) * t175 + (t101 * t174 - t170 * t94) * t171 + t159 - pkin(1) * t317 + t209 * t11 + (t324 * qJD(5) + t3 * t168 + t4 * t173) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248, t276 * t178, -t261, t248, -t260, qJDD(2), t171 * t239 + t315, -g(3) * t171 + t175 * t239, 0, 0, t287, t58, t45, -t287, t46, t164, (-t127 * t273 + t174 * t164 + t170 * t226) * pkin(2) + t330, (-t128 * t273 + (-qJDD(2) - t164) * t170 + t174 * t226) * pkin(2) + t193, (-t196 * t170 - t45 * t174) * pkin(2), (t315 + (t153 * t178 + t137) * t171 + t245 * qJDD(2) * pkin(2)) * pkin(2), t339, t22, t18, -t339, t19, t224, -t88 * t236 + t125 * t224 + (t123 * t80 + t304) * qJD(1) + t352, -t87 * t236 - t124 * t224 + (t123 * t201 - t303) * qJD(1) + t194, -t124 * t28 - t125 * t198 - t201 * t88 + t80 * t87 + t210, t47 * t124 + t102 * t87 - t192 * t125 - t104 * t88 - t123 * t288 - (pkin(3) * t316 - t137 * t151) * t171 + (pkin(3) * t284 + g(3) * t151) * t175, t7, t1, t6, t8, t5, t313, t118 * t17 + t88 * t55 + t195 * t168 + (-t328 + t357) * t173 + t243, -t118 * t16 + t88 * t57 + t195 * t173 + (t328 - t349) * t168 + t220, (-t117 * t17 + t38 * t57 - t55 * t87 + (t117 * t57 - t29) * qJD(5)) * t173 + (-t117 * t16 + t38 * t55 + t57 * t87 - t4 + (t117 * t55 - t30) * qJD(5)) * t168 + t359, t44 * t118 + t98 * t88 - (-pkin(2) * t137 - t333) * t171 + (g(3) * pkin(2) + t335) * t175 + (-t29 * t38 + t30 * t87) * t173 + (-t29 * t87 - t30 * t38) * t168 + t185 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, t58, t45, -t287, t46, t164, pkin(2) * t170 * t227 + t330, (-t170 * qJDD(2) + t174 * t227) * pkin(2) + t193, 0, 0, t339, t22, t18, -t339, t19, t224, t236 * t271 + (-t131 * t80 + t304) * qJD(1) + (t323 * t224 + (-qJD(2) - 0.2e1 * qJD(3) - qJD(4)) * t269) * pkin(3) - t251 + t349, t236 * t272 + (-t131 * t201 - t303) * qJD(1) + (-t169 * t224 - t236 * t242) * pkin(3) + t194, (-t129 * t80 + t130 * t201) * qJD(2) + (-t323 * t198 - t169 * t28 + (-t169 * t201 + t323 * t80) * qJD(4)) * pkin(3) + t210, t131 * t288 + (-t102 * t129 + t104 * t130) * qJD(2) + (t47 * t169 - t192 * t323 + (t102 * t323 - t104 * t169) * qJD(4) + t341) * pkin(3), t7, t1, t6, t8, t5, t313, t150 * t17 - t31 * t337 + t203 * t55 + (-t44 - t358) * t173 + t186 * t168 + t243, -t150 * t16 + t358 * t168 + t186 * t173 + t203 * t57 + t32 * t337 + t220, t31 * t57 + t32 * t55 + (-t55 * t223 - t149 * t17 + (t149 * t57 - t29) * qJD(5)) * t173 + (t57 * t223 - t149 * t16 - t4 + (t149 * t55 - t30) * qJD(5)) * t168 + t359, t44 * t150 - t30 * t32 - t29 * t31 - t98 * t271 + t335 * t175 + t333 * t171 + (t169 * t98 + t324 * t323) * qJD(4) * pkin(3) + t185 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, t22, t18, -t339, t19, t224, t102 * t236 + t201 * t275 + t352, t104 * t236 - t275 * t80 + t194, 0, 0, t7, t1, t6, t300 * t337 - t15, t5, t313, -pkin(4) * t17 - t102 * t55 - t33 * t337 + (-pkin(6) * t26 - t309) * t168 + (-pkin(6) * t289 + t357) * t173 + t243, pkin(4) * t16 + t34 * t337 - t102 * t57 - t98 * t340 - t349 * t168 + (t268 * t337 - t24) * pkin(6) + t220, t34 * t55 + t33 * t57 + (-t348 + (-t17 + t290) * pkin(6)) * t173 + (-t4 - t347 + (qJD(5) * t55 - t16) * pkin(6)) * t168 + t354, -t44 * pkin(4) - t30 * t34 - t29 * t33 - t98 * t102 + (t76 * t170 + t174 * t187) * t175 + (-t170 * t187 + t76 * t174) * t171 + ((t171 * t247 - t175 * t246) * t318 + t185) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, -t55 ^ 2 + t57 ^ 2, t337 * t55 - t16, -t314, t337 * t57 - t17, t26, t355 * t168 + t231 * t173 - t98 * t57 + t347 + t9, t348 + t98 * t55 + t355 * t173 + (-t10 - t231) * t168, 0, 0;];
tau_reg = t12;
