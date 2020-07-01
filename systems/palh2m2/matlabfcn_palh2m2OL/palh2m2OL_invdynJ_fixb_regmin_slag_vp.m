% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% palh2m2OL
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
% tau_reg [6x38]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m2OL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:38
% EndTime: 2020-06-30 18:07:11
% DurationCPUTime: 5.97s
% Computational Cost: add. (9322->406), mult. (22435->547), div. (0->0), fcn. (19272->18), ass. (0->236)
t182 = sin(qJ(6));
t188 = cos(qJ(6));
t183 = sin(qJ(5));
t185 = sin(qJ(3));
t190 = cos(qJ(3));
t191 = cos(qJ(2));
t286 = qJD(1) * t191;
t186 = sin(qJ(2));
t287 = qJD(1) * t186;
t122 = t185 * t287 - t190 * t286;
t124 = -t185 * t286 - t190 * t287;
t184 = sin(qJ(4));
t189 = cos(qJ(4));
t239 = -t122 * t189 + t184 * t124;
t328 = cos(qJ(5));
t340 = t184 * t122 + t189 * t124;
t53 = t183 * t340 + t239 * t328;
t382 = qJD(6) + t53;
t178 = qJD(2) + qJD(3);
t173 = qJD(4) + t178;
t159 = qJD(5) + t173;
t349 = t183 * t239 - t328 * t340;
t43 = t159 * t188 + t182 * t349;
t361 = t382 * t43;
t280 = qJD(6) * t182;
t307 = t188 * t349;
t45 = -t182 * t159 + t307;
t42 = t45 * t280;
t128 = t185 * t186 - t190 * t191;
t129 = t185 * t191 + t186 * t190;
t236 = t128 * t184 - t129 * t189;
t237 = t189 * t128 + t129 * t184;
t355 = t178 * qJD(1);
t195 = t340 * qJD(4) - t237 * qJDD(1) + t236 * t355;
t266 = qJD(5) * t328;
t281 = qJD(5) * t183;
t204 = t129 * t355;
t197 = -qJDD(1) * t128 - t204;
t225 = t129 * qJDD(1);
t96 = t178 * t128;
t198 = -t96 * qJD(1) + t225;
t282 = qJD(4) * t189;
t283 = qJD(4) * t184;
t33 = -t122 * t282 + t124 * t283 + t184 * t197 + t189 * t198;
t13 = t183 * t195 + t239 * t266 + t281 * t340 + t328 * t33;
t177 = qJDD(2) + qJDD(3);
t172 = qJDD(4) + t177;
t157 = qJDD(5) + t172;
t6 = -qJD(6) * t43 + t188 * t13 - t157 * t182;
t7 = t188 * (qJD(6) * t349 + t157) + t13 * t182 - t159 * t280;
t388 = (-t6 + t361) * t188 + t42 + (t45 * t53 + t7) * t182;
t14 = qJD(5) * t349 + t183 * t33 - t328 * t195;
t12 = -qJDD(6) + t14;
t312 = t182 * t12;
t386 = -t188 * t382 ^ 2 + t349 * t45 + t312;
t309 = t188 * t45;
t322 = t182 * t6;
t385 = -t309 * t382 - t322;
t384 = t349 * t53;
t327 = pkin(4) * t190;
t168 = qJDD(2) * t327;
t313 = pkin(4) * qJD(2);
t273 = t185 * t313;
t111 = pkin(2) * t177 - qJD(3) * t273 + t168;
t103 = t189 * t111;
t267 = t185 * t282;
t130 = pkin(2) * t178 + t190 * t313;
t268 = t130 * t283;
t276 = qJDD(2) * t185;
t285 = qJD(3) * t190;
t199 = t103 + (-t184 * t276 + (-t184 * t285 - t267) * qJD(2)) * pkin(4) - t268;
t56 = pkin(5) * t172 + t199;
t248 = t184 * t273;
t252 = -qJD(4) * t248 + t111 * t184;
t264 = qJD(2) * t285;
t220 = (t264 + t276) * pkin(4);
t374 = (qJD(4) * t130 + t220) * t189;
t63 = t252 + t374;
t101 = t184 * t130 + t189 * t273;
t270 = t328 * t101;
t100 = t189 * t130 - t248;
t95 = pkin(5) * t173 + t100;
t67 = t183 * t95 + t270;
t214 = -qJD(5) * t67 - t183 * t63 + t328 * t56;
t19 = t157 * pkin(3) + t214;
t181 = qJ(2) + qJ(3);
t176 = qJ(4) + t181;
t169 = qJ(5) + t176;
t155 = sin(t169);
t187 = sin(qJ(1));
t192 = cos(qJ(1));
t244 = g(1) * t192 + g(2) * t187;
t359 = t155 * t244;
t383 = t19 + t359;
t381 = t349 ^ 2 - t53 ^ 2;
t380 = -t159 * t53 + t13;
t156 = cos(t169);
t150 = g(3) * t155;
t23 = -t101 * t281 + t183 * t56 + t95 * t266 + t328 * t63;
t343 = -t23 + t150;
t163 = -t191 * pkin(4) - pkin(1);
t141 = qJD(1) * t163;
t99 = pkin(2) * t122 + t141;
t64 = -pkin(5) * t239 + t99;
t379 = t244 * t156 - t64 * t53 + t343;
t46 = t382 * t280;
t262 = t188 * t12 + t46;
t378 = t182 * t382 * t53 - t349 * t43 + t262;
t373 = t382 * t349;
t166 = sin(t176);
t167 = cos(t176);
t303 = t167 * t192;
t304 = t167 * t187;
t372 = g(1) * t303 + g(2) * t304 + g(3) * t166 - t239 * t99 - t252;
t371 = pkin(3) * t349;
t28 = -pkin(3) * t53 + t64;
t21 = -t182 * t28 + t188 * t67;
t324 = g(3) * t156;
t370 = t182 * t324 - t21 * t349;
t369 = t159 * t349 - t14;
t20 = -t182 * t67 - t188 * t28;
t368 = t383 * t188 + t20 * t349;
t367 = -t349 * t64 + t214 - t324 + t359;
t363 = pkin(5) * t340;
t358 = t372 - t374;
t27 = -t239 ^ 2 + t340 ^ 2;
t350 = -g(3) * t167 + t244 * t166 + t340 * t99;
t25 = -t173 * t340 + t195;
t161 = pkin(2) * t189 + pkin(5);
t269 = t328 * t184;
t116 = pkin(2) * t269 + t183 * t161;
t279 = t124 * pkin(2);
t68 = -t279 - t363;
t34 = t68 + t371;
t345 = (qJD(6) * t116 - t34) * t382;
t170 = pkin(4) * t287;
t162 = pkin(2) + t327;
t145 = t189 * t162;
t297 = t184 * t185;
t115 = -pkin(4) * t297 + pkin(5) + t145;
t296 = t185 * t189;
t120 = pkin(4) * t296 + t162 * t184;
t77 = t183 * t115 + t328 * t120;
t344 = (qJD(6) * t77 - t170 - t34) * t382;
t316 = t340 * t239;
t253 = t328 * t115 - t183 * t120;
t105 = pkin(2) * t128 + t163;
t24 = -t173 * t239 + t33;
t97 = t178 * t129;
t331 = qJD(4) * t236 + t184 * t96 - t189 * t97;
t321 = (-t363 + t371) * t382;
t235 = -t184 * t190 - t296;
t118 = t235 * t313;
t234 = t189 * t190 - t297;
t119 = t234 * t313;
t298 = t183 * t184;
t315 = t183 * t118 + t328 * t119 - t161 * t266 - (-t184 * t281 + (t328 * t189 - t298) * qJD(4)) * pkin(2);
t314 = -t328 * t118 + t183 * t119 - t161 * t281 + (-t184 * t266 + (-t183 * t189 - t269) * qJD(4)) * pkin(2);
t306 = qJD(6) * t382;
t305 = t124 * t122;
t302 = t182 * t187;
t301 = t182 * t192;
t300 = t183 * t101;
t295 = t187 * t188;
t294 = t188 * t192;
t179 = t186 ^ 2;
t289 = -t191 ^ 2 + t179;
t278 = qJD(1) * qJD(2);
t277 = qJDD(1) * t163;
t275 = t191 * qJDD(1);
t274 = qJDD(1) * pkin(1);
t171 = t186 * t313;
t271 = t188 * t306;
t265 = t186 * t278;
t158 = pkin(4) * t265;
t22 = -pkin(2) * t197 - pkin(4) * t275 - pkin(5) * t195 + t158 - t274;
t263 = t14 * pkin(3) + qJD(6) * t67 + t22;
t66 = t328 * t95 - t300;
t59 = pkin(3) * t159 + t66;
t259 = t382 * t59;
t258 = -qJD(6) * t28 + t23;
t83 = t97 * pkin(2) + t171;
t250 = qJD(2) * (-qJD(3) + t178);
t249 = qJD(3) * (-qJD(2) - t178);
t247 = -0.2e1 * pkin(1) * t278;
t246 = t382 * t266;
t245 = qJD(2) * t267;
t243 = g(1) * t187 - g(2) * t192;
t230 = t183 * t236 - t237 * t328;
t41 = -qJD(4) * t237 - t184 * t97 - t189 * t96;
t16 = t230 * qJD(5) + t183 * t331 + t328 * t41;
t58 = -t183 * t237 - t236 * t328;
t242 = t12 * t58 - t16 * t382;
t241 = -pkin(2) * t298 + t328 * t161;
t174 = sin(t181);
t175 = cos(t181);
t232 = g(3) * t174 + t122 * t141 + t244 * t175;
t229 = (-t59 + t66) * t382;
t70 = t328 * t100 - t300;
t228 = (-t59 + t70) * t382;
t194 = qJD(1) ^ 2;
t224 = pkin(1) * t194 + t244;
t222 = t243 + 0.2e1 * t274;
t221 = -g(3) * t175 + t124 * t141 + t244 * t174 + t168;
t71 = pkin(5) * t237 + t105;
t36 = -pkin(3) * t230 + t71;
t219 = t59 * t16 + t19 * t58 + t36 * t306;
t89 = t162 * t282 + (qJD(3) * t234 - t185 * t283) * pkin(4);
t90 = -t162 * t283 + (qJD(3) * t235 - t267) * pkin(4);
t38 = t253 * qJD(5) + t183 * t90 + t328 * t89;
t218 = t12 * t77 - t38 * t382 - t259;
t217 = t116 * t12 + t315 * t382 - t259;
t216 = t103 + t350;
t17 = t58 * qJD(5) + t183 * t41 - t328 * t331;
t35 = -pkin(5) * t331 + t83;
t212 = t59 * qJD(6) * t58 - (t17 * pkin(3) + t35) * t382 + t36 * t12 - t263 * t230;
t208 = (-pkin(2) * t173 - t130) * qJD(4) - t220;
t193 = qJD(2) ^ 2;
t160 = t328 * pkin(5) + pkin(3);
t117 = t158 + t277;
t112 = pkin(3) + t241;
t110 = t156 * t294 - t302;
t109 = -t156 * t301 - t295;
t108 = -t156 * t295 - t301;
t107 = t156 * t302 - t294;
t102 = t170 - t279;
t75 = pkin(3) + t253;
t74 = -t122 ^ 2 + t124 ^ 2;
t69 = -t183 * t100 - t270;
t65 = t170 + t68;
t62 = pkin(2) * t204 + t105 * qJDD(1) + t158;
t61 = -t124 * t178 + t197;
t60 = t122 * t178 + t198;
t39 = -t77 * qJD(5) - t183 * t89 + t328 * t90;
t26 = t28 * t280;
t1 = [qJDD(1), t243, t244, qJDD(1) * t179 + 0.2e1 * t191 * t265, 0.2e1 * t186 * t275 - 0.2e1 * t278 * t289, qJDD(2) * t186 + t191 * t193, qJDD(2) * t191 - t186 * t193, 0, t186 * t247 + t191 * t222, -t186 * t222 + t191 * t247, t124 * t96 + t198 * t129, -0.2e1 * t128 * t225 + t96 * t122 + t124 * t97 + (t128 ^ 2 - t129 ^ 2) * t355, t129 * t177 - t178 * t96, -t128 * t177 - t178 * t97, 0, t122 * t171 + 0.2e1 * t141 * t97 + t175 * t243 + (t117 + t277) * t128, t117 * t129 - t124 * t171 - t141 * t96 + t163 * t198 - t174 * t243, -t236 * t33 - t340 * t41, -t195 * t236 - t33 * t237 + t41 * t239 - t331 * t340, -t172 * t236 + t173 * t41, -t172 * t237 + t173 * t331, 0, g(1) * t304 - g(2) * t303 - t105 * t195 + t62 * t237 - t83 * t239 - t99 * t331, t105 * t33 - t166 * t243 - t236 * t62 - t340 * t83 + t41 * t99, t13 * t58 + t16 * t349, t13 * t230 - t14 * t58 + t16 * t53 - t17 * t349, t157 * t58 + t159 * t16, t157 * t230 - t159 * t17, 0, t14 * t71 + t156 * t243 + t17 * t64 - t22 * t230 - t35 * t53, t13 * t71 - t155 * t243 + t16 * t64 + t22 * t58 + t349 * t35, -t58 * t42 + (t16 * t45 + t58 * t6) * t188, (-t182 * t45 - t188 * t43) * t16 + (-t322 - t188 * t7 + (t182 * t43 - t309) * qJD(6)) * t58, -t17 * t45 - t188 * t242 + t230 * t6 - t46 * t58, t17 * t43 + t182 * t242 - t230 * t7 - t271 * t58, -t12 * t230 - t17 * t382, -g(1) * t108 - g(2) * t110 - t20 * t17 + t26 * t230 + (-t23 * t230 + t219) * t182 + t212 * t188, -g(1) * t107 - g(2) * t109 + t21 * t17 + (-t230 * t258 + t219) * t188 - t212 * t182; 0, 0, 0, -t186 * t194 * t191, t289 * t194, t186 * qJDD(1), t275, qJDD(2), -g(3) * t191 + t186 * t224, g(3) * t186 + t191 * t224, -t305, t74, t60, t61, t177, (-t122 * t287 + t177 * t190 + t185 * t249) * pkin(4) + t221, (t124 * t287 + (-qJDD(2) - t177) * t185 + t190 * t249) * pkin(4) + t232, t316, t27, t24, t25, t172, t90 * t173 + t145 * t172 - t268 + t102 * t239 + (-t245 + (-t264 + (-qJDD(2) - t172) * t185) * t184) * pkin(4) + t216, t102 * t340 - t120 * t172 - t173 * t89 + t358, -t384, t381, t380, t369, t157, t157 * t253 + t39 * t159 + t53 * t65 + t367, -t77 * t157 - t38 * t159 - t349 * t65 + t379, t385, t388, t386, t378, t373, t39 * t43 + t7 * t75 + (-t324 - t344) * t188 + t218 * t182 + t368, t39 * t45 + t6 * t75 + t218 * t188 + (-t383 + t344) * t182 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t74, t60, t61, t177, pkin(4) * t185 * t250 + t221, (t190 * t250 - t276) * pkin(4) + t232, t316, t27, t24, t25, t172, -pkin(4) * t245 - t118 * t173 + (-t124 * t239 + t189 * t172) * pkin(2) + t208 * t184 + t216, t119 * t173 + (-t124 * t340 - t172 * t184) * pkin(2) + t208 * t189 + t372, -t384, t381, t380, t369, t157, t241 * t157 + t314 * t159 + t53 * t68 + t367, -t116 * t157 + t315 * t159 - t349 * t68 + t379, t385, t388, t386, t378, t373, t112 * t7 + t314 * t43 + (-t324 - t345) * t188 + t217 * t182 + t368, t112 * t6 + t314 * t45 + t217 * t188 + (-t383 + t345) * t182 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t27, t24, t25, t172, t101 * t173 + t199 + t350, t100 * t173 + t358, -t384, t381, t380, t369, t157, -t69 * t159 + (t328 * t157 - t159 * t281 - t340 * t53) * pkin(5) + t367, t70 * t159 + (-t157 * t183 - t159 * t266 + t340 * t349) * pkin(5) + t379, t385, t388, t386, t378, t373, t160 * t7 - t69 * t43 + (t321 - t324) * t188 + t228 * t182 + (-t182 * t246 + (-qJD(5) * t43 - t271 + t312) * t183) * pkin(5) + t368, t160 * t6 - t69 * t45 + t228 * t188 + (-t383 - t321) * t182 + (-t188 * t246 + (-qJD(5) * t45 + t262) * t183) * pkin(5) + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, t381, t380, t369, t157, t67 * t159 + t367, t66 * t159 + t379, t385, t388, t386, t378, t373, -t188 * t324 + t67 * t43 + (t307 * t382 + t7) * pkin(3) + t229 * t182 + t368, pkin(3) * t6 + t67 * t45 + t229 * t188 + (-pkin(3) * t373 - t383) * t182 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t6 + t361, t382 * t45 - t7, -t12, -g(1) * t109 + g(2) * t107 + t343 * t182 - t263 * t188 + t21 * t382 - t45 * t59 + t26, g(1) * t110 - g(2) * t108 + t20 * t382 + t43 * t59 + t263 * t182 + (-t258 + t150) * t188;];
tau_reg = t1;
