% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbar1turnTE
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
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnTE_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_regmin_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:33
% EndTime: 2020-06-27 16:23:08
% DurationCPUTime: 7.52s
% Computational Cost: add. (71790->307), mult. (101353->727), div. (3060->15), fcn. (27231->6), ass. (0->276)
t110 = 0.1e1 / pkin(4);
t102 = sin(qJ(2));
t240 = qJD(2) * t102;
t209 = pkin(2) * t240;
t186 = pkin(1) * t209;
t114 = pkin(2) ^ 2;
t115 = pkin(1) ^ 2;
t104 = cos(qJ(2));
t288 = pkin(2) * t104;
t254 = -0.2e1 * pkin(1) * t288 + t115;
t96 = t114 + t254;
t94 = 0.1e1 / t96 ^ 2;
t164 = t94 * t186;
t93 = 0.1e1 / t96;
t311 = -t93 / 0.2e1;
t308 = -pkin(3) - pkin(4);
t89 = (pkin(2) - t308) * (pkin(2) + t308) + t254;
t307 = pkin(4) - pkin(3);
t90 = (pkin(2) - t307) * (pkin(2) + t307) + t254;
t183 = (-t89 - t90) * pkin(1) * pkin(2);
t79 = t102 * t183;
t78 = qJD(2) * t79;
t272 = t89 * t90;
t116 = sqrt(-t272);
t83 = 0.1e1 / t116;
t276 = t78 * t83;
t213 = t102 * t276;
t97 = pkin(1) - t288;
t233 = 0.2e1 * t97 * pkin(1);
t247 = t104 * t116;
t331 = pkin(4) ^ 2;
t241 = pkin(3) ^ 2 - t331;
t92 = t96 - t241;
t36 = (-t213 + (-t247 + (t92 + t233) * t102) * qJD(2)) * pkin(2);
t249 = t102 * t116;
t73 = -pkin(2) * t249 + t92 * t97;
t133 = t73 * t164 + t36 * t311;
t27 = t133 * t110;
t68 = 0.1e1 / t73 ^ 2;
t289 = pkin(2) * t102;
t87 = t92 * t289;
t74 = t116 * t97 + t87;
t280 = t68 * t74;
t310 = t93 / 0.2e1;
t100 = t102 ^ 2;
t252 = t100 * t114;
t212 = pkin(1) * t252;
t180 = qJD(2) * t212;
t196 = t116 * t240;
t239 = qJD(2) * t104;
t208 = pkin(2) * t239;
t268 = pkin(2) * t196 + t92 * t208;
t274 = t83 * t97;
t37 = t78 * t274 + 0.2e1 * t180 + t268;
t132 = -t74 * t164 + t37 * t310;
t29 = t132 * t110;
t67 = 0.1e1 / t73;
t151 = t27 * t280 + t29 * t67;
t291 = pkin(1) * t102;
t91 = t96 + t241;
t88 = t91 * t291;
t98 = pkin(1) * t104 - pkin(2);
t75 = -t116 * t98 + t88;
t259 = t104 * t75;
t72 = -pkin(1) * t249 - t91 * t98;
t263 = t102 * t72;
t334 = t259 + t263;
t113 = 0.1e1 / pkin(3);
t153 = -0.2e1 * t164;
t234 = -0.2e1 * pkin(2) * t98;
t198 = t91 + t234;
t35 = (-t213 + (t198 * t102 - t247) * qJD(2)) * pkin(1);
t28 = (t72 * t153 + t35 * t93) * t113;
t65 = 0.1e1 / t72 ^ 2;
t282 = t65 * t75;
t251 = t100 * t115;
t197 = qJD(2) * t251;
t179 = pkin(2) * t197;
t267 = (t239 * t91 + t196) * pkin(1);
t273 = t83 * t98;
t38 = -t78 * t273 + 0.2e1 * t179 + t267;
t30 = (t75 * t153 + t38 * t93) * t113;
t64 = 0.1e1 / t72;
t285 = t30 * t64;
t333 = -t28 * t282 + t285;
t275 = t83 * t79;
t39 = t88 + (-t247 + (t234 - t275) * t102) * pkin(1);
t205 = -t75 / 0.2e1 + t39 / 0.2e1;
t312 = -t72 / 0.2e1;
t43 = -t79 * t273 + 0.2e1 * pkin(2) * t251 + (t104 * t91 + t249) * pkin(1);
t206 = t312 - t43 / 0.2e1;
t301 = pkin(2) * t94;
t229 = pkin(1) * t301;
t250 = t102 * t104;
t326 = (t100 * t72 + t75 * t250) * t229;
t332 = (t205 * t102 - t206 * t104) * t93 - t326;
t330 = -0.4e1 * t73;
t329 = 0.6e1 * t102;
t290 = pkin(1) * t113;
t69 = t67 * t68;
t70 = t74 ^ 2;
t279 = t69 * t70;
t21 = -t36 * t279 + t37 * t280;
t51 = t68 * t70 + 0.1e1;
t48 = 0.1e1 / t51 ^ 2;
t328 = t21 * t48;
t138 = (t259 / 0.2e1 + t263 / 0.2e1) * t93;
t299 = pkin(4) * t96;
t47 = 0.1e1 / t51;
t224 = t47 * t299;
t189 = t68 * t224;
t327 = t110 * t74 * t189;
t232 = qJD(1) * qJD(2);
t194 = t102 * t232;
t173 = pkin(2) * t194;
t230 = t104 * qJDD(1);
t325 = pkin(2) * t230 - t173;
t107 = qJD(2) ^ 2;
t248 = t104 * t107;
t324 = qJDD(2) * t102 + t248;
t318 = 0.2e1 * t96;
t178 = 0.2e1 * t48 * t318;
t277 = t74 * t96;
t195 = 0.4e1 * t69 * t277;
t226 = pkin(1) * t289;
t199 = -0.4e1 * t226;
t227 = t68 * t318;
t258 = 0.4e1 * t83 / t272;
t185 = t78 * t79 * t258;
t157 = -t185 / 0.4e1;
t298 = -t102 / 0.2e1;
t184 = -0.4e1 * t114 * t251;
t76 = (t104 * t183 + t184) * qJD(2);
t128 = (t102 * t157 + 0.2e1 * (-t104 * t78 + t76 * t298) * t83) * t93;
t95 = t93 * t94;
t255 = t114 * t95;
t154 = t197 * t255;
t159 = t73 * t94 - t93 * t97;
t214 = t114 * t291;
t165 = 0.6e1 * t104 * t214;
t223 = t67 * t299;
t166 = t110 * t47 * t223;
t217 = t93 * t276;
t261 = t102 * t94;
t284 = t37 * t94;
t303 = pkin(1) * t94;
t321 = 0.4e1 * t74;
t40 = t87 + (-t247 + (t233 - t275) * t102) * pkin(2);
t44 = t79 * t274 + 0.2e1 * t212 + (t104 * t92 + t249) * pkin(2);
t305 = -0.2e1 * ((0.4e1 * t180 + t268) * t311 + t154 * t330 + (-t128 / 0.2e1 + (t36 * t261 + (t104 * t159 + t40 * t261) * qJD(2)) * pkin(1)) * pkin(2)) * t327 - 0.2e1 * ((t97 * t185 / 0.4e1 + t76 * t274 + qJD(2) * t165) * t310 + t154 * t321 + ((t217 / 0.2e1 - pkin(1) * t284) * t102 + ((t247 + (-t92 + t275) * t102) * t310 + (-t102 * t44 - t104 * t74) * t303) * qJD(2)) * pkin(2)) * t166;
t319 = -0.2e1 * t96;
t1 = (t151 * (-t40 * t279 + t44 * t280) * t178 + ((t199 * t67 + t227 * t40) * t29 + (t40 * t195 + (t199 * t74 + t44 * t319) * t68) * t27) * t47) * pkin(4) + t305;
t172 = 0.2e1 * t224;
t12 = t151 * t172;
t221 = pkin(2) * t261;
t193 = pkin(1) * t221;
t31 = (t73 * t193 + t40 * t311) * t110;
t33 = (-t74 * t193 + t44 * t310) * t110;
t149 = t31 * t280 + t33 * t67;
t13 = t149 * t172;
t167 = -0.4e1 * t186;
t3 = (t149 * t21 * t178 + ((t167 * t67 + t227 * t36) * t33 + (t36 * t195 + (t167 * t74 + t37 * t319) * t68) * t31) * t47) * pkin(4) + t305;
t323 = (-qJD(2) * t13 + t12) * t193 - (t3 / 0.2e1 - t1 / 0.2e1) * t93;
t320 = 0.2e1 * t75;
t316 = -t36 / 0.2e1;
t315 = t37 / 0.2e1;
t314 = t40 / 0.2e1;
t313 = -t44 / 0.2e1;
t309 = t94 / 0.4e1;
t174 = 0.2e1 * t226;
t66 = t64 * t65;
t228 = t66 * t320;
t71 = t75 ^ 2;
t281 = t66 * t71;
t52 = t65 * t71 + 0.1e1;
t50 = 0.1e1 / t52 ^ 2;
t283 = t50 * t96;
t147 = 0.8e1 * t154;
t235 = -0.2e1 * t301;
t264 = t102 * t43;
t287 = ((t115 * t208 * t329 + t157 * t98 - t76 * t273) * t93 + t75 * t147 + ((t38 * t235 + t217) * t102 + ((t247 + (-t91 + t275) * t102) * t93 + 0.2e1 * (-t259 - t264) * t301) * qJD(2)) * pkin(1)) * t113 * t64;
t49 = 0.1e1 / t52;
t271 = t93 * t98;
t139 = -0.2e1 * t72 * t94 - 0.2e1 * t271;
t300 = pkin(3) * t96;
t225 = t49 * t300;
t190 = t65 * t225;
t170 = t75 * t190;
t222 = -0.2e1 * t261;
t8 = (((0.4e1 * t179 + t267) * t93 + t72 * t147) * t113 + (t128 + (t35 * t222 + (t104 * t139 + t222 * t39) * qJD(2)) * pkin(2)) * t290) * t170;
t2 = -t8 + (-0.2e1 * t333 * (-t39 * t281 + t43 * t282) * t283 + (t333 * t174 + (-t30 * t39 * t65 + t287 + (t228 * t39 - t43 * t65) * t28) * t96) * t49) * pkin(3);
t177 = -0.2e1 * t193;
t32 = (t72 * t177 + t39 * t93) * t113;
t34 = (t75 * t177 + t43 * t93) * t113;
t148 = t32 * t282 - t34 * t64;
t22 = -t35 * t281 + t38 * t282;
t220 = t22 * t283;
t4 = -t8 + (0.2e1 * t148 * t220 + (-t148 * qJD(2) * t174 + (-t34 * t35 * t65 + t287 + (t228 * t35 - t38 * t65) * t32) * t96) * t49) * pkin(3);
t306 = -t2 + t4;
t304 = pkin(1) * t93;
t302 = pkin(2) * t93;
t297 = t102 / 0.2e1;
t296 = -t104 / 0.2e1;
t103 = sin(qJ(1));
t295 = g(1) * t103;
t105 = cos(qJ(1));
t294 = g(1) * t105;
t293 = g(2) * t103;
t292 = g(2) * t105;
t278 = t74 * t94;
t256 = t113 * t93;
t200 = t256 / 0.2e1;
t201 = -t256 / 0.2e1;
t260 = t104 * t72;
t262 = t102 * t75;
t270 = (t200 * t260 + t201 * t262) * t103;
t269 = t334 * t105 * t200;
t265 = pkin(2) * qJD(1);
t191 = t64 * t225;
t14 = -t170 * t32 + t191 * t34 + 0.1e1;
t257 = t113 * t14;
t253 = qJD(2) * t94;
t246 = t105 * t113;
t245 = t107 * t115;
t244 = t107 * t116;
t108 = qJD(1) ^ 2;
t111 = 0.1e1 / t331;
t243 = t108 * t111;
t242 = -t104 ^ 2 + t100;
t238 = qJDD(2) * t91;
t237 = qJDD(2) * t92;
t236 = qJDD(2) * t93;
t215 = qJD(1) * t304;
t211 = t102 * t265;
t210 = t104 * t265;
t207 = t74 * t310;
t203 = t107 * t252;
t202 = t108 * t310;
t192 = t95 * t226;
t188 = -qJDD(1) * t73 / 0.2e1;
t182 = pkin(2) * t201;
t181 = pkin(2) * t200;
t163 = t95 * t186;
t162 = t108 * t115 * t221;
t161 = t293 + t294;
t160 = -t292 + t295;
t158 = pkin(2) * t236 * t257;
t156 = qJD(2) * t182;
t155 = qJD(2) * t181;
t145 = t294 / 0.2e1 + t293 / 0.2e1;
t142 = t78 ^ 2 * t258 / 0.4e1 + t83 * (t107 * t184 + t183 * t324);
t141 = 0.2e1 * qJD(2) * t276 + qJDD(2) * t116;
t140 = t94 * t107 * t214 * t257;
t137 = (t262 / 0.2e1 - t260 / 0.2e1) * t93;
t45 = t113 * t138;
t46 = t113 * t137;
t134 = (-t100 * t75 + t72 * t250) * t229;
t131 = -t142 + t244;
t130 = -t107 * t92 + t141;
t129 = 0.2e1 * t115 * t94 * t173 + (-qJDD(1) * pkin(1) - t295 / 0.2e1 + t292 / 0.2e1) * t93;
t127 = t113 * t332;
t126 = t113 * (t134 + (-t206 * t102 - t205 * t104) * t93);
t125 = (t38 * t296 + t35 * t298) * t93 + (t137 + t326) * qJD(2);
t124 = (t35 * t296 + t38 * t297) * t93 + (t134 + t138) * qJD(2);
t42 = qJD(1) * t46;
t41 = t334 * qJD(1) * t201;
t20 = qJD(1) * t126;
t19 = qJD(1) * t127;
t18 = t124 * t113;
t17 = t125 * t113;
t16 = (t124 * qJD(1) + qJDD(1) * t137) * t113;
t15 = (qJD(1) * t125 - qJDD(1) * t138) * t113;
t11 = -t170 * t28 + t191 * t30 + qJD(2);
t6 = qJDD(2) + ((-t142 * t271 + (0.8e1 * t75 * t95 * t203 + (0.2e1 * qJDD(2) * t100 + t248 * t329) * t302) * t115) * t113 + (((t238 + t244) * t93 + t75 * t107 * t235) * t104 + ((-t107 * t91 + t141) * t93 + (-0.4e1 * t38 * qJD(2) - 0.2e1 * t75 * qJDD(2)) * t301) * t102) * t290) * t191 - ((0.8e1 * t72 * t255 + 0.4e1 * t302) * t113 * t100 * t245 + ((-t141 * t93 + (t198 * t93 + t72 * t235) * t107) * t104 + ((t131 + t238) * t93 + (qJDD(2) * t139 - 0.4e1 * t35 * t253) * pkin(2)) * t102) * t290) * t170 + (-t28 * t38 - t30 * t35) * t190 + (t65 * t50 * t22 + t66 * t49 * t35) * t28 * t300 * t320 + (0.2e1 * t186 * t333 * t49 - 0.2e1 * t220 * t285) * pkin(3);
t5 = -0.2e1 * ((t107 * t165 + t142 * t97) * t310 + (t95 * t245 * t321 + pkin(1) * t236) * t252 + (((t237 + t244) * t104 + t130 * t102) * t310 + (-t74 * t248 + (-0.2e1 * qJD(2) * t37 - qJDD(2) * t74) * t102) * t303) * pkin(2)) * t166 + 0.4e1 * t29 * t223 * t328 - 0.2e1 * ((t115 * t330 * t95 - 0.2e1 * t304) * t203 + ((t159 * t107 * pkin(1) + t130 * t310) * t104 + ((t131 + t237) * t311 + (qJDD(2) * t159 + 0.2e1 * t36 * t253) * pkin(1)) * t102) * pkin(2)) * t327 + 0.2e1 * (-t27 * t37 + t29 * t36) * t189 + 0.2e1 * (0.2e1 * (t69 * t47 * t36 + t68 * t328) * t27 * t277 - 0.2e1 * t151 * t47 * t186) * pkin(4);
t7 = [qJDD(1), t160, t161, qJDD(1) * t100 + 0.2e1 * t104 * t194, 0.2e1 * t102 * t230 - 0.2e1 * t242 * t232, t324, qJDD(2) * t104 - t102 * t107, 0, t160 * t104, -t160 * t102, -t15 * t45 + t17 * t41, t15 * t46 - t16 * t45 + t17 * t42 + t18 * t41, t11 * t17 - t45 * t6, t16 * t46 + t18 * t42, t11 * t18 + t46 * t6, 0, -g(1) * t270 + t16 * t288 + t18 * t210 - t42 * t209 + (-t292 + t325) * t46, -g(2) * t269 - t15 * t288 - t17 * t210 + t41 * t209 - (-t295 - t325) * t45, (t70 * qJDD(1) * t309 + (-t70 * t163 + t278 * t315) * qJD(1)) * t111, (t188 * t278 + (t278 * t316 + (-t284 / 0.2e1 + 0.2e1 * t74 * t163) * t73) * qJD(1)) * t111, (-t12 * t132 + t207 * t5) * t110, (t73 * t5 * t311 - t133 * t12) * t110, 0, (t129 * t73 - t215 * t36) * t110, (t129 * t74 - t215 * t37) * t110; 0, 0, 0, -t108 * t250, t242 * t108, t102 * qJDD(1), t230, qJDD(2), -g(3) * t104 + t102 * t161, g(3) * t102 + t104 * t161, t41 * t19, t19 * t42 - t20 * t41, t11 * t19 + t14 * t15 + t306 * t41, -t42 * t20, -t11 * t20 + t14 * t16 + t306 * t42, t306 * t11 + t14 * t6, t35 * t14 * t156 + t158 * t312 + t42 * t211 - t20 * t210 - g(1) * (((t39 * t296 + t264 / 0.2e1) * t93 + t134) * t246 + t269) - t126 * t293 + g(3) * t127 + (t39 * t155 + t35 * t182) * t11 + (t2 * t155 + t4 * t156 + t6 * t182 + t140) * t72, t38 * t14 * t155 - t41 * t211 - t19 * t210 - g(1) * t332 * t246 - g(2) * (((t104 * t43 / 0.2e1 + t39 * t297) * t93 - t326) * t113 * t103 + t270) - g(3) * t126 + (t43 * t156 + t38 * t181) * t11 + (t6 * t181 - t140 + t158 / 0.2e1 + t4 * t155 + t2 * t156) * t75, (-t44 * t278 / 0.4e1 + t70 * t192 / 0.2e1) * t243, (t40 * t278 / 0.4e1 + (-t74 * t192 + t44 * t309) * t73) * t243, (-qJDD(1) * t13 * t207 + ((-t12 * t313 - t13 * t315) * t93 - t323 * t74) * qJD(1)) * t110, (-t93 * t13 * t188 + ((-t12 * t314 - t13 * t316) * t93 + t323 * t73) * qJD(1)) * t110, -t13 * t5 - (-t1 + t3) * t12, (-t73 * t162 + (g(3) * t313 + t145 * t40) * t93 + (t40 * t202 + (g(3) * t74 - t161 * t73) * t221) * pkin(1)) * t110, (-t74 * t162 + (g(3) * t314 + t145 * t44) * t93 + (t44 * t202 + (-g(3) * t73 - t161 * t74) * t221) * pkin(1)) * t110;];
tau_reg = t7;
