% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% fourbar1turnOL
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = fourbar1turnOL_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:07
% EndTime: 2020-04-12 19:41:15
% DurationCPUTime: 2.42s
% Computational Cost: add. (1976->269), mult. (4667->423), div. (0->0), fcn. (3406->8), ass. (0->228)
t275 = sin(qJ(3));
t297 = qJDD(2) + qJDD(3);
t279 = cos(qJ(3));
t280 = cos(qJ(2));
t310 = t279 * t280;
t276 = sin(qJ(2));
t322 = t275 * t276;
t210 = (-t310 + t322) * qJD(1);
t321 = t275 * t280;
t211 = (-t276 * t279 - t321) * qJD(1);
t331 = t210 * t211;
t335 = t297 + t331;
t339 = t275 * t335;
t338 = t279 * t335;
t277 = sin(qJ(1));
t281 = cos(qJ(1));
t246 = t281 * g(1) + t277 * g(2);
t205 = t280 * g(3) - t276 * t246;
t284 = qJD(1) ^ 2;
t257 = t280 * t284 * t276;
t243 = qJDD(2) + t257;
t173 = t243 * pkin(2) - t205;
t206 = -t276 * g(3) - t280 * t246;
t273 = t280 ^ 2;
t265 = t273 * t284;
t283 = qJD(2) ^ 2;
t254 = -t265 - t283;
t178 = t254 * pkin(2) + t206;
t139 = -t279 * t173 + t275 * t178;
t140 = -t275 * t173 - t279 * t178;
t113 = t279 * t139 + t275 * t140;
t337 = t276 * t113;
t336 = t280 * t113;
t245 = t277 * g(1) - t281 * g(2);
t221 = t281 * t245;
t179 = -t277 * t246 + t221;
t262 = t280 * qJDD(1);
t303 = qJD(1) * qJD(2);
t293 = t276 * t303;
t231 = t262 - t293;
t334 = t231 - t293;
t269 = qJD(2) + qJD(3);
t202 = t269 * t210;
t291 = t280 * t303;
t300 = t276 * qJDD(1);
t228 = t291 + t300;
t285 = -t210 * qJD(3) + t279 * t228 + t275 * t231;
t128 = t285 - t202;
t208 = t210 ^ 2;
t209 = t211 ^ 2;
t267 = t269 ^ 2;
t224 = -t284 * pkin(1) - t246;
t274 = sin(qJ(4));
t278 = cos(qJ(4));
t200 = t278 * g(3) + t274 * t224;
t201 = -t274 * g(3) + t278 * t224;
t147 = t278 * t200 - t274 * t201;
t333 = pkin(1) * t147;
t143 = -t211 * qJD(3) + t275 * t228 - t279 * t231;
t330 = t269 * t211;
t125 = t143 + t330;
t129 = t202 + t285;
t102 = -t275 * t125 - t279 * t129;
t332 = pkin(2) * t102;
t329 = t269 * t275;
t328 = t269 * t279;
t270 = t274 ^ 2;
t327 = t270 * t284;
t271 = t276 ^ 2;
t326 = t271 * t284;
t255 = t278 * t284 * t274;
t241 = qJDD(4) + t255;
t325 = t274 * t241;
t242 = qJDD(4) - t255;
t324 = t274 * t242;
t165 = -t331 + t297;
t323 = t275 * t165;
t320 = t276 * t243;
t244 = qJDD(2) - t257;
t319 = t276 * t244;
t172 = t334 * pkin(2) + t245;
t318 = t277 * t172;
t223 = qJDD(1) * pkin(1) + t245;
t317 = t277 * t223;
t316 = t277 * t245;
t314 = t278 * t223;
t313 = t278 * t242;
t312 = t279 * t165;
t311 = t279 * t172;
t309 = t280 * t244;
t308 = t280 * t245;
t307 = t281 * t172;
t306 = t281 * t223;
t272 = t278 ^ 2;
t305 = t270 + t272;
t304 = t271 + t273;
t302 = qJD(1) * qJD(4);
t301 = t274 * qJDD(1);
t299 = t277 * qJDD(1);
t261 = t278 * qJDD(1);
t298 = t281 * qJDD(1);
t296 = t277 * t331;
t295 = t281 * t331;
t294 = t274 * t302;
t292 = t278 * t302;
t148 = t274 * t200 + t278 * t201;
t162 = t276 * t205 + t280 * t206;
t180 = -t281 * t246 - t316;
t290 = t277 * t255;
t289 = t277 * t257;
t288 = t281 * t255;
t287 = t281 * t257;
t236 = -t277 * t284 + t298;
t286 = -t277 * g(3) - pkin(5) * t236;
t114 = t275 * t139 - t279 * t140;
t161 = t280 * t205 - t276 * t206;
t282 = qJD(4) ^ 2;
t264 = t272 * t284;
t253 = t265 - t283;
t252 = -t264 - t282;
t251 = t264 - t282;
t250 = -t283 - t326;
t249 = t283 - t326;
t248 = -t282 - t327;
t247 = t282 - t327;
t240 = t265 - t326;
t239 = t265 + t326;
t238 = t264 - t327;
t237 = t264 + t327;
t235 = t281 * t284 + t299;
t234 = t304 * qJDD(1);
t233 = t305 * qJDD(1);
t232 = t262 - 0.2e1 * t293;
t230 = t261 - 0.2e1 * t294;
t229 = t261 - t294;
t227 = 0.2e1 * t291 + t300;
t226 = t292 + t301;
t225 = 0.2e1 * t292 + t301;
t219 = t280 * t243;
t218 = t278 * t241;
t217 = t304 * t303;
t216 = t305 * t302;
t207 = t281 * g(3) - pkin(5) * t235;
t198 = -t209 + t267;
t197 = t208 - t267;
t196 = t280 * t228 - t271 * t303;
t195 = t278 * t226 - t270 * t302;
t194 = -t276 * t231 - t273 * t303;
t193 = -t274 * t229 - t272 * t302;
t191 = -t209 - t267;
t190 = -t276 * t250 - t309;
t189 = -t276 * t249 + t219;
t188 = t280 * t254 - t320;
t187 = t280 * t253 - t319;
t186 = -t274 * t248 - t313;
t185 = -t274 * t247 + t218;
t184 = t278 * t252 - t325;
t183 = t278 * t251 - t324;
t182 = t278 * t248 - t324;
t181 = t274 * t252 + t218;
t177 = t281 * t234 - t277 * t239;
t176 = t281 * t233 - t277 * t237;
t175 = t277 * t234 + t281 * t239;
t174 = t277 * t233 + t281 * t237;
t171 = -t276 * t227 + t280 * t232;
t170 = -t274 * t225 + t278 * t230;
t167 = -t209 + t208;
t163 = -t267 - t208;
t160 = t281 * t190 + t277 * t227;
t159 = t281 * t188 - t277 * t232;
t158 = t281 * t186 + t277 * t225;
t157 = t281 * t184 - t277 * t230;
t156 = t277 * t190 - t281 * t227;
t155 = t277 * t188 + t281 * t232;
t154 = t277 * t186 - t281 * t225;
t153 = t277 * t184 + t281 * t230;
t152 = (-t210 * t279 - t211 * t275) * t269;
t151 = (-t210 * t275 + t211 * t279) * t269;
t150 = -pkin(1) * t182 + t201;
t149 = -pkin(1) * t181 + t200;
t145 = -t208 - t209;
t142 = t281 * t162 - t316;
t141 = t277 * t162 + t221;
t138 = -t279 * t197 + t323;
t137 = t275 * t198 - t338;
t136 = -t275 * t197 - t312;
t135 = -t279 * t198 - t339;
t134 = t275 * t191 + t312;
t133 = -t279 * t191 + t323;
t132 = t281 * t148 - t317;
t131 = t277 * t148 + t306;
t126 = t143 - t330;
t124 = t211 * t329 + t279 * t285;
t123 = -t211 * t328 + t275 * t285;
t122 = t275 * t143 + t210 * t328;
t121 = -t279 * t143 + t210 * t329;
t120 = -t279 * t163 + t339;
t119 = -t275 * t163 - t338;
t118 = -t276 * t151 + t280 * t152;
t117 = pkin(2) * t126 - t311;
t116 = pkin(2) * t128 + t275 * t172;
t115 = -pkin(2) * t133 + t140;
t111 = -t276 * t136 + t280 * t138;
t110 = -t276 * t135 + t280 * t137;
t109 = -t276 * t133 + t280 * t134;
t108 = -t276 * t116 + t172 * t310;
t107 = -t276 * t117 + t172 * t321;
t106 = -pkin(2) * t119 - t139;
t105 = -t279 * t126 - t128 * t275;
t104 = -t279 * t125 + t275 * t129;
t103 = -t275 * t126 + t128 * t279;
t101 = -t276 * t123 + t280 * t124;
t100 = -t276 * t121 + t280 * t122;
t99 = -t276 * t119 + t280 * t120;
t98 = -pkin(2) * t145 + t114;
t97 = t281 * t109 - t277 * t128;
t96 = t277 * t109 + t281 * t128;
t95 = -t277 * t126 + t281 * t99;
t94 = t281 * t126 + t277 * t99;
t93 = t280 * t114 + t337;
t92 = -t276 * t103 + t280 * t105;
t91 = -t276 * t102 + t280 * t104;
t90 = t281 * t93 - t318;
t89 = t277 * t93 + t307;
t88 = -t276 * t98 + t336;
t87 = t277 * t145 + t281 * t91;
t86 = -t281 * t145 + t277 * t91;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t235, -t236, 0, t180, 0, 0, 0, 0, 0, 0, t159, t160, t177, t142, 0, 0, 0, 0, 0, 0, t95, t97, t87, t90, 0, 0, 0, 0, 0, 0, t157, t158, t176, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t236, -t235, 0, t179, 0, 0, 0, 0, 0, 0, t155, t156, t175, t141, 0, 0, 0, 0, 0, 0, t94, t96, t86, t89, 0, 0, 0, 0, 0, 0, t153, t154, t174, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t276 * t254 + t219, t280 * t250 - t319, 0, -t161, 0, 0, 0, 0, 0, 0, t280 * t119 + t276 * t120, t280 * t133 + t276 * t134, t280 * t102 + t276 * t104, t276 * t114 - t336, 0, 0, 0, 0, 0, 0, t181, t182, 0, -t147; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t236, 0, -t235, 0, t286, -t207, -t179, -pkin(5) * t179, t281 * t196 - t289, t281 * t171 - t277 * t240, t281 * t189 + t276 * t299, t281 * t194 + t289, t281 * t187 + t262 * t277, t277 * qJDD(2) + t281 * t217, -pkin(5) * t155 - t277 * t205 - t221 * t276, -pkin(5) * t156 - t277 * t206 - t221 * t280, -pkin(5) * t175 + t281 * t161, -pkin(5) * t141, t281 * t101 - t296, -t277 * t167 + t281 * t92, t281 * t110 - t277 * t129, t281 * t100 + t296, t281 * t111 + t125 * t277, t281 * t118 + t277 * t297, -pkin(5) * t94 - t277 * t106 + t281 * t107, -pkin(5) * t96 + t281 * t108 - t277 * t115, -pkin(5) * t86 + t277 * t332 + t281 * t88, -pkin(5) * t89 + (-t113 * t277 - t276 * t307) * pkin(2), t281 * t195 - t290, t281 * t170 - t277 * t238, t281 * t185 + t274 * t299, t281 * t193 + t290, t281 * t183 + t261 * t277, t277 * qJDD(4) + t281 * t216, -pkin(5) * t153 - t277 * t149 - t274 * t306, -pkin(5) * t154 - t277 * t150 - t278 * t306, -pkin(5) * t174 + t281 * t147, -pkin(5) * t131 - t277 * t333; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t235, 0, t236, 0, t207, t286, t180, pkin(5) * t180, t277 * t196 + t287, t277 * t171 + t281 * t240, t277 * t189 - t276 * t298, t277 * t194 - t287, t277 * t187 - t280 * t298, -t281 * qJDD(2) + t277 * t217, pkin(5) * t159 + t281 * t205 - t276 * t316, pkin(5) * t160 + t281 * t206 - t277 * t308, pkin(5) * t177 + t277 * t161, pkin(5) * t142, t277 * t101 + t295, t281 * t167 + t277 * t92, t277 * t110 + t281 * t129, t277 * t100 - t295, t277 * t111 - t125 * t281, t277 * t118 - t281 * t297, pkin(5) * t95 + t281 * t106 + t277 * t107, pkin(5) * t97 + t277 * t108 + t281 * t115, pkin(5) * t87 + t277 * t88 - t281 * t332, pkin(5) * t90 + (t113 * t281 - t276 * t318) * pkin(2), t277 * t195 + t288, t277 * t170 + t281 * t238, t277 * t185 - t274 * t298, t277 * t193 - t288, t277 * t183 - t278 * t298, -t281 * qJDD(4) + t277 * t216, pkin(5) * t157 + t281 * t149 - t274 * t317, pkin(5) * t158 + t281 * t150 - t277 * t314, pkin(5) * t176 + t277 * t147, pkin(5) * t132 + t281 * t333; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t245, t246, 0, 0, (t228 + t291) * t276, t280 * t227 + t276 * t232, t280 * t249 + t320, t334 * t280, t276 * t253 + t309, 0, t308, -t276 * t245, t162, 0, t280 * t123 + t276 * t124, t280 * t103 + t276 * t105, t280 * t135 + t276 * t137, t280 * t121 + t276 * t122, t280 * t136 + t276 * t138, t280 * t151 + t276 * t152, t280 * t117 + t172 * t322, t280 * t116 + t276 * t311, t280 * t98 + t337, t280 * pkin(2) * t172, (t226 + t292) * t274, t278 * t225 + t274 * t230, t278 * t247 + t325, (t229 - t294) * t278, t274 * t251 + t313, 0, pkin(1) * t230 + t314, -pkin(1) * t225 - t274 * t223, pkin(1) * t237 + t148, pkin(1) * t223;];
tauB_reg = t1;
