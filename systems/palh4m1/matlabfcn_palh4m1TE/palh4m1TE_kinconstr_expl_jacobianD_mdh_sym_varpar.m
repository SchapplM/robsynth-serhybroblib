% Jacobian time derivative of explicit kinematic constraints of
% palh4m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% WD [9x5]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = palh4m1TE_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh4m1TE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 21:35:20
% EndTime: 2020-04-11 21:35:23
% DurationCPUTime: 3.28s
% Computational Cost: add. (13254->188), mult. (13357->424), div. (411->20), fcn. (3321->4), ass. (0->191)
t188 = cos(qJ(5));
t312 = pkin(2) * t188;
t269 = pkin(1) * t312;
t242 = qJD(5) * t269;
t340 = 0.2e1 * t242;
t187 = sin(qJ(5));
t313 = pkin(2) * t187;
t268 = pkin(1) * t313;
t174 = -0.2e1 * t268;
t186 = qJ(2) + pkin(6);
t181 = t186 ^ 2;
t189 = pkin(3) ^ 2;
t191 = pkin(2) ^ 2;
t192 = pkin(1) ^ 2;
t278 = t191 + t192;
t256 = -t189 + t278;
t163 = t174 + t181 + t256;
t172 = pkin(1) * t187 - pkin(2);
t180 = -pkin(3) - t186;
t175 = pkin(2) - t180;
t178 = pkin(2) + t180;
t279 = t174 + t192;
t157 = t175 * t178 + t279;
t179 = pkin(3) - t186;
t176 = pkin(2) - t179;
t177 = pkin(2) + t179;
t158 = t176 * t177 + t279;
t292 = t157 * t158;
t193 = sqrt(-t292);
t283 = t188 * t193;
t140 = -pkin(1) * t283 - t163 * t172;
t133 = 0.1e1 / t140 ^ 2;
t142 = pkin(1) * t188 * t163 - t172 * t193;
t138 = t142 ^ 2;
t121 = t133 * t138 + 0.1e1;
t117 = 0.1e1 / t121;
t170 = t174 + t278;
t339 = (-qJD(2) * t170 + t186 * t340) * t117;
t171 = -0.2e1 * t242;
t155 = t171 + (-t176 + t177) * qJD(2);
t156 = t171 + (-t175 + t178) * qJD(2);
t131 = -t155 * t157 - t156 * t158;
t152 = 0.1e1 / t193;
t294 = t152 * t188;
t246 = -t294 / 0.2e1;
t271 = qJD(5) * t193;
t251 = t187 * t271;
t338 = t131 * t246 + t251;
t272 = qJD(5) * t188;
t250 = t163 * t272;
t282 = t338 * pkin(1);
t275 = qJD(2) * t186;
t166 = -t242 + t275;
t327 = -0.2e1 * t166;
t105 = -pkin(1) * t250 + t172 * t327 + t282;
t132 = 0.1e1 / t140;
t287 = t170 * t186;
t258 = t133 * t287;
t233 = t117 * t258;
t259 = t132 * t287;
t218 = -t187 * t163 - t283;
t296 = t152 * t172;
t248 = -t296 / 0.2e1;
t322 = 0.2e1 * t188;
t106 = t131 * t248 + (t218 * qJD(5) + t166 * t322) * pkin(1);
t302 = t133 * t142;
t306 = t105 * t132 * t133;
t334 = 0.1e1 / t121 ^ 2 * (t106 * t302 - t138 * t306);
t337 = -t105 * t233 - t339 * t132 - 0.2e1 * t259 * t334;
t232 = t142 * t258;
t336 = 0.2e1 * t117 * t142 * t287 * t306 - t106 * t233 + 0.2e1 * t232 * t334 + t339 * t302;
t162 = t181 - t256 + 0.2e1 * t268;
t159 = 0.1e1 / t162;
t160 = 0.1e1 / t162 ^ 2;
t335 = 0.4e1 * t159 * t160;
t182 = 0.1e1 / t186;
t325 = 0.2e1 * t182;
t164 = -t181 + t189 + t170;
t173 = pkin(1) - t313;
t265 = pkin(2) * t283;
t141 = t164 * t173 - t265;
t136 = 0.1e1 / t141 ^ 2;
t165 = t242 + t275;
t295 = t152 * t173;
t247 = t295 / 0.2e1;
t217 = t187 * t164 + t283;
t330 = t217 * qJD(5);
t108 = t131 * t247 + (-0.2e1 * t188 * t165 - t330) * pkin(2);
t143 = t164 * t312 + t173 * t193;
t139 = t143 ^ 2;
t122 = t136 * t139 + 0.1e1;
t301 = t136 * t143;
t249 = t164 * t272;
t281 = t338 * pkin(2);
t326 = -0.2e1 * t173;
t107 = -pkin(2) * t249 + t165 * t326 + t281;
t135 = 0.1e1 / t141;
t305 = t107 * t135 * t136;
t307 = (t108 * t301 - t139 * t305) / t122 ^ 2;
t333 = t136 * t307;
t298 = t152 * t131;
t167 = 0.1e1 / t170;
t183 = 0.1e1 / t186 ^ 2;
t168 = 0.1e1 / t170 ^ 2;
t224 = t168 * t242;
t276 = qJD(2) * t182 / t181;
t262 = 0.2e1 * t276;
t332 = t167 * t262 - 0.2e1 * t183 * t224;
t277 = qJD(2) * t183;
t254 = t167 * t277;
t331 = t224 * t325 - t254;
t329 = 0.2e1 * pkin(2);
t328 = -0.2e1 * t140;
t185 = t188 ^ 2;
t324 = -0.2e1 * t185;
t323 = -0.2e1 * t186;
t321 = -0.2e1 * qJD(2);
t320 = pkin(1) * pkin(2);
t319 = t152 / 0.2e1;
t318 = -t167 / 0.2e1;
t317 = t167 / 0.2e1;
t316 = t182 / 0.2e1;
t315 = pkin(1) * t191;
t314 = pkin(2) * t168;
t311 = pkin(2) * t192;
t291 = t158 * t160;
t146 = -t157 * t291 + 0.1e1;
t144 = 0.1e1 / t146;
t310 = pkin(3) * t144;
t309 = pkin(3) * t170;
t308 = pkin(3) * t186;
t266 = t165 * t335;
t304 = (-t156 * t291 + (-t155 * t160 + t158 * t266) * t157) / t146 ^ 2;
t119 = 0.1e1 / t122;
t303 = t119 * t136;
t300 = t142 * t183;
t280 = t157 + t158;
t129 = ((t155 + t156) * t188 - t280 * t187 * qJD(5)) * t320;
t299 = t152 * t129;
t148 = t280 * t269;
t297 = t152 * t148;
t293 = 0.1e1 / t292 * t298;
t290 = t167 * t182;
t289 = t167 * t183;
t169 = t167 * t168;
t288 = t169 * t191;
t286 = t186 * t188;
t285 = t187 * t186;
t284 = t187 * t193;
t274 = qJD(5) * t185;
t273 = qJD(5) * t186;
t270 = t188 * qJD(2);
t147 = -0.2e1 * t157 * t179 - 0.2e1 * t158 * t180;
t226 = t147 * t246;
t126 = pkin(2) * t226 + t173 * t323;
t267 = -0.2e1 * t126 * t143;
t264 = t119 * t309;
t263 = t144 * t308;
t261 = -0.2e1 * t270;
t260 = t159 * t304;
t257 = t148 * t293;
t253 = qJD(5) * t297;
t252 = t192 * t274;
t245 = t293 / 0.4e1;
t244 = t187 * t319;
t243 = t168 * t269;
t241 = 0.2e1 * t135 * t307;
t240 = t136 * t264;
t238 = t165 * t263;
t237 = t160 * t193 * t308;
t235 = 0.4e1 * t187 * t272;
t234 = t117 * t259;
t231 = 0.2e1 * t243;
t229 = t188 * t253;
t227 = t147 * t245;
t221 = t117 * t232;
t219 = t252 * t288;
t216 = qJD(5) * t147 * t244;
t213 = pkin(3) * t119 * t340;
t212 = t172 * t329 - t163 - t297;
t130 = 0.2e1 * qJD(2) * t280 - 0.2e1 * t155 * t180 - 0.2e1 * t156 * t179;
t206 = t130 * t319 + t227;
t205 = -t257 / 0.4e1 - t299 / 0.2e1;
t190 = 0.1e1 / pkin(3);
t127 = 0.2e1 * pkin(1) * t286 + t147 * t248;
t125 = pkin(1) * t226 + t172 * t323;
t115 = -t217 * pkin(2) + t148 * t295 + t315 * t324;
t114 = t218 * pkin(1) - t148 * t296 + t311 * t324;
t113 = (t284 + (pkin(1) * t326 - t164 - t297) * t188) * pkin(2);
t112 = (t212 * t188 + t284) * pkin(1);
t1 = [0, 0, 0, 0, 0; 0, ((-t206 * t172 + (0.2e1 * t270 + (t226 - 0.2e1 * t285) * qJD(5)) * pkin(1)) * t290 - t106 * t289 + t332 * t142 + t331 * t127) * t234 - ((t172 * t321 + (t216 + (-t206 - 0.2e1 * t273) * t188) * pkin(1)) * t290 - t105 * t289 + t332 * t140 + t331 * t125) * t221 + (t337 * (t127 * t182 - t300) + t336 * (t125 * t182 - t140 * t183)) * t167, 0, 0, (0.8e1 * t142 * t182 * t219 + ((t235 * t311 + t282) * t182 - t114 * t277 + t205 * t172 * t325) * t167 + ((t187 * t327 - t229 - t250) * t290 + (t261 * t300 + (t106 * t322 + 0.2e1 * (t114 * t188 - t142 * t187) * qJD(5)) * t182) * t314) * pkin(1)) * t234 - (-t112 * t254 + (0.8e1 * t140 * t288 + t167 * t329) * t182 * t252 + ((t298 * t317 + (-t212 * t167 + t314 * t328) * qJD(5)) * t182 * t187 + ((t271 - t257 / 0.2e1 - t299 + t327) * t290 + (t277 * t328 + (qJD(5) * t112 + t105) * t325) * t314) * t188) * pkin(1)) * t221 + t337 * (t114 * t167 + t142 * t231) * t182 + t336 * (t112 * t167 + t140 * t231) * t182; 0, 0, 0, 0, 0; 0, ((t193 * t262 + t182 * t227 + (t130 * t316 + (-t131 / 0.2e1 - t147 * qJD(2) / 0.2e1) * t183) * t152) * t159 * t263 - 0.2e1 * (t162 * t276 - t165 * t183) * t144 * t237 + (t159 * qJD(2) * t310 - 0.2e1 * t160 * t238 - t260 * t308) * (t152 * t147 * t316 - t183 * t193) + ((-qJD(2) * t193 - t186 * t298 / 0.2e1) * t160 * t310 + t193 * t238 * t335 + t237 * t304) * (-t162 * t183 + 0.2e1)) * t190, 0, 0, (0.2e1 * t283 * t304 + (-t131 * t294 + 0.2e1 * t251) * t144) * t160 * t320 - t260 * t297 + 0.2e1 * (t159 * t299 / 0.2e1 + pkin(1) * t265 * t266 + (-t152 * t160 * t165 + t159 * t245) * t148) * t144; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t126 * t108 * t303 + t267 * t333 + (t107 * t303 + t241) * (-0.2e1 * pkin(2) * t286 + t147 * t247) + (-(t206 * t173 + (t261 + (t226 + 0.2e1 * t285) * qJD(5)) * pkin(2)) * t135 + (t173 * t321 + (t216 + (-t206 + 0.2e1 * t273) * t188) * pkin(2)) * t301 + t267 * t305) * t119, 0, 0, 0.2e1 * (-((t173 * t257 / 0.2e1 + t129 * t295 + t235 * t315 + t281) * t317 + 0.4e1 * t143 * t219 + ((0.2e1 * t187 * t165 - t229 - t249) * t317 + (t108 * t188 + (t115 * t188 - t143 * t187) * qJD(5)) * t168 * pkin(1)) * pkin(2)) * t135 * t264 - ((-0.4e1 * t141 * t169 * t192 - pkin(1) * t167) * t191 * t274 + ((t131 * t244 + t330 + t187 * t253 + (t165 + t205) * t322) * t318 + (-t107 * t168 * t188 + (-t167 * t173 * t187 + (-t113 * t188 + t141 * t187) * t168) * qJD(5)) * pkin(1)) * pkin(2)) * t143 * t240 + (t107 * t240 + t135 * t213 + t241 * t309) * (t115 * t317 + t143 * t243) + (0.2e1 * (t305 * t119 + t333) * t143 * t309 - t108 * t240 + t213 * t301) * (t113 * t318 - t141 * t243)) * t190; 0, 0, 0, 0, 0;];
WD = t1;
