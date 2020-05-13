% Jacobian of explicit kinematic constraints of
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% W [16x4]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = palh1m1TE_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:08:07
% EndTime: 2020-04-12 20:08:12
% DurationCPUTime: 5.24s
% Computational Cost: add. (75588->188), mult. (114208->365), div. (4912->39), fcn. (72404->22), ass. (0->216)
t228 = sin(pkin(20));
t231 = cos(pkin(20));
t232 = sin(qJ(3));
t235 = cos(qJ(3));
t217 = t228 * t235 + t231 * t232;
t334 = t217 * pkin(6);
t295 = pkin(1) * t334;
t212 = 0.2e1 * t295;
t247 = pkin(6) ^ 2;
t254 = pkin(1) ^ 2;
t300 = t247 + t254;
t206 = t212 + t300;
t364 = 0.2e1 * t206;
t233 = sin(qJ(2));
t236 = cos(pkin(19));
t343 = sin(pkin(19));
t344 = cos(qJ(2));
t220 = t233 * t236 - t343 * t344;
t336 = pkin(7) * t220;
t294 = pkin(1) * t336;
t215 = -0.2e1 * t294;
t246 = pkin(7) ^ 2;
t301 = t246 + t254;
t209 = t215 + t301;
t244 = pkin(8) ^ 2;
t250 = pkin(3) ^ 2;
t203 = -t244 + t250 + t209;
t213 = pkin(1) - t336;
t302 = t215 + t254;
t356 = -pkin(8) - pkin(3);
t194 = (pkin(7) - t356) * (pkin(7) + t356) + t302;
t355 = -pkin(8) + pkin(3);
t195 = (pkin(7) - t355) * (pkin(7) + t355) + t302;
t319 = t195 * t194;
t257 = sqrt(-t319);
t221 = t233 * t343 + t236 * t344;
t335 = pkin(7) * t221;
t182 = t203 * t335 + t213 * t257;
t305 = t235 * t182;
t309 = t221 * t257;
t292 = pkin(7) * t309;
t181 = t203 * t213 - t292;
t308 = t232 * t181;
t207 = 0.1e1 / t209;
t251 = 0.1e1 / pkin(3);
t313 = t207 * t251;
t166 = (t308 / 0.2e1 + t305 / 0.2e1) * t313;
t306 = t235 * t181;
t307 = t232 * t182;
t167 = (-t306 / 0.2e1 + t307 / 0.2e1) * t313;
t225 = pkin(23) + pkin(22);
t223 = sin(t225);
t224 = cos(t225);
t151 = t166 * t224 - t167 * t223;
t339 = t151 * pkin(5);
t296 = pkin(4) * t339;
t150 = 0.2e1 * t296;
t248 = pkin(5) ^ 2;
t304 = t150 + t248;
t354 = -pkin(9) - pkin(11);
t138 = (pkin(4) - t354) * (pkin(4) + t354) + t304;
t353 = -pkin(9) + pkin(11);
t139 = (pkin(4) - t353) * (pkin(4) + t353) + t304;
t327 = t139 * t138;
t255 = sqrt(-t327);
t277 = t166 * t223 + t224 * t167;
t363 = t277 * t255;
t240 = pkin(11) ^ 2;
t242 = pkin(9) ^ 2;
t249 = pkin(4) ^ 2;
t299 = t248 + t249;
t289 = -t242 + t299;
t144 = t150 + t240 + t289;
t149 = pkin(4) * t151 + pkin(5);
t127 = -pkin(4) * t363 + t144 * t149;
t129 = pkin(4) * t144 * t277 + t149 * t255;
t230 = cos(pkin(21));
t147 = t150 + t299;
t145 = 0.1e1 / t147;
t241 = 0.1e1 / pkin(11);
t325 = t145 * t241;
t227 = sin(pkin(21));
t347 = t227 / 0.2e1;
t121 = (-t127 * t230 / 0.2e1 + t129 * t347) * t325;
t118 = 0.1e1 / t121;
t146 = 0.1e1 / t147 ^ 2;
t340 = pkin(5) * t146;
t297 = pkin(4) * t340;
t278 = t230 * t297;
t279 = t227 * t297;
t119 = 0.1e1 / t121 ^ 2;
t120 = (t127 * t347 + t129 * t230 / 0.2e1) * t325;
t331 = t119 * t120;
t362 = (t127 * t279 + t129 * t278) * t118 - (-t127 * t278 + t129 * t279) * t331;
t361 = -0.2e1 * pkin(1);
t360 = -0.2e1 * t221 ^ 2;
t359 = pkin(4) * pkin(5);
t358 = -pkin(2) - pkin(13);
t357 = -pkin(2) + pkin(13);
t352 = t145 / 0.2e1;
t310 = t220 * t257;
t298 = pkin(1) * t335;
t321 = 0.1e1 / t257 * (t194 + t195) * t298;
t158 = (t310 + (t213 * t361 - t203 - t321) * t221) * pkin(7);
t351 = -t158 / 0.2e1;
t160 = t213 * t321 + t246 * pkin(1) * t360 + (-t220 * t203 - t309) * pkin(7);
t350 = t160 / 0.2e1;
t204 = 0.1e1 / t206;
t349 = t204 / 0.2e1;
t226 = sin(pkin(23));
t348 = t226 / 0.2e1;
t346 = -t235 / 0.2e1;
t237 = cos(pkin(18));
t345 = t237 / 0.2e1;
t238 = pkin(13) ^ 2;
t252 = pkin(2) ^ 2;
t290 = -t252 + t300;
t199 = t212 + t238 + t290;
t211 = pkin(1) * t217 + pkin(6);
t303 = t212 + t247;
t192 = (pkin(1) - t358) * (pkin(1) + t358) + t303;
t193 = (pkin(1) - t357) * (pkin(1) + t357) + t303;
t320 = t193 * t192;
t256 = sqrt(-t320);
t218 = t228 * t232 - t231 * t235;
t341 = pkin(1) * t218;
t179 = t199 * t341 + t211 * t256;
t342 = pkin(1) * t179;
t338 = pkin(5) * t277;
t198 = -t238 + t252 + t206;
t210 = -pkin(1) - t334;
t333 = t218 * pkin(6);
t178 = t198 * t333 - t210 * t256;
t337 = pkin(6) * t178;
t332 = 0.1e1 / (t119 * t120 ^ 2 + 0.1e1) * t241;
t273 = (t138 + t139) * t359;
t130 = t277 * t273;
t135 = 0.1e1 / t255;
t330 = t130 * t135;
t280 = 0.1e1 / t209 ^ 2 * t298;
t136 = ((t160 * t346 + t232 * t351) * t207 + (-t305 - t308) * t280) * t251;
t137 = ((t158 * t346 + t232 * t350) * t207 + (-t306 + t307) * t280) * t251;
t131 = t136 * t224 + t137 * t223;
t329 = t131 * t255;
t328 = t135 * t149;
t326 = t145 * t230;
t322 = 0.1e1 / t256 * (t192 + t193) * pkin(1) * t333;
t318 = t204 * t218 ^ 2;
t205 = 0.1e1 / t206 ^ 2;
t317 = t205 * t218;
t229 = cos(pkin(23));
t316 = t207 * t229;
t234 = sin(pkin(18));
t315 = t207 * t234;
t245 = 0.1e1 / pkin(8);
t314 = t207 * t245;
t312 = t217 * t256;
t311 = t218 * t256;
t293 = pkin(6) * t311;
t291 = -t250 + t301;
t123 = 0.2e1 * t131 * t273;
t288 = -t123 * t135 / 0.2e1;
t287 = t145 * t347;
t286 = -t326 / 0.2e1;
t285 = t326 / 0.2e1;
t284 = t207 * t348;
t283 = t207 * t345;
t282 = -0.2e1 * t249 * t338;
t281 = -0.2e1 * t149 * pkin(5) - t144;
t143 = -t240 + t242 + t147;
t148 = -pkin(4) - t339;
t126 = -pkin(5) * t363 - t143 * t148;
t125 = 0.1e1 / t126 ^ 2;
t128 = t143 * t338 - t148 * t255;
t275 = 0.1e1 / (t125 * t128 ^ 2 + 0.1e1) * t147;
t142 = t240 - t289 - 0.2e1 * t296;
t141 = 0.1e1 / t142 ^ 2;
t272 = 0.2e1 * t141 * t255 * t359;
t202 = t215 + t244 + t291;
t214 = pkin(1) * t220 - pkin(7);
t180 = -pkin(1) * t309 - t202 * t214;
t267 = t180 * t280;
t266 = t181 * t280;
t265 = t182 * t280;
t183 = pkin(1) * t221 * t202 - t214 * t257;
t264 = t183 * t280;
t263 = 0.1e1 / t126 * t275;
t262 = -t151 * t255 - t277 * t330;
t261 = pkin(4) * (t126 * t146 + t145 * t148);
t260 = pkin(5) * t125 * t128 * t275;
t132 = -t136 * t223 + t137 * t224;
t259 = -t132 * t255 + t277 * t288;
t258 = pkin(4) * (-t145 * t248 * t277 + t128 * t340);
t201 = t244 - t291 + 0.2e1 * t294;
t200 = 0.1e1 / t201 ^ 2;
t197 = t238 - t290 - 0.2e1 * t295;
t196 = 0.1e1 / t197 ^ 2;
t177 = -pkin(1) * t311 + t199 * t211;
t176 = -t198 * t210 - t293;
t175 = 0.1e1 / t177 ^ 2;
t174 = 0.1e1 / t176 ^ 2;
t169 = (t183 * t345 + t180 * t234 / 0.2e1) * t314;
t168 = (t180 * t345 - t234 * t183 / 0.2e1) * t314;
t165 = 0.1e1 / t168 ^ 2;
t164 = (t229 * t182 / 0.2e1 + t181 * t348) * t313;
t163 = (-t229 * t181 / 0.2e1 + t182 * t348) * t313;
t162 = 0.1e1 / t163 ^ 2;
t159 = -t214 * t321 + t254 * pkin(7) * t360 + (-t220 * t202 - t309) * pkin(1);
t157 = (t310 + (0.2e1 * t214 * pkin(7) - t202 - t321) * t221) * pkin(1);
t140 = 0.1e1 / t142;
t133 = 0.1e1 / (-t141 * t327 + 0.1e1);
t117 = t130 * t328 + t277 * t282 + (t151 * t144 - t363) * pkin(4);
t116 = (t277 * t281 + t262) * pkin(4);
t114 = t123 * t328 / 0.2e1 + t131 * t282 + (t132 * t144 - t329) * pkin(4);
t113 = (t131 * t281 + t259) * pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, ((t113 * t287 + t114 * t285) * t118 - (t113 * t286 + t114 * t287) * t331 + t362 * t131) * t332, ((t116 * t287 + t117 * t285) * t118 - (t116 * t286 + t117 * t287) * t331 + t362 * t277) * t332, 0; 0, 0, 0, 1; 0, ((t159 * t283 + t237 * t264 + t157 * t315 / 0.2e1 + t234 * t267) / t168 - (t157 * t283 + t237 * t267 - t159 * t315 / 0.2e1 - t234 * t264) * t169 * t165) / (t165 * t169 ^ 2 + 0.1e1) * t245, 0, 0; 0, ((t158 * t284 + t226 * t266 + t229 * t265 + t316 * t350) / t163 - (t160 * t284 + t226 * t265 - t229 * t266 + t316 * t351) * t164 * t162) / (t162 * t164 ^ 2 + 0.1e1) * t251, 0, 0; 0, 0, (((-t210 * t322 + (t217 * t198 - t311) * pkin(6)) * t349 + (-t247 * t318 + t317 * t337) * pkin(1)) / t176 - ((-t312 + (-t198 - t322) * t218) * t349 + (t176 * t205 + t204 * t210) * t341) * t174 * t337) / (t174 * t178 ^ 2 + 0.1e1) * t364, 0; 0, 0, (0.1e1 / t197 * t322 + t196 * t293 * t361) / (-t196 * t320 + 0.1e1), 0; 0, 0.2e1 * ((t148 * t288 + (t132 * t143 - t329) * pkin(5)) * t352 + t131 * t258) * t263 - 0.2e1 * ((-t131 * t143 + t259) * t352 + t131 * t261) * t260, 0.2e1 * ((-t148 * t330 + (t151 * t143 - t363) * pkin(5)) * t352 + t277 * t258) * t263 - 0.2e1 * ((-t143 * t277 + t262) * t352 + t277 * t261) * t260, 0; 0, (-0.1e1 / t201 * t321 + 0.2e1 * pkin(1) * t200 * t292) / (-t200 * t319 + 0.1e1), 0, 0; 0, 0, (-((t211 * t322 + (t217 * t199 - t311) * pkin(1)) * t349 + (-t254 * t318 + t317 * t342) * pkin(6)) / t177 - (-(-t312 + (-t199 - t322) * t218) * t204 / 0.2e1 + (-t177 * t205 + t204 * t211) * t333) * t175 * t342) / (t175 * t179 ^ 2 + 0.1e1) * t364, 0; 0, (t131 * t272 + t140 * t288) * t133, (-t140 * t330 + t272 * t277) * t133, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
W = t1;
