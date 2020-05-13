% Calculate potential energy for
% picker2Dm1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 19:54
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 09:11:00
% EndTime: 2020-05-10 09:11:19
% DurationCPUTime: 11.74s
% Computational Cost: add. (173294->526), mult. (465796->660), div. (7180->13), fcn. (123736->38), ass. (0->287)
t183 = pkin(1) ^ 2;
t181 = t183 ^ 2;
t186 = pkin(7) ^ 2;
t141 = cos(qJ(1));
t107 = t141 ^ 2;
t190 = pkin(3) ^ 2;
t342 = sin(qJ(2));
t279 = pkin(3) * t342;
t252 = pkin(7) * t279;
t244 = 0.2e1 * t252;
t228 = t186 + t244;
t201 = t342 ^ 2;
t318 = t190 * t201;
t287 = 0.2e1 * t318;
t217 = -t190 + t228 + t287;
t215 = t107 * t217;
t175 = pkin(4) ^ 2;
t223 = t183 + t228;
t221 = -t175 + t223;
t240 = t279 + pkin(7);
t225 = t240 * t141;
t138 = sin(qJ(1));
t140 = cos(qJ(2));
t328 = t140 * t138;
t283 = pkin(3) * t328;
t251 = pkin(1) * t283;
t347 = 0.2e1 * t190;
t300 = t186 + t347;
t260 = t175 - t300;
t301 = -t186 + t175;
t262 = t183 - t301;
t285 = pkin(1) * t328;
t286 = t342 * pkin(7);
t304 = t183 - t186;
t352 = 0.1e1 / pkin(3);
t334 = t352 / 0.2e1;
t359 = 0.1e1 / t334;
t213 = sqrt(-0.4e1 * pkin(1) * (0.2e1 * (-t285 + t286) * pkin(3) + t262) * t225 + 0.4e1 * t221 * t251 + 0.4e1 * t304 * t318 - 0.4e1 * t262 * t252 - t181 - (t186 - (t359 + pkin(4)) * pkin(4)) * (t186 + (t359 - pkin(4)) * pkin(4)) + (-0.4e1 * t215 + 0.2e1 * t260) * t183);
t176 = 0.1e1 / pkin(4);
t276 = t176 * t334;
t100 = t183 + t186;
t102 = -0.3e1 * t183 + t186;
t106 = t141 * t107;
t108 = 0.10e2 / 0.3e1 * t183;
t109 = -0.20e2 / 0.3e1 * t183;
t117 = -t175 / 0.6e1;
t118 = -t175 / 0.4e1;
t119 = -t175 / 0.3e1;
t120 = -t175 / 0.2e1;
t122 = -0.3e1 / 0.2e1 * t175;
t126 = 0.2e1 / 0.3e1 * t190;
t129 = -t190 / 0.3e1;
t131 = 0.4e1 / 0.3e1 * t183;
t133 = t183 / 0.2e1;
t143 = 0.15e2 * t181;
t145 = 0.10e2 * t183;
t150 = 0.18e2 * t186;
t151 = -0.2e1 * t175;
t152 = -0.5e1 * t175;
t153 = -0.6e1 * t175;
t192 = t190 ^ 2;
t154 = 0.5e1 * t192;
t157 = 0.7e1 * t181;
t158 = 0.5e1 * t181;
t159 = 0.6e1 * t183;
t160 = 0.2e1 * t183;
t185 = t186 ^ 2;
t162 = 0.3e1 * t185;
t163 = 0.8e1 * t186;
t164 = 0.6e1 * t186;
t165 = 0.4e1 * t186;
t166 = 0.3e1 * t186;
t167 = 0.2e1 * t186;
t174 = t175 ^ 2;
t191 = pkin(3) * t190;
t177 = t191 ^ 2;
t187 = pkin(1) * t183;
t196 = pkin(7) * t186;
t305 = t181 + t185;
t308 = t167 - t175;
t321 = t186 * t175;
t222 = t308 * t183 + t174 / 0.6e1 + t305 - t321;
t219 = 0.5e1 / 0.6e1 * t192 + t222;
t96 = pkin(3) * t140;
t284 = t183 * t96;
t220 = t107 * (-t138 * t187 + t284);
t307 = t174 - t192;
t224 = 0.6e1 * t185 + t307 - 0.6e1 * t321;
t227 = t186 - t251;
t242 = -0.4e1 * t251;
t298 = t190 + t186;
t259 = t183 + t298;
t84 = t120 + t259;
t229 = t84 * t242;
t317 = t174 / 0.2e1 - t192 / 0.2e1;
t239 = -0.3e1 * t321 + t162 + t317;
t243 = -0.6e1 * t251;
t312 = 0.15e2 * t183 + t166;
t331 = t177 + t100 * ((t122 + t167) * t183 - 0.3e1 / 0.2e1 * t321 + t305 + t317);
t234 = ((t108 + t308) * t190 + t219) * t243 + (t143 + (t150 - 0.9e1 * t175) * t183 + t239) * t190 + (t122 + t312) * t192 + t331;
t349 = 0.3e1 * t183;
t258 = t349 + t298;
t341 = pkin(1) * t138;
t236 = t258 * t96 - t341 * (0.3e1 * t190 + t100);
t121 = -0.2e1 / 0.3e1 * t175;
t130 = -0.2e1 / 0.3e1 * t190;
t310 = t151 - 0.2e1 * t190;
t263 = t164 + t310;
t264 = t121 + t100;
t311 = t145 + t167;
t316 = t130 + t186;
t237 = -t341 * (t154 + (t145 + t263) * t190 + (t130 + t264) * t100) + (t192 + (t121 + t130 + t311) * t190 + t158 + t263 * t183 + t186 * (t121 + t316)) * t96;
t241 = -0.2e1 * t251;
t319 = t187 * t106;
t246 = t319 * t96;
t323 = t181 * t107 ^ 2;
t247 = t323 * t96;
t97 = pkin(1) * t141;
t88 = t97 + pkin(7);
t248 = t88 * t279;
t265 = t121 + t126 + t167;
t330 = t192 + (t126 + t264) * t100;
t125 = 0.4e1 / 0.3e1 * t190;
t266 = t119 + t100;
t82 = t125 + t266;
t249 = t82 * t242 + t330 + (t159 + t265) * t190;
t278 = t342 * t201 * t191 * t88;
t250 = -0.8e1 * t278;
t280 = 0.16e2 * t319;
t255 = pkin(7) * t280;
t290 = pkin(7) * t319;
t257 = 0.8e1 * t290;
t299 = -t190 + t186;
t261 = -t175 + t299;
t127 = t190 / 0.3e1;
t267 = t117 + t127 + t186;
t268 = t175 / 0.3e1 + t127 + t167;
t269 = 0.2e1 / 0.3e1 * t175 + t126 + t165;
t270 = 0.4e1 / 0.3e1 * t175 + t125 - 0.2e1 * t186;
t326 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t273 = t138 * t326;
t292 = 0.6e1 * t97;
t274 = pkin(7) * t292;
t293 = 0.4e1 * t97;
t275 = pkin(7) * t293;
t277 = -t341 / 0.2e1;
t322 = t183 * t107;
t281 = 0.12e2 * t322;
t327 = t140 * t141;
t282 = pkin(3) * t327;
t288 = 0.4e1 * t322;
t289 = 0.8e1 * t323;
t295 = 0.2e1 * t341;
t296 = pkin(7) * t97;
t297 = 0.4e1 * pkin(7);
t302 = t185 + t192;
t303 = t185 - t181;
t309 = t152 - 0.5e1 * t190;
t313 = 0.4e1 / 0.7e1 * t186 - t175 / 0.7e1;
t314 = t133 + t186;
t315 = t183 / 0.3e1 + t186;
t320 = t186 * t183;
t325 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t329 = t201 ^ 2 * t192;
t95 = t100 ^ 2;
t338 = t95 * (t183 + t261);
t346 = 0.4e1 * t190;
t348 = -0.6e1 * t190;
t350 = 0.4e1 * t181;
t353 = t166 - t175 - t190;
t354 = t118 + t190 / 0.2e1;
t71 = -t192 / 0.6e1 + t222;
t93 = t129 + t186;
t74 = t93 * t241;
t78 = t97 + t240;
t99 = t186 - 0.3e1 * t190;
t81 = t99 * t257;
t85 = t353 * t145;
t86 = -t175 + t259;
t294 = 0.2e1 * t97;
t89 = pkin(7) * t294;
t91 = (t165 + t175) * t183;
t94 = -t183 / 0.3e1 + t186;
t98 = -0.30e2 * t175 + 0.60e2 * t186;
t340 = ((-0.24e2 * (0.4e1 / 0.3e1 * t322 + t89 + t94) * t329 * t341 + (0.4e1 * t220 + (t96 + t295) * t89 + t102 * t96 + (t120 + t258) * t295) * t250 - 0.12e2 * (-0.8e1 / 0.3e1 * t247 + ((t131 + t267) * t96 - (0.7e1 / 0.6e1 * t190 + t117 + t314) * t341) * t288 + (-t190 * t304 - 0.5e1 / 0.3e1 * t181 + t268 * t183 + t186 * (t119 + t93)) * t96 + (-t192 + (t109 + t269) * t190 - 0.3e1 * t181 + t270 * t183 + t185) * t277 + (-t138 * t181 * t106 + ((t183 + t267) * t96 + (-t183 + t300) * t277) * t97) * t297) * t318 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t190 + t133 + t117) * t140 * t359 + pkin(1) * t273) * t322 + (-0.8e1 * t246 + ((t119 + t126 + t349 + t186) * t96 - (0.8e1 / 0.3e1 * t190 + t266) * t341) * t293) * pkin(7) + t237) * t248 + 0.24e2 * t93 * t247 + ((t186 + 0.5e1 / 0.2e1 * t190 + 0.3e1 / 0.2e1 * t183 + t120) * t96 + t99 * t341 / 0.2e1) * t255 - 0.6e1 * ((-0.3e1 * t192 + (t109 + t270) * t190 + t269 * t183 + t303) * t96 - 0.2e1 * (-0.5e1 / 0.3e1 * t192 + (-t183 + t268) * t190 + t186 * (t129 + t266)) * t341) * t322 - 0.6e1 * t237 * t296 - (t177 + (0.21e2 * t183 + t353) * t192 + (t186 * t310 + t162 + 0.35e2 * t181 + t85) * t190 + (t157 + (t163 + t309) * t183 + t186 * t261) * t100) * t96 + (0.7e1 * t177 + (0.35e2 * t183 + 0.15e2 * t186 + t309) * t192 + (0.21e2 * t181 + t85 + 0.9e1 * t185 + (t153 + t348) * t186) * t190 + t338) * t341) * t213 + (0.16e2 * (t289 + t255 + (-0.8e1 * t181 + 0.12e2 * t320) * t107 + (0.4e1 * pkin(1) * t196 - 0.12e2 * pkin(7) * t187) * t141 - 0.6e1 * t320 + t305) * t329 + 0.32e2 * (t257 + (-0.4e1 * t187 * t283 + t350 + (t346 + t151 + t163) * t183) * t107 + (-t183 + t227 + t354) * t275 + t241 * t325 + t102 * t84) * t278 + 0.24e2 * (t316 * t289 + 0.14e2 * (-0.32e2 / 0.21e2 * (t186 + t190 / 0.4e1 + t183 / 0.4e1 - t175 / 0.8e1) * t251 + 0.5e1 / 0.42e2 * t192 + (0.16e2 / 0.21e2 * t183 + t313) * t190 + t181 / 0.7e1 + t313 * t183 + t185 - 0.3e1 / 0.7e1 * t321 + t174 / 0.42e2) * t322 + t94 * t229 - t304 * t192 + (-0.10e2 / 0.3e1 * t181 + 0.2e1 * t185 - t321 + t91) * t190 + t71 * t325 + ((-0.2e1 / 0.3e1 * t251 + t118 + t314) * t280 + (-0.8e1 / 0.3e1 * (t315 + t354) * t251 + 0.5e1 / 0.18e2 * t192 + (0.4e1 / 0.3e1 * t186 + t131 + t119) * t190 + t185 + 0.2e1 / 0.3e1 * t320 - 0.2e1 / 0.3e1 * t321 - t181 / 0.3e1 + t174 / 0.18e2) * t292) * pkin(7)) * t318 + 0.8e1 * (t81 + (t326 * t84 + t74) * t281 + (t229 + (t159 + t308) * t190 + t219) * t274 + t234) * t248 + 0.16e2 * (t186 * t348 + t302) * t323 + 0.32e2 * (t241 * t326 + t84 * t99) * t290 + 0.24e2 * (t93 * t229 - t177 + (-t108 + t301) * t192 + (t91 + t192 / 0.6e1 - t174 / 0.6e1 + t303) * t190 + t71 * t186) * t322 + 0.8e1 * t234 * t296 - 0.8e1 * ((t122 + t166 + 0.7e1 * t183) * t192 + (t157 + (t152 + 0.10e2 * t186) * t183 + t239) * t190 + t331) * t251 + t192 ^ 2 + (t151 + t165 + 0.28e2 * t183) * t177 + (t183 * t98 + 0.70e2 * t181 + t224) * t192 + (t185 * t153 + t224 * t159 + t307 * t167 + t98 * t181 + 0.28e2 * t187 ^ 2 + 0.4e1 * t196 ^ 2) * t190 + t86 * t338) * t78) / ((t250 * t341 - 0.4e1 * (0.2e1 * t220 - t304 * t96 + (0.2e1 * pkin(7) * t282 + t138 * (t160 + t190)) * pkin(1)) * t318 - 0.4e1 * (-0.2e1 * t107 * t284 + (t96 - t341) * t89 + t236) * t248 + 0.8e1 * pkin(7) * t246 + ((pkin(3) * t350 + 0.8e1 * t183 * t191) * t140 + 0.4e1 * t187 * t273) * t107 - 0.4e1 * t236 * t296 - (t190 * t311 + t158 + t302 + 0.6e1 * t320) * t96 + (t154 + (t145 + t164) * t190 + t95) * t341) * t213 + (0.8e1 * (t89 + t288 + t102) * t278 + 0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t251 + 0.4e1 / 0.9e1 * t190 - t175 / 0.9e1 + t315) * t322 + t94 * t241 + t82 * t325 + (t319 + (t117 + t126 + t227) * t97) * t297) * t318 + 0.6e1 * (t299 * t288 + (t82 + t241) * t275 + t249) * t248 + t81 + (t326 * t82 + t74) * t281 + t249 * t274 + ((t108 + t265) * t190 + t330) * t243 + t177 + (-t175 + t190 + t312) * t192 + (t143 + (t150 + t153 + 0.6e1 * t190) * t183 + t162 + (t151 + t347) * t186) * t190 + t95 * t86) * t78);
t27 = atan2(t213 * t276, t276 * t340);
t25 = sin(t27);
t26 = cos(t27);
t169 = 0.2e1 * pkin(1);
t216 = t225 - t283;
t214 = 0.1e1 / (t169 * t216 + t190 + t223);
t238 = t349 - t260;
t351 = -0.2e1 * pkin(1);
t209 = t214 * ((t138 * t240 + t282) * t213 + t215 * t351 - ((-0.4e1 * t285 + 0.2e1 * t286) * pkin(3) + t238) * t225 + (t244 + t238) * t283 + (t287 - t244 - t346 - t262) * pkin(1));
t208 = t352 * t209;
t218 = t347 + t221;
t212 = t214 * ((pkin(1) + t216) * t213 + (t217 * t294 + t218 * t240) * t138 + (t141 * t218 + (0.4e1 * t107 - 0.2e1) * pkin(1) * t240) * t96);
t206 = atan2(t212 * t334, t208 / 0.2e1);
t204 = sin(t206);
t205 = cos(t206);
t39 = -t138 * t205 - t141 * t204;
t40 = t138 * t204 - t141 * t205;
t12 = t25 * t40 + t26 * t39;
t256 = -t25 * t39 + t40 * t26;
t207 = t209 / 0.4e1;
t324 = t176 / t190;
t23 = (t207 * t340 - t213 * t212 / 0.4e1) * t324;
t211 = t212 / 0.4e1;
t24 = (t207 * t213 + t211 * t340) * t324;
t75 = -t138 * t342 - t327;
t76 = -t141 * t342 + t328;
t10 = atan2(t23 * t75 + t24 * t76, -t23 * t76 + t24 * t75);
t8 = sin(t10);
t9 = cos(t10);
t363 = t12 * t8 - t256 * t9;
t362 = t12 * t9 + t256 * t8;
t171 = 0.1e1 / pkin(6);
t335 = 0.1e1 / pkin(2) / 0.2e1;
t168 = 0.1e1 / t335;
t170 = pkin(6) ^ 2;
t139 = sin(pkin(9));
t142 = cos(pkin(9));
t210 = t352 * t211;
t172 = pkin(5) ^ 2;
t135 = sin(pkin(8));
t136 = cos(pkin(8));
t73 = -t135 * t138 - t136 * t141;
t345 = pkin(5) * t73;
t70 = t345 * t351;
t332 = t172 + t70;
t63 = 0.1e1 / (t183 + t332);
t337 = 0.1e1 / pkin(5) * t63;
t58 = sqrt(-(-(t169 + pkin(5)) * pkin(5) + t332) * (pkin(5) * (t169 - pkin(5)) + t332));
t72 = t135 * t141 - t136 * t138;
t339 = t58 * t72;
t64 = t70 + 0.2e1 * t172;
t68 = -pkin(1) + t345;
t56 = -pkin(5) * t339 - t64 * t68;
t57 = pkin(5) * t64 * t72 - t58 * t68;
t41 = (t207 * t352 * t56 + t57 * t210) * t337;
t42 = (-t57 * t208 / 0.4e1 + t56 * t210) * t337;
t33 = t139 * t42 + t142 * t41;
t344 = pkin(6) * t33;
t32 = t168 * t344;
t333 = t170 + t32;
t19 = sqrt(-(-(t168 + pkin(6)) * pkin(6) + t333) * (pkin(6) * (t168 - pkin(6)) + t333));
t17 = atan2(t171 * t19 * t335, -t33);
t15 = sin(t17);
t16 = cos(t17);
t272 = t171 / (pkin(2) ^ 2 + t333) / 0.2e1;
t30 = t32 + 0.2e1 * t170;
t31 = -pkin(2) - t344;
t34 = t139 * t41 - t142 * t42;
t343 = pkin(6) * t34;
t7 = atan2((-t19 * t31 + t30 * t343) * t272, (-t19 * t343 - t30 * t31) * t272);
t5 = sin(t7);
t6 = cos(t7);
t235 = t39 * t5 - t40 * t6;
t3 = t39 * t6 + t40 * t5;
t361 = -t15 * t3 - t16 * t235;
t360 = -t15 * t235 + t16 * t3;
t336 = 0.1e1 / pkin(1) * t63;
t291 = r_base(1) - t97;
t271 = t337 / 0.2e1;
t254 = t40 * pkin(3) + t291;
t253 = t40 * pkin(2) + t291;
t245 = r_base(2) - t341;
t22 = atan2(t34, t33);
t20 = sin(t22);
t21 = cos(t22);
t233 = t20 * t40 + t21 * t39;
t232 = t20 * t39 - t21 * t40;
t231 = t39 * pkin(3) + t245;
t230 = t39 * pkin(2) + t245;
t69 = -pkin(1) * t73 + pkin(5);
t65 = t70 + t160;
t55 = atan2((pkin(1) * t65 * t72 + t58 * t69) * t336 / 0.2e1, -(-pkin(1) * t339 + t65 * t69) * t336 / 0.2e1);
t54 = atan2(t57 * t271, t56 * t271);
t53 = cos(t55);
t52 = cos(t54);
t51 = sin(t55);
t50 = sin(t54);
t47 = -t138 * t50 + t141 * t52;
t46 = t138 * t52 + t141 * t50;
t45 = -t135 * t51 + t136 * t53;
t44 = t135 * t53 + t136 * t51;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t141 + rSges(2,2) * t138 + r_base(1)) + g(2) * (-rSges(2,1) * t138 - rSges(2,2) * t141 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t40 - rSges(3,2) * t39 + t291) + g(2) * (rSges(3,1) * t39 + rSges(3,2) * t40 + t245) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t235 + rSges(4,2) * t3 + t253) + g(2) * (-rSges(4,1) * t3 + rSges(4,2) * t235 + t230) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t256 - rSges(5,2) * t12 + t254) + g(2) * (rSges(5,1) * t12 + rSges(5,2) * t256 + t231) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t136 + rSges(6,1) * t45 - rSges(6,2) * t44 + r_base(1)) + g(2) * (pkin(5) * t135 + rSges(6,1) * t44 + rSges(6,2) * t45 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t232 + rSges(7,2) * t233 + t291) + g(2) * (-rSges(7,1) * t233 + rSges(7,2) * t232 + t245) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t342 + t140 * rSges(8,2) + pkin(7) + r_base(1)) + g(2) * (-t140 * rSges(8,1) + rSges(8,2) * t342 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t47 - rSges(9,2) * t46 + t291) + g(2) * (rSges(9,1) * t46 + rSges(9,2) * t47 + t245) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (t235 * pkin(6) + rSges(10,1) * t361 - rSges(10,2) * t360 + t253) + g(2) * (-t3 * pkin(6) + rSges(10,1) * t360 + rSges(10,2) * t361 + t230) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (t256 * pkin(4) + rSges(11,1) * t363 + t362 * rSges(11,2) + t254) + g(2) * (t12 * pkin(4) - t362 * rSges(11,1) + rSges(11,2) * t363 + t231) + g(3) * (r_base(3) + rSges(11,3)));
U = t1;
