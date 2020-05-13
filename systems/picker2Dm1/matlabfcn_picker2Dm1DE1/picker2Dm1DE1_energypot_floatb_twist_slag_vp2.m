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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 19:54
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1DE1_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 09:10:56
% EndTime: 2020-05-10 09:11:14
% DurationCPUTime: 11.64s
% Computational Cost: add. (173294->509), mult. (465818->634), div. (7180->13), fcn. (123736->38), ass. (0->288)
t185 = pkin(1) ^ 2;
t183 = t185 ^ 2;
t188 = pkin(7) ^ 2;
t143 = cos(qJ(1));
t109 = t143 ^ 2;
t192 = pkin(3) ^ 2;
t344 = sin(qJ(2));
t279 = pkin(3) * t344;
t252 = pkin(7) * t279;
t245 = 0.2e1 * t252;
t231 = t188 + t245;
t203 = t344 ^ 2;
t318 = t192 * t203;
t287 = 0.2e1 * t318;
t219 = -t192 + t231 + t287;
t217 = t109 * t219;
t177 = pkin(4) ^ 2;
t225 = t185 + t231;
t223 = -t177 + t225;
t241 = t279 + pkin(7);
t228 = t241 * t143;
t140 = sin(qJ(1));
t142 = cos(qJ(2));
t328 = t142 * t140;
t283 = pkin(3) * t328;
t251 = pkin(1) * t283;
t350 = 0.2e1 * t192;
t300 = t188 + t350;
t260 = t177 - t300;
t301 = -t188 + t177;
t262 = t185 - t301;
t285 = pkin(1) * t328;
t286 = t344 * pkin(7);
t304 = t185 - t188;
t355 = 0.1e1 / pkin(3);
t334 = t355 / 0.2e1;
t362 = 0.1e1 / t334;
t215 = sqrt(-0.4e1 * pkin(1) * (0.2e1 * (-t285 + t286) * pkin(3) + t262) * t228 + 0.4e1 * t223 * t251 + 0.4e1 * t304 * t318 - 0.4e1 * t262 * t252 - t183 - (t188 - (t362 + pkin(4)) * pkin(4)) * (t188 + (t362 - pkin(4)) * pkin(4)) + (-0.4e1 * t217 + 0.2e1 * t260) * t185);
t178 = 0.1e1 / pkin(4);
t276 = t178 * t334;
t100 = -0.30e2 * t177 + 0.60e2 * t188;
t101 = t188 - 0.3e1 * t192;
t102 = t185 + t188;
t104 = -0.3e1 * t185 + t188;
t108 = t143 * t109;
t110 = 0.10e2 / 0.3e1 * t185;
t111 = -0.20e2 / 0.3e1 * t185;
t119 = -t177 / 0.6e1;
t120 = -t177 / 0.4e1;
t121 = -t177 / 0.3e1;
t122 = -t177 / 0.2e1;
t124 = -0.3e1 / 0.2e1 * t177;
t128 = 0.2e1 / 0.3e1 * t192;
t131 = -t192 / 0.3e1;
t133 = 0.4e1 / 0.3e1 * t185;
t135 = t185 / 0.2e1;
t145 = 0.15e2 * t183;
t147 = 0.10e2 * t185;
t152 = 0.18e2 * t188;
t153 = -0.2e1 * t177;
t154 = -0.5e1 * t177;
t155 = -0.6e1 * t177;
t194 = t192 ^ 2;
t156 = 0.5e1 * t194;
t159 = 0.7e1 * t183;
t160 = 0.5e1 * t183;
t161 = 0.6e1 * t185;
t162 = 0.2e1 * t185;
t187 = t188 ^ 2;
t164 = 0.3e1 * t187;
t165 = 0.8e1 * t188;
t166 = 0.6e1 * t188;
t167 = 0.4e1 * t188;
t168 = 0.3e1 * t188;
t169 = 0.2e1 * t188;
t176 = t177 ^ 2;
t193 = pkin(3) * t192;
t179 = t193 ^ 2;
t189 = pkin(1) * t185;
t198 = pkin(7) * t188;
t305 = t183 + t187;
t308 = t169 - t177;
t321 = t188 * t177;
t224 = t308 * t185 + t176 / 0.6e1 + t305 - t321;
t221 = 0.5e1 / 0.6e1 * t194 + t224;
t98 = pkin(3) * t142;
t284 = t185 * t98;
t222 = t109 * (-t140 * t189 + t284);
t307 = t176 - t194;
t227 = 0.6e1 * t187 + t307 - 0.6e1 * t321;
t230 = t188 - t251;
t243 = -0.4e1 * t251;
t298 = t192 + t188;
t259 = t185 + t298;
t84 = t122 + t259;
t232 = t84 * t243;
t317 = t176 / 0.2e1 - t194 / 0.2e1;
t240 = -0.3e1 * t321 + t164 + t317;
t244 = -0.6e1 * t251;
t312 = 0.15e2 * t185 + t168;
t331 = t179 + t102 * ((t124 + t169) * t185 - 0.3e1 / 0.2e1 * t321 + t305 + t317);
t235 = ((t110 + t308) * t192 + t221) * t244 + (t145 + (t152 - 0.9e1 * t177) * t185 + t240) * t192 + (t124 + t312) * t194 + t331;
t352 = 0.3e1 * t185;
t258 = t352 + t298;
t343 = pkin(1) * t140;
t237 = t258 * t98 - t343 * (0.3e1 * t192 + t102);
t123 = -0.2e1 / 0.3e1 * t177;
t132 = -0.2e1 / 0.3e1 * t192;
t310 = t153 - 0.2e1 * t192;
t263 = t166 + t310;
t264 = t123 + t102;
t311 = t147 + t169;
t316 = t132 + t188;
t238 = -t343 * (t156 + (t147 + t263) * t192 + (t132 + t264) * t102) + (t194 + (t123 + t132 + t311) * t192 + t160 + t263 * t185 + t188 * (t123 + t316)) * t98;
t242 = -0.2e1 * t251;
t319 = t189 * t108;
t246 = t319 * t98;
t323 = t183 * t109 ^ 2;
t247 = t323 * t98;
t99 = pkin(1) * t143;
t90 = t99 + pkin(7);
t248 = t90 * t279;
t265 = t123 + t128 + t169;
t330 = t194 + (t128 + t264) * t102;
t127 = 0.4e1 / 0.3e1 * t192;
t266 = t121 + t102;
t82 = t127 + t266;
t249 = t82 * t243 + t330 + (t161 + t265) * t192;
t278 = t344 * t203 * t193 * t90;
t250 = -0.8e1 * t278;
t280 = 0.16e2 * t319;
t253 = pkin(7) * t280;
t290 = pkin(7) * t319;
t257 = 0.8e1 * t290;
t299 = -t192 + t188;
t261 = -t177 + t299;
t129 = t192 / 0.3e1;
t267 = t119 + t129 + t188;
t268 = t177 / 0.3e1 + t129 + t169;
t269 = 0.2e1 / 0.3e1 * t177 + t128 + t167;
t270 = 0.4e1 / 0.3e1 * t177 + t127 - 0.2e1 * t188;
t326 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t273 = t140 * t326;
t291 = 0.6e1 * t99;
t274 = pkin(7) * t291;
t292 = 0.4e1 * t99;
t275 = pkin(7) * t292;
t277 = -t343 / 0.2e1;
t322 = t185 * t109;
t281 = 0.12e2 * t322;
t327 = t142 * t143;
t282 = pkin(3) * t327;
t288 = 0.4e1 * t322;
t289 = 0.8e1 * t323;
t294 = 0.2e1 * t343;
t295 = pkin(7) * t99;
t297 = 0.4e1 * pkin(7);
t302 = t187 + t194;
t303 = t187 - t183;
t309 = t154 - 0.5e1 * t192;
t313 = 0.4e1 / 0.7e1 * t188 - t177 / 0.7e1;
t314 = t135 + t188;
t315 = t185 / 0.3e1 + t188;
t320 = t188 * t185;
t325 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t329 = t203 ^ 2 * t194;
t97 = t102 ^ 2;
t340 = t97 * (t185 + t261);
t349 = 0.4e1 * t192;
t351 = -0.6e1 * t192;
t353 = 0.4e1 * t183;
t356 = t168 - t177 - t192;
t357 = t120 + t192 / 0.2e1;
t71 = -t194 / 0.6e1 + t224;
t95 = t131 + t188;
t74 = t95 * t242;
t78 = t99 + t241;
t81 = t101 * t257;
t85 = t356 * t147;
t86 = -t177 + t259;
t293 = 0.2e1 * t99;
t91 = pkin(7) * t293;
t93 = (t167 + t177) * t185;
t96 = -t185 / 0.3e1 + t188;
t342 = ((-0.24e2 * (0.4e1 / 0.3e1 * t322 + t91 + t96) * t329 * t343 + (0.4e1 * t222 + (t98 + t294) * t91 + t104 * t98 + (t122 + t258) * t294) * t250 - 0.12e2 * (-0.8e1 / 0.3e1 * t247 + ((t133 + t267) * t98 - (0.7e1 / 0.6e1 * t192 + t119 + t314) * t343) * t288 + (-t192 * t304 - 0.5e1 / 0.3e1 * t183 + t268 * t185 + t188 * (t121 + t95)) * t98 + (-t194 + (t111 + t269) * t192 - 0.3e1 * t183 + t270 * t185 + t187) * t277 + (-t140 * t183 * t108 + ((t185 + t267) * t98 + (-t185 + t300) * t277) * t99) * t297) * t318 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t192 + t135 + t119) * t142 * t362 + pkin(1) * t273) * t322 + (-0.8e1 * t246 + ((t121 + t128 + t352 + t188) * t98 - (0.8e1 / 0.3e1 * t192 + t266) * t343) * t292) * pkin(7) + t238) * t248 + 0.24e2 * t95 * t247 + ((t188 + 0.5e1 / 0.2e1 * t192 + 0.3e1 / 0.2e1 * t185 + t122) * t98 + t101 * t343 / 0.2e1) * t253 - 0.6e1 * ((-0.3e1 * t194 + (t111 + t270) * t192 + t269 * t185 + t303) * t98 - 0.2e1 * (-0.5e1 / 0.3e1 * t194 + (-t185 + t268) * t192 + t188 * (t131 + t266)) * t343) * t322 - 0.6e1 * t238 * t295 - (t179 + (0.21e2 * t185 + t356) * t194 + (t188 * t310 + t164 + 0.35e2 * t183 + t85) * t192 + (t159 + (t165 + t309) * t185 + t188 * t261) * t102) * t98 + (0.7e1 * t179 + (0.35e2 * t185 + 0.15e2 * t188 + t309) * t194 + (0.21e2 * t183 + t85 + 0.9e1 * t187 + (t155 + t351) * t188) * t192 + t340) * t343) * t215 + (0.16e2 * (t289 + t253 + (-0.8e1 * t183 + 0.12e2 * t320) * t109 + (0.4e1 * pkin(1) * t198 - 0.12e2 * pkin(7) * t189) * t143 - 0.6e1 * t320 + t305) * t329 + 0.32e2 * (t257 + (-0.4e1 * t189 * t283 + t353 + (t349 + t153 + t165) * t185) * t109 + (-t185 + t230 + t357) * t275 + t242 * t325 + t104 * t84) * t278 + 0.24e2 * (t316 * t289 + 0.14e2 * (-0.32e2 / 0.21e2 * (t188 + t192 / 0.4e1 + t185 / 0.4e1 - t177 / 0.8e1) * t251 + 0.5e1 / 0.42e2 * t194 + (0.16e2 / 0.21e2 * t185 + t313) * t192 + t183 / 0.7e1 + t313 * t185 + t187 - 0.3e1 / 0.7e1 * t321 + t176 / 0.42e2) * t322 + t96 * t232 - t304 * t194 + (-0.10e2 / 0.3e1 * t183 + 0.2e1 * t187 - t321 + t93) * t192 + t71 * t325 + ((-0.2e1 / 0.3e1 * t251 + t120 + t314) * t280 + (-0.8e1 / 0.3e1 * (t315 + t357) * t251 + 0.5e1 / 0.18e2 * t194 + (0.4e1 / 0.3e1 * t188 + t133 + t121) * t192 + t187 + 0.2e1 / 0.3e1 * t320 - 0.2e1 / 0.3e1 * t321 - t183 / 0.3e1 + t176 / 0.18e2) * t291) * pkin(7)) * t318 + 0.8e1 * (t81 + (t326 * t84 + t74) * t281 + (t232 + (t161 + t308) * t192 + t221) * t274 + t235) * t248 + 0.16e2 * (t188 * t351 + t302) * t323 + 0.32e2 * (t101 * t84 + t242 * t326) * t290 + 0.24e2 * (t95 * t232 - t179 + (-t110 + t301) * t194 + (t93 + t194 / 0.6e1 - t176 / 0.6e1 + t303) * t192 + t71 * t188) * t322 + 0.8e1 * t235 * t295 - 0.8e1 * ((t124 + t168 + 0.7e1 * t185) * t194 + (t159 + (t154 + 0.10e2 * t188) * t185 + t240) * t192 + t331) * t251 + t194 ^ 2 + (t153 + t167 + 0.28e2 * t185) * t179 + (t100 * t185 + 0.70e2 * t183 + t227) * t194 + (t100 * t183 + t155 * t187 + t227 * t161 + t307 * t169 + 0.28e2 * t189 ^ 2 + 0.4e1 * t198 ^ 2) * t192 + t86 * t340) * t78) / ((t250 * t343 - 0.4e1 * (0.2e1 * t222 - t304 * t98 + (0.2e1 * pkin(7) * t282 + t140 * (t162 + t192)) * pkin(1)) * t318 - 0.4e1 * (-0.2e1 * t109 * t284 + (t98 - t343) * t91 + t237) * t248 + 0.8e1 * pkin(7) * t246 + ((pkin(3) * t353 + 0.8e1 * t185 * t193) * t142 + 0.4e1 * t189 * t273) * t109 - 0.4e1 * t237 * t295 - (t192 * t311 + t160 + t302 + 0.6e1 * t320) * t98 + (t156 + (t147 + t166) * t192 + t97) * t343) * t215 + (0.8e1 * (t91 + t288 + t104) * t278 + 0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t251 + 0.4e1 / 0.9e1 * t192 - t177 / 0.9e1 + t315) * t322 + t96 * t242 + t82 * t325 + (t319 + (t119 + t128 + t230) * t99) * t297) * t318 + 0.6e1 * (t299 * t288 + (t82 + t242) * t275 + t249) * t248 + t81 + (t326 * t82 + t74) * t281 + t249 * t274 + ((t110 + t265) * t192 + t330) * t244 + t179 + (-t177 + t192 + t312) * t194 + (t145 + (t152 + t155 + 0.6e1 * t192) * t185 + t164 + (t153 + t350) * t188) * t192 + t97 * t86) * t78);
t27 = atan2(t215 * t276, t276 * t342);
t25 = sin(t27);
t26 = cos(t27);
t171 = 0.2e1 * pkin(1);
t218 = t228 - t283;
t216 = 0.1e1 / (t171 * t218 + t192 + t225);
t239 = t352 - t260;
t354 = -0.2e1 * pkin(1);
t211 = t216 * ((t140 * t241 + t282) * t215 + t217 * t354 - ((-0.4e1 * t285 + 0.2e1 * t286) * pkin(3) + t239) * t228 + (t245 + t239) * t283 + (t287 - t245 - t349 - t262) * pkin(1));
t210 = t355 * t211;
t220 = t350 + t223;
t214 = t216 * ((pkin(1) + t218) * t215 + (t219 * t293 + t220 * t241) * t140 + (t143 * t220 + (0.4e1 * t109 - 0.2e1) * pkin(1) * t241) * t98);
t208 = atan2(t214 * t334, t210 / 0.2e1);
t206 = sin(t208);
t207 = cos(t208);
t39 = -t140 * t207 - t143 * t206;
t40 = t140 * t206 - t143 * t207;
t12 = t25 * t40 + t26 * t39;
t254 = -t25 * t39 + t40 * t26;
t209 = t211 / 0.4e1;
t324 = t178 / t192;
t23 = (t209 * t342 - t215 * t214 / 0.4e1) * t324;
t213 = t214 / 0.4e1;
t24 = (t209 * t215 + t213 * t342) * t324;
t75 = -t140 * t344 - t327;
t76 = -t143 * t344 + t328;
t10 = atan2(t23 * t75 + t24 * t76, -t23 * t76 + t24 * t75);
t8 = sin(t10);
t9 = cos(t10);
t366 = t12 * t8 - t254 * t9;
t365 = t12 * t9 + t254 * t8;
t173 = 0.1e1 / pkin(6);
t335 = 0.1e1 / pkin(2) / 0.2e1;
t170 = 0.1e1 / t335;
t172 = pkin(6) ^ 2;
t141 = sin(pkin(9));
t144 = cos(pkin(9));
t212 = t355 * t213;
t174 = pkin(5) ^ 2;
t137 = sin(pkin(8));
t138 = cos(pkin(8));
t73 = -t137 * t140 - t138 * t143;
t347 = pkin(5) * t73;
t70 = t347 * t354;
t332 = t174 + t70;
t63 = 0.1e1 / (t185 + t332);
t337 = 0.1e1 / pkin(5) * t63;
t58 = sqrt(-(-(t171 + pkin(5)) * pkin(5) + t332) * (pkin(5) * (t171 - pkin(5)) + t332));
t72 = t137 * t143 - t138 * t140;
t341 = t58 * t72;
t64 = t70 + 0.2e1 * t174;
t68 = -pkin(1) + t347;
t56 = -pkin(5) * t341 - t64 * t68;
t57 = pkin(5) * t64 * t72 - t58 * t68;
t41 = (t209 * t355 * t56 + t57 * t212) * t337;
t42 = (-t57 * t210 / 0.4e1 + t56 * t212) * t337;
t33 = t141 * t42 + t144 * t41;
t346 = pkin(6) * t33;
t32 = t170 * t346;
t333 = t172 + t32;
t19 = sqrt(-(-(t170 + pkin(6)) * pkin(6) + t333) * (pkin(6) * (t170 - pkin(6)) + t333));
t17 = atan2(t173 * t19 * t335, -t33);
t15 = sin(t17);
t16 = cos(t17);
t272 = t173 / (pkin(2) ^ 2 + t333) / 0.2e1;
t30 = t32 + 0.2e1 * t172;
t31 = -pkin(2) - t346;
t34 = t141 * t41 - t144 * t42;
t345 = pkin(6) * t34;
t7 = atan2((-t19 * t31 + t30 * t345) * t272, (-t19 * t345 - t30 * t31) * t272);
t5 = sin(t7);
t6 = cos(t7);
t236 = t39 * t5 - t40 * t6;
t3 = t39 * t6 + t40 * t5;
t364 = t15 * t3 + t16 * t236;
t363 = -t15 * t236 + t16 * t3;
t348 = pkin(5) * m(6);
t339 = t40 * pkin(3) - t99;
t338 = t40 * pkin(2) - t99;
t336 = 0.1e1 / pkin(1) * t63;
t296 = -m(3) - m(9) - m(7);
t271 = t337 / 0.2e1;
t256 = t39 * pkin(3) - t343;
t255 = t39 * pkin(2) - t343;
t22 = atan2(t34, t33);
t20 = sin(t22);
t21 = cos(t22);
t234 = t20 * t40 + t21 * t39;
t233 = t20 * t39 - t21 * t40;
t226 = -m(2) - m(8) - m(4) - m(10) - m(5) - m(11) - m(6) - m(1);
t69 = -pkin(1) * t73 + pkin(5);
t65 = t70 + t162;
t55 = atan2((pkin(1) * t65 * t72 + t58 * t69) * t336 / 0.2e1, -(-pkin(1) * t341 + t65 * t69) * t336 / 0.2e1);
t54 = atan2(t57 * t271, t56 * t271);
t53 = cos(t55);
t52 = cos(t54);
t51 = sin(t55);
t50 = sin(t54);
t47 = -t140 * t50 + t143 * t52;
t46 = t140 * t52 + t143 * t50;
t45 = -t137 * t51 + t138 * t53;
t44 = t137 * t53 + t138 * t51;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t226 + t296) * r_base(3)) * g(3) + (-t137 * t348 - t363 * mrSges(10,1) + t364 * mrSges(10,2) + t365 * mrSges(11,1) - t366 * mrSges(11,2) - mrSges(3,1) * t39 + mrSges(2,1) * t140 + mrSges(2,2) * t143 + t296 * (r_base(2) - t343) - t344 * mrSges(8,2) - mrSges(3,2) * t40 - mrSges(9,1) * t46 - mrSges(9,2) * t47 + t3 * mrSges(4,1) - m(10) * (-pkin(6) * t3 + t255) + t226 * r_base(2) - t233 * mrSges(7,2) + t234 * mrSges(7,1) - t236 * mrSges(4,2) - t254 * mrSges(5,2) - m(4) * t255 - m(11) * (pkin(4) * t12 + t256) - m(5) * t256 - t12 * mrSges(5,1) + t142 * mrSges(8,1) - t44 * mrSges(6,1) - t45 * mrSges(6,2) - mrSges(1,2)) * g(2) + (-t3 * mrSges(4,2) - mrSges(9,1) * t47 + mrSges(9,2) * t46 - t138 * t348 + t363 * mrSges(10,2) + t364 * mrSges(10,1) - t365 * mrSges(11,2) - t366 * mrSges(11,1) + mrSges(2,1) * t143 - mrSges(2,2) * t140 - m(4) * t338 - m(5) * t339 - m(8) * pkin(7) - t344 * mrSges(8,1) + t296 * (-t99 + r_base(1)) - mrSges(3,1) * t40 + mrSges(3,2) * t39 + t12 * mrSges(5,2) + t226 * r_base(1) - t233 * mrSges(7,1) - t234 * mrSges(7,2) - m(10) * (pkin(6) * t236 + t338) - t236 * mrSges(4,1) - m(11) * (pkin(4) * t254 + t339) - t254 * mrSges(5,1) - t142 * mrSges(8,2) + t44 * mrSges(6,2) - t45 * mrSges(6,1) - mrSges(1,1)) * g(1);
U = t1;
