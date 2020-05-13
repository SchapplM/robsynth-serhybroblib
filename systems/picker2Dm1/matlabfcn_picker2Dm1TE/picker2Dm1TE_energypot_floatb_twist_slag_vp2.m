% Calculate potential energy for
% picker2Dm1TE
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
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:11:19
% EndTime: 2020-05-10 00:11:37
% DurationCPUTime: 9.92s
% Computational Cost: add. (86766->502), mult. (233100->653), div. (3590->12), fcn. (61716->14), ass. (0->274)
t134 = cos(pkin(9));
t131 = sin(qJ(1));
t159 = 2 * pkin(1);
t170 = pkin(3) ^ 2;
t133 = cos(qJ(1));
t335 = sin(qJ(2));
t267 = pkin(3) * t335;
t230 = t267 + pkin(7);
t215 = t230 * t133;
t132 = cos(qJ(2));
t314 = t132 * t131;
t271 = pkin(3) * t314;
t206 = t215 - t271;
t175 = pkin(1) ^ 2;
t178 = pkin(7) ^ 2;
t274 = t335 * pkin(7);
t244 = 0.2e1 * t274;
t236 = pkin(3) * t244;
t221 = t178 + t236;
t213 = t175 + t221;
t204 = 0.1e1 / (t159 * t206 + t170 + t213);
t101 = t133 ^ 2;
t192 = t335 ^ 2;
t309 = t170 * t192;
t276 = 0.2e1 * t309;
t207 = -t170 + t221 + t276;
t205 = t101 * t207;
t165 = pkin(4) ^ 2;
t294 = t165 - t178;
t348 = 0.2e1 * t170;
t250 = t348 - t294;
t345 = 3 * t175;
t227 = t345 + t250;
t94 = t175 + t178;
t253 = -t165 + t94;
t313 = t132 * t133;
t270 = pkin(3) * t313;
t273 = pkin(1) * t314;
t347 = 0.4e1 * t170;
t350 = 0.1e1 / pkin(3);
t319 = t350 / 0.2e1;
t157 = 0.1e1 / t319;
t173 = t175 ^ 2;
t211 = -t165 + t213;
t241 = -0.4e1 * t267;
t242 = pkin(1) * t271;
t291 = t175 - t178;
t41 = sqrt(-0.4e1 * pkin(1) * (0.2e1 * (-t273 + t274) * pkin(3) + t253) * t215 + 0.4e1 * t211 * t242 + 0.4e1 * t291 * t309 + pkin(7) * t253 * t241 - t173 - (t178 - (t157 + pkin(4)) * pkin(4)) * (t178 + (t157 - pkin(4)) * pkin(4)) + (-0.4e1 * t205 - 0.2e1 * t250) * t175);
t199 = t204 * ((t131 * t230 + t270) * t41 - ((t244 - 0.4e1 * t273) * pkin(3) + t227) * t215 + (t236 + t227) * t271 + (-0.2e1 * t205 + t276 - t236 - t347 - t253) * pkin(1));
t197 = t199 / 0.4e1;
t208 = t348 + t211;
t90 = pkin(1) * t133;
t282 = 0.2e1 * t90;
t89 = pkin(3) * t132;
t203 = t204 * ((pkin(1) + t206) * t41 + (t207 * t282 + t208 * t230) * t131 + (t133 * t208 + (0.4e1 * t101 - 0.2e1) * pkin(1) * t230) * t89);
t201 = t203 / 0.4e1;
t200 = t350 * t201;
t162 = pkin(5) ^ 2;
t128 = sin(pkin(8));
t129 = cos(pkin(8));
t61 = -t128 * t131 - t129 * t133;
t342 = pkin(5) * t61;
t284 = pkin(1) * t342;
t317 = t162 - 0.2e1 * t284;
t51 = 0.1e1 / (t175 + t317);
t327 = 0.1e1 / pkin(5) * t51;
t46 = sqrt(-(-(t159 + pkin(5)) * pkin(5) + t317) * (pkin(5) * (t159 - pkin(5)) + t317));
t60 = t128 * t133 - t129 * t131;
t332 = t46 * t60;
t52 = t162 - t284;
t56 = -pkin(1) + t342;
t42 = -pkin(5) * t332 - 0.2e1 * t52 * t56;
t349 = 0.2e1 * t60;
t44 = pkin(5) * t349 * t52 - t46 * t56;
t195 = (t197 * t350 * t42 + t44 * t200) * t327;
t198 = t350 * t199;
t27 = (-t44 * t198 / 0.4e1 + t42 * t200) * t327;
t336 = sin(pkin(9));
t24 = t134 * t195 + t27 * t336;
t161 = 0.1e1 / pkin(6);
t160 = pkin(6) ^ 2;
t158 = 2 * pkin(2);
t341 = pkin(6) * t24;
t23 = t341 * t158;
t318 = t160 + t23;
t328 = t161 / ((pkin(2) ^ 2) + t318);
t196 = -t198 / 0.2e1;
t202 = t350 * t203;
t338 = t131 / 0.2e1;
t34 = t133 * t196 + t202 * t338;
t33 = t131 * t196 - t133 * t202 / 0.2e1;
t343 = t33 / 0.2e1;
t14 = sqrt(-(-(t158 + pkin(6)) * pkin(6) + t318) * (pkin(6) * (t158 - pkin(6)) + t318));
t21 = t23 + 0.2e1 * t160;
t22 = -pkin(2) - t341;
t25 = t134 * t27 - t195 * t336;
t340 = pkin(6) * t25;
t6 = t14 * t340 - t21 * t22;
t7 = -t14 * t22 - t21 * t340;
t352 = (t34 * t7 / 0.2e1 + t6 * t343) * t328;
t359 = t24 * t352;
t166 = 0.1e1 / pkin(4);
t100 = t133 * t101;
t110 = -t165 / 0.6e1;
t111 = -t165 / 0.4e1;
t112 = -t165 / 0.3e1;
t113 = -t165 / 0.2e1;
t119 = 0.2e1 / 0.3e1 * t170;
t122 = -t170 / 0.3e1;
t124 = 0.4e1 / 0.3e1 * t175;
t126 = t175 / 0.2e1;
t135 = 15 * t173;
t136 = 15 * t175;
t137 = 10 * t175;
t142 = -0.2e1 * t165;
t143 = -0.5e1 * t165;
t183 = t170 ^ 2;
t144 = 0.5e1 * t183;
t145 = 7 * t173;
t146 = 5 * t173;
t147 = 7 * t175;
t148 = 6 * t175;
t177 = t178 ^ 2;
t151 = 0.3e1 * t177;
t152 = 0.8e1 * t178;
t153 = 0.4e1 * t178;
t155 = 0.2e1 * t178;
t164 = t165 ^ 2;
t182 = pkin(3) * t170;
t167 = t182 ^ 2;
t179 = pkin(1) * t175;
t187 = pkin(7) * t178;
t292 = t173 + t177;
t297 = t155 - t165;
t307 = t178 * t165;
t212 = t297 * t175 + t164 / 0.6e1 + t292 - t307;
t210 = 0.5e1 / 0.6e1 * t183 + t212;
t220 = t178 - t242;
t232 = -0.4e1 * t242;
t252 = t170 + t94;
t74 = t113 + t252;
t222 = t74 * t232;
t304 = t164 / 0.2e1 - t183 / 0.2e1;
t229 = -0.3e1 * t307 + t151 + t304;
t233 = -0.6e1 * t242;
t115 = -0.3e1 / 0.2e1 * t165;
t154 = 0.3e1 * t178;
t303 = t115 + t154;
t316 = t167 + t94 * ((t115 + t155) * t175 - 0.3e1 / 0.2e1 * t307 + t292 + t304);
t87 = 0.10e2 / 0.3e1 * t175;
t224 = ((t87 + t297) * t170 + t210) * t233 + (t135 + (-0.9e1 * t165 + 0.18e2 * t178) * t175 + t229) * t170 + (t136 + t303) * t183 + t316;
t288 = t178 + t345;
t251 = t170 + t288;
t334 = pkin(1) * t131;
t225 = t251 * t89 - t334 * (0.3e1 * t170 + t94);
t114 = -0.2e1 / 0.3e1 * t165;
t123 = -0.2e1 / 0.3e1 * t170;
t255 = t114 + t94;
t298 = t137 + t155;
t302 = t123 + t178;
t92 = -t165 - t170;
t81 = t154 + t92;
t320 = t81 * t175;
t226 = -t334 * (t144 + ((5 * t175) + t81) * t348 + (t123 + t255) * t94) + (t183 + (t114 + t123 + t298) * t170 + t146 + 0.2e1 * t320 + t178 * (t114 + t302)) * t89;
t296 = t164 - t183;
t228 = -0.6e1 * t307 + 0.6e1 * t177 + t296;
t231 = -0.2e1 * t242;
t305 = t179 * t100;
t237 = t305 * t89;
t256 = t114 + t119 + t155;
t315 = t183 + (t119 + t255) * t94;
t118 = 0.4e1 / 0.3e1 * t170;
t257 = t112 + t94;
t72 = t118 + t257;
t238 = t72 * t232 + t315 + (t148 + t256) * t170;
t326 = t173 * t101 ^ 2;
t240 = t326 * t89;
t268 = 0.16e2 * t305;
t243 = pkin(7) * t268;
t245 = 0.20e2 / 0.3e1 * t175;
t277 = pkin(7) * t305;
t249 = 0.8e1 * t277;
t295 = -t165 + t170;
t254 = t154 + t295;
t120 = t170 / 0.3e1;
t258 = t110 + t120 + t178;
t259 = t165 / 0.3e1 + t120 + t155;
t260 = 0.2e1 / 0.3e1 * t165 + t119 + t153;
t261 = 0.4e1 / 0.3e1 * t165 + t118 - 0.2e1 * t178;
t312 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t262 = t131 * t312;
t280 = 0.6e1 * t90;
t263 = pkin(7) * t280;
t281 = 0.4e1 * t90;
t264 = pkin(7) * t281;
t265 = -t334 / 0.2e1;
t308 = t175 * t101;
t269 = 0.12e2 * t308;
t272 = t175 * t89;
t275 = 0.4e1 * t308;
t323 = t182 * t335 * t192;
t278 = -0.8e1 * t323;
t279 = 0.8e1 * t326;
t283 = 0.2e1 * t334;
t285 = pkin(7) * t90;
t287 = 0.4e1 * pkin(7);
t289 = t177 + t183;
t290 = t177 - t173;
t293 = -t170 + t178;
t299 = 0.4e1 / 0.7e1 * t178 - t165 / 0.7e1;
t300 = t126 + t178;
t301 = t175 / 0.3e1 + t178;
t306 = t178 * t175;
t311 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t321 = (-t131 * t179 + t272) * t101;
t322 = t183 * t192 ^ 2;
t324 = t178 * t92;
t88 = t94 ^ 2;
t331 = t88 * (-t170 + t253);
t346 = 4 * t173;
t351 = t111 + t170 / 0.2e1;
t59 = -t183 / 0.6e1 + t212;
t85 = t122 + t178;
t62 = t85 * t231;
t67 = t90 + t230;
t93 = -0.3e1 * t170 + t178;
t70 = t93 * t249;
t71 = 0.10e2 * t320;
t75 = -t165 + t252;
t79 = t90 + pkin(7);
t80 = pkin(7) * t282;
t83 = (t153 + t165) * t175;
t86 = -t175 / 0.3e1 + t178;
t91 = -0.30e2 * t165 + 0.60e2 * t178;
t96 = -(3 * t175) + t178;
t333 = ((-0.24e2 * (0.4e1 / 0.3e1 * t308 + t80 + t86) * t322 * t334 - 0.12e2 * (-0.8e1 / 0.3e1 * t240 + ((t124 + t258) * t89 - (0.7e1 / 0.6e1 * t170 + t110 + t300) * t334) * t275 + (-t170 * t291 - 0.5e1 / 0.3e1 * t173 + t259 * t175 + t178 * (t112 + t85)) * t89 + (-t183 + (-t245 + t260) * t170 - (3 * t173) + t261 * t175 + t177) * t265 + (-t131 * t173 * t100 + ((t175 + t258) * t89 + (t348 - t291) * t265) * t90) * t287) * t309 + 0.24e2 * t85 * t240 + ((t178 + 0.5e1 / 0.2e1 * t170 + 0.3e1 / 0.2e1 * t175 + t113) * t89 + t93 * t334 / 0.2e1) * t243 - 0.6e1 * ((-0.3e1 * t183 + (-t245 + t261) * t170 + t260 * t175 + t290) * t89 - 0.2e1 * (-0.5e1 / 0.3e1 * t183 + (-t175 + t259) * t170 + t178 * (t122 + t257)) * t334) * t308 - 0.6e1 * t226 * t285 - (t167 + ((21 * t175) + t81) * t183 + (t151 + (35 * t173) + t71 + 0.2e1 * t324) * t170 + (t145 + (t143 + t152 - 0.5e1 * t170) * t175 + t178 * (-t165 + t293)) * t94) * t89 + (0.7e1 * t167 + (t147 + t81) * t144 + ((21 * t173) + 0.9e1 * t177 + t71 + 0.6e1 * t324) * t170 + t331) * t334) * t41 + (0.16e2 * (t279 + t243 + (-(8 * t173) + 0.12e2 * t306) * t101 + (0.4e1 * pkin(1) * t187 - 0.12e2 * pkin(7) * t179) * t133 - 0.6e1 * t306 + t292) * t322 + 0.24e2 * (t302 * t279 + 0.14e2 * (-0.32e2 / 0.21e2 * (t178 + t170 / 0.4e1 + t175 / 0.4e1 - t165 / 0.8e1) * t242 + 0.5e1 / 0.42e2 * t183 + (0.16e2 / 0.21e2 * t175 + t299) * t170 + t173 / 0.7e1 + t299 * t175 + t177 - 0.3e1 / 0.7e1 * t307 + t164 / 0.42e2) * t308 + t86 * t222 - t291 * t183 + (-0.10e2 / 0.3e1 * t173 + 0.2e1 * t177 - t307 + t83) * t170 + t59 * t311 + ((-0.2e1 / 0.3e1 * t242 + t111 + t300) * t268 + (-0.8e1 / 0.3e1 * (t301 + t351) * t242 + 0.5e1 / 0.18e2 * t183 + (0.4e1 / 0.3e1 * t178 + t124 + t112) * t170 + t177 + 0.2e1 / 0.3e1 * t306 - 0.2e1 / 0.3e1 * t307 - t173 / 0.3e1 + t164 / 0.18e2) * t280) * pkin(7)) * t309 + 0.16e2 * (-0.6e1 * t178 * t170 + t289) * t326 + 0.32e2 * (t231 * t312 + t74 * t93) * t277 + 0.24e2 * (t85 * t222 - t167 + (-t87 + t294) * t183 + (t83 + t183 / 0.6e1 - t164 / 0.6e1 + t290) * t170 + t59 * t178) * t308 + 0.8e1 * t224 * t285 - 0.8e1 * ((t147 + t303) * t183 + (t145 + (t143 + 0.10e2 * t178) * t175 + t229) * t170 + t316) * t242 + t183 ^ 2 + (t142 + t153 + (28 * t175)) * t167 + (t175 * t91 + (70 * t173) + t228) * t183 + (t228 * t148 + t296 * t155 - 0.6e1 * t177 * t165 + t91 * t173 + (28 * t179 ^ 2) + 0.4e1 * t187 ^ 2) * t170 + t75 * t331) * t67 + (((0.4e1 * t321 + (t89 + t283) * t80 + t96 * t89 + (t113 + t251) * t283) * t278 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t170 + t126 + t110) * t132 * t157 + pkin(1) * t262) * t308 + (-0.8e1 * t237 + ((t112 + t119 + t288) * t89 - (0.8e1 / 0.3e1 * t170 + t257) * t334) * t281) * pkin(7) + t226) * t267) * t41 + (0.32e2 * (t249 + (-0.4e1 * t179 * t271 + t346 + (t347 + t142 + t152) * t175) * t101 + (-t175 + t220 + t351) * t264 + t231 * t311 + t96 * t74) * t323 + 0.8e1 * (t70 + (t312 * t74 + t62) * t269 + (t222 + (t148 + t297) * t170 + t210) * t263 + t224) * t267) * t67) * t79) / ((-0.4e1 * (-t291 * t89 + 0.2e1 * t321 + (0.2e1 * pkin(7) * t270 + t131 * (t170 + (2 * t175))) * pkin(1)) * t309 + 0.8e1 * pkin(7) * t237 + ((pkin(3) * t346 + 0.8e1 * t175 * t182) * t132 + 0.4e1 * t179 * t262) * t101 - 0.4e1 * t225 * t285 - (t170 * t298 + t146 + t289 + 0.6e1 * t306) * t89 + (t144 + (t137 + 0.6e1 * t178) * t170 + t88) * t334) * t41 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t242 + 0.4e1 / 0.9e1 * t170 - t165 / 0.9e1 + t301) * t308 + t86 * t231 + t72 * t311 + (t305 + (t110 + t119 + t220) * t90) * t287) * t309 + t70 + (t312 * t72 + t62) * t269 + t238 * t263 + ((t87 + t256) * t170 + t315) * t233 + t167 + (t136 + t254) * t183 + (t254 * t148 + t295 * t155 + t135 + t151) * t170 + t88 * t75) * t67 + ((t278 * t334 + (-0.2e1 * t101 * t272 + (t89 - t334) * t80 + t225) * t241) * t41 + (0.8e1 * (t80 + t275 + t96) * t323 + 0.6e1 * (t293 * t275 + (t72 + t231) * t264 + t238) * t267) * t67) * t79);
t239 = t319 * t333;
t216 = (-t33 * t350 * t41 / 0.2e1 + t34 * t239) * t166;
t353 = (t319 * t34 * t41 + t239 * t33) * t166;
t310 = t166 / pkin(3) ^ 2;
t15 = (t197 * t333 - t41 * t203 / 0.4e1) * t310;
t16 = (t197 * t41 + t201 * t333) * t310;
t63 = -t131 * t335 - t313;
t64 = -t133 * t335 + t314;
t8 = t15 * t63 + t16 * t64;
t9 = t15 * t64 - t16 * t63;
t358 = t216 * t9 + t353 * t8;
t357 = -t216 * t8 + t353 * t9;
t218 = (-t34 * t6 / 0.2e1 + t343 * t7) * t328;
t354 = t218 * t24;
t344 = pkin(5) * m(6);
t339 = -t128 / 0.2e1;
t337 = t133 / 0.2e1;
t330 = t34 * pkin(3) - t90;
t329 = t34 * pkin(2) - t90;
t325 = 0.1e1 / pkin(1) * t51;
t286 = -m(7) - m(9) - m(3);
t266 = t14 * t161 / pkin(2);
t248 = t33 * pkin(3) - t334;
t247 = t33 * pkin(2) - t334;
t246 = t34 * t24 + t25 * t33;
t235 = -t266 / 0.2e1;
t234 = t266 / 0.2e1;
t223 = t24 * t33 - t25 * t34;
t214 = -m(6) - m(11) - m(5) - m(10) - m(4) - m(8) - m(2) - m(1);
t57 = -pkin(1) * t61 + pkin(5);
t53 = t175 - t284;
t45 = pkin(1) * t349 * t53 + t46 * t57;
t43 = -pkin(1) * t332 + 0.2e1 * t53 * t57;
t39 = (t42 * t337 - t131 * t44 / 0.2e1) * t327;
t38 = (t337 * t44 + t338 * t42) * t327;
t37 = (-t129 * t43 / 0.2e1 + t45 * t339) * t325;
t36 = (t43 * t339 + t129 * t45 / 0.2e1) * t325;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t214 + t286) * r_base(3)) * g(3) + (-(-t234 * t352 + t354) * mrSges(10,2) - m(10) * (-pkin(6) * t352 + t247) + t352 * mrSges(4,1) - m(11) * (pkin(4) * t353 + t248) - t353 * mrSges(5,1) - t357 * mrSges(11,1) - t358 * mrSges(11,2) + t246 * mrSges(7,2) - m(4) * t247 - m(5) * t248 - mrSges(9,1) * t38 - mrSges(9,2) * t39 - mrSges(3,1) * t33 - mrSges(3,2) * t34 - t128 * t344 - t36 * mrSges(6,1) + t132 * mrSges(8,1) + mrSges(2,1) * t131 + mrSges(2,2) * t133 - t37 * mrSges(6,2) + t286 * (r_base(2) - t334) - t335 * mrSges(8,2) + t223 * mrSges(7,1) - t218 * mrSges(4,2) - t216 * mrSges(5,2) - (t218 * t235 - t359) * mrSges(10,1) + t214 * r_base(2) - mrSges(1,2)) * g(2) + (-(t235 * t352 + t354) * mrSges(10,1) - t352 * mrSges(4,2) + t353 * mrSges(5,2) - m(4) * t329 - m(5) * t330 + t357 * mrSges(11,2) - t358 * mrSges(11,1) + t246 * mrSges(7,1) - mrSges(3,1) * t34 + mrSges(3,2) * t33 - mrSges(9,1) * t39 + mrSges(9,2) * t38 - t129 * t344 + t36 * mrSges(6,2) - m(8) * pkin(7) - t132 * mrSges(8,2) + mrSges(2,1) * t133 - mrSges(2,2) * t131 + t286 * (-t90 + r_base(1)) - t335 * mrSges(8,1) - t37 * mrSges(6,1) - t223 * mrSges(7,2) - m(10) * (pkin(6) * t218 + t329) - t218 * mrSges(4,1) - m(11) * (pkin(4) * t216 + t330) - t216 * mrSges(5,1) - (t218 * t234 + t359) * mrSges(10,2) + t214 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
