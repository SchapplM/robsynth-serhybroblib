% Calculate potential energy for
% picker2Dm2DE1
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
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 14:28:06
% EndTime: 2020-05-09 14:28:20
% DurationCPUTime: 8.65s
% Computational Cost: add. (79926->488), mult. (235156->602), div. (1676->6), fcn. (42008->31), ass. (0->246)
t130 = sin(qJ(1));
t165 = pkin(3) ^ 2;
t169 = (pkin(1) ^ 2);
t171 = pkin(7) ^ 2;
t91 = t169 + t171;
t227 = t165 + t91;
t133 = cos(qJ(1));
t129 = sin(qJ(2));
t85 = pkin(3) * t129;
t304 = t85 + pkin(7);
t241 = t304 * t133;
t132 = cos(qJ(2));
t288 = t130 * t132;
t247 = pkin(3) * t288;
t221 = pkin(1) * t247;
t68 = -0.2e1 * t221;
t309 = 0.1e1 / pkin(3);
t291 = t309 / 0.2e1;
t158 = 0.1e1 / t291;
t75 = pkin(7) * t129 * t158;
t196 = 0.1e1 / (0.2e1 * pkin(1) * t241 + t227 + t68 + t75);
t96 = t129 ^ 2;
t299 = t165 * t96;
t256 = 0.2e1 * t299;
t269 = -t165 + t171;
t200 = t75 + t256 + t269;
t99 = t133 ^ 2;
t198 = t200 * t99;
t210 = -pkin(1) + t247;
t149 = 3 * t169;
t160 = pkin(4) ^ 2;
t270 = t160 - t171;
t228 = -0.2e1 * t165 + t270;
t211 = t149 - t228;
t214 = -0.4e1 * t221;
t229 = -t160 + t91;
t287 = t132 * t133;
t246 = pkin(3) * t287;
t307 = 0.4e1 * t165;
t167 = t169 ^ 2;
t217 = t75 + t229;
t258 = -0.4e1 * t85;
t267 = t169 - t171;
t305 = 2 * t169;
t315 = 0.4e1 * pkin(1);
t38 = sqrt(-0.4e1 * t169 * t198 + 0.4e1 * t267 * t299 + pkin(7) * t229 * t258 - t167 + t228 * t305 - (t171 - (t158 + pkin(4)) * pkin(4)) * (t171 + (t158 - pkin(4)) * pkin(4)) + (-(t68 + t217) * t241 + t217 * t247) * t315);
t193 = t196 * ((t130 * t304 + t246) * t38 - (t211 + t75 + t214) * t241 + t210 * t75 + t211 * t247 + (-0.2e1 * t198 + t256 - t307 - t229) * pkin(1));
t191 = t193 / 0.4e1;
t308 = 0.2e1 * t165;
t202 = t308 + t217;
t87 = pkin(1) * t133;
t261 = 0.2e1 * t87;
t86 = pkin(3) * t132;
t195 = t196 * ((t241 - t210) * t38 + (t200 * t261 + t202 * t304) * t130 + (t133 * t202 + (0.4e1 * t99 - 0.2e1) * pkin(1) * t304) * t86);
t161 = 0.1e1 / pkin(4);
t284 = t161 / pkin(3) ^ 2;
t108 = -t160 / 0.6e1;
t109 = -t160 / 0.4e1;
t110 = -t160 / 0.3e1;
t111 = -t160 / 0.2e1;
t117 = 0.2e1 / 0.3e1 * t165;
t120 = -t165 / 0.3e1;
t122 = 0.4e1 / 0.3e1 * t169;
t124 = t169 / 0.2e1;
t135 = 15 * t167;
t136 = 15 * t169;
t137 = 10 * t169;
t142 = -0.2e1 * t160;
t143 = -0.5e1 * t160;
t176 = t165 ^ 2;
t144 = 0.5e1 * t176;
t145 = 7 * t167;
t146 = 5 * t167;
t147 = 7 * t169;
t148 = 6 * t169;
t170 = t171 ^ 2;
t152 = 0.3e1 * t170;
t153 = 0.8e1 * t171;
t154 = 0.4e1 * t171;
t156 = 0.2e1 * t171;
t159 = t160 ^ 2;
t175 = pkin(3) * t165;
t162 = t175 ^ 2;
t172 = pkin(1) * t169;
t180 = pkin(7) * t171;
t268 = t167 + t170;
t273 = t156 - t160;
t283 = t171 * t160;
t199 = t273 * t169 + t159 / 0.6e1 + t268 - t283;
t197 = 0.5e1 / 0.6e1 * t176 + t199;
t201 = t171 - t221;
t71 = t111 + t227;
t203 = t71 * t214;
t281 = t159 / 0.2e1 - t176 / 0.2e1;
t213 = -0.3e1 * t283 + t152 + t281;
t215 = -0.6e1 * t221;
t113 = -0.3e1 / 0.2e1 * t160;
t155 = 0.3e1 * t171;
t280 = t113 + t155;
t290 = t162 + t91 * ((t113 + t156) * t169 - 0.3e1 / 0.2e1 * t283 + t268 + t281);
t83 = 0.10e2 / 0.3e1 * t169;
t207 = ((t83 + t273) * t165 + t197) * t215 + (t135 + (-0.9e1 * t160 + 0.18e2 * t171) * t169 + t213) * t165 + (t136 + t280) * t176 + t290;
t274 = t149 + t171;
t231 = t165 + t274;
t303 = pkin(1) * t130;
t208 = t231 * t86 - t303 * (0.3e1 * t165 + t91);
t112 = -0.2e1 / 0.3e1 * t160;
t121 = -0.2e1 / 0.3e1 * t165;
t232 = t112 + t91;
t275 = t137 + t156;
t279 = t121 + t171;
t89 = -t160 - t165;
t77 = t155 + t89;
t292 = t77 * t169;
t209 = -t303 * (t144 + ((5 * t169) + t77) * t308 + (t121 + t232) * t91) + (t176 + (t112 + t121 + t275) * t165 + t146 + 0.2e1 * t292 + t171 * (t112 + t279)) * t86;
t272 = t159 - t176;
t212 = -0.6e1 * t283 + 0.6e1 * t170 + t272;
t233 = t112 + t117 + t156;
t289 = t176 + (t117 + t232) * t91;
t116 = 0.4e1 / 0.3e1 * t165;
t234 = t110 + t91;
t69 = t116 + t234;
t218 = t69 * t214 + t289 + (t148 + t233) * t165;
t98 = t133 * t99;
t295 = t172 * t98;
t219 = t295 * t86;
t298 = t167 * t99 ^ 2;
t220 = t298 * t86;
t224 = 0.20e2 / 0.3e1 * t169;
t250 = 0.16e2 * t295;
t226 = pkin(7) * t250;
t271 = -t160 + t165;
t230 = t155 + t271;
t118 = t165 / 0.3e1;
t235 = t108 + t118 + t171;
t236 = t160 / 0.3e1 + t118 + t156;
t237 = 0.2e1 / 0.3e1 * t160 + t117 + t154;
t238 = 0.4e1 / 0.3e1 * t160 + t116 - 0.2e1 * t171;
t257 = pkin(7) * t295;
t239 = 0.8e1 * t257;
t286 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t240 = t130 * t286;
t259 = 0.6e1 * t87;
t242 = pkin(7) * t259;
t260 = 0.4e1 * t87;
t243 = pkin(7) * t260;
t245 = -t303 / 0.2e1;
t248 = t169 * t86;
t297 = t169 * t99;
t251 = 0.12e2 * t297;
t294 = t175 * t129 * t96;
t253 = -0.8e1 * t294;
t254 = 0.4e1 * t297;
t255 = 0.8e1 * t298;
t262 = 0.2e1 * t303;
t263 = pkin(7) * t87;
t264 = 0.4e1 * pkin(7);
t265 = t170 + t176;
t266 = t170 - t167;
t276 = 0.4e1 / 0.7e1 * t171 - t160 / 0.7e1;
t277 = t124 + t171;
t278 = t169 / 0.3e1 + t171;
t282 = t171 * t169;
t285 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t293 = t176 * t96 ^ 2;
t296 = t171 * t89;
t84 = t91 ^ 2;
t300 = t84 * (-t165 + t229);
t301 = (-t130 * t172 + t248) * t99;
t306 = 4 * t167;
t310 = t109 + t165 / 0.2e1;
t55 = -t176 / 0.6e1 + t199;
t81 = t120 + t171;
t58 = t81 * t68;
t74 = t87 + pkin(7);
t63 = t85 + t74;
t90 = -0.3e1 * t165 + t171;
t66 = t90 * t239;
t67 = 0.10e2 * t292;
t72 = -t160 + t227;
t76 = pkin(7) * t261;
t79 = (t154 + t160) * t169;
t82 = -t169 / 0.3e1 + t171;
t88 = -0.30e2 * t160 + 0.60e2 * t171;
t93 = -(3 * t169) + t171;
t302 = ((-0.24e2 * (0.4e1 / 0.3e1 * t297 + t76 + t82) * t293 * t303 - 0.12e2 * (-0.8e1 / 0.3e1 * t220 + ((t122 + t235) * t86 - (0.7e1 / 0.6e1 * t165 + t108 + t277) * t303) * t254 + (-t165 * t267 - 0.5e1 / 0.3e1 * t167 + t236 * t169 + t171 * (t110 + t81)) * t86 + (-t176 + (-t224 + t237) * t165 - (3 * t167) + t238 * t169 + t170) * t245 + (-t130 * t167 * t98 + ((t169 + t235) * t86 + (t308 - t267) * t245) * t87) * t264) * t299 + 0.24e2 * t81 * t220 + ((t171 + 0.5e1 / 0.2e1 * t165 + 0.3e1 / 0.2e1 * t169 + t111) * t86 + t90 * t303 / 0.2e1) * t226 - 0.6e1 * ((-0.3e1 * t176 + (-t224 + t238) * t165 + t237 * t169 + t266) * t86 - 0.2e1 * (-0.5e1 / 0.3e1 * t176 + (-t169 + t236) * t165 + t171 * (t120 + t234)) * t303) * t297 - 0.6e1 * t209 * t263 - (t162 + ((21 * t169) + t77) * t176 + (t152 + (35 * t167) + t67 + 0.2e1 * t296) * t165 + (t145 + (t143 + t153 - 0.5e1 * t165) * t169 + t171 * (-t160 + t269)) * t91) * t86 + (0.7e1 * t162 + (t147 + t77) * t144 + ((21 * t167) + 0.9e1 * t170 + t67 + 0.6e1 * t296) * t165 + t300) * t303) * t38 + (0.16e2 * (t255 + t226 + (-(8 * t167) + 0.12e2 * t282) * t99 + (-0.12e2 * pkin(7) * t172 + t180 * t315) * t133 - 0.6e1 * t282 + t268) * t293 + 0.24e2 * (t279 * t255 + 0.14e2 * (-0.32e2 / 0.21e2 * (t171 + t165 / 0.4e1 + t169 / 0.4e1 - t160 / 0.8e1) * t221 + 0.5e1 / 0.42e2 * t176 + (0.16e2 / 0.21e2 * t169 + t276) * t165 + t167 / 0.7e1 + t276 * t169 + t170 - 0.3e1 / 0.7e1 * t283 + t159 / 0.42e2) * t297 + t82 * t203 - t267 * t176 + (-0.10e2 / 0.3e1 * t167 + 0.2e1 * t170 - t283 + t79) * t165 + t55 * t285 + ((-0.2e1 / 0.3e1 * t221 + t109 + t277) * t250 + (-0.8e1 / 0.3e1 * (t278 + t310) * t221 + 0.5e1 / 0.18e2 * t176 + (0.4e1 / 0.3e1 * t171 + t122 + t110) * t165 + t170 + 0.2e1 / 0.3e1 * t282 - 0.2e1 / 0.3e1 * t283 - t167 / 0.3e1 + t159 / 0.18e2) * t259) * pkin(7)) * t299 + 0.16e2 * (-0.6e1 * t171 * t165 + t265) * t298 + 0.32e2 * (t286 * t68 + t71 * t90) * t257 + 0.24e2 * (t81 * t203 - t162 + (-t83 + t270) * t176 + (t79 + t176 / 0.6e1 - t159 / 0.6e1 + t266) * t165 + t55 * t171) * t297 + 0.8e1 * t207 * t263 - 0.8e1 * ((t147 + t280) * t176 + (t145 + (t143 + 0.10e2 * t171) * t169 + t213) * t165 + t290) * t221 + t176 ^ 2 + (t142 + t154 + (28 * t169)) * t162 + (t169 * t88 + (70 * t167) + t212) * t176 + (t148 * t212 + t156 * t272 - 0.6e1 * t170 * t160 + t88 * t167 + 0.28e2 * t172 ^ 2 + 0.4e1 * t180 ^ 2) * t165 + t72 * t300) * t63 + (((0.4e1 * t301 + (t86 + t262) * t76 + t93 * t86 + (t111 + t231) * t262) * t253 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t165 + t124 + t108) * t132 * t158 + pkin(1) * t240) * t297 + (-0.8e1 * t219 + ((t110 + t117 + t274) * t86 - (0.8e1 / 0.3e1 * t165 + t234) * t303) * t260) * pkin(7) + t209) * t85) * t38 + (0.32e2 * (t239 + (-0.4e1 * t172 * t247 + t306 + (t307 + t142 + t153) * t169) * t99 + (-t169 + t201 + t310) * t243 + t68 * t285 + t93 * t71) * t294 + 0.8e1 * (t66 + (t286 * t71 + t58) * t251 + (t203 + (t148 + t273) * t165 + t197) * t242 + t207) * t85) * t63) * t74) / ((-0.4e1 * (-t267 * t86 + 0.2e1 * t301 + (0.2e1 * pkin(7) * t246 + t130 * (t165 + t305)) * pkin(1)) * t299 + 0.8e1 * pkin(7) * t219 + ((pkin(3) * t306 + 0.8e1 * t169 * t175) * t132 + 0.4e1 * t172 * t240) * t99 - 0.4e1 * t208 * t263 - (t165 * t275 + t146 + t265 + 0.6e1 * t282) * t86 + (t144 + (t137 + 0.6e1 * t171) * t165 + t84) * t303) * t38 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t221 + 0.4e1 / 0.9e1 * t165 - t160 / 0.9e1 + t278) * t297 + t82 * t68 + t69 * t285 + (t295 + (t108 + t117 + t201) * t87) * t264) * t299 + t66 + (t286 * t69 + t58) * t251 + t218 * t242 + ((t83 + t233) * t165 + t289) * t215 + t162 + (t136 + t230) * t176 + (t148 * t230 + t156 * t271 + t135 + t152) * t165 + t84 * t72) * t63 + ((t253 * t303 + (-0.2e1 * t99 * t248 + (t86 - t303) * t76 + t208) * t258) * t38 + (0.8e1 * (t76 + t254 + t93) * t294 + 0.6e1 * (t269 * t254 + (t68 + t69) * t243 + t218) * t85) * t63) * t74);
t13 = (t191 * t302 - t38 * t195 / 0.4e1) * t284;
t14 = (t38 * t191 + t195 * t302 / 0.4e1) * t284;
t59 = -t129 * t130 - t287;
t60 = -t129 * t133 + t288;
t3 = atan2(t13 * t59 + t14 * t60, -t13 * t60 + t14 * t59);
t1 = sin(t3);
t2 = cos(t3);
t244 = t161 * t291;
t17 = atan2(t38 * t244, t244 * t302);
t15 = sin(t17);
t16 = cos(t17);
t192 = t309 * t193;
t190 = t192 / 0.2e1;
t194 = t195 * t291;
t189 = atan2(t194, t190);
t187 = sin(t189);
t188 = cos(t189);
t31 = -t130 * t188 - t133 * t187;
t32 = t130 * t187 - t133 * t188;
t225 = -t15 * t31 + t32 * t16;
t5 = t15 * t32 + t16 * t31;
t319 = t1 * t5 - t2 * t225;
t318 = t1 * t225 + t2 * t5;
t131 = sin(pkin(9));
t134 = cos(pkin(9));
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t56 = t126 * t133 - t127 * t130;
t57 = t126 * t130 + t127 * t133;
t34 = -t56 * t192 / 0.2e1 + t57 * t194;
t35 = t190 * t57 + t194 * t56;
t25 = t131 * t35 - t134 * t34;
t26 = t131 * t34 + t134 * t35;
t22 = atan2(t25, t26);
t18 = sin(t22);
t20 = cos(t22);
t10 = t18 * t32 + t20 * t31;
t23 = atan2(t25, -t26);
t19 = sin(t23);
t206 = t18 * t31 - t20 * t32;
t21 = cos(t23);
t317 = -t10 * t19 - t206 * t21;
t316 = t10 * t21 - t19 * t206;
t252 = r_base(1) - t87;
t223 = t32 * pkin(3) + t252;
t222 = t32 * pkin(2) + t252;
t216 = r_base(2) - t303;
t205 = t31 * pkin(3) + t216;
t204 = t31 * pkin(2) + t216;
t50 = atan2(t56, -t57);
t49 = atan2(t56, t57);
t48 = cos(t50);
t47 = cos(t49);
t46 = sin(t50);
t45 = sin(t49);
t42 = -t130 * t45 + t133 * t47;
t41 = t130 * t47 + t133 * t45;
t40 = -t126 * t46 + t127 * t48;
t39 = t126 * t48 + t127 * t46;
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t133 + rSges(2,2) * t130 + r_base(1)) + g(2) * (-rSges(2,1) * t130 - rSges(2,2) * t133 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t32 - rSges(3,2) * t31 + t252) + g(2) * (rSges(3,1) * t31 + rSges(3,2) * t32 + t216) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t206 + rSges(4,2) * t10 + t222) + g(2) * (-rSges(4,1) * t10 + rSges(4,2) * t206 + t204) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t225 - rSges(5,2) * t5 + t223) + g(2) * (rSges(5,1) * t5 + rSges(5,2) * t225 + t205) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t127 + rSges(6,1) * t40 - rSges(6,2) * t39 + r_base(1)) + g(2) * (pkin(5) * t126 + rSges(6,1) * t39 + rSges(6,2) * t40 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t206 + rSges(7,2) * t10 + t252) + g(2) * (-rSges(7,1) * t10 + rSges(7,2) * t206 + t216) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t129 + rSges(8,2) * t132 + pkin(7) + r_base(1)) + g(2) * (-rSges(8,1) * t132 + rSges(8,2) * t129 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t42 - rSges(9,2) * t41 + t252) + g(2) * (rSges(9,1) * t41 + rSges(9,2) * t42 + t216) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (t206 * pkin(6) + rSges(10,1) * t317 - rSges(10,2) * t316 + t222) + g(2) * (-t10 * pkin(6) + rSges(10,1) * t316 + rSges(10,2) * t317 + t204) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (t225 * pkin(4) + rSges(11,1) * t319 + t318 * rSges(11,2) + t223) + g(2) * (t5 * pkin(4) - t318 * rSges(11,1) + rSges(11,2) * t319 + t205) + g(3) * (r_base(3) + rSges(11,3)));
U = t4;
