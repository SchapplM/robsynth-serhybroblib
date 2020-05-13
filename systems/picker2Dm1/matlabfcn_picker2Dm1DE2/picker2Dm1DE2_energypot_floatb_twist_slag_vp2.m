% Calculate potential energy for
% picker2Dm1DE2
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
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:21:16
% EndTime: 2020-05-10 20:21:30
% DurationCPUTime: 7.72s
% Computational Cost: add. (65214->490), mult. (176728->601), div. (2626->12), fcn. (45644->38), ass. (0->268)
t332 = 4 * pkin(1);
t308 = 0.1e1 / pkin(2) / 0.2e1;
t170 = 0.1e1 / t308;
t177 = pkin(4) ^ 2;
t119 = -t177 / 0.4e1;
t182 = pkin(3) ^ 2;
t331 = t119 + t182 / 0.2e1;
t330 = 2 * pkin(7);
t136 = sin(pkin(8));
t137 = cos(pkin(8));
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t61 = t136 * t143 - t137 * t140;
t329 = 0.2e1 * t61;
t109 = t143 ^ 2;
t328 = -0.2e1 * t109;
t155 = 0.2e1 * t182;
t327 = 0.4e1 * t182;
t187 = pkin(1) ^ 2;
t185 = t187 ^ 2;
t326 = 4 * t185;
t325 = 2 * t187;
t159 = 6 * t187;
t190 = pkin(7) ^ 2;
t167 = 2 * t190;
t195 = t182 ^ 2;
t154 = 0.5e1 * t195;
t324 = pkin(5) * m(6);
t171 = 2 * pkin(1);
t174 = pkin(5) ^ 2;
t62 = -t136 * t140 - t137 * t143;
t322 = pkin(5) * t62;
t271 = pkin(1) * t322;
t306 = t174 - 0.2e1 * t271;
t46 = sqrt(-(-(t171 + pkin(5)) * pkin(5) + t306) * (pkin(5) * (t171 - pkin(5)) + t306));
t317 = t46 * t61;
t53 = t174 - t271;
t57 = -pkin(1) + t322;
t44 = -pkin(5) * t317 - 0.2e1 * t53 * t57;
t323 = t44 / 0.4e1;
t97 = pkin(1) * t143;
t84 = t97 + pkin(7);
t141 = sin(pkin(9));
t144 = cos(pkin(9));
t183 = 0.1e1 / pkin(3);
t101 = t187 + t190;
t234 = t182 + t101;
t269 = 0.2e1 * t97;
t142 = cos(qJ(2));
t300 = t140 * t142;
t259 = pkin(3) * t300;
t227 = pkin(1) * t259;
t75 = -0.2e1 * t227;
t139 = sin(qJ(2));
t95 = pkin(3) * t139;
t83 = t95 + pkin(7);
t85 = t95 * t330;
t49 = 0.1e1 / (t269 * t83 + t234 + t75 + t85);
t313 = t183 * t49;
t52 = 0.1e1 / (t187 + t306);
t314 = 0.1e1 / pkin(5) * t52;
t223 = t313 * t314;
t216 = -pkin(1) + t259;
t160 = 3 * t187;
t280 = t177 - t190;
t217 = t155 + t160 - t280;
t220 = -0.4e1 * t227;
t235 = -t177 + t101;
t299 = t142 * t143;
t258 = pkin(3) * t299;
t106 = t139 ^ 2;
t301 = t106 * t182;
t263 = 0.2e1 * t301;
t315 = t143 * t83;
t169 = 0.2e1 * pkin(3);
t222 = t85 + t235;
t295 = t187 * t109;
t261 = -0.4e1 * t295;
t266 = -0.4e1 * t95;
t277 = t187 - t190;
t279 = -t182 + t190;
t63 = t85 + t263 + t279;
t43 = sqrt(t63 * t261 + 0.4e1 * t277 * t301 + pkin(7) * t235 * t266 - t185 + (-0.2e1 * t182 + t280) * t325 - (t190 - (t169 + pkin(4)) * pkin(4)) * (t190 + (t169 - pkin(4)) * pkin(4)) + (-(t75 + t222) * t315 + t222 * t259) * t332);
t41 = (t140 * t83 + t258) * t43 - (t217 + t85 + t220) * t315 + t216 * t85 + t217 * t259 + (t63 * t328 - t235 + t263 - t327) * pkin(1);
t68 = t155 + t222;
t96 = pkin(3) * t142;
t42 = (-t216 + t315) * t43 + (t269 * t63 + t68 * t83) * t140 + (t143 * t68 + (0.4e1 * t109 - 0.2e1) * t83 * pkin(1)) * t96;
t45 = pkin(5) * t329 * t53 - t46 * t57;
t28 = (t41 * t323 + t42 * t45 / 0.4e1) * t223;
t29 = (-t41 * t45 / 0.4e1 + t42 * t323) * t223;
t26 = t141 * t29 + t144 * t28;
t321 = pkin(6) * t26;
t27 = t141 * t28 - t144 * t29;
t320 = pkin(6) * t27;
t319 = pkin(1) * t140;
t100 = -0.3e1 * t182 + t190;
t103 = -3 * t187 + t190;
t108 = t143 * t109;
t118 = -t177 / 0.6e1;
t120 = -t177 / 0.3e1;
t121 = -t177 / 0.2e1;
t127 = 0.2e1 / 0.3e1 * t182;
t130 = -t182 / 0.3e1;
t132 = 0.4e1 / 0.3e1 * t187;
t134 = t187 / 0.2e1;
t145 = 15 * t185;
t146 = 15 * t187;
t147 = 10 * t187;
t152 = -0.2e1 * t177;
t153 = -0.5e1 * t177;
t156 = 7 * t185;
t157 = 5 * t185;
t158 = 7 * t187;
t189 = t190 ^ 2;
t163 = 3 * t189;
t164 = 8 * t190;
t165 = 4 * t190;
t176 = t177 ^ 2;
t194 = pkin(3) * t182;
t179 = t194 ^ 2;
t191 = pkin(1) * t187;
t199 = pkin(7) * t190;
t278 = t185 + t189;
t283 = t167 - t177;
t294 = t190 * t177;
t207 = t283 * t187 + t176 / 0.6e1 + t278 - t294;
t206 = 0.5e1 / 0.6e1 * t195 + t207;
t209 = t190 - t227;
t78 = t121 + t234;
t210 = t78 * t220;
t291 = t176 / 0.2e1 - t195 / 0.2e1;
t219 = -0.3e1 * t294 + t163 + t291;
t221 = -0.6e1 * t227;
t123 = -0.3e1 / 0.2e1 * t177;
t166 = 3 * t190;
t290 = t123 + t166;
t305 = t179 + t101 * ((t123 + t167) * t187 - 0.3e1 / 0.2e1 * t294 + t278 + t291);
t93 = 0.10e2 / 0.3e1 * t187;
t213 = ((t93 + t283) * t182 + t206) * t221 + (t145 + (-0.9e1 * t177 + (18 * t190)) * t187 + t219) * t182 + (t146 + t290) * t195 + t305;
t284 = t160 + t190;
t237 = t182 + t284;
t214 = t237 * t96 - t319 * (0.3e1 * t182 + t101);
t122 = -0.2e1 / 0.3e1 * t177;
t131 = -0.2e1 / 0.3e1 * t182;
t238 = t122 + t101;
t285 = t147 + t167;
t289 = t131 + t190;
t99 = -t177 - t182;
t87 = t166 + t99;
t309 = t87 * t187;
t215 = -t319 * (t154 + ((5 * t187) + t87) * t155 + (t131 + t238) * t101) + (t195 + (t122 + t131 + t285) * t182 + t157 + 0.2e1 * t309 + t190 * (t122 + t289)) * t96;
t282 = t176 - t195;
t218 = -0.6e1 * t294 + (6 * t189) + t282;
t256 = t191 * t96;
t224 = t108 * t256;
t296 = t185 * t109 ^ 2;
t225 = t296 * t96;
t239 = t122 + t127 + t167;
t304 = t195 + (t127 + t238) * t101;
t126 = 0.4e1 / 0.3e1 * t182;
t240 = t120 + t101;
t76 = t126 + t240;
t226 = t76 * t220 + t304 + (t159 + t239) * t182;
t292 = t191 * t108;
t254 = 0.16e2 * t292;
t228 = pkin(7) * t254;
t229 = 0.20e2 / 0.3e1 * t187;
t265 = pkin(7) * t292;
t230 = 0.8e1 * t265;
t281 = -t177 + t182;
t236 = t166 + t281;
t128 = t182 / 0.3e1;
t241 = t118 + t128 + t190;
t242 = t177 / 0.3e1 + t128 + t167;
t243 = 0.2e1 / 0.3e1 * t177 + t127 + t165;
t244 = 0.4e1 / 0.3e1 * t177 + t126 - (2 * t190);
t298 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t248 = t140 * t298;
t267 = 0.6e1 * t97;
t249 = pkin(7) * t267;
t268 = 0.4e1 * t97;
t250 = pkin(7) * t268;
t252 = -t319 / 0.2e1;
t255 = 0.12e2 * t295;
t257 = t187 * t96;
t260 = 0.4e1 * t295;
t262 = 0.8e1 * t296;
t302 = t139 * t106 * t194;
t264 = -0.8e1 * t302;
t270 = 0.2e1 * t319;
t272 = pkin(7) * t97;
t274 = 4 * pkin(7);
t275 = t189 + t195;
t276 = t189 - t185;
t286 = 0.4e1 / 0.7e1 * t190 - t177 / 0.7e1;
t287 = t134 + t190;
t288 = t187 / 0.3e1 + t190;
t293 = t190 * t187;
t297 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t303 = t106 ^ 2 * t195;
t310 = (-t140 * t191 + t257) * t109;
t311 = t190 * t99;
t94 = t101 ^ 2;
t316 = t94 * (-t182 + t235);
t60 = -t195 / 0.6e1 + t207;
t91 = t130 + t190;
t64 = t91 * t75;
t70 = t95 + t84;
t73 = t100 * t230;
t74 = 0.10e2 * t309;
t79 = -t177 + t234;
t86 = pkin(7) * t269;
t89 = (t165 + t177) * t187;
t92 = -t187 / 0.3e1 + t190;
t98 = -0.30e2 * t177 + (60 * t190);
t318 = ((-0.24e2 * (0.4e1 / 0.3e1 * t295 + t86 + t92) * t303 * t319 - 0.12e2 * (-0.8e1 / 0.3e1 * t225 + ((t132 + t241) * t96 - (0.7e1 / 0.6e1 * t182 + t118 + t287) * t319) * t260 + (-t182 * t277 - 0.5e1 / 0.3e1 * t185 + t242 * t187 + t190 * (t120 + t91)) * t96 + (-t195 + (-t229 + t243) * t182 - (3 * t185) + t244 * t187 + t189) * t252 + (-t140 * t185 * t108 + ((t187 + t241) * t96 + (t155 - t277) * t252) * t97) * t274) * t301 + 0.24e2 * t91 * t225 + ((t190 + 0.5e1 / 0.2e1 * t182 + 0.3e1 / 0.2e1 * t187 + t121) * t96 + t100 * t319 / 0.2e1) * t228 - 0.6e1 * ((-0.3e1 * t195 + (-t229 + t244) * t182 + t243 * t187 + t276) * t96 - 0.2e1 * (-0.5e1 / 0.3e1 * t195 + (-t187 + t242) * t182 + t190 * (t130 + t240)) * t319) * t295 - 0.6e1 * t215 * t272 - (t179 + ((21 * t187) + t87) * t195 + (t163 + (35 * t185) + t74 + 0.2e1 * t311) * t182 + (t156 + (t153 + t164 - 0.5e1 * t182) * t187 + t190 * (-t177 + t279)) * t101) * t96 + (0.7e1 * t179 + (t158 + t87) * t154 + ((21 * t185) + (9 * t189) + t74 + 0.6e1 * t311) * t182 + t316) * t319) * t43 + (0.16e2 * (t262 + t228 + (-8 * t185 + 12 * t293) * t109 + (-12 * pkin(7) * t191 + t199 * t332) * t143 - (6 * t293) + t278) * t303 + 0.24e2 * (t289 * t262 + 0.14e2 * (-0.32e2 / 0.21e2 * (t190 + t182 / 0.4e1 + t187 / 0.4e1 - t177 / 0.8e1) * t227 + 0.5e1 / 0.42e2 * t195 + (0.16e2 / 0.21e2 * t187 + t286) * t182 + t185 / 0.7e1 + t286 * t187 + t189 - 0.3e1 / 0.7e1 * t294 + t176 / 0.42e2) * t295 + t92 * t210 - t277 * t195 + (-0.10e2 / 0.3e1 * t185 + (2 * t189) - t294 + t89) * t182 + t60 * t297 + ((-0.2e1 / 0.3e1 * t227 + t119 + t287) * t254 + (-0.8e1 / 0.3e1 * (t288 + t331) * t227 + 0.5e1 / 0.18e2 * t195 + (0.4e1 / 0.3e1 * t190 + t132 + t120) * t182 + t189 + 0.2e1 / 0.3e1 * t293 - 0.2e1 / 0.3e1 * t294 - t185 / 0.3e1 + t176 / 0.18e2) * t267) * pkin(7)) * t301 + 0.16e2 * (-0.6e1 * t190 * t182 + t275) * t296 + 0.32e2 * (t100 * t78 + t298 * t75) * t265 + 0.24e2 * (t91 * t210 - t179 + (-t93 + t280) * t195 + (t89 + t195 / 0.6e1 - t176 / 0.6e1 + t276) * t182 + t60 * t190) * t295 + 0.8e1 * t213 * t272 - 0.8e1 * ((t158 + t290) * t195 + (t156 + (t153 + (10 * t190)) * t187 + t219) * t182 + t305) * t227 + t195 ^ 2 + (t152 + t165 + (28 * t187)) * t179 + (t187 * t98 + (70 * t185) + t218) * t195 + (t159 * t218 + t167 * t282 - 0.6e1 * t189 * t177 + t98 * t185 + (28 * t191 ^ 2) + (4 * t199 ^ 2)) * t182 + t79 * t316) * t70 + (((0.4e1 * t310 + (t96 + t270) * t86 + t103 * t96 + (t121 + t237) * t270) * t264 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t182 + t134 + t118) * t96 + pkin(1) * t248) * t261 + (-0.8e1 * t224 + ((t120 + t127 + t284) * t96 - (0.8e1 / 0.3e1 * t182 + t240) * t319) * t268) * pkin(7) + t215) * t95) * t43 + (0.32e2 * (t230 + (-0.4e1 * t140 * t256 + t326 + (t327 + t152 + t164) * t187) * t109 + (-t187 + t209 + t331) * t250 + t75 * t297 + t103 * t78) * t302 + 0.8e1 * (t73 + (t298 * t78 + t64) * t255 + (t210 + (t159 + t283) * t182 + t206) * t249 + t213) * t95) * t70) * t84) / ((-0.4e1 * (-t277 * t96 + 0.2e1 * t310 + (t258 * t330 + t140 * (t182 + t325)) * pkin(1)) * t301 + 0.8e1 * pkin(7) * t224 + ((pkin(3) * t326 + 0.8e1 * t187 * t194) * t142 + 0.4e1 * t191 * t248) * t109 - 0.4e1 * t214 * t272 - (t182 * t285 + t157 + t275 + (6 * t293)) * t96 + (t154 + (t147 + 6 * t190) * t182 + t94) * t319) * t43 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t227 + 0.4e1 / 0.9e1 * t182 - t177 / 0.9e1 + t288) * t295 + t92 * t75 + t76 * t297 + (t292 + (t118 + t127 + t209) * t97) * t274) * t301 + t73 + (t298 * t76 + t64) * t255 + t226 * t249 + ((t93 + t239) * t182 + t304) * t221 + t179 + (t146 + t236) * t195 + (t159 * t236 + t167 * t281 + t145 + t163) * t182 + t94 * t79) * t70 + ((t264 * t319 + (t257 * t328 + (t96 - t319) * t86 + t214) * t266) * t43 + (0.8e1 * (t86 + t260 + t103) * t302 + 0.6e1 * (t279 * t260 + (t75 + t76) * t250 + t226) * t95) * t70) * t84);
t312 = 0.1e1 / pkin(1) * t52;
t245 = t313 / 0.2e1;
t33 = qJ(1) + atan2(t42 * t245, t41 * t245);
t172 = pkin(6) ^ 2;
t25 = t170 * t321;
t307 = t172 + t25;
t273 = -m(9) - m(7) - m(3);
t14 = sqrt(-(-(t170 + pkin(6)) * pkin(6) + t307) * (pkin(6) * (t170 - pkin(6)) + t307));
t23 = t25 + 0.2e1 * t172;
t24 = -pkin(2) - t321;
t173 = 0.1e1 / pkin(6);
t247 = t173 / (pkin(2) ^ 2 + t307) / 0.2e1;
t6 = atan2((-t14 * t24 + t23 * t320) * t247, (-t14 * t320 - t23 * t24) * t247) + t33;
t178 = 0.1e1 / pkin(4);
t231 = t178 * t183 / 0.2e1;
t17 = atan2(t43 * t231, t231 * t318) + t33;
t253 = t178 / pkin(3) ^ 2 * t49;
t251 = t318 / 0.4e1;
t246 = t314 / 0.2e1;
t32 = cos(t33);
t233 = -pkin(2) * t32 - t97;
t232 = -pkin(3) * t32 - t97;
t31 = sin(t33);
t212 = -pkin(2) * t31 - t319;
t211 = -pkin(3) * t31 - t319;
t208 = -m(8) - m(2) - m(10) - m(5) - m(11) - m(4) - m(6) - m(1);
t66 = -t139 * t143 + t300;
t65 = -t139 * t140 - t299;
t58 = -pkin(1) * t62 + pkin(5);
t54 = t187 - t271;
t40 = qJ(1) + atan2(t45 * t246, t44 * t246);
t39 = pkin(8) + atan2((pkin(1) * t329 * t54 + t46 * t58) * t312 / 0.2e1, -(-pkin(1) * t317 + 0.2e1 * t54 * t58) * t312 / 0.2e1);
t38 = cos(t40);
t37 = sin(t40);
t36 = cos(t39);
t35 = sin(t39);
t19 = (t41 * t43 / 0.4e1 + t42 * t251) * t253;
t18 = (t41 * t251 - t42 * t43 / 0.4e1) * t253;
t16 = cos(t17);
t15 = sin(t17);
t13 = atan2(t27, t26) + t33;
t12 = cos(t13);
t11 = sin(t13);
t10 = atan2(t18 * t65 + t19 * t66, -t18 * t66 + t19 * t65) + t17;
t9 = cos(t10);
t8 = sin(t10);
t5 = cos(t6);
t4 = sin(t6);
t3 = atan2(t173 * t14 * t308, -t26) + t6;
t2 = cos(t3);
t1 = sin(t3);
t7 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t208 + t273) * r_base(3)) * g(3) + (-t8 * mrSges(11,1) - t4 * mrSges(4,1) + t16 * mrSges(5,2) + t15 * mrSges(5,1) - t9 * mrSges(11,2) - t136 * t324 + mrSges(8,1) * t142 - mrSges(8,2) * t139 - t5 * mrSges(4,2) + t208 * r_base(2) - m(11) * (-pkin(4) * t15 + t211) - m(5) * t211 - m(10) * (pkin(6) * t4 + t212) - m(4) * t212 - mrSges(9,1) * t37 - mrSges(9,2) * t38 - mrSges(7,1) * t11 - mrSges(7,2) * t12 + mrSges(3,1) * t31 + mrSges(3,2) * t32 + mrSges(2,1) * t140 + mrSges(2,2) * t143 + t273 * (r_base(2) - t319) + t2 * mrSges(10,2) - t36 * mrSges(6,2) + t1 * mrSges(10,1) - t35 * mrSges(6,1) - mrSges(1,2)) * g(2) + (t8 * mrSges(11,2) - mrSges(7,1) * t12 + mrSges(7,2) * t11 + mrSges(3,1) * t32 - mrSges(3,2) * t31 + t4 * mrSges(4,2) + t16 * mrSges(5,1) - t15 * mrSges(5,2) - t9 * mrSges(11,1) - t137 * t324 + mrSges(2,1) * t143 - t5 * mrSges(4,1) - m(11) * (-pkin(4) * t16 + t232) - m(5) * t232 - m(10) * (pkin(6) * t5 + t233) - m(4) * t233 + t208 * r_base(1) - (m(8) * pkin(7)) + t273 * (-t97 + r_base(1)) - mrSges(9,1) * t38 + mrSges(9,2) * t37 - mrSges(2,2) * t140 + t2 * mrSges(10,1) - t142 * mrSges(8,2) - t139 * mrSges(8,1) - t1 * mrSges(10,2) + t35 * mrSges(6,2) - t36 * mrSges(6,1) - mrSges(1,1)) * g(1);
U = t7;
