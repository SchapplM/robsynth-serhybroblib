% Calculate potential energy for
% picker2Dm2TE
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
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:34:08
% EndTime: 2020-05-09 10:34:18
% DurationCPUTime: 6.43s
% Computational Cost: add. (40082->486), mult. (117758->609), div. (838->5), fcn. (20852->10), ass. (0->229)
t117 = cos(pkin(9));
t149 = 0.1e1 / pkin(3);
t114 = sin(qJ(1));
t148 = pkin(3) ^ 2;
t153 = (pkin(1) ^ 2);
t155 = pkin(7) ^ 2;
t75 = t153 + t155;
t213 = t148 + t75;
t116 = cos(qJ(1));
t113 = sin(qJ(2));
t69 = pkin(3) * t113;
t291 = t69 + pkin(7);
t220 = t291 * t116;
t115 = cos(qJ(2));
t272 = t114 * t115;
t233 = pkin(3) * t272;
t207 = pkin(1) * t233;
t51 = -0.2e1 * t207;
t296 = 0.2e1 * pkin(7);
t59 = t69 * t296;
t179 = 0.1e1 / (0.2e1 * pkin(1) * t220 + t213 + t51 + t59);
t80 = t113 ^ 2;
t284 = t148 * t80;
t242 = 0.2e1 * t284;
t255 = -t148 + t155;
t185 = t59 + t242 + t255;
t83 = t116 ^ 2;
t181 = t185 * t83;
t196 = -pkin(1) + t233;
t132 = 3 * t153;
t143 = pkin(4) ^ 2;
t256 = t143 - t155;
t214 = -0.2e1 * t148 + t256;
t197 = t132 - t214;
t198 = -0.4e1 * t207;
t215 = -t143 + t75;
t271 = t115 * t116;
t232 = pkin(3) * t271;
t141 = 0.2e1 * pkin(3);
t151 = t153 ^ 2;
t201 = t59 + t215;
t244 = -0.4e1 * t69;
t253 = t153 - t155;
t292 = 2 * t153;
t303 = 0.4e1 * pkin(1);
t29 = sqrt(-0.4e1 * t153 * t181 + 0.4e1 * t253 * t284 + pkin(7) * t215 * t244 - t151 + t214 * t292 - (t155 - (t141 + pkin(4)) * pkin(4)) * (t155 + (t141 - pkin(4)) * pkin(4)) + (-(t51 + t201) * t220 + t201 * t233) * t303);
t294 = 0.4e1 * t148;
t175 = t179 * ((t114 * t291 + t232) * t29 - (t197 + t59 + t198) * t220 + t196 * t59 + t197 * t233 + (-0.2e1 * t181 + t242 - t294 - t215) * pkin(1));
t173 = t175 / 0.2e1;
t295 = 0.2e1 * t148;
t189 = t295 + t201;
t71 = pkin(1) * t116;
t247 = 0.2e1 * t71;
t70 = pkin(3) * t115;
t178 = t179 * ((t220 - t196) * t29 + (t185 * t247 + t189 * t291) * t114 + (t116 * t189 + (0.4e1 * t83 - 0.2e1) * pkin(1) * t291) * t70);
t177 = -t178 / 0.2e1;
t110 = sin(pkin(8));
t111 = cos(pkin(8));
t39 = -t110 * t116 + t111 * t114;
t40 = t110 * t114 + t111 * t116;
t171 = t149 * (t173 * t40 + t177 * t39);
t176 = t178 / 0.2e1;
t20 = (t173 * t39 + t176 * t40) * t149;
t290 = sin(pkin(9));
t16 = t117 * t20 - t171 * t290;
t17 = t117 * t171 + t20 * t290;
t174 = -t175 / 0.2e1;
t25 = (t114 * t174 + t116 * t177) * t149;
t26 = (t114 * t176 + t116 * t174) * t149;
t211 = -t16 * t25 - t26 * t17;
t8 = t16 * t26 - t17 * t25;
t307 = -t16 * t8 + t17 * t211;
t306 = t16 * t211 + t17 * t8;
t172 = t175 / 0.4e1;
t144 = 0.1e1 / pkin(4);
t267 = t144 / pkin(3) ^ 2;
t101 = 0.2e1 / 0.3e1 * t148;
t104 = -t148 / 0.3e1;
t106 = 0.4e1 / 0.3e1 * t153;
t108 = t153 / 0.2e1;
t118 = 15 * t151;
t119 = 15 * t153;
t120 = 10 * t153;
t125 = -0.2e1 * t143;
t126 = -0.5e1 * t143;
t160 = t148 ^ 2;
t127 = 0.5e1 * t160;
t128 = 7 * t151;
t129 = 5 * t151;
t130 = 7 * t153;
t131 = 6 * t153;
t154 = t155 ^ 2;
t135 = 0.3e1 * t154;
t136 = 0.8e1 * t155;
t137 = 0.4e1 * t155;
t139 = 0.2e1 * t155;
t142 = t143 ^ 2;
t159 = pkin(3) * t148;
t145 = t159 ^ 2;
t156 = pkin(1) * t153;
t164 = pkin(7) * t155;
t254 = t151 + t154;
t259 = t139 - t143;
t266 = t155 * t143;
t183 = t259 * t153 + t142 / 0.6e1 + t254 - t266;
t180 = t183 + 0.5e1 / 0.6e1 * t160;
t186 = t155 - t207;
t95 = -t143 / 0.2e1;
t54 = t95 + t213;
t190 = t54 * t198;
t199 = -0.6e1 * t207;
t285 = t142 / 0.2e1 - t160 / 0.2e1;
t204 = t135 - 0.3e1 * t266 + t285;
t138 = 0.3e1 * t155;
t97 = -0.3e1 / 0.2e1 * t143;
t273 = t97 + t138;
t275 = t145 + t75 * ((t97 + t139) * t153 - 0.3e1 / 0.2e1 * t266 + t254 + t285);
t67 = 0.10e2 / 0.3e1 * t153;
t193 = ((t67 + t259) * t148 + t180) * t199 + (t118 + (-0.9e1 * t143 + 0.18e2 * t155) * t153 + t204) * t148 + (t119 + t273) * t160 + t275;
t260 = t132 + t155;
t217 = t148 + t260;
t289 = pkin(1) * t114;
t194 = t217 * t70 - t289 * (0.3e1 * t148 + t75);
t105 = -0.2e1 / 0.3e1 * t148;
t96 = -0.2e1 / 0.3e1 * t143;
t222 = t96 + t75;
t261 = t120 + t139;
t264 = t105 + t155;
t73 = -t143 - t148;
t61 = t138 + t73;
t277 = t61 * t153;
t195 = -t289 * (t127 + ((5 * t153) + t61) * t295 + (t105 + t222) * t75) + (t160 + (t105 + t96 + t261) * t148 + t129 + 0.2e1 * t277 + t155 * (t96 + t264)) * t70;
t258 = t142 - t160;
t202 = 0.6e1 * t154 - 0.6e1 * t266 + t258;
t226 = t101 + t139 + t96;
t274 = t160 + (t101 + t222) * t75;
t100 = 0.4e1 / 0.3e1 * t148;
t94 = -t143 / 0.3e1;
t221 = t94 + t75;
t52 = t100 + t221;
t203 = t52 * t198 + t274 + (t131 + t226) * t148;
t82 = t116 * t83;
t280 = t156 * t82;
t205 = t280 * t70;
t283 = t151 * t83 ^ 2;
t206 = t283 * t70;
t210 = 0.20e2 / 0.3e1 * t153;
t236 = 0.16e2 * t280;
t212 = pkin(7) * t236;
t257 = -t143 + t148;
t216 = t138 + t257;
t243 = pkin(7) * t280;
t218 = 0.8e1 * t243;
t270 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t219 = t114 * t270;
t102 = t148 / 0.3e1;
t92 = -t143 / 0.6e1;
t223 = t102 + t92 + t155;
t224 = t102 + t143 / 0.3e1 + t139;
t225 = t101 + 0.2e1 / 0.3e1 * t143 + t137;
t227 = t100 + 0.4e1 / 0.3e1 * t143 - 0.2e1 * t155;
t245 = 0.6e1 * t71;
t228 = pkin(7) * t245;
t246 = 0.4e1 * t71;
t229 = pkin(7) * t246;
t231 = -t289 / 0.2e1;
t234 = t153 * t70;
t282 = t153 * t83;
t237 = 0.12e2 * t282;
t279 = t159 * t113 * t80;
t239 = -0.8e1 * t279;
t240 = 0.4e1 * t282;
t241 = 0.8e1 * t283;
t248 = 0.2e1 * t289;
t249 = pkin(7) * t71;
t250 = 0.4e1 * pkin(7);
t251 = t154 + t160;
t252 = t154 - t151;
t262 = t108 + t155;
t263 = t153 / 0.3e1 + t155;
t265 = t155 * t153;
t269 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t276 = 0.4e1 / 0.7e1 * t155 - t143 / 0.7e1;
t278 = t160 * t80 ^ 2;
t281 = t155 * t73;
t68 = t75 ^ 2;
t286 = t68 * (-t148 + t215);
t287 = (-t114 * t156 + t234) * t83;
t293 = 4 * t151;
t93 = -t143 / 0.4e1;
t298 = t148 / 0.2e1 + t93;
t38 = -t160 / 0.6e1 + t183;
t65 = t104 + t155;
t41 = t65 * t51;
t58 = t71 + pkin(7);
t46 = t69 + t58;
t74 = -0.3e1 * t148 + t155;
t49 = t74 * t218;
t50 = 0.10e2 * t277;
t55 = -t143 + t213;
t60 = pkin(7) * t247;
t63 = (t137 + t143) * t153;
t66 = -t153 / 0.3e1 + t155;
t72 = -0.30e2 * t143 + 0.60e2 * t155;
t77 = -(3 * t153) + t155;
t288 = ((-0.24e2 * (0.4e1 / 0.3e1 * t282 + t60 + t66) * t278 * t289 - 0.12e2 * (-0.8e1 / 0.3e1 * t206 + ((t106 + t223) * t70 - (0.7e1 / 0.6e1 * t148 + t92 + t262) * t289) * t240 + (-t148 * t253 - 0.5e1 / 0.3e1 * t151 + t224 * t153 + t155 * (t94 + t65)) * t70 + (-t160 + (-t210 + t225) * t148 - (3 * t151) + t227 * t153 + t154) * t231 + (-t114 * t151 * t82 + ((t153 + t223) * t70 + (t295 - t253) * t231) * t71) * t250) * t284 + 0.24e2 * t65 * t206 + ((t155 + 0.5e1 / 0.2e1 * t148 + 0.3e1 / 0.2e1 * t153 + t95) * t70 + t74 * t289 / 0.2e1) * t212 - 0.6e1 * ((-0.3e1 * t160 + (-t210 + t227) * t148 + t225 * t153 + t252) * t70 - 0.2e1 * (-0.5e1 / 0.3e1 * t160 + (-t153 + t224) * t148 + t155 * (t104 + t221)) * t289) * t282 - 0.6e1 * t195 * t249 - (t145 + ((21 * t153) + t61) * t160 + (t135 + (35 * t151) + t50 + 0.2e1 * t281) * t148 + (t128 + (t126 + t136 - 0.5e1 * t148) * t153 + t155 * (-t143 + t255)) * t75) * t70 + (0.7e1 * t145 + (t130 + t61) * t127 + ((21 * t151) + 0.9e1 * t154 + t50 + 0.6e1 * t281) * t148 + t286) * t289) * t29 + (0.16e2 * (t241 + t212 + (-(8 * t151) + 0.12e2 * t265) * t83 + (-0.12e2 * pkin(7) * t156 + t164 * t303) * t116 - 0.6e1 * t265 + t254) * t278 + 0.24e2 * (t264 * t241 + 0.14e2 * (-0.32e2 / 0.21e2 * (t155 + t148 / 0.4e1 + t153 / 0.4e1 - t143 / 0.8e1) * t207 + 0.5e1 / 0.42e2 * t160 + (0.16e2 / 0.21e2 * t153 + t276) * t148 + t151 / 0.7e1 + t276 * t153 + t154 - 0.3e1 / 0.7e1 * t266 + t142 / 0.42e2) * t282 + t66 * t190 - t253 * t160 + (-0.10e2 / 0.3e1 * t151 + 0.2e1 * t154 - t266 + t63) * t148 + t38 * t269 + ((-0.2e1 / 0.3e1 * t207 + t93 + t262) * t236 + (-0.8e1 / 0.3e1 * (t263 + t298) * t207 + 0.5e1 / 0.18e2 * t160 + (0.4e1 / 0.3e1 * t155 + t106 + t94) * t148 + t154 + 0.2e1 / 0.3e1 * t265 - 0.2e1 / 0.3e1 * t266 - t151 / 0.3e1 + t142 / 0.18e2) * t245) * pkin(7)) * t284 + 0.16e2 * (-0.6e1 * t155 * t148 + t251) * t283 + 0.32e2 * (t270 * t51 + t54 * t74) * t243 + 0.24e2 * (t65 * t190 - t145 + (-t67 + t256) * t160 + (t63 + t160 / 0.6e1 - t142 / 0.6e1 + t252) * t148 + t38 * t155) * t282 + 0.8e1 * t193 * t249 - 0.8e1 * ((t130 + t273) * t160 + (t128 + (t126 + 0.10e2 * t155) * t153 + t204) * t148 + t275) * t207 + t160 ^ 2 + (t125 + t137 + (28 * t153)) * t145 + (t153 * t72 + (70 * t151) + t202) * t160 + (t202 * t131 + t258 * t139 - 0.6e1 * t154 * t143 + t72 * t151 + 0.28e2 * t156 ^ 2 + 0.4e1 * t164 ^ 2) * t148 + t55 * t286) * t46 + (((0.4e1 * t287 + (t70 + t248) * t60 + t77 * t70 + (t95 + t217) * t248) * t239 - 0.6e1 * (-0.4e1 * (0.2e1 * (0.5e1 / 0.6e1 * t148 + t108 + t92) * t70 + pkin(1) * t219) * t282 + (-0.8e1 * t205 + ((t101 + t94 + t260) * t70 - (0.8e1 / 0.3e1 * t148 + t221) * t289) * t246) * pkin(7) + t195) * t69) * t29 + (0.32e2 * (t218 + (-0.4e1 * t156 * t233 + t293 + (t294 + t125 + t136) * t153) * t83 + (-t153 + t186 + t298) * t229 + t51 * t269 + t77 * t54) * t279 + 0.8e1 * (t49 + (t270 * t54 + t41) * t237 + (t190 + (t131 + t259) * t148 + t180) * t228 + t193) * t69) * t46) * t58) / ((-0.4e1 * (-t253 * t70 + 0.2e1 * t287 + (t232 * t296 + t114 * (t148 + t292)) * pkin(1)) * t284 + 0.8e1 * pkin(7) * t205 + ((pkin(3) * t293 + 0.8e1 * t153 * t159) * t115 + 0.4e1 * t156 * t219) * t83 - 0.4e1 * t194 * t249 - (t148 * t261 + t129 + t251 + 0.6e1 * t265) * t70 + (t127 + (t120 + 0.6e1 * t155) * t148 + t68) * t289) * t29 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t207 + 0.4e1 / 0.9e1 * t148 - t143 / 0.9e1 + t263) * t282 + t66 * t51 + t52 * t269 + (t280 + (t101 + t186 + t92) * t71) * t250) * t284 + t49 + (t270 * t52 + t41) * t237 + t203 * t228 + ((t67 + t226) * t148 + t274) * t199 + t145 + (t119 + t216) * t160 + (t216 * t131 + t257 * t139 + t118 + t135) * t148 + t68 * t55) * t46 + ((t239 * t289 + (-0.2e1 * t83 * t234 + (t70 - t289) * t60 + t194) * t244) * t29 + (0.8e1 * (t60 + t240 + t77) * t279 + 0.6e1 * (t255 * t240 + (t51 + t52) * t229 + t203) * t69) * t46) * t58);
t11 = (t172 * t288 - t29 * t178 / 0.4e1) * t267;
t12 = (t29 * t172 + t178 * t288 / 0.4e1) * t267;
t42 = -t113 * t114 - t271;
t43 = -t113 * t116 + t272;
t1 = t11 * t42 + t12 * t43;
t230 = t288 / 0.2e1;
t268 = t144 * t149;
t184 = (-t25 * t29 / 0.2e1 + t26 * t230) * t268;
t2 = t11 * t43 - t12 * t42;
t297 = (t25 * t230 + t26 * t29 / 0.2e1) * t268;
t305 = t1 * t297 + t184 * t2;
t304 = -t1 * t184 + t2 * t297;
t238 = r_base(1) - t71;
t209 = t26 * pkin(3) + t238;
t208 = t26 * pkin(2) + t238;
t200 = r_base(2) - t289;
t192 = t25 * pkin(3) + t200;
t191 = t25 * pkin(2) + t200;
t188 = t114 * t40 - t116 * t39;
t187 = t114 * t39 + t116 * t40;
t33 = t110 * t39 - t111 * t40;
t32 = -t110 * t40 - t111 * t39;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t116 + rSges(2,2) * t114 + r_base(1)) + g(2) * (-rSges(2,1) * t114 - rSges(2,2) * t116 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t26 - rSges(3,2) * t25 + t238) + g(2) * (rSges(3,1) * t25 + rSges(3,2) * t26 + t200) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t211 - rSges(4,2) * t8 + t208) + g(2) * (rSges(4,1) * t8 + rSges(4,2) * t211 + t191) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t184 - rSges(5,2) * t297 + t209) + g(2) * (rSges(5,1) * t297 + rSges(5,2) * t184 + t192) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t111 + rSges(6,1) * t33 - rSges(6,2) * t32 + r_base(1)) + g(2) * (pkin(5) * t110 + rSges(6,1) * t32 + rSges(6,2) * t33 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t211 - rSges(7,2) * t8 + t238) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t211 + t200) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t113 + rSges(8,2) * t115 + pkin(7) + r_base(1)) + g(2) * (-rSges(8,1) * t115 + rSges(8,2) * t113 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t187 - rSges(9,2) * t188 + t238) + g(2) * (rSges(9,1) * t188 + rSges(9,2) * t187 + t200) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (t211 * pkin(6) + rSges(10,1) * t307 - t306 * rSges(10,2) + t208) + g(2) * (t8 * pkin(6) + t306 * rSges(10,1) + rSges(10,2) * t307 + t191) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (t184 * pkin(4) + rSges(11,1) * t305 - t304 * rSges(11,2) + t209) + g(2) * (t297 * pkin(4) + rSges(11,1) * t304 + rSges(11,2) * t305 + t192) + g(3) * (r_base(3) + rSges(11,3)));
U = t3;
