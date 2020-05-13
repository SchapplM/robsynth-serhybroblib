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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2DE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2DE1_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 14:28:05
% EndTime: 2020-05-09 14:28:18
% DurationCPUTime: 8.44s
% Computational Cost: add. (79926->465), mult. (235178->565), div. (1676->6), fcn. (42008->31), ass. (0->246)
t132 = sin(qJ(1));
t167 = pkin(3) ^ 2;
t171 = (pkin(1) ^ 2);
t173 = pkin(7) ^ 2;
t93 = t171 + t173;
t231 = t167 + t93;
t135 = cos(qJ(1));
t131 = sin(qJ(2));
t87 = pkin(3) * t131;
t309 = t87 + pkin(7);
t244 = t309 * t135;
t134 = cos(qJ(2));
t293 = t132 * t134;
t252 = pkin(3) * t293;
t224 = pkin(1) * t252;
t68 = -0.2e1 * t224;
t315 = 0.1e1 / pkin(3);
t296 = t315 / 0.2e1;
t160 = 0.1e1 / t296;
t77 = pkin(7) * t131 * t160;
t198 = 0.1e1 / (0.2e1 * pkin(1) * t244 + t231 + t68 + t77);
t101 = t135 ^ 2;
t98 = t131 ^ 2;
t303 = t167 * t98;
t259 = 0.2e1 * t303;
t272 = -t167 + t173;
t205 = t77 + t259 + t272;
t201 = t205 * t101;
t214 = -pkin(1) + t252;
t151 = 3 * t171;
t162 = pkin(4) ^ 2;
t273 = t162 - t173;
t232 = -0.2e1 * t167 + t273;
t215 = t151 - t232;
t218 = -0.4e1 * t224;
t233 = -t162 + t93;
t292 = t134 * t135;
t251 = pkin(3) * t292;
t313 = 0.4e1 * t167;
t169 = t171 ^ 2;
t220 = t77 + t233;
t260 = -0.4e1 * t87;
t270 = t171 - t173;
t311 = 2 * t171;
t320 = 0.4e1 * pkin(1);
t38 = sqrt(-0.4e1 * t171 * t201 + 0.4e1 * t270 * t303 + pkin(7) * t233 * t260 - t169 + t232 * t311 - (t173 - (t160 + pkin(4)) * pkin(4)) * (t173 + (t160 - pkin(4)) * pkin(4)) + (-(t68 + t220) * t244 + t220 * t252) * t320);
t195 = t198 * ((t132 * t309 + t251) * t38 - (t215 + t77 + t218) * t244 + t214 * t77 + t215 * t252 + (-0.2e1 * t201 + t259 - t313 - t233) * pkin(1));
t193 = t195 / 0.4e1;
t314 = 0.2e1 * t167;
t207 = t314 + t220;
t89 = pkin(1) * t135;
t263 = 0.2e1 * t89;
t88 = pkin(3) * t134;
t197 = t198 * ((t244 - t214) * t38 + (t205 * t263 + t207 * t309) * t132 + (t135 * t207 + (0.4e1 * t101 - 0.2e1) * pkin(1) * t309) * t88);
t163 = 0.1e1 / pkin(4);
t289 = t163 / pkin(3) ^ 2;
t100 = t135 * t101;
t110 = -t162 / 0.6e1;
t111 = -t162 / 0.4e1;
t112 = -t162 / 0.3e1;
t113 = -t162 / 0.2e1;
t119 = 0.2e1 / 0.3e1 * t167;
t122 = -t167 / 0.3e1;
t124 = 0.4e1 / 0.3e1 * t171;
t126 = t171 / 0.2e1;
t137 = 15 * t169;
t138 = 15 * t171;
t139 = 10 * t171;
t144 = -0.2e1 * t162;
t145 = -0.5e1 * t162;
t178 = t167 ^ 2;
t146 = 0.5e1 * t178;
t147 = 7 * t169;
t148 = 5 * t169;
t149 = 7 * t171;
t150 = 6 * t171;
t172 = t173 ^ 2;
t154 = 0.3e1 * t172;
t155 = 0.8e1 * t173;
t156 = 0.4e1 * t173;
t158 = 0.2e1 * t173;
t161 = t162 ^ 2;
t177 = pkin(3) * t167;
t164 = t177 ^ 2;
t174 = pkin(1) * t171;
t182 = pkin(7) * t173;
t271 = t169 + t172;
t276 = t158 - t162;
t287 = t173 * t162;
t202 = t276 * t171 + t161 / 0.6e1 + t271 - t287;
t200 = 0.5e1 / 0.6e1 * t178 + t202;
t206 = t173 - t224;
t71 = t113 + t231;
t208 = t71 * t218;
t284 = t161 / 0.2e1 - t178 / 0.2e1;
t217 = -0.3e1 * t287 + t154 + t284;
t219 = -0.6e1 * t224;
t115 = -0.3e1 / 0.2e1 * t162;
t157 = 0.3e1 * t173;
t283 = t115 + t157;
t295 = t164 + t93 * ((t115 + t158) * t171 - 0.3e1 / 0.2e1 * t287 + t271 + t284);
t85 = 0.10e2 / 0.3e1 * t171;
t211 = ((t85 + t276) * t167 + t200) * t219 + (t137 + (-0.9e1 * t162 + 0.18e2 * t173) * t171 + t217) * t167 + (t138 + t283) * t178 + t295;
t277 = t151 + t173;
t235 = t167 + t277;
t308 = pkin(1) * t132;
t212 = t235 * t88 - t308 * (0.3e1 * t167 + t93);
t114 = -0.2e1 / 0.3e1 * t162;
t123 = -0.2e1 / 0.3e1 * t167;
t236 = t114 + t93;
t278 = t139 + t158;
t282 = t123 + t173;
t91 = -t162 - t167;
t79 = t157 + t91;
t297 = t79 * t171;
t213 = -t308 * (t146 + ((5 * t171) + t79) * t314 + (t123 + t236) * t93) + (t178 + (t114 + t123 + t278) * t167 + t148 + 0.2e1 * t297 + t173 * (t114 + t282)) * t88;
t275 = t161 - t178;
t216 = -0.6e1 * t287 + 0.6e1 * t172 + t275;
t285 = t174 * t100;
t221 = t285 * t88;
t237 = t114 + t119 + t158;
t294 = t178 + (t119 + t236) * t93;
t118 = 0.4e1 / 0.3e1 * t167;
t238 = t112 + t93;
t69 = t118 + t238;
t222 = t69 * t218 + t294 + (t150 + t237) * t167;
t302 = t169 * t101 ^ 2;
t223 = t302 * t88;
t249 = 0.16e2 * t285;
t225 = pkin(7) * t249;
t226 = 0.20e2 / 0.3e1 * t171;
t256 = pkin(7) * t285;
t230 = 0.8e1 * t256;
t274 = -t162 + t167;
t234 = t157 + t274;
t120 = t167 / 0.3e1;
t239 = t110 + t120 + t173;
t240 = t162 / 0.3e1 + t120 + t158;
t241 = 0.2e1 / 0.3e1 * t162 + t119 + t156;
t242 = 0.4e1 / 0.3e1 * t162 + t118 - 0.2e1 * t173;
t291 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t243 = t132 * t291;
t261 = 0.6e1 * t89;
t245 = pkin(7) * t261;
t262 = 0.4e1 * t89;
t246 = pkin(7) * t262;
t248 = -t308 / 0.2e1;
t288 = t171 * t101;
t250 = 0.12e2 * t288;
t253 = t171 * t88;
t255 = 0.4e1 * t288;
t300 = t177 * t131 * t98;
t257 = -0.8e1 * t300;
t258 = 0.8e1 * t302;
t264 = 0.2e1 * t308;
t265 = pkin(7) * t89;
t267 = 0.4e1 * pkin(7);
t268 = t172 + t178;
t269 = t172 - t169;
t279 = 0.4e1 / 0.7e1 * t173 - t162 / 0.7e1;
t280 = t126 + t173;
t281 = t171 / 0.3e1 + t173;
t286 = t173 * t171;
t290 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t298 = (-t132 * t174 + t253) * t101;
t299 = t178 * t98 ^ 2;
t301 = t173 * t91;
t86 = t93 ^ 2;
t306 = t86 * (-t167 + t233);
t312 = 4 * t169;
t316 = t111 + t167 / 0.2e1;
t55 = -t178 / 0.6e1 + t202;
t83 = t122 + t173;
t58 = t83 * t68;
t76 = t89 + pkin(7);
t63 = t87 + t76;
t92 = -0.3e1 * t167 + t173;
t66 = t92 * t230;
t67 = 0.10e2 * t297;
t72 = -t162 + t231;
t78 = pkin(7) * t263;
t81 = (t156 + t162) * t171;
t84 = -t171 / 0.3e1 + t173;
t90 = -0.30e2 * t162 + 0.60e2 * t173;
t95 = -(3 * t171) + t173;
t307 = ((-0.24e2 * (0.4e1 / 0.3e1 * t288 + t78 + t84) * t299 * t308 - 0.12e2 * (-0.8e1 / 0.3e1 * t223 + ((t124 + t239) * t88 - (0.7e1 / 0.6e1 * t167 + t110 + t280) * t308) * t255 + (-t167 * t270 - 0.5e1 / 0.3e1 * t169 + t240 * t171 + t173 * (t112 + t83)) * t88 + (-t178 + (-t226 + t241) * t167 - (3 * t169) + t242 * t171 + t172) * t248 + (-t132 * t169 * t100 + ((t171 + t239) * t88 + (t314 - t270) * t248) * t89) * t267) * t303 + 0.24e2 * t83 * t223 + ((t173 + 0.5e1 / 0.2e1 * t167 + 0.3e1 / 0.2e1 * t171 + t113) * t88 + t92 * t308 / 0.2e1) * t225 - 0.6e1 * ((-0.3e1 * t178 + (-t226 + t242) * t167 + t241 * t171 + t269) * t88 - 0.2e1 * (-0.5e1 / 0.3e1 * t178 + (-t171 + t240) * t167 + t173 * (t122 + t238)) * t308) * t288 - 0.6e1 * t213 * t265 - (t164 + ((21 * t171) + t79) * t178 + (t154 + (35 * t169) + t67 + 0.2e1 * t301) * t167 + (t147 + (t145 + t155 - 0.5e1 * t167) * t171 + t173 * (-t162 + t272)) * t93) * t88 + (0.7e1 * t164 + (t149 + t79) * t146 + ((21 * t169) + 0.9e1 * t172 + t67 + 0.6e1 * t301) * t167 + t306) * t308) * t38 + (0.16e2 * (t258 + t225 + (-(8 * t169) + 0.12e2 * t286) * t101 + (-0.12e2 * pkin(7) * t174 + t182 * t320) * t135 - 0.6e1 * t286 + t271) * t299 + 0.24e2 * (t282 * t258 + 0.14e2 * (-0.32e2 / 0.21e2 * (t173 + t167 / 0.4e1 + t171 / 0.4e1 - t162 / 0.8e1) * t224 + 0.5e1 / 0.42e2 * t178 + (0.16e2 / 0.21e2 * t171 + t279) * t167 + t169 / 0.7e1 + t279 * t171 + t172 - 0.3e1 / 0.7e1 * t287 + t161 / 0.42e2) * t288 + t84 * t208 - t270 * t178 + (-0.10e2 / 0.3e1 * t169 + 0.2e1 * t172 - t287 + t81) * t167 + t55 * t290 + ((-0.2e1 / 0.3e1 * t224 + t111 + t280) * t249 + (-0.8e1 / 0.3e1 * (t281 + t316) * t224 + 0.5e1 / 0.18e2 * t178 + (0.4e1 / 0.3e1 * t173 + t124 + t112) * t167 + t172 + 0.2e1 / 0.3e1 * t286 - 0.2e1 / 0.3e1 * t287 - t169 / 0.3e1 + t161 / 0.18e2) * t261) * pkin(7)) * t303 + 0.16e2 * (-0.6e1 * t173 * t167 + t268) * t302 + 0.32e2 * (t291 * t68 + t71 * t92) * t256 + 0.24e2 * (t83 * t208 - t164 + (-t85 + t273) * t178 + (t81 + t178 / 0.6e1 - t161 / 0.6e1 + t269) * t167 + t55 * t173) * t288 + 0.8e1 * t211 * t265 - 0.8e1 * ((t149 + t283) * t178 + (t147 + (t145 + 0.10e2 * t173) * t171 + t217) * t167 + t295) * t224 + t178 ^ 2 + (t144 + t156 + (28 * t171)) * t164 + (t171 * t90 + (70 * t169) + t216) * t178 + (t216 * t150 + t275 * t158 - 0.6e1 * t172 * t162 + t90 * t169 + 0.28e2 * t174 ^ 2 + 0.4e1 * t182 ^ 2) * t167 + t72 * t306) * t63 + (((0.4e1 * t298 + (t88 + t264) * t78 + t95 * t88 + (t113 + t235) * t264) * t257 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t167 + t126 + t110) * t134 * t160 + pkin(1) * t243) * t288 + (-0.8e1 * t221 + ((t112 + t119 + t277) * t88 - (0.8e1 / 0.3e1 * t167 + t238) * t308) * t262) * pkin(7) + t213) * t87) * t38 + (0.32e2 * (t230 + (-0.4e1 * t174 * t252 + t312 + (t313 + t144 + t155) * t171) * t101 + (-t171 + t206 + t316) * t246 + t68 * t290 + t95 * t71) * t300 + 0.8e1 * (t66 + (t291 * t71 + t58) * t250 + (t208 + (t150 + t276) * t167 + t200) * t245 + t211) * t87) * t63) * t76) / ((-0.4e1 * (-t270 * t88 + 0.2e1 * t298 + (0.2e1 * pkin(7) * t251 + t132 * (t167 + t311)) * pkin(1)) * t303 + 0.8e1 * pkin(7) * t221 + ((pkin(3) * t312 + 0.8e1 * t171 * t177) * t134 + 0.4e1 * t174 * t243) * t101 - 0.4e1 * t212 * t265 - (t167 * t278 + t148 + t268 + 0.6e1 * t286) * t88 + (t146 + (t139 + 0.6e1 * t173) * t167 + t86) * t308) * t38 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t224 + 0.4e1 / 0.9e1 * t167 - t162 / 0.9e1 + t281) * t288 + t84 * t68 + t69 * t290 + (t285 + (t110 + t119 + t206) * t89) * t267) * t303 + t66 + (t291 * t69 + t58) * t250 + t222 * t245 + ((t85 + t237) * t167 + t294) * t219 + t164 + (t138 + t234) * t178 + (t234 * t150 + t274 * t158 + t137 + t154) * t167 + t86 * t72) * t63 + ((t257 * t308 + (-0.2e1 * t101 * t253 + (t88 - t308) * t78 + t212) * t260) * t38 + (0.8e1 * (t78 + t255 + t95) * t300 + 0.6e1 * (t272 * t255 + (t68 + t69) * t246 + t222) * t87) * t63) * t76);
t13 = (t193 * t307 - t38 * t197 / 0.4e1) * t289;
t14 = (t38 * t193 + t197 * t307 / 0.4e1) * t289;
t59 = -t131 * t132 - t292;
t60 = -t131 * t135 + t293;
t3 = atan2(t13 * t59 + t14 * t60, -t13 * t60 + t14 * t59);
t1 = sin(t3);
t2 = cos(t3);
t247 = t163 * t296;
t17 = atan2(t38 * t247, t247 * t307);
t15 = sin(t17);
t16 = cos(t17);
t194 = t315 * t195;
t192 = t194 / 0.2e1;
t196 = t197 * t296;
t191 = atan2(t196, t192);
t189 = sin(t191);
t190 = cos(t191);
t31 = -t132 * t190 - t135 * t189;
t32 = t132 * t189 - t135 * t190;
t227 = -t15 * t31 + t32 * t16;
t5 = t15 * t32 + t16 * t31;
t322 = t1 * t5 - t2 * t227;
t321 = t1 * t227 + t2 * t5;
t319 = m(4) + m(10);
t310 = pkin(5) * m(6);
t305 = t32 * pkin(3) - t89;
t266 = -m(7) - m(9) - m(3);
t229 = t31 * pkin(3) - t308;
t133 = sin(pkin(9));
t136 = cos(pkin(9));
t128 = sin(pkin(8));
t129 = cos(pkin(8));
t56 = t128 * t135 - t129 * t132;
t57 = t128 * t132 + t129 * t135;
t34 = -t56 * t194 / 0.2e1 + t57 * t196;
t35 = t192 * t57 + t196 * t56;
t25 = t133 * t35 - t136 * t34;
t26 = t133 * t34 + t136 * t35;
t22 = atan2(t25, t26);
t18 = sin(t22);
t20 = cos(t22);
t210 = t18 * t32 + t20 * t31;
t209 = t18 * t31 - t20 * t32;
t204 = -m(6) - m(5) - m(11) - m(8) - m(2) - m(1) - t319;
t23 = atan2(t25, -t26);
t19 = sin(t23);
t21 = cos(t23);
t203 = t19 * mrSges(10,1) + t21 * mrSges(10,2) - mrSges(4,2) - mrSges(7,2);
t199 = -m(10) * pkin(6) + t21 * mrSges(10,1) - t19 * mrSges(10,2) - mrSges(4,1) - mrSges(7,1);
t50 = atan2(t56, -t57);
t49 = atan2(t56, t57);
t48 = cos(t50);
t47 = cos(t49);
t46 = sin(t50);
t45 = sin(t49);
t42 = -t132 * t45 + t135 * t47;
t41 = t132 * t47 + t135 * t45;
t40 = -t128 * t46 + t129 * t48;
t39 = t128 * t48 + t129 * t46;
t4 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t204 + t266) * r_base(3)) * g(3) + (-m(5) * t229 - m(11) * (pkin(4) * t5 + t229) - t5 * mrSges(5,1) + t321 * mrSges(11,1) - t322 * mrSges(11,2) + mrSges(8,1) * t134 - mrSges(8,2) * t131 + mrSges(2,1) * t132 + mrSges(2,2) * t135 - mrSges(9,1) * t41 - mrSges(9,2) * t42 - mrSges(3,1) * t31 - mrSges(3,2) * t32 - t128 * t310 - t227 * mrSges(5,2) - t39 * mrSges(6,1) - t40 * mrSges(6,2) - mrSges(1,2) - t199 * t210 + t203 * t209 + t204 * r_base(2) + t266 * (r_base(2) - t308) + t319 * (-t31 * pkin(2) + t308)) * g(2) + (-t129 * t310 + mrSges(2,1) * t135 - mrSges(2,2) * t132 + t5 * mrSges(5,2) - t131 * mrSges(8,1) + t199 * t209 + t203 * t210 - t134 * mrSges(8,2) - t227 * mrSges(5,1) - t322 * mrSges(11,1) - t321 * mrSges(11,2) - mrSges(9,1) * t42 + mrSges(9,2) * t41 - mrSges(3,1) * t32 + mrSges(3,2) * t31 + t204 * r_base(1) + t266 * (-t89 + r_base(1)) - m(8) * pkin(7) - m(11) * (pkin(4) * t227 + t305) - m(5) * t305 + t39 * mrSges(6,2) - t40 * mrSges(6,1) - mrSges(1,1) - t319 * (t32 * pkin(2) - t89)) * g(1);
U = t4;
