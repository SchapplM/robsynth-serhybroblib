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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2TE_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:34:07
% EndTime: 2020-05-09 10:34:17
% DurationCPUTime: 6.26s
% Computational Cost: add. (40082->463), mult. (117780->572), div. (838->5), fcn. (20852->10), ass. (0->229)
t116 = sin(qJ(1));
t150 = pkin(3) ^ 2;
t155 = (pkin(1) ^ 2);
t157 = pkin(7) ^ 2;
t77 = t155 + t157;
t216 = t150 + t77;
t118 = cos(qJ(1));
t115 = sin(qJ(2));
t71 = pkin(3) * t115;
t296 = t71 + pkin(7);
t223 = t296 * t118;
t117 = cos(qJ(2));
t275 = t116 * t117;
t236 = pkin(3) * t275;
t210 = pkin(1) * t236;
t51 = -0.2e1 * t210;
t302 = 0.2e1 * pkin(7);
t61 = t71 * t302;
t181 = 0.1e1 / (0.2e1 * pkin(1) * t223 + t216 + t51 + t61);
t82 = t115 ^ 2;
t288 = t150 * t82;
t244 = 0.2e1 * t288;
t258 = -t150 + t157;
t190 = t61 + t244 + t258;
t85 = t118 ^ 2;
t184 = t190 * t85;
t200 = -pkin(1) + t236;
t134 = 3 * t155;
t145 = pkin(4) ^ 2;
t259 = t145 - t157;
t217 = -0.2e1 * t150 + t259;
t201 = t134 - t217;
t202 = -0.4e1 * t210;
t218 = -t145 + t77;
t274 = t117 * t118;
t235 = pkin(3) * t274;
t143 = 0.2e1 * pkin(3);
t153 = t155 ^ 2;
t204 = t61 + t218;
t246 = -0.4e1 * t71;
t256 = t155 - t157;
t298 = 2 * t155;
t308 = 0.4e1 * pkin(1);
t29 = sqrt(-0.4e1 * t155 * t184 + 0.4e1 * t256 * t288 + pkin(7) * t218 * t246 - t153 + t217 * t298 - (t157 - (t143 + pkin(4)) * pkin(4)) * (t157 + (t143 - pkin(4)) * pkin(4)) + (-(t51 + t204) * t223 + t204 * t236) * t308);
t300 = 0.4e1 * t150;
t177 = t181 * ((t116 * t296 + t235) * t29 - (t201 + t61 + t202) * t223 + t200 * t61 + t201 * t236 + (-0.2e1 * t184 + t244 - t300 - t218) * pkin(1));
t174 = t177 / 0.4e1;
t301 = 0.2e1 * t150;
t194 = t301 + t204;
t73 = pkin(1) * t118;
t249 = 0.2e1 * t73;
t72 = pkin(3) * t117;
t180 = t181 * ((t223 - t200) * t29 + (t190 * t249 + t194 * t296) * t116 + (t118 * t194 + (0.4e1 * t85 - 0.2e1) * pkin(1) * t296) * t72);
t146 = 0.1e1 / pkin(4);
t270 = t146 / pkin(3) ^ 2;
t103 = 0.2e1 / 0.3e1 * t150;
t106 = -t150 / 0.3e1;
t108 = 0.4e1 / 0.3e1 * t155;
t110 = t155 / 0.2e1;
t120 = 15 * t153;
t121 = 15 * t155;
t122 = 10 * t155;
t127 = -0.2e1 * t145;
t128 = -0.5e1 * t145;
t162 = t150 ^ 2;
t129 = 0.5e1 * t162;
t130 = 7 * t153;
t131 = 5 * t153;
t132 = 7 * t155;
t133 = 6 * t155;
t156 = t157 ^ 2;
t137 = 0.3e1 * t156;
t138 = 0.8e1 * t157;
t139 = 0.4e1 * t157;
t141 = 0.2e1 * t157;
t144 = t145 ^ 2;
t161 = pkin(3) * t150;
t147 = t161 ^ 2;
t158 = pkin(1) * t155;
t166 = pkin(7) * t157;
t257 = t153 + t156;
t262 = t141 - t145;
t269 = t157 * t145;
t186 = t262 * t155 + t144 / 0.6e1 + t257 - t269;
t183 = 0.5e1 / 0.6e1 * t162 + t186;
t191 = t157 - t210;
t97 = -t145 / 0.2e1;
t54 = t97 + t216;
t195 = t54 * t202;
t203 = -0.6e1 * t210;
t280 = -t162 / 0.2e1 + t144 / 0.2e1;
t206 = t137 - 0.3e1 * t269 + t280;
t99 = -0.3e1 / 0.2e1 * t145;
t277 = t147 + t77 * ((t99 + t141) * t155 - 0.3e1 / 0.2e1 * t269 + t257 + t280);
t140 = 0.3e1 * t157;
t278 = t140 + t99;
t69 = 0.10e2 / 0.3e1 * t155;
t197 = ((t69 + t262) * t150 + t183) * t203 + (t120 + (-0.9e1 * t145 + 0.18e2 * t157) * t155 + t206) * t150 + (t121 + t278) * t162 + t277;
t263 = t134 + t157;
t220 = t150 + t263;
t294 = pkin(1) * t116;
t198 = t220 * t72 - t294 * (0.3e1 * t150 + t77);
t107 = -0.2e1 / 0.3e1 * t150;
t98 = -0.2e1 / 0.3e1 * t145;
t225 = t98 + t77;
t264 = t122 + t141;
t267 = t107 + t157;
t75 = -t145 - t150;
t63 = t140 + t75;
t281 = t63 * t155;
t199 = -t294 * (t129 + ((5 * t155) + t63) * t301 + (t107 + t225) * t77) + (t162 + (t107 + t98 + t264) * t150 + t131 + 0.2e1 * t281 + t157 * (t98 + t267)) * t72;
t261 = t144 - t162;
t205 = 0.6e1 * t156 - 0.6e1 * t269 + t261;
t229 = t103 + t141 + t98;
t276 = t162 + (t103 + t225) * t77;
t102 = 0.4e1 / 0.3e1 * t150;
t96 = -t145 / 0.3e1;
t224 = t96 + t77;
t52 = t102 + t224;
t207 = t52 * t202 + t276 + (t133 + t229) * t150;
t84 = t118 * t85;
t284 = t158 * t84;
t208 = t284 * t72;
t287 = t153 * t85 ^ 2;
t209 = t287 * t72;
t211 = 0.20e2 / 0.3e1 * t155;
t239 = 0.16e2 * t284;
t215 = pkin(7) * t239;
t260 = -t145 + t150;
t219 = t140 + t260;
t245 = pkin(7) * t284;
t221 = 0.8e1 * t245;
t273 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t222 = t116 * t273;
t104 = t150 / 0.3e1;
t94 = -t145 / 0.6e1;
t226 = t104 + t94 + t157;
t227 = t104 + t145 / 0.3e1 + t141;
t228 = t103 + 0.2e1 / 0.3e1 * t145 + t139;
t230 = t102 + 0.4e1 / 0.3e1 * t145 - 0.2e1 * t157;
t247 = 0.6e1 * t73;
t231 = pkin(7) * t247;
t248 = 0.4e1 * t73;
t232 = pkin(7) * t248;
t234 = -t294 / 0.2e1;
t237 = t155 * t72;
t286 = t155 * t85;
t240 = 0.12e2 * t286;
t283 = t161 * t115 * t82;
t241 = -0.8e1 * t283;
t242 = 0.4e1 * t286;
t243 = 0.8e1 * t287;
t250 = 0.2e1 * t294;
t251 = pkin(7) * t73;
t253 = 0.4e1 * pkin(7);
t254 = t156 + t162;
t255 = t156 - t153;
t265 = t110 + t157;
t266 = t155 / 0.3e1 + t157;
t268 = t157 * t155;
t272 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t279 = 0.4e1 / 0.7e1 * t157 - t145 / 0.7e1;
t282 = t162 * t82 ^ 2;
t285 = t157 * t75;
t70 = t77 ^ 2;
t291 = t70 * (-t150 + t218);
t292 = (-t116 * t158 + t237) * t85;
t299 = 4 * t153;
t95 = -t145 / 0.4e1;
t304 = t150 / 0.2e1 + t95;
t38 = -t162 / 0.6e1 + t186;
t67 = t106 + t157;
t41 = t67 * t51;
t60 = t73 + pkin(7);
t46 = t71 + t60;
t76 = -0.3e1 * t150 + t157;
t49 = t76 * t221;
t50 = 0.10e2 * t281;
t55 = -t145 + t216;
t62 = pkin(7) * t249;
t65 = (t139 + t145) * t155;
t68 = -t155 / 0.3e1 + t157;
t74 = -0.30e2 * t145 + 0.60e2 * t157;
t79 = -(3 * t155) + t157;
t293 = ((-0.24e2 * (0.4e1 / 0.3e1 * t286 + t62 + t68) * t282 * t294 - 0.12e2 * (-0.8e1 / 0.3e1 * t209 + ((t108 + t226) * t72 - (0.7e1 / 0.6e1 * t150 + t94 + t265) * t294) * t242 + (-t150 * t256 - 0.5e1 / 0.3e1 * t153 + t227 * t155 + t157 * (t96 + t67)) * t72 + (-t162 + (-t211 + t228) * t150 - (3 * t153) + t230 * t155 + t156) * t234 + (-t116 * t153 * t84 + ((t155 + t226) * t72 + (t301 - t256) * t234) * t73) * t253) * t288 + 0.24e2 * t67 * t209 + ((t157 + 0.5e1 / 0.2e1 * t150 + 0.3e1 / 0.2e1 * t155 + t97) * t72 + t76 * t294 / 0.2e1) * t215 - 0.6e1 * ((-0.3e1 * t162 + (-t211 + t230) * t150 + t228 * t155 + t255) * t72 - 0.2e1 * (-0.5e1 / 0.3e1 * t162 + (-t155 + t227) * t150 + t157 * (t106 + t224)) * t294) * t286 - 0.6e1 * t199 * t251 - (t147 + ((21 * t155) + t63) * t162 + (t137 + (35 * t153) + t50 + 0.2e1 * t285) * t150 + (t130 + (t128 + t138 - 0.5e1 * t150) * t155 + t157 * (-t145 + t258)) * t77) * t72 + (0.7e1 * t147 + (t132 + t63) * t129 + ((21 * t153) + 0.9e1 * t156 + t50 + 0.6e1 * t285) * t150 + t291) * t294) * t29 + (0.16e2 * (t243 + t215 + (-(8 * t153) + 0.12e2 * t268) * t85 + (-0.12e2 * pkin(7) * t158 + t166 * t308) * t118 - 0.6e1 * t268 + t257) * t282 + 0.24e2 * (t267 * t243 + 0.14e2 * (-0.32e2 / 0.21e2 * (t157 + t150 / 0.4e1 + t155 / 0.4e1 - t145 / 0.8e1) * t210 + 0.5e1 / 0.42e2 * t162 + (0.16e2 / 0.21e2 * t155 + t279) * t150 + t153 / 0.7e1 + t279 * t155 + t156 - 0.3e1 / 0.7e1 * t269 + t144 / 0.42e2) * t286 + t68 * t195 - t256 * t162 + (-0.10e2 / 0.3e1 * t153 + 0.2e1 * t156 - t269 + t65) * t150 + t38 * t272 + ((-0.2e1 / 0.3e1 * t210 + t95 + t265) * t239 + (-0.8e1 / 0.3e1 * (t266 + t304) * t210 + 0.5e1 / 0.18e2 * t162 + (0.4e1 / 0.3e1 * t157 + t108 + t96) * t150 + t156 + 0.2e1 / 0.3e1 * t268 - 0.2e1 / 0.3e1 * t269 - t153 / 0.3e1 + t144 / 0.18e2) * t247) * pkin(7)) * t288 + 0.16e2 * (-0.6e1 * t157 * t150 + t254) * t287 + 0.32e2 * (t273 * t51 + t54 * t76) * t245 + 0.24e2 * (t67 * t195 - t147 + (-t69 + t259) * t162 + (t65 + t162 / 0.6e1 - t144 / 0.6e1 + t255) * t150 + t38 * t157) * t286 + 0.8e1 * t197 * t251 - 0.8e1 * ((t132 + t278) * t162 + (t130 + (t128 + 0.10e2 * t157) * t155 + t206) * t150 + t277) * t210 + t162 ^ 2 + (t127 + t139 + (28 * t155)) * t147 + (t155 * t74 + (70 * t153) + t205) * t162 + (t205 * t133 + t261 * t141 - 0.6e1 * t156 * t145 + t74 * t153 + 0.28e2 * t158 ^ 2 + 0.4e1 * t166 ^ 2) * t150 + t55 * t291) * t46 + (((0.4e1 * t292 + (t72 + t250) * t62 + t79 * t72 + (t97 + t220) * t250) * t241 - 0.6e1 * (-0.4e1 * (0.2e1 * (0.5e1 / 0.6e1 * t150 + t110 + t94) * t72 + pkin(1) * t222) * t286 + (-0.8e1 * t208 + ((t103 + t96 + t263) * t72 - (0.8e1 / 0.3e1 * t150 + t224) * t294) * t248) * pkin(7) + t199) * t71) * t29 + (0.32e2 * (t221 + (-0.4e1 * t158 * t236 + t299 + (t300 + t127 + t138) * t155) * t85 + (-t155 + t191 + t304) * t232 + t51 * t272 + t79 * t54) * t283 + 0.8e1 * (t49 + (t273 * t54 + t41) * t240 + (t195 + (t133 + t262) * t150 + t183) * t231 + t197) * t71) * t46) * t60) / ((-0.4e1 * (-t256 * t72 + 0.2e1 * t292 + (t235 * t302 + t116 * (t150 + t298)) * pkin(1)) * t288 + 0.8e1 * pkin(7) * t208 + ((pkin(3) * t299 + 0.8e1 * t155 * t161) * t117 + 0.4e1 * t158 * t222) * t85 - 0.4e1 * t198 * t251 - (t150 * t264 + t131 + t254 + 0.6e1 * t268) * t72 + (t129 + (t122 + 0.6e1 * t157) * t150 + t70) * t294) * t29 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t210 + 0.4e1 / 0.9e1 * t150 - t145 / 0.9e1 + t266) * t286 + t68 * t51 + t52 * t272 + (t284 + (t103 + t191 + t94) * t73) * t253) * t288 + t49 + (t273 * t52 + t41) * t240 + t207 * t231 + ((t69 + t229) * t150 + t276) * t203 + t147 + (t121 + t219) * t162 + (t219 * t133 + t260 * t141 + t120 + t137) * t150 + t70 * t55) * t46 + ((t241 * t294 + (-0.2e1 * t85 * t237 + (t72 - t294) * t62 + t198) * t246) * t29 + (0.8e1 * (t62 + t242 + t79) * t283 + 0.6e1 * (t258 * t242 + (t51 + t52) * t232 + t207) * t71) * t46) * t60);
t11 = (t174 * t293 - t29 * t180 / 0.4e1) * t270;
t12 = (t29 * t174 + t180 * t293 / 0.4e1) * t270;
t42 = -t115 * t116 - t274;
t43 = -t115 * t118 + t275;
t1 = t11 * t42 + t12 * t43;
t233 = t293 / 0.2e1;
t151 = 0.1e1 / pkin(3);
t176 = -t177 / 0.2e1;
t179 = -t180 / 0.2e1;
t25 = (t116 * t176 + t118 * t179) * t151;
t178 = t180 / 0.2e1;
t26 = (t116 * t178 + t118 * t176) * t151;
t271 = t146 * t151;
t189 = (-t25 * t29 / 0.2e1 + t26 * t233) * t271;
t2 = t11 * t43 - t12 * t42;
t303 = (t25 * t233 + t26 * t29 / 0.2e1) * t271;
t310 = t1 * t303 + t189 * t2;
t309 = -t1 * t189 + t2 * t303;
t307 = m(10) + m(4);
t297 = pkin(5) * m(6);
t295 = sin(pkin(9));
t290 = t26 * pkin(3) - t73;
t252 = -m(9) - m(3) - m(7);
t214 = t25 * pkin(3) - t294;
t119 = cos(pkin(9));
t175 = t177 / 0.2e1;
t112 = sin(pkin(8));
t113 = cos(pkin(8));
t39 = -t112 * t118 + t113 * t116;
t40 = t112 * t116 + t113 * t118;
t173 = t151 * (t175 * t40 + t179 * t39);
t20 = (t175 * t39 + t178 * t40) * t151;
t16 = t119 * t20 - t173 * t295;
t17 = t119 * t173 + t20 * t295;
t212 = -t16 * t25 - t26 * t17;
t196 = t16 * t26 - t17 * t25;
t193 = t116 * t40 - t118 * t39;
t192 = t116 * t39 + t118 * t40;
t188 = -m(5) - m(11) - m(6) - m(8) - m(2) - m(1) - t307;
t187 = -t16 * mrSges(10,1) - t17 * mrSges(10,2) - mrSges(4,2) - mrSges(7,2);
t182 = -m(10) * pkin(6) - t17 * mrSges(10,1) + t16 * mrSges(10,2) - mrSges(4,1) - mrSges(7,1);
t33 = t112 * t39 - t113 * t40;
t32 = -t112 * t40 - t113 * t39;
t3 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t188 + t252) * r_base(3)) * g(3) + (-m(5) * t214 - m(11) * (pkin(4) * t303 + t214) - t112 * t297 - t32 * mrSges(6,1) - t33 * mrSges(6,2) + mrSges(2,1) * t116 + mrSges(2,2) * t118 - t193 * mrSges(9,1) - t192 * mrSges(9,2) - mrSges(3,1) * t25 - mrSges(3,2) * t26 + mrSges(8,1) * t117 - mrSges(8,2) * t115 - t189 * mrSges(5,2) - t309 * mrSges(11,1) - t310 * mrSges(11,2) - t303 * mrSges(5,1) - mrSges(1,2) + t182 * t196 + t187 * t212 + t188 * r_base(2) + t252 * (r_base(2) - t294) + t307 * (-t25 * pkin(2) + t294)) * g(2) + (t252 * (-t73 + r_base(1)) + t32 * mrSges(6,2) - t33 * mrSges(6,1) + t188 * r_base(1) - m(8) * pkin(7) - t192 * mrSges(9,1) + t193 * mrSges(9,2) - m(11) * (pkin(4) * t189 + t290) - m(5) * t290 + mrSges(2,1) * t118 - mrSges(2,2) * t116 - t187 * t196 - t117 * mrSges(8,2) - t115 * mrSges(8,1) - mrSges(3,1) * t26 + mrSges(3,2) * t25 + t182 * t212 + t303 * mrSges(5,2) - t310 * mrSges(11,1) + t309 * mrSges(11,2) - t189 * mrSges(5,1) - t113 * t297 - mrSges(1,1) - t307 * (t26 * pkin(2) - t73)) * g(1);
U = t3;
