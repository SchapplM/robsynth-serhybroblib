% Calculate potential energy for
% picker2Dm2DE2
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
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 19:15:43
% EndTime: 2020-05-09 19:15:51
% DurationCPUTime: 5.30s
% Computational Cost: add. (31130->460), mult. (92616->543), div. (594->5), fcn. (15660->31), ass. (0->235)
t297 = 4 * pkin(1);
t296 = m(10) + m(4);
t154 = pkin(4) ^ 2;
t102 = -t154 / 0.4e1;
t159 = pkin(3) ^ 2;
t295 = t102 + t159 / 0.2e1;
t294 = 2 * pkin(7);
t126 = cos(qJ(1));
t92 = t126 ^ 2;
t293 = -0.2e1 * t92;
t138 = 0.2e1 * t159;
t292 = 0.4e1 * t159;
t164 = pkin(1) ^ 2;
t162 = t164 ^ 2;
t291 = 4 * t162;
t290 = 2 * t164;
t142 = 6 * t164;
t166 = pkin(7) ^ 2;
t150 = 2 * t166;
t171 = t159 ^ 2;
t137 = 0.5e1 * t171;
t289 = pkin(5) * m(6);
t119 = sin(pkin(8));
t120 = cos(pkin(8));
t123 = sin(qJ(1));
t45 = t119 * t123 + t120 * t126;
t288 = t45 / 0.2e1;
t80 = pkin(1) * t126;
t67 = t80 + pkin(7);
t287 = pkin(1) * t123;
t122 = sin(qJ(2));
t78 = pkin(3) * t122;
t125 = cos(qJ(2));
t79 = pkin(3) * t125;
t101 = -t154 / 0.6e1;
t103 = -t154 / 0.3e1;
t104 = -t154 / 0.2e1;
t110 = 0.2e1 / 0.3e1 * t159;
t113 = -t159 / 0.3e1;
t115 = 0.4e1 / 0.3e1 * t164;
t117 = t164 / 0.2e1;
t128 = 15 * t162;
t129 = 15 * t164;
t130 = 10 * t164;
t135 = -0.2e1 * t154;
t136 = -0.5e1 * t154;
t139 = 7 * t162;
t140 = 5 * t162;
t141 = 7 * t164;
t165 = t166 ^ 2;
t146 = 3 * t165;
t147 = 8 * t166;
t148 = 4 * t166;
t153 = t154 ^ 2;
t170 = pkin(3) * t159;
t156 = t170 ^ 2;
t167 = pkin(1) * t164;
t175 = pkin(7) * t166;
t251 = t162 + t165;
t256 = t150 - t154;
t265 = t166 * t154;
t183 = t256 * t164 + t153 / 0.6e1 + t251 - t265;
t182 = 0.5e1 / 0.6e1 * t171 + t183;
t269 = t123 * t125;
t231 = pkin(3) * t269;
t203 = pkin(1) * t231;
t185 = t166 - t203;
t194 = -0.4e1 * t203;
t84 = t164 + t166;
t209 = t159 + t84;
t61 = t104 + t209;
t186 = t61 * t194;
t195 = -0.6e1 * t203;
t272 = -t171 / 0.2e1 + t153 / 0.2e1;
t199 = t146 - 0.3e1 * t265 + t272;
t106 = -0.3e1 / 0.2e1 * t154;
t149 = 3 * t166;
t263 = t106 + t149;
t271 = t156 + t84 * ((t106 + t150) * t164 - 0.3e1 / 0.2e1 * t265 + t251 + t272);
t76 = 0.10e2 / 0.3e1 * t164;
t189 = ((t76 + t256) * t159 + t182) * t195 + (t128 + (-0.9e1 * t154 + (18 * t166)) * t164 + t199) * t159 + (t129 + t263) * t171 + t271;
t143 = 3 * t164;
t257 = t143 + t166;
t212 = t159 + t257;
t190 = t212 * t79 - t287 * (0.3e1 * t159 + t84);
t105 = -0.2e1 / 0.3e1 * t154;
t114 = -0.2e1 / 0.3e1 * t159;
t213 = t105 + t84;
t258 = t130 + t150;
t262 = t114 + t166;
t82 = -t154 - t159;
t70 = t149 + t82;
t273 = t70 * t164;
t191 = -t287 * (t137 + ((5 * t164) + t70) * t138 + (t114 + t213) * t84) + (t171 + (t105 + t114 + t258) * t159 + t140 + 0.2e1 * t273 + t166 * (t105 + t262)) * t79;
t255 = t153 - t171;
t197 = (6 * t165) - 0.6e1 * t265 + t255;
t214 = t105 + t110 + t150;
t270 = t171 + (t110 + t213) * t84;
t109 = 0.4e1 / 0.3e1 * t159;
t215 = t103 + t84;
t59 = t109 + t215;
t200 = t59 * t194 + t270 + (t142 + t214) * t159;
t228 = t167 * t79;
t91 = t126 * t92;
t201 = t91 * t228;
t279 = t162 * t92 ^ 2;
t202 = t279 * t79;
t204 = 0.20e2 / 0.3e1 * t164;
t276 = t167 * t91;
t232 = 0.16e2 * t276;
t205 = pkin(7) * t232;
t254 = -t154 + t159;
t211 = t149 + t254;
t111 = t159 / 0.3e1;
t216 = t101 + t111 + t166;
t239 = pkin(7) * t276;
t217 = 0.8e1 * t239;
t267 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t219 = t123 * t267;
t220 = t111 + t154 / 0.3e1 + t150;
t221 = t110 + 0.2e1 / 0.3e1 * t154 + t148;
t222 = t109 + 0.4e1 / 0.3e1 * t154 - (2 * t166);
t241 = 0.6e1 * t80;
t223 = pkin(7) * t241;
t242 = 0.4e1 * t80;
t224 = pkin(7) * t242;
t226 = -t287 / 0.2e1;
t229 = t164 * t79;
t268 = t125 * t126;
t230 = pkin(3) * t268;
t278 = t164 * t92;
t233 = 0.12e2 * t278;
t89 = t122 ^ 2;
t275 = t170 * t122 * t89;
t234 = -0.8e1 * t275;
t235 = 0.4e1 * t278;
t236 = -0.4e1 * t278;
t237 = 0.8e1 * t279;
t240 = -0.4e1 * t78;
t244 = 0.2e1 * t287;
t245 = pkin(7) * t80;
t247 = 4 * pkin(7);
t248 = t165 + t171;
t249 = t165 - t162;
t250 = t164 - t166;
t252 = -t159 + t166;
t253 = t154 - t166;
t259 = 0.4e1 / 0.7e1 * t166 - t154 / 0.7e1;
t260 = t117 + t166;
t261 = t164 / 0.3e1 + t166;
t264 = t166 * t164;
t266 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t274 = t171 * t89 ^ 2;
t277 = t166 * t82;
t281 = t159 * t89;
t210 = -t154 + t84;
t77 = t84 ^ 2;
t284 = t77 * (-t159 + t210);
t285 = (-t123 * t167 + t229) * t92;
t152 = 0.2e1 * pkin(3);
t68 = t78 * t294;
t196 = t68 + t210;
t66 = t78 + pkin(7);
t282 = t126 * t66;
t238 = 0.2e1 * t281;
t46 = t68 + t238 + t252;
t58 = -0.2e1 * t203;
t29 = sqrt(t46 * t236 + 0.4e1 * t250 * t281 + pkin(7) * t210 * t240 - t162 + (-0.2e1 * t159 + t253) * t290 - (t166 - (t152 + pkin(4)) * pkin(4)) * (t166 + (t152 - pkin(4)) * pkin(4)) + (-(t58 + t196) * t282 + t196 * t231) * t297);
t43 = -t171 / 0.6e1 + t183;
t74 = t113 + t166;
t47 = t74 * t58;
t53 = t78 + t67;
t83 = -0.3e1 * t159 + t166;
t56 = t83 * t217;
t57 = 0.10e2 * t273;
t62 = -t154 + t209;
t243 = 0.2e1 * t80;
t69 = pkin(7) * t243;
t72 = (t148 + t154) * t164;
t75 = -t164 / 0.3e1 + t166;
t81 = -0.30e2 * t154 + (60 * t166);
t86 = -3 * t164 + t166;
t286 = ((-0.24e2 * (0.4e1 / 0.3e1 * t278 + t69 + t75) * t274 * t287 - 0.12e2 * (-0.8e1 / 0.3e1 * t202 + ((t115 + t216) * t79 - (0.7e1 / 0.6e1 * t159 + t101 + t260) * t287) * t235 + (-t159 * t250 - 0.5e1 / 0.3e1 * t162 + t220 * t164 + t166 * (t103 + t74)) * t79 + (-t171 + (-t204 + t221) * t159 - (3 * t162) + t222 * t164 + t165) * t226 + (-t123 * t162 * t91 + ((t164 + t216) * t79 + (t138 - t250) * t226) * t80) * t247) * t281 + 0.24e2 * t74 * t202 + ((t166 + 0.5e1 / 0.2e1 * t159 + 0.3e1 / 0.2e1 * t164 + t104) * t79 + t83 * t287 / 0.2e1) * t205 - 0.6e1 * ((-0.3e1 * t171 + (-t204 + t222) * t159 + t221 * t164 + t249) * t79 - 0.2e1 * (-0.5e1 / 0.3e1 * t171 + (-t164 + t220) * t159 + t166 * (t113 + t215)) * t287) * t278 - 0.6e1 * t191 * t245 - (t156 + ((21 * t164) + t70) * t171 + (t146 + (35 * t162) + t57 + 0.2e1 * t277) * t159 + (t139 + (t136 + t147 - 0.5e1 * t159) * t164 + t166 * (-t154 + t252)) * t84) * t79 + (0.7e1 * t156 + (t141 + t70) * t137 + ((21 * t162) + (9 * t165) + t57 + 0.6e1 * t277) * t159 + t284) * t287) * t29 + (0.16e2 * (t237 + t205 + (-8 * t162 + 12 * t264) * t92 + (-12 * pkin(7) * t167 + t175 * t297) * t126 - (6 * t264) + t251) * t274 + 0.24e2 * (t262 * t237 + 0.14e2 * (-0.32e2 / 0.21e2 * (t166 + t159 / 0.4e1 + t164 / 0.4e1 - t154 / 0.8e1) * t203 + 0.5e1 / 0.42e2 * t171 + (0.16e2 / 0.21e2 * t164 + t259) * t159 + t162 / 0.7e1 + t259 * t164 + t165 - 0.3e1 / 0.7e1 * t265 + t153 / 0.42e2) * t278 + t75 * t186 - t250 * t171 + (-0.10e2 / 0.3e1 * t162 + (2 * t165) - t265 + t72) * t159 + t43 * t266 + ((-0.2e1 / 0.3e1 * t203 + t102 + t260) * t232 + (-0.8e1 / 0.3e1 * (t261 + t295) * t203 + 0.5e1 / 0.18e2 * t171 + (0.4e1 / 0.3e1 * t166 + t115 + t103) * t159 + t165 + 0.2e1 / 0.3e1 * t264 - 0.2e1 / 0.3e1 * t265 - t162 / 0.3e1 + t153 / 0.18e2) * t241) * pkin(7)) * t281 + 0.16e2 * (-0.6e1 * t166 * t159 + t248) * t279 + 0.32e2 * (t267 * t58 + t61 * t83) * t239 + 0.24e2 * (t74 * t186 - t156 + (-t76 + t253) * t171 + (t72 + t171 / 0.6e1 - t153 / 0.6e1 + t249) * t159 + t43 * t166) * t278 + 0.8e1 * t189 * t245 - 0.8e1 * ((t141 + t263) * t171 + (t139 + (t136 + (10 * t166)) * t164 + t199) * t159 + t271) * t203 + t171 ^ 2 + (t135 + t148 + (28 * t164)) * t156 + (t164 * t81 + (70 * t162) + t197) * t171 + (t197 * t142 + t255 * t150 - 0.6e1 * t165 * t154 + t81 * t162 + (28 * t167 ^ 2) + (4 * t175 ^ 2)) * t159 + t62 * t284) * t53 + (((0.4e1 * t285 + (t79 + t244) * t69 + t86 * t79 + (t104 + t212) * t244) * t234 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t159 + t117 + t101) * t79 + pkin(1) * t219) * t236 + (-0.8e1 * t201 + ((t103 + t110 + t257) * t79 - (0.8e1 / 0.3e1 * t159 + t215) * t287) * t242) * pkin(7) + t191) * t78) * t29 + (0.32e2 * (t217 + (-0.4e1 * t123 * t228 + t291 + (t292 + t135 + t147) * t164) * t92 + (-t164 + t185 + t295) * t224 + t58 * t266 + t86 * t61) * t275 + 0.8e1 * (t56 + (t267 * t61 + t47) * t233 + (t186 + (t142 + t256) * t159 + t182) * t223 + t189) * t78) * t53) * t67) / ((-0.4e1 * (-t250 * t79 + 0.2e1 * t285 + (t230 * t294 + t123 * (t159 + t290)) * pkin(1)) * t281 + 0.8e1 * pkin(7) * t201 + ((pkin(3) * t291 + 0.8e1 * t164 * t170) * t125 + 0.4e1 * t167 * t219) * t92 - 0.4e1 * t190 * t245 - (t159 * t258 + t140 + t248 + (6 * t264)) * t79 + (t137 + (t130 + 6 * t166) * t159 + t77) * t287) * t29 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t203 + 0.4e1 / 0.9e1 * t159 - t154 / 0.9e1 + t261) * t278 + t75 * t58 + t59 * t266 + (t276 + (t101 + t110 + t185) * t80) * t247) * t281 + t56 + (t267 * t59 + t47) * t233 + t200 * t223 + ((t76 + t214) * t159 + t270) * t195 + t156 + (t129 + t211) * t171 + (t211 * t142 + t254 * t150 + t128 + t146) * t159 + t77 * t62) * t53 + ((t234 * t287 + (t229 * t293 + (t79 - t287) * t69 + t190) * t240) * t29 + (0.8e1 * (t69 + t235 + t86) * t275 + 0.6e1 * (t252 * t235 + (t58 + t59) * t224 + t200) * t78) * t53) * t67);
t283 = mrSges(7,2) + mrSges(4,2);
t160 = 0.1e1 / pkin(3);
t38 = 0.1e1 / (t243 * t66 + t209 + t58 + t68);
t280 = t160 * t38;
t218 = t280 / 0.2e1;
t192 = -pkin(1) + t231;
t193 = t138 + t143 - t253;
t27 = (t123 * t66 + t230) * t29 - (t193 + t68 + t194) * t282 + t192 * t68 + t193 * t231 + (t46 * t293 - t210 + t238 - t292) * pkin(1);
t51 = t138 + t196;
t28 = (-t192 + t282) * t29 + (t243 * t46 + t51 * t66) * t123 + (t126 * t51 + (0.4e1 * t92 - 0.2e1) * t66 * pkin(1)) * t79;
t25 = qJ(1) + atan2(t28 * t218, t27 * t218);
t246 = -m(9) - m(3) - m(7);
t155 = 0.1e1 / pkin(4);
t206 = t155 * t160 / 0.2e1;
t9 = atan2(t29 * t206, t206 * t286) + t25;
t124 = sin(pkin(9));
t127 = cos(pkin(9));
t44 = t119 * t126 - t120 * t123;
t21 = (-t27 * t44 / 0.2e1 + t28 * t288) * t280;
t22 = (t27 * t288 + t28 * t44 / 0.2e1) * t280;
t18 = t124 * t22 - t127 * t21;
t19 = t124 * t21 + t127 * t22;
t14 = atan2(t18, t19) + t25;
t227 = t155 / pkin(3) ^ 2 * t38;
t225 = t286 / 0.4e1;
t24 = cos(t25);
t207 = -pkin(3) * t24 - t80;
t198 = -m(10) * pkin(6) - mrSges(4,1) - mrSges(7,1);
t23 = sin(t25);
t187 = -pkin(3) * t23 - t287;
t184 = -m(6) - m(5) - m(11) - m(8) - m(2) - m(1) - t296;
t49 = -t122 * t126 + t269;
t48 = -t122 * t123 - t268;
t37 = qJ(1) + atan2(t44, t45);
t36 = pkin(8) + atan2(t44, -t45);
t35 = cos(t37);
t34 = sin(t37);
t33 = cos(t36);
t32 = sin(t36);
t13 = cos(t14);
t12 = sin(t14);
t11 = (t27 * t29 / 0.4e1 + t28 * t225) * t227;
t10 = (t27 * t225 - t28 * t29 / 0.4e1) * t227;
t8 = cos(t9);
t7 = sin(t9);
t6 = atan2(t18, -t19) + t14;
t5 = cos(t6);
t4 = sin(t6);
t3 = atan2(t10 * t48 + t11 * t49, -t10 * t49 + t11 * t48) + t9;
t2 = cos(t3);
t1 = sin(t3);
t15 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t184 + t246) * r_base(3)) * g(3) + (-t119 * t289 - mrSges(9,1) * t34 - mrSges(9,2) * t35 + mrSges(3,1) * t23 + mrSges(3,2) * t24 + mrSges(2,1) * t123 + mrSges(2,2) * t126 + mrSges(8,1) * t125 - mrSges(8,2) * t122 - m(11) * (-pkin(4) * t7 + t187) - m(5) * t187 + t246 * (r_base(2) - t287) + t198 * t12 - t283 * t13 - t1 * mrSges(11,1) - t2 * mrSges(11,2) - t32 * mrSges(6,1) - t33 * mrSges(6,2) + t296 * (pkin(2) * t23 + t287) + t8 * mrSges(5,2) + t7 * mrSges(5,1) + t184 * r_base(2) + t4 * mrSges(10,1) + t5 * mrSges(10,2) - mrSges(1,2)) * g(2) + (-t120 * t289 + mrSges(3,1) * t24 - mrSges(3,2) * t23 - mrSges(9,1) * t35 + mrSges(9,2) * t34 + mrSges(2,1) * t126 - mrSges(2,2) * t123 + t246 * (-t80 + r_base(1)) + t198 * t13 + t283 * t12 - t125 * mrSges(8,2) - t122 * mrSges(8,1) - m(11) * (-pkin(4) * t8 + t207) - m(5) * t207 + t1 * mrSges(11,2) - t2 * mrSges(11,1) - (m(8) * pkin(7)) + t32 * mrSges(6,2) - t33 * mrSges(6,1) + t296 * (pkin(2) * t24 + t80) + t8 * mrSges(5,1) - t7 * mrSges(5,2) + t184 * r_base(1) - t4 * mrSges(10,2) + t5 * mrSges(10,1) - mrSges(1,1)) * g(1);
U = t15;
