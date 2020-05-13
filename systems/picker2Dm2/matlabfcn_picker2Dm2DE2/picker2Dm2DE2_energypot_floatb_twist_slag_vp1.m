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
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2DE2_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 19:15:43
% EndTime: 2020-05-09 19:15:52
% DurationCPUTime: 5.40s
% Computational Cost: add. (31130->481), mult. (92594->576), div. (594->5), fcn. (15660->31), ass. (0->233)
t291 = 4 * pkin(1);
t152 = pkin(4) ^ 2;
t100 = -t152 / 0.4e1;
t157 = pkin(3) ^ 2;
t290 = t100 + t157 / 0.2e1;
t289 = 2 * pkin(7);
t124 = cos(qJ(1));
t90 = t124 ^ 2;
t288 = -0.2e1 * t90;
t136 = 0.2e1 * t157;
t287 = 0.4e1 * t157;
t162 = pkin(1) ^ 2;
t160 = t162 ^ 2;
t286 = 4 * t160;
t285 = 2 * t162;
t140 = 6 * t162;
t164 = pkin(7) ^ 2;
t148 = 2 * t164;
t169 = t157 ^ 2;
t135 = 0.5e1 * t169;
t117 = sin(pkin(8));
t118 = cos(pkin(8));
t121 = sin(qJ(1));
t45 = t117 * t121 + t118 * t124;
t284 = t45 / 0.2e1;
t78 = pkin(1) * t124;
t65 = t78 + pkin(7);
t283 = pkin(1) * t121;
t120 = sin(qJ(2));
t76 = pkin(3) * t120;
t123 = cos(qJ(2));
t77 = pkin(3) * t123;
t101 = -t152 / 0.3e1;
t102 = -t152 / 0.2e1;
t108 = 0.2e1 / 0.3e1 * t157;
t111 = -t157 / 0.3e1;
t113 = 0.4e1 / 0.3e1 * t162;
t115 = t162 / 0.2e1;
t126 = 15 * t160;
t127 = 15 * t162;
t128 = 10 * t162;
t133 = -0.2e1 * t152;
t134 = -0.5e1 * t152;
t137 = 7 * t160;
t138 = 5 * t160;
t139 = 7 * t162;
t163 = t164 ^ 2;
t144 = 3 * t163;
t145 = 8 * t164;
t146 = 4 * t164;
t151 = t152 ^ 2;
t168 = pkin(3) * t157;
t154 = t168 ^ 2;
t165 = pkin(1) * t162;
t173 = pkin(7) * t164;
t248 = t160 + t163;
t253 = t148 - t152;
t261 = t164 * t152;
t181 = t253 * t162 + t151 / 0.6e1 + t248 - t261;
t180 = 0.5e1 / 0.6e1 * t169 + t181;
t265 = t121 * t123;
t228 = pkin(3) * t265;
t202 = pkin(1) * t228;
t184 = t164 - t202;
t193 = -0.4e1 * t202;
t82 = t162 + t164;
t206 = t157 + t82;
t61 = t102 + t206;
t185 = t61 * t193;
t194 = -0.6e1 * t202;
t269 = -t169 / 0.2e1 + t151 / 0.2e1;
t198 = t144 - 0.3e1 * t261 + t269;
t104 = -0.3e1 / 0.2e1 * t152;
t147 = 3 * t164;
t259 = t104 + t147;
t267 = t154 + t82 * ((t104 + t148) * t162 - 0.3e1 / 0.2e1 * t261 + t248 + t269);
t74 = 0.10e2 / 0.3e1 * t162;
t188 = ((t74 + t253) * t157 + t180) * t194 + (t126 + (-0.9e1 * t152 + (18 * t164)) * t162 + t198) * t157 + (t127 + t259) * t169 + t267;
t141 = 3 * t162;
t254 = t141 + t164;
t209 = t157 + t254;
t189 = t209 * t77 - t283 * (0.3e1 * t157 + t82);
t103 = -0.2e1 / 0.3e1 * t152;
t112 = -0.2e1 / 0.3e1 * t157;
t210 = t103 + t82;
t255 = t128 + t148;
t258 = t112 + t164;
t80 = -t152 - t157;
t68 = t147 + t80;
t270 = t68 * t162;
t190 = -t283 * (t135 + ((5 * t162) + t68) * t136 + (t112 + t210) * t82) + (t169 + (t103 + t112 + t255) * t157 + t138 + 0.2e1 * t270 + t164 * (t103 + t258)) * t77;
t252 = t151 - t169;
t197 = (6 * t163) - 0.6e1 * t261 + t252;
t211 = t103 + t108 + t148;
t266 = t169 + (t108 + t210) * t82;
t107 = 0.4e1 / 0.3e1 * t157;
t212 = t101 + t82;
t59 = t107 + t212;
t199 = t59 * t193 + t266 + (t140 + t211) * t157;
t225 = t165 * t77;
t89 = t124 * t90;
t200 = t89 * t225;
t276 = t160 * t90 ^ 2;
t201 = t276 * t77;
t203 = 0.20e2 / 0.3e1 * t162;
t273 = t165 * t89;
t229 = 0.16e2 * t273;
t204 = pkin(7) * t229;
t251 = -t152 + t157;
t208 = t147 + t251;
t237 = pkin(7) * t273;
t213 = 0.8e1 * t237;
t263 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t215 = t121 * t263;
t109 = t157 / 0.3e1;
t99 = -t152 / 0.6e1;
t216 = t109 + t99 + t164;
t217 = t109 + t152 / 0.3e1 + t148;
t218 = t108 + 0.2e1 / 0.3e1 * t152 + t146;
t219 = t107 + 0.4e1 / 0.3e1 * t152 - (2 * t164);
t239 = 0.6e1 * t78;
t220 = pkin(7) * t239;
t240 = 0.4e1 * t78;
t221 = pkin(7) * t240;
t223 = -t283 / 0.2e1;
t226 = t162 * t77;
t264 = t123 * t124;
t227 = pkin(3) * t264;
t275 = t162 * t90;
t230 = 0.12e2 * t275;
t87 = t120 ^ 2;
t272 = t168 * t120 * t87;
t232 = -0.8e1 * t272;
t233 = 0.4e1 * t275;
t234 = -0.4e1 * t275;
t235 = 0.8e1 * t276;
t238 = -0.4e1 * t76;
t242 = 0.2e1 * t283;
t243 = pkin(7) * t78;
t244 = 4 * pkin(7);
t245 = t163 + t169;
t246 = t163 - t160;
t247 = t162 - t164;
t249 = -t157 + t164;
t250 = t152 - t164;
t256 = t115 + t164;
t257 = t162 / 0.3e1 + t164;
t260 = t164 * t162;
t262 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t268 = 0.4e1 / 0.7e1 * t164 - t152 / 0.7e1;
t271 = t169 * t87 ^ 2;
t274 = t164 * t80;
t278 = t157 * t87;
t207 = -t152 + t82;
t75 = t82 ^ 2;
t280 = t75 * (-t157 + t207);
t281 = (-t121 * t165 + t226) * t90;
t150 = 0.2e1 * pkin(3);
t66 = t76 * t289;
t196 = t66 + t207;
t64 = t76 + pkin(7);
t279 = t124 * t64;
t236 = 0.2e1 * t278;
t46 = t66 + t236 + t249;
t58 = -0.2e1 * t202;
t29 = sqrt(t46 * t234 + 0.4e1 * t247 * t278 + pkin(7) * t207 * t238 - t160 + (-0.2e1 * t157 + t250) * t285 - (t164 - (t150 + pkin(4)) * pkin(4)) * (t164 + (t150 - pkin(4)) * pkin(4)) + (-(t58 + t196) * t279 + t196 * t228) * t291);
t43 = -t169 / 0.6e1 + t181;
t72 = t111 + t164;
t47 = t72 * t58;
t53 = t76 + t65;
t81 = -0.3e1 * t157 + t164;
t56 = t81 * t213;
t57 = 0.10e2 * t270;
t62 = -t152 + t206;
t241 = 0.2e1 * t78;
t67 = pkin(7) * t241;
t70 = (t146 + t152) * t162;
t73 = -t162 / 0.3e1 + t164;
t79 = -0.30e2 * t152 + (60 * t164);
t84 = -3 * t162 + t164;
t282 = ((-0.24e2 * (0.4e1 / 0.3e1 * t275 + t67 + t73) * t271 * t283 - 0.12e2 * (-0.8e1 / 0.3e1 * t201 + ((t113 + t216) * t77 - (0.7e1 / 0.6e1 * t157 + t99 + t256) * t283) * t233 + (-t157 * t247 - 0.5e1 / 0.3e1 * t160 + t217 * t162 + t164 * (t101 + t72)) * t77 + (-t169 + (-t203 + t218) * t157 - (3 * t160) + t219 * t162 + t163) * t223 + (-t121 * t160 * t89 + ((t162 + t216) * t77 + (t136 - t247) * t223) * t78) * t244) * t278 + 0.24e2 * t72 * t201 + ((t164 + 0.5e1 / 0.2e1 * t157 + 0.3e1 / 0.2e1 * t162 + t102) * t77 + t81 * t283 / 0.2e1) * t204 - 0.6e1 * ((-0.3e1 * t169 + (-t203 + t219) * t157 + t218 * t162 + t246) * t77 - 0.2e1 * (-0.5e1 / 0.3e1 * t169 + (-t162 + t217) * t157 + t164 * (t111 + t212)) * t283) * t275 - 0.6e1 * t190 * t243 - (t154 + ((21 * t162) + t68) * t169 + (t144 + (35 * t160) + t57 + 0.2e1 * t274) * t157 + (t137 + (t134 + t145 - 0.5e1 * t157) * t162 + t164 * (-t152 + t249)) * t82) * t77 + (0.7e1 * t154 + (t139 + t68) * t135 + ((21 * t160) + (9 * t163) + t57 + 0.6e1 * t274) * t157 + t280) * t283) * t29 + (0.16e2 * (t235 + t204 + (-8 * t160 + 12 * t260) * t90 + (-12 * pkin(7) * t165 + t173 * t291) * t124 - (6 * t260) + t248) * t271 + 0.24e2 * (t258 * t235 + 0.14e2 * (-0.32e2 / 0.21e2 * (t164 + t157 / 0.4e1 + t162 / 0.4e1 - t152 / 0.8e1) * t202 + 0.5e1 / 0.42e2 * t169 + (0.16e2 / 0.21e2 * t162 + t268) * t157 + t160 / 0.7e1 + t268 * t162 + t163 - 0.3e1 / 0.7e1 * t261 + t151 / 0.42e2) * t275 + t73 * t185 - t247 * t169 + (-0.10e2 / 0.3e1 * t160 + (2 * t163) - t261 + t70) * t157 + t43 * t262 + ((-0.2e1 / 0.3e1 * t202 + t100 + t256) * t229 + (-0.8e1 / 0.3e1 * (t257 + t290) * t202 + 0.5e1 / 0.18e2 * t169 + (0.4e1 / 0.3e1 * t164 + t113 + t101) * t157 + t163 + 0.2e1 / 0.3e1 * t260 - 0.2e1 / 0.3e1 * t261 - t160 / 0.3e1 + t151 / 0.18e2) * t239) * pkin(7)) * t278 + 0.16e2 * (-0.6e1 * t164 * t157 + t245) * t276 + 0.32e2 * (t263 * t58 + t61 * t81) * t237 + 0.24e2 * (t72 * t185 - t154 + (-t74 + t250) * t169 + (t70 + t169 / 0.6e1 - t151 / 0.6e1 + t246) * t157 + t43 * t164) * t275 + 0.8e1 * t188 * t243 - 0.8e1 * ((t139 + t259) * t169 + (t137 + (t134 + (10 * t164)) * t162 + t198) * t157 + t267) * t202 + t169 ^ 2 + (t133 + t146 + (28 * t162)) * t154 + (t162 * t79 + (70 * t160) + t197) * t169 + (t197 * t140 + t252 * t148 - 0.6e1 * t163 * t152 + t79 * t160 + (28 * t165 ^ 2) + (4 * t173 ^ 2)) * t157 + t62 * t280) * t53 + (((0.4e1 * t281 + (t77 + t242) * t67 + t84 * t77 + (t102 + t209) * t242) * t232 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t157 + t115 + t99) * t77 + pkin(1) * t215) * t234 + (-0.8e1 * t200 + ((t101 + t108 + t254) * t77 - (0.8e1 / 0.3e1 * t157 + t212) * t283) * t240) * pkin(7) + t190) * t76) * t29 + (0.32e2 * (t213 + (-0.4e1 * t121 * t225 + t286 + (t287 + t133 + t145) * t162) * t90 + (-t162 + t184 + t290) * t221 + t58 * t262 + t84 * t61) * t272 + 0.8e1 * (t56 + (t263 * t61 + t47) * t230 + (t185 + (t140 + t253) * t157 + t180) * t220 + t188) * t76) * t53) * t65) / ((-0.4e1 * (-t247 * t77 + 0.2e1 * t281 + (t227 * t289 + t121 * (t157 + t285)) * pkin(1)) * t278 + 0.8e1 * pkin(7) * t200 + ((pkin(3) * t286 + 0.8e1 * t162 * t168) * t123 + 0.4e1 * t165 * t215) * t90 - 0.4e1 * t189 * t243 - (t157 * t255 + t138 + t245 + (6 * t260)) * t77 + (t135 + (t128 + 6 * t164) * t157 + t75) * t283) * t29 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t202 + 0.4e1 / 0.9e1 * t157 - t152 / 0.9e1 + t257) * t275 + t73 * t58 + t59 * t262 + (t273 + (t108 + t184 + t99) * t78) * t244) * t278 + t56 + (t263 * t59 + t47) * t230 + t199 * t220 + ((t74 + t211) * t157 + t266) * t194 + t154 + (t127 + t208) * t169 + (t208 * t140 + t251 * t148 + t126 + t144) * t157 + t75 * t62) * t53 + ((t232 * t283 + (t226 * t288 + (t77 - t283) * t67 + t189) * t238) * t29 + (0.8e1 * (t67 + t233 + t84) * t272 + 0.6e1 * (t249 * t233 + (t58 + t59) * t221 + t199) * t76) * t53) * t65);
t158 = 0.1e1 / pkin(3);
t38 = 0.1e1 / (t241 * t64 + t206 + t58 + t66);
t277 = t158 * t38;
t214 = t277 / 0.2e1;
t191 = -pkin(1) + t228;
t192 = t136 + t141 - t250;
t27 = (t121 * t64 + t227) * t29 - (t192 + t66 + t193) * t279 + t191 * t66 + t192 * t228 + (t46 * t288 - t207 + t236 - t287) * pkin(1);
t51 = t136 + t196;
t28 = (-t191 + t279) * t29 + (t241 * t46 + t51 * t64) * t121 + (t124 * t51 + (0.4e1 * t90 - 0.2e1) * t64 * pkin(1)) * t77;
t25 = qJ(1) + atan2(t28 * t214, t27 * t214);
t231 = r_base(1) - t78;
t153 = 0.1e1 / pkin(4);
t205 = t153 * t158 / 0.2e1;
t9 = atan2(t29 * t205, t205 * t282) + t25;
t122 = sin(pkin(9));
t125 = cos(pkin(9));
t44 = t117 * t124 - t118 * t121;
t21 = (-t27 * t44 / 0.2e1 + t28 * t284) * t277;
t22 = (t27 * t284 + t28 * t44 / 0.2e1) * t277;
t18 = t122 * t22 - t125 * t21;
t19 = t122 * t21 + t125 * t22;
t14 = atan2(t18, t19) + t25;
t224 = t153 / pkin(3) ^ 2 * t38;
t222 = t282 / 0.4e1;
t195 = r_base(2) - t283;
t24 = cos(t25);
t187 = -pkin(2) * t24 + t231;
t186 = -pkin(3) * t24 + t231;
t23 = sin(t25);
t183 = -pkin(2) * t23 + t195;
t182 = -pkin(3) * t23 + t195;
t49 = -t120 * t124 + t265;
t48 = -t120 * t121 - t264;
t37 = qJ(1) + atan2(t44, t45);
t36 = pkin(8) + atan2(t44, -t45);
t35 = cos(t37);
t34 = sin(t37);
t33 = cos(t36);
t32 = sin(t36);
t13 = cos(t14);
t12 = sin(t14);
t11 = (t27 * t29 / 0.4e1 + t28 * t222) * t224;
t10 = (t27 * t222 - t28 * t29 / 0.4e1) * t224;
t8 = cos(t9);
t7 = sin(t9);
t6 = atan2(t18, -t19) + t14;
t5 = cos(t6);
t4 = sin(t6);
t3 = atan2(t10 * t48 + t11 * t49, -t10 * t49 + t11 * t48) + t9;
t2 = cos(t3);
t1 = sin(t3);
t15 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t124 + rSges(2,2) * t121 + r_base(1)) + g(2) * (-rSges(2,1) * t121 - rSges(2,2) * t124 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,1) * t24 + rSges(3,2) * t23 + t231) + g(2) * (-rSges(3,1) * t23 - rSges(3,2) * t24 + t195) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t13 - rSges(4,2) * t12 + t187) + g(2) * (rSges(4,1) * t12 + rSges(4,2) * t13 + t183) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (-rSges(5,1) * t8 + rSges(5,2) * t7 + t186) + g(2) * (-rSges(5,1) * t7 - rSges(5,2) * t8 + t182) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t118 + rSges(6,1) * t33 - rSges(6,2) * t32 + r_base(1)) + g(2) * (pkin(5) * t117 + rSges(6,1) * t32 + rSges(6,2) * t33 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t13 - rSges(7,2) * t12 + t231) + g(2) * (rSges(7,1) * t12 + rSges(7,2) * t13 + t195) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t120 + rSges(8,2) * t123 + pkin(7) + r_base(1)) + g(2) * (-rSges(8,1) * t123 + rSges(8,2) * t120 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t35 - rSges(9,2) * t34 + t231) + g(2) * (rSges(9,1) * t34 + rSges(9,2) * t35 + t195) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (pkin(6) * t13 - rSges(10,1) * t5 + rSges(10,2) * t4 + t187) + g(2) * (pkin(6) * t12 - rSges(10,1) * t4 - rSges(10,2) * t5 + t183) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (-pkin(4) * t8 + rSges(11,1) * t2 - rSges(11,2) * t1 + t186) + g(2) * (-pkin(4) * t7 + rSges(11,1) * t1 + rSges(11,2) * t2 + t182) + g(3) * (r_base(3) + rSges(11,3)));
U = t15;
