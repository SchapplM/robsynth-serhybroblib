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
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1TE_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:11:22
% EndTime: 2020-05-10 00:11:40
% DurationCPUTime: 9.89s
% Computational Cost: add. (86766->519), mult. (233078->680), div. (3590->12), fcn. (61716->14), ass. (0->274)
t132 = cos(pkin(9));
t129 = sin(qJ(1));
t157 = 2 * pkin(1);
t168 = pkin(3) ^ 2;
t131 = cos(qJ(1));
t333 = sin(qJ(2));
t267 = pkin(3) * t333;
t229 = t267 + pkin(7);
t212 = t229 * t131;
t130 = cos(qJ(2));
t312 = t130 * t129;
t269 = pkin(3) * t312;
t204 = t212 - t269;
t173 = pkin(1) ^ 2;
t176 = pkin(7) ^ 2;
t272 = t333 * pkin(7);
t245 = 0.2e1 * t272;
t235 = pkin(3) * t245;
t218 = t176 + t235;
t211 = t173 + t218;
t202 = 0.1e1 / (t157 * t204 + t168 + t211);
t190 = t333 ^ 2;
t307 = t168 * t190;
t275 = 0.2e1 * t307;
t205 = -t168 + t218 + t275;
t99 = t131 ^ 2;
t203 = t99 * t205;
t163 = pkin(4) ^ 2;
t294 = t163 - t176;
t251 = -0.2e1 * t168 + t294;
t342 = 3 * t173;
t226 = t342 - t251;
t92 = t173 + t176;
t252 = -t163 + t92;
t311 = t130 * t131;
t268 = pkin(3) * t311;
t271 = pkin(1) * t312;
t345 = 0.4e1 * t168;
t348 = 0.1e1 / pkin(3);
t317 = t348 / 0.2e1;
t155 = 0.1e1 / t317;
t171 = t173 ^ 2;
t209 = -t163 + t211;
t241 = -0.4e1 * t267;
t242 = pkin(1) * t269;
t291 = t173 - t176;
t343 = 2 * t173;
t41 = sqrt(-0.4e1 * t173 * t203 - 0.4e1 * pkin(1) * (0.2e1 * (-t271 + t272) * pkin(3) + t252) * t212 + 0.4e1 * t209 * t242 + 0.4e1 * t291 * t307 + pkin(7) * t252 * t241 - t171 + t251 * t343 - (t176 - (t155 + pkin(4)) * pkin(4)) * (t176 + (t155 - pkin(4)) * pkin(4)));
t197 = t202 * ((t129 * t229 + t268) * t41 - ((t245 - 0.4e1 * t271) * pkin(3) + t226) * t212 + (t235 + t226) * t269 + (-0.2e1 * t203 + t275 - t235 - t345 - t252) * pkin(1));
t195 = t197 / 0.4e1;
t346 = 0.2e1 * t168;
t206 = t346 + t209;
t88 = pkin(1) * t131;
t283 = 0.2e1 * t88;
t87 = pkin(3) * t130;
t201 = t202 * ((pkin(1) + t204) * t41 + (t205 * t283 + t206 * t229) * t129 + (t131 * t206 + (0.4e1 * t99 - 0.2e1) * pkin(1) * t229) * t87);
t199 = t201 / 0.4e1;
t198 = t348 * t199;
t160 = pkin(5) ^ 2;
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t61 = -t126 * t129 - t127 * t131;
t340 = pkin(5) * t61;
t285 = pkin(1) * t340;
t315 = t160 - 0.2e1 * t285;
t51 = 0.1e1 / (t173 + t315);
t326 = 0.1e1 / pkin(5) * t51;
t46 = sqrt(-(-(t157 + pkin(5)) * pkin(5) + t315) * (pkin(5) * (t157 - pkin(5)) + t315));
t60 = t126 * t131 - t127 * t129;
t330 = t46 * t60;
t52 = t160 - t285;
t56 = -pkin(1) + t340;
t42 = -pkin(5) * t330 - 0.2e1 * t52 * t56;
t347 = 0.2e1 * t60;
t44 = pkin(5) * t347 * t52 - t46 * t56;
t193 = (t195 * t348 * t42 + t44 * t198) * t326;
t196 = t348 * t197;
t27 = (-t44 * t196 / 0.4e1 + t42 * t198) * t326;
t334 = sin(pkin(9));
t24 = t132 * t193 + t27 * t334;
t159 = 0.1e1 / pkin(6);
t158 = pkin(6) ^ 2;
t156 = 2 * pkin(2);
t339 = pkin(6) * t24;
t23 = t339 * t156;
t316 = t158 + t23;
t327 = t159 / ((pkin(2) ^ 2) + t316);
t194 = -t196 / 0.2e1;
t200 = t348 * t201;
t336 = t129 / 0.2e1;
t34 = t131 * t194 + t200 * t336;
t33 = t129 * t194 - t131 * t200 / 0.2e1;
t341 = t33 / 0.2e1;
t14 = sqrt(-(-(t156 + pkin(6)) * pkin(6) + t316) * (pkin(6) * (t156 - pkin(6)) + t316));
t21 = t23 + 0.2e1 * t158;
t22 = -pkin(2) - t339;
t25 = t132 * t27 - t193 * t334;
t338 = pkin(6) * t25;
t6 = t14 * t338 - t21 * t22;
t7 = -t14 * t22 - t21 * t338;
t350 = (t34 * t7 / 0.2e1 + t6 * t341) * t327;
t357 = t24 * t350;
t164 = 0.1e1 / pkin(4);
t108 = -t163 / 0.6e1;
t109 = -t163 / 0.4e1;
t110 = -t163 / 0.3e1;
t111 = -t163 / 0.2e1;
t117 = 0.2e1 / 0.3e1 * t168;
t120 = -t168 / 0.3e1;
t122 = 0.4e1 / 0.3e1 * t173;
t124 = t173 / 0.2e1;
t133 = 15 * t171;
t134 = 15 * t173;
t135 = 10 * t173;
t140 = -0.2e1 * t163;
t141 = -0.5e1 * t163;
t181 = t168 ^ 2;
t142 = 0.5e1 * t181;
t143 = 7 * t171;
t144 = 5 * t171;
t145 = 7 * t173;
t146 = 6 * t173;
t175 = t176 ^ 2;
t149 = 0.3e1 * t175;
t150 = 0.8e1 * t176;
t151 = 0.4e1 * t176;
t153 = 0.2e1 * t176;
t162 = t163 ^ 2;
t180 = pkin(3) * t168;
t165 = t180 ^ 2;
t177 = pkin(1) * t173;
t185 = pkin(7) * t176;
t292 = t171 + t175;
t297 = t153 - t163;
t306 = t176 * t163;
t210 = t297 * t173 + t162 / 0.6e1 + t292 - t306;
t208 = 0.5e1 / 0.6e1 * t181 + t210;
t217 = t176 - t242;
t231 = -0.4e1 * t242;
t250 = t168 + t92;
t74 = t111 + t250;
t219 = t74 * t231;
t304 = t162 / 0.2e1 - t181 / 0.2e1;
t228 = -0.3e1 * t306 + t149 + t304;
t232 = -0.6e1 * t242;
t113 = -0.3e1 / 0.2e1 * t163;
t152 = 0.3e1 * t176;
t303 = t113 + t152;
t314 = t165 + t92 * ((t113 + t153) * t173 - 0.3e1 / 0.2e1 * t306 + t292 + t304);
t85 = 0.10e2 / 0.3e1 * t173;
t223 = ((t85 + t297) * t168 + t208) * t232 + (t133 + (-0.9e1 * t163 + 0.18e2 * t176) * t173 + t228) * t168 + (t134 + t303) * t181 + t314;
t288 = t176 + t342;
t249 = t168 + t288;
t332 = pkin(1) * t129;
t224 = t249 * t87 - t332 * (0.3e1 * t168 + t92);
t112 = -0.2e1 / 0.3e1 * t163;
t121 = -0.2e1 / 0.3e1 * t168;
t254 = t112 + t92;
t298 = t135 + t153;
t302 = t121 + t176;
t90 = -t163 - t168;
t79 = t152 + t90;
t318 = t79 * t173;
t225 = -t332 * (t142 + ((5 * t173) + t79) * t346 + (t121 + t254) * t92) + (t181 + (t112 + t121 + t298) * t168 + t144 + 0.2e1 * t318 + t176 * (t112 + t302)) * t87;
t296 = t162 - t181;
t227 = -0.6e1 * t306 + 0.6e1 * t175 + t296;
t230 = -0.2e1 * t242;
t255 = t112 + t117 + t153;
t313 = t181 + (t117 + t254) * t92;
t116 = 0.4e1 / 0.3e1 * t168;
t256 = t110 + t92;
t72 = t116 + t256;
t237 = t72 * t231 + t313 + (t146 + t255) * t168;
t98 = t131 * t99;
t321 = t177 * t98;
t239 = t321 * t87;
t325 = t171 * t99 ^ 2;
t240 = t325 * t87;
t246 = 0.20e2 / 0.3e1 * t173;
t273 = 0.16e2 * t321;
t248 = pkin(7) * t273;
t295 = -t163 + t168;
t253 = t152 + t295;
t118 = t168 / 0.3e1;
t257 = t108 + t118 + t176;
t258 = t163 / 0.3e1 + t118 + t153;
t259 = 0.2e1 / 0.3e1 * t163 + t117 + t151;
t260 = 0.4e1 / 0.3e1 * t163 + t116 - 0.2e1 * t176;
t280 = pkin(7) * t321;
t261 = 0.8e1 * t280;
t310 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t262 = t129 * t310;
t281 = 0.6e1 * t88;
t263 = pkin(7) * t281;
t282 = 0.4e1 * t88;
t264 = pkin(7) * t282;
t265 = -t332 / 0.2e1;
t270 = t173 * t87;
t324 = t173 * t99;
t274 = 0.12e2 * t324;
t320 = t180 * t333 * t190;
t277 = -0.8e1 * t320;
t278 = 0.4e1 * t324;
t279 = 0.8e1 * t325;
t284 = 0.2e1 * t332;
t286 = pkin(7) * t88;
t287 = 0.4e1 * pkin(7);
t289 = t175 + t181;
t290 = t175 - t171;
t293 = -t168 + t176;
t299 = 0.4e1 / 0.7e1 * t176 - t163 / 0.7e1;
t300 = t124 + t176;
t301 = t173 / 0.3e1 + t176;
t305 = t176 * t173;
t309 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t319 = t181 * t190 ^ 2;
t322 = t176 * t90;
t86 = t92 ^ 2;
t328 = t86 * (-t168 + t252);
t329 = (-t129 * t177 + t270) * t99;
t344 = 4 * t171;
t349 = t109 + t168 / 0.2e1;
t59 = -t181 / 0.6e1 + t210;
t83 = t120 + t176;
t62 = t83 * t230;
t67 = t88 + t229;
t91 = -0.3e1 * t168 + t176;
t70 = t91 * t261;
t71 = 0.10e2 * t318;
t75 = -t163 + t250;
t77 = t88 + pkin(7);
t78 = pkin(7) * t283;
t81 = (t151 + t163) * t173;
t84 = -t173 / 0.3e1 + t176;
t89 = -0.30e2 * t163 + 0.60e2 * t176;
t94 = -(3 * t173) + t176;
t331 = ((-0.24e2 * (0.4e1 / 0.3e1 * t324 + t78 + t84) * t319 * t332 - 0.12e2 * (-0.8e1 / 0.3e1 * t240 + ((t122 + t257) * t87 - (0.7e1 / 0.6e1 * t168 + t108 + t300) * t332) * t278 + (-t168 * t291 - 0.5e1 / 0.3e1 * t171 + t258 * t173 + t176 * (t110 + t83)) * t87 + (-t181 + (-t246 + t259) * t168 - (3 * t171) + t260 * t173 + t175) * t265 + (-t129 * t171 * t98 + ((t173 + t257) * t87 + (t346 - t291) * t265) * t88) * t287) * t307 + 0.24e2 * t83 * t240 + ((t176 + 0.5e1 / 0.2e1 * t168 + 0.3e1 / 0.2e1 * t173 + t111) * t87 + t91 * t332 / 0.2e1) * t248 - 0.6e1 * ((-0.3e1 * t181 + (-t246 + t260) * t168 + t259 * t173 + t290) * t87 - 0.2e1 * (-0.5e1 / 0.3e1 * t181 + (-t173 + t258) * t168 + t176 * (t120 + t256)) * t332) * t324 - 0.6e1 * t225 * t286 - (t165 + ((21 * t173) + t79) * t181 + (t149 + (35 * t171) + t71 + 0.2e1 * t322) * t168 + (t143 + (t141 + t150 - 0.5e1 * t168) * t173 + t176 * (-t163 + t293)) * t92) * t87 + (0.7e1 * t165 + (t145 + t79) * t142 + ((21 * t171) + 0.9e1 * t175 + t71 + 0.6e1 * t322) * t168 + t328) * t332) * t41 + (0.16e2 * (t279 + t248 + (-(8 * t171) + 0.12e2 * t305) * t99 + (0.4e1 * pkin(1) * t185 - 0.12e2 * pkin(7) * t177) * t131 - 0.6e1 * t305 + t292) * t319 + 0.24e2 * (t302 * t279 + 0.14e2 * (-0.32e2 / 0.21e2 * (t176 + t168 / 0.4e1 + t173 / 0.4e1 - t163 / 0.8e1) * t242 + 0.5e1 / 0.42e2 * t181 + (0.16e2 / 0.21e2 * t173 + t299) * t168 + t171 / 0.7e1 + t299 * t173 + t175 - 0.3e1 / 0.7e1 * t306 + t162 / 0.42e2) * t324 + t84 * t219 - t291 * t181 + (-0.10e2 / 0.3e1 * t171 + 0.2e1 * t175 - t306 + t81) * t168 + t59 * t309 + ((-0.2e1 / 0.3e1 * t242 + t109 + t300) * t273 + (-0.8e1 / 0.3e1 * (t301 + t349) * t242 + 0.5e1 / 0.18e2 * t181 + (0.4e1 / 0.3e1 * t176 + t122 + t110) * t168 + t175 + 0.2e1 / 0.3e1 * t305 - 0.2e1 / 0.3e1 * t306 - t171 / 0.3e1 + t162 / 0.18e2) * t281) * pkin(7)) * t307 + 0.16e2 * (-0.6e1 * t176 * t168 + t289) * t325 + 0.32e2 * (t230 * t310 + t74 * t91) * t280 + 0.24e2 * (t83 * t219 - t165 + (-t85 + t294) * t181 + (t81 + t181 / 0.6e1 - t162 / 0.6e1 + t290) * t168 + t59 * t176) * t324 + 0.8e1 * t223 * t286 - 0.8e1 * ((t145 + t303) * t181 + (t143 + (t141 + 0.10e2 * t176) * t173 + t228) * t168 + t314) * t242 + t181 ^ 2 + (t140 + t151 + (28 * t173)) * t165 + (t173 * t89 + (70 * t171) + t227) * t181 + (t227 * t146 + t296 * t153 - 0.6e1 * t175 * t163 + t89 * t171 + (28 * t177 ^ 2) + 0.4e1 * t185 ^ 2) * t168 + t75 * t328) * t67 + (((0.4e1 * t329 + (t87 + t284) * t78 + t94 * t87 + (t111 + t249) * t284) * t277 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t168 + t124 + t108) * t130 * t155 + pkin(1) * t262) * t324 + (-0.8e1 * t239 + ((t110 + t117 + t288) * t87 - (0.8e1 / 0.3e1 * t168 + t256) * t332) * t282) * pkin(7) + t225) * t267) * t41 + (0.32e2 * (t261 + (-0.4e1 * t177 * t269 + t344 + (t345 + t140 + t150) * t173) * t99 + (-t173 + t217 + t349) * t264 + t230 * t309 + t94 * t74) * t320 + 0.8e1 * (t70 + (t310 * t74 + t62) * t274 + (t219 + (t146 + t297) * t168 + t208) * t263 + t223) * t267) * t67) * t77) / ((-0.4e1 * (-t291 * t87 + 0.2e1 * t329 + (0.2e1 * pkin(7) * t268 + t129 * (t168 + t343)) * pkin(1)) * t307 + 0.8e1 * pkin(7) * t239 + ((pkin(3) * t344 + 0.8e1 * t173 * t180) * t130 + 0.4e1 * t177 * t262) * t99 - 0.4e1 * t224 * t286 - (t168 * t298 + t144 + t289 + 0.6e1 * t305) * t87 + (t142 + (t135 + 0.6e1 * t176) * t168 + t86) * t332) * t41 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t242 + 0.4e1 / 0.9e1 * t168 - t163 / 0.9e1 + t301) * t324 + t84 * t230 + t72 * t309 + (t321 + (t108 + t117 + t217) * t88) * t287) * t307 + t70 + (t310 * t72 + t62) * t274 + t237 * t263 + ((t85 + t255) * t168 + t313) * t232 + t165 + (t134 + t253) * t181 + (t253 * t146 + t295 * t153 + t133 + t149) * t168 + t86 * t75) * t67 + ((t277 * t332 + (-0.2e1 * t99 * t270 + (t87 - t332) * t78 + t224) * t241) * t41 + (0.8e1 * (t78 + t278 + t94) * t320 + 0.6e1 * (t293 * t278 + (t72 + t230) * t264 + t237) * t267) * t67) * t77);
t238 = t317 * t331;
t213 = (-t33 * t348 * t41 / 0.2e1 + t34 * t238) * t164;
t351 = (t317 * t34 * t41 + t238 * t33) * t164;
t308 = t164 / pkin(3) ^ 2;
t15 = (t195 * t331 - t41 * t201 / 0.4e1) * t308;
t16 = (t195 * t41 + t199 * t331) * t308;
t63 = -t129 * t333 - t311;
t64 = -t131 * t333 + t312;
t8 = t15 * t63 + t16 * t64;
t9 = t15 * t64 - t16 * t63;
t356 = t213 * t9 + t351 * t8;
t355 = -t213 * t8 + t351 * t9;
t215 = (-t34 * t6 / 0.2e1 + t341 * t7) * t327;
t352 = t215 * t24;
t337 = -t126 / 0.2e1;
t335 = t131 / 0.2e1;
t323 = 0.1e1 / pkin(1) * t51;
t276 = r_base(1) - t88;
t266 = t14 * t159 / pkin(2);
t247 = -t34 * t24 - t25 * t33;
t244 = t34 * pkin(3) + t276;
t243 = t34 * pkin(2) + t276;
t236 = r_base(2) - t332;
t234 = -t266 / 0.2e1;
t233 = t266 / 0.2e1;
t222 = t24 * t33 - t25 * t34;
t221 = t33 * pkin(3) + t236;
t220 = t33 * pkin(2) + t236;
t57 = -pkin(1) * t61 + pkin(5);
t53 = t173 - t285;
t45 = pkin(1) * t347 * t53 + t46 * t57;
t43 = -pkin(1) * t330 + 0.2e1 * t53 * t57;
t39 = (t42 * t335 - t129 * t44 / 0.2e1) * t326;
t38 = (t335 * t44 + t336 * t42) * t326;
t37 = (-t127 * t43 / 0.2e1 + t45 * t337) * t323;
t36 = (t43 * t337 + t127 * t45 / 0.2e1) * t323;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-rSges(2,1) * t131 + rSges(2,2) * t129 + r_base(1)) + g(2) * (-rSges(2,1) * t129 - rSges(2,2) * t131 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t34 - rSges(3,2) * t33 + t276) + g(2) * (rSges(3,1) * t33 + rSges(3,2) * t34 + t236) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t215 + rSges(4,2) * t350 + t243) + g(2) * (-rSges(4,1) * t350 + rSges(4,2) * t215 + t220) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t213 - rSges(5,2) * t351 + t244) + g(2) * (rSges(5,1) * t351 + rSges(5,2) * t213 + t221) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (pkin(5) * t127 + rSges(6,1) * t37 - rSges(6,2) * t36 + r_base(1)) + g(2) * (pkin(5) * t126 + rSges(6,1) * t36 + rSges(6,2) * t37 + r_base(2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (rSges(7,1) * t247 + rSges(7,2) * t222 + t276) + g(2) * (-rSges(7,1) * t222 + rSges(7,2) * t247 + t236) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (rSges(8,1) * t333 + t130 * rSges(8,2) + pkin(7) + r_base(1)) + g(2) * (-t130 * rSges(8,1) + rSges(8,2) * t333 + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (rSges(9,1) * t39 - rSges(9,2) * t38 + t276) + g(2) * (rSges(9,1) * t38 + rSges(9,2) * t39 + t236) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (t215 * pkin(6) + (t234 * t350 + t352) * rSges(10,1) + (t215 * t233 + t357) * rSges(10,2) + t243) + g(2) * (-t350 * pkin(6) + (t215 * t234 - t357) * rSges(10,1) + (-t233 * t350 + t352) * rSges(10,2) + t220) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (t213 * pkin(4) + rSges(11,1) * t356 - t355 * rSges(11,2) + t244) + g(2) * (t351 * pkin(4) + rSges(11,1) * t355 + rSges(11,2) * t356 + t221) + g(3) * (r_base(3) + rSges(11,3)));
U = t1;
