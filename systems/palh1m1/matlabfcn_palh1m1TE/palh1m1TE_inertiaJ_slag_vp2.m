% Calculate joint inertia matrix for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_inertiaJ_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1TE_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:25:32
% EndTime: 2020-04-12 20:32:47
% DurationCPUTime: 119.38s
% Computational Cost: add. (3208644->399), mult. (4845948->724), div. (212810->34), fcn. (3072050->26), ass. (0->317)
t210 = pkin(7) ^ 2;
t218 = pkin(1) ^ 2;
t193 = sin(qJ(2));
t197 = cos(qJ(2));
t198 = cos(pkin(19));
t352 = sin(pkin(19));
t169 = t193 * t198 - t197 * t352;
t342 = pkin(7) * t169;
t386 = -2 * pkin(1);
t288 = t342 * t386 + t218;
t158 = t210 + t288;
t284 = pkin(3) ^ 2 - pkin(8) ^ 2;
t151 = t158 + t284;
t161 = pkin(1) - t342;
t375 = pkin(7) + pkin(8);
t376 = pkin(7) - pkin(8);
t144 = (pkin(3) + t375) * (-pkin(3) + t376) + t288;
t145 = (-pkin(3) + t375) * (pkin(3) + t376) + t288;
t220 = sqrt(-t145 * t144);
t172 = t193 * t352 + t197 * t198;
t341 = pkin(7) * t172;
t133 = t151 * t341 + t161 * t220;
t196 = cos(qJ(3));
t291 = t196 * t133;
t298 = t172 * t220;
t132 = -pkin(7) * t298 + t151 * t161;
t192 = sin(qJ(3));
t294 = t192 * t132;
t156 = 0.1e1 / t158;
t215 = 0.1e1 / pkin(3);
t301 = t156 * t215;
t121 = (t294 / 0.2e1 + t291 / 0.2e1) * t301;
t292 = t196 * t132;
t293 = t192 * t133;
t122 = (-t292 / 0.2e1 + t293 / 0.2e1) * t301;
t180 = pkin(23) + pkin(22);
t177 = sin(t180);
t178 = cos(t180);
t251 = t121 * t177 + t178 * t122;
t189 = cos(pkin(21));
t213 = pkin(4) ^ 2;
t212 = pkin(5) ^ 2;
t98 = t121 * t178 - t122 * t177;
t367 = t98 * pkin(5);
t399 = 2 * pkin(4);
t311 = t367 * t399 + t212;
t94 = t213 + t311;
t93 = 0.1e1 / t94 ^ 2;
t368 = pkin(5) * t93;
t282 = pkin(4) * t368;
t256 = t189 * t282;
t242 = t251 * t256;
t185 = sin(pkin(21));
t257 = t185 * t282;
t243 = t251 * t257;
t92 = 0.1e1 / t94;
t319 = t189 * t92;
t265 = t319 / 0.2e1;
t266 = -t319 / 0.2e1;
t360 = t185 / 0.2e1;
t267 = t92 * t360;
t205 = 0.1e1 / pkin(11);
t316 = t205 * t92;
t374 = (-pkin(9) - pkin(11));
t88 = ((pkin(4) - t374) * (pkin(4) + t374)) + t311;
t373 = (pkin(11) - pkin(9));
t89 = ((pkin(4) - t373) * (pkin(4) + t373)) + t311;
t221 = sqrt(-t89 * t88);
t392 = t251 * t221;
t286 = pkin(9) ^ 2 - pkin(11) ^ 2;
t91 = t94 - t286;
t96 = pkin(4) * t98 + pkin(5);
t66 = -pkin(4) * t392 + t91 * t96;
t68 = pkin(4) * t251 * t91 + t221 * t96;
t51 = (t66 * t360 + t68 * t189 / 0.2e1) * t316;
t46 = t51 ^ 2;
t52 = (-t66 * t189 / 0.2e1 + t68 * t360) * t316;
t223 = t52 ^ 2;
t48 = 0.1e1 / t223;
t317 = t205 / (t46 * t48 + 0.1e1);
t332 = t48 * t51;
t77 = 0.1e1 / t221;
t379 = -t77 / 0.2e1;
t268 = t251 * t379;
t233 = pkin(5) * (t88 + t89) * t399;
t69 = t251 * t233;
t229 = -t98 * t221 + t268 * t69;
t262 = -0.2e1 * pkin(5) * t96 - t91;
t43 = (t251 * t262 + t229) * pkin(4);
t346 = pkin(5) * t251;
t255 = -0.2e1 * t213 * t346;
t269 = t77 * t96 / 0.2e1;
t44 = t69 * t269 + t251 * t255 + (t91 * t98 - t392) * pkin(4);
t47 = 0.1e1 / t52;
t18 = 0.1e1 + ((t242 * t68 + t243 * t66 + t265 * t44 + t267 * t43) * t47 - (-t242 * t66 + t243 * t68 + t266 * t43 + t267 * t44) * t332) * t317;
t370 = pkin(5) * t51;
t14 = pkin(12) * t18 + t370;
t191 = sin(qJ(4));
t195 = cos(qJ(4));
t391 = t191 ^ 2 + t195 ^ 2;
t246 = 0.2e1 * t391;
t402 = t14 * t246;
t216 = pkin(2) ^ 2;
t211 = pkin(6) ^ 2;
t285 = t211 + t218;
t400 = pkin(13) ^ 2;
t263 = t285 - t400;
t186 = sin(pkin(20));
t190 = cos(pkin(20));
t166 = t186 * t196 + t190 * t192;
t344 = pkin(6) * t166;
t279 = pkin(1) * t344;
t148 = t216 - t263 - 0.2e1 * t279;
t401 = t148 ^ 2;
t160 = 0.2e1 * t279;
t289 = t160 + t211;
t377 = (-pkin(2) - pkin(13));
t142 = ((pkin(1) - t377) * (pkin(1) + t377)) + t289;
t372 = (pkin(13) - pkin(2));
t143 = ((pkin(1) - t372) * (pkin(1) + t372)) + t289;
t305 = t143 * t142;
t219 = sqrt(-t305);
t398 = t219 / 0.2e1;
t397 = Ifges(5,3) * t18;
t149 = t160 + t216 + t263;
t159 = -pkin(1) - t344;
t167 = t186 * t192 - t190 * t196;
t343 = pkin(6) * t167;
t130 = t149 * t343 - t159 * t219;
t155 = t160 + t285;
t152 = 0.1e1 / t155;
t217 = 0.1e1 / pkin(2);
t304 = t152 * t217;
t357 = -t193 / 0.2e1;
t300 = t167 * t219;
t273 = pkin(6) * t300;
t129 = -t149 * t159 - t273;
t364 = t129 / 0.2e1;
t113 = (t130 * t357 + t197 * t364) * t304;
t114 = (-t197 * t130 / 0.2e1 + t129 * t357) * t304;
t202 = 0.1e1 / pkin(13);
t290 = t202 * t217;
t353 = -t219 / 0.2e1;
t363 = -t148 / 0.2e1;
t86 = (t113 * t363 + t114 * t353) * t290;
t87 = (t113 * t398 + t114 * t363) * t290;
t396 = Ifges(10,5) * t86 + Ifges(10,6) * t87;
t90 = t94 + t286;
t95 = -pkin(4) - t367;
t67 = -t221 * t95 + t346 * t90;
t227 = pkin(4) * (-t212 * t251 * t92 + t368 * t67);
t281 = pkin(1) * t341;
t306 = 0.2e1 / t220 * (t144 + t145) * t281;
t383 = -0.2e1 * t172 ^ 2;
t112 = t161 * t306 / 0.2e1 + t210 * pkin(1) * t383 + (-t151 * t169 - t298) * pkin(7);
t254 = 0.1e1 / t158 ^ 2 * t281;
t355 = -t196 / 0.2e1;
t260 = -t306 / 0.2e1;
t299 = t169 * t220;
t107 = (t299 + (t161 * t386 - t151 + t260) * t172) * pkin(7);
t366 = -t107 / 0.2e1;
t81 = ((t112 * t355 + t192 * t366) * t156 + (-t291 - t294) * t254) * t215;
t365 = t112 / 0.2e1;
t82 = ((t107 * t355 + t192 * t365) * t156 + (-t292 + t293) * t254) * t215;
t71 = t177 * t82 + t178 * t81;
t54 = t71 * t233;
t72 = -t177 * t81 + t178 * t82;
t230 = -t72 * t221 + t268 * t54;
t207 = 0.1e1 / pkin(9);
t65 = -pkin(5) * t392 - t90 * t95;
t64 = 0.1e1 / t65 ^ 2;
t253 = pkin(9) * t207 / (t64 * t67 ^ 2 + 0.1e1) * t94;
t231 = pkin(5) * t64 * t67 * t253;
t232 = pkin(4) * (t65 * t93 + t92 * t95);
t241 = 0.1e1 / t65 * t253;
t270 = t95 * t379;
t314 = t221 * t71;
t378 = t92 / 0.2e1;
t187 = cos(pkin(23));
t296 = t187 * t132;
t183 = sin(pkin(23));
t297 = t183 * t133;
t118 = (-t296 / 0.2e1 + t297 / 0.2e1) * t301;
t222 = t118 ^ 2;
t116 = 0.1e1 / t222;
t295 = t187 * t133;
t308 = t132 * t183;
t119 = (t295 / 0.2e1 + t308 / 0.2e1) * t301;
t117 = t119 ^ 2;
t361 = t183 / 0.2e1;
t57 = 0.1e1 + (-((t112 * t361 + t187 * t366) * t156 + (-t296 + t297) * t254) * t119 * t116 + ((t107 * t361 + t187 * t365) * t156 + (t295 + t308) * t254) / t118) / (t116 * t117 + 0.1e1) * t215;
t19 = 0.2e1 * ((t54 * t270 + (t72 * t90 - t314) * pkin(5)) * t378 + t71 * t227) * t241 - 0.2e1 * ((-t71 * t90 + t230) * t378 + t71 * t232) * t231 + t57;
t395 = Ifges(11,3) * t19;
t394 = Ifges(9,3) + m(10) * (-t305 / 0.2e1 + t401 / 0.2e1) / t400 / 0.2e1;
t271 = t202 * t219 * mrSges(10,2);
t272 = t148 * t202 * mrSges(10,1);
t390 = t271 / 0.2e1 - t272 / 0.2e1;
t389 = t271 - t272;
t387 = Ifges(9,5) * t113 + Ifges(9,6) * t114 + (t87 * t353 + t148 * t86 / 0.2e1) * t202 * mrSges(10,3);
t369 = pkin(5) * t52;
t15 = -pkin(10) * t18 - t369;
t384 = 0.2e1 * t15;
t381 = m(6) / 0.2e1;
t170 = t192 * t197 + t193 * t196;
t171 = -t192 * t193 + t196 * t197;
t38 = -t170 * t51 + t171 * t52;
t380 = -t38 / 0.2e1;
t371 = pkin(4) * t57;
t362 = t152 / 0.2e1;
t188 = cos(pkin(22));
t359 = -t188 / 0.2e1;
t358 = -t191 / 0.2e1;
t356 = t195 / 0.2e1;
t199 = cos(pkin(18));
t354 = t199 / 0.2e1;
t351 = pkin(1) * t118;
t350 = pkin(1) * t119;
t349 = pkin(1) * t167;
t348 = pkin(1) * t192;
t347 = pkin(1) * t196;
t345 = pkin(6) * t130;
t340 = mrSges(5,3) * t38;
t244 = t71 * t256;
t245 = t71 * t257;
t34 = (t262 * t71 + t230) * pkin(4);
t35 = t54 * t269 + t71 * t255 + (t72 * t91 - t314) * pkin(4);
t12 = 0.1e1 + ((t244 * t68 + t245 * t66 + t265 * t35 + t267 * t34) * t47 - (-t244 * t66 + t245 * t68 + t266 * t34 + t267 * t35) * t332) * t317;
t339 = mrSges(6,3) * t12;
t338 = mrSges(6,3) * t18;
t337 = Ifges(6,6) * t38;
t176 = pkin(5) + t348;
t40 = t176 * t52 + t347 * t51;
t334 = t40 * mrSges(5,1);
t41 = t51 * t176 - t347 * t52;
t333 = t41 * mrSges(5,2);
t322 = t12 * t195;
t323 = t12 * t191;
t331 = Ifges(6,5) * t323 + Ifges(6,6) * t322;
t320 = t18 * t195;
t321 = t18 * t191;
t330 = Ifges(6,5) * t321 + Ifges(6,6) * t320;
t39 = t170 * t52 + t171 * t51;
t318 = t195 * t39;
t329 = Ifges(6,5) * t318 - Ifges(6,3) * t38;
t128 = 0.1e1 / t129 ^ 2;
t153 = 0.1e1 / t155 ^ 2;
t307 = (t142 + t143) * pkin(1) * t343 / t398;
t261 = -t307 / 0.2e1;
t74 = (((t159 * t261 + (t149 * t166 - t300) * pkin(6)) * t362 + (-t152 * t167 * t211 + t153 * t345) * t349) / t364 - 0.2e1 * ((-t166 * t219 + (t261 - t149) * t167) * t362 + (t129 * t153 + t152 * t159) * t349) * t128 * t345) * pkin(2) / (t128 * t130 ^ 2 + 0.1e1) * t155 * t217;
t184 = sin(pkin(22));
t315 = t207 * t92;
t49 = (t184 * t65 / 0.2e1 + t67 * t359) * t315;
t50 = (t65 * t359 - t184 * t67 / 0.2e1) * t315;
t55 = t188 * t371 + t351;
t56 = t184 * t371 + t350;
t31 = -t49 * t56 + t50 * t55;
t328 = mrSges(11,1) * t31;
t327 = mrSges(6,1) * t195;
t30 = t49 * t55 + t50 * t56;
t326 = mrSges(11,2) * t30;
t325 = Ifges(6,4) * t191;
t324 = Ifges(6,4) * t195;
t174 = t193 * pkin(1) - pkin(16);
t194 = sin(pkin(18));
t303 = t156 * t194;
t209 = 0.1e1 / pkin(8);
t302 = t156 * t209;
t287 = Ifges(4,5) * t170 + Ifges(4,6) * t171;
t278 = mrSges(5,2) * t370;
t277 = mrSges(5,1) * t369;
t275 = mrSges(4,1) * t348;
t274 = mrSges(4,2) * t347;
t259 = t156 * t354;
t154 = -pkin(5) * t171 + t174;
t240 = t194 * t254;
t239 = t199 * t254;
t238 = mrSges(6,2) * t191 - t327;
t237 = mrSges(6,1) * t191 + mrSges(6,2) * t195;
t236 = Ifges(6,1) * t191 + t324;
t235 = Ifges(6,2) * t195 + t325;
t27 = -mrSges(6,3) * t191 * t39 + mrSges(6,2) * t38;
t28 = -mrSges(6,1) * t38 - mrSges(6,3) * t318;
t234 = -t191 * t28 + t195 * t27;
t102 = t118 * t197 - t119 * t193;
t103 = -t118 * t193 - t119 * t197;
t32 = t102 * t50 + t103 * t49;
t33 = -t102 * t49 + t103 * t50;
t228 = t32 * Ifges(11,5) + t33 * Ifges(11,6);
t21 = -t337 + (-Ifges(6,2) * t191 + t324) * t39;
t22 = -Ifges(6,5) * t38 + (Ifges(6,1) * t195 - t325) * t39;
t224 = t21 * t356 + t191 * t22 / 0.2e1 + t39 * Ifges(5,5) + t38 * Ifges(5,6);
t162 = pkin(1) * t169 - pkin(7);
t150 = t158 - t284;
t147 = 0.1e1 / t401;
t134 = pkin(1) * t172 * t150 - t162 * t220;
t131 = -pkin(1) * t298 - t150 * t162;
t124 = (t134 * t354 + t131 * t194 / 0.2e1) * t302;
t123 = (t131 * t354 - t194 * t134 / 0.2e1) * t302;
t120 = 0.1e1 / t123 ^ 2;
t111 = t162 * t260 + t218 * pkin(7) * t383 + (-t150 * t169 - t298) * pkin(1);
t108 = -pkin(2) * t114 - pkin(16);
t106 = (t299 + (0.2e1 * t162 * pkin(7) - t150 + t260) * t172) * pkin(1);
t78 = (-t102 * t184 - t103 * t188) * pkin(4) + t174;
t70 = (0.1e1 / t148 * t307 / 0.2e1 + t147 * t273 * t386) / (-t147 * t305 + 0.1e1) + t74;
t58 = ((t111 * t259 + t134 * t239 + t106 * t303 / 0.2e1 + t131 * t240) / t123 - (t106 * t259 + t131 * t239 - t111 * t303 / 0.2e1 - t134 * t240) * t124 * t120) / (t120 * t124 ^ 2 + 0.1e1) * t209;
t26 = t237 * t39;
t25 = -pkin(10) * t38 - pkin(12) * t39 + t154;
t20 = 0.2e1 * ((t69 * t270 + (t90 * t98 - t392) * pkin(5)) * t378 + t251 * t227) * t241 - 0.2e1 * ((-t251 * t90 + t229) * t378 + t251 * t232) * t231;
t9 = pkin(12) * t12 + t41;
t8 = -pkin(10) * t12 - t40;
t6 = t236 * t18;
t5 = t235 * t18;
t4 = t238 * t18;
t3 = t236 * t12;
t2 = t235 * t12;
t1 = t238 * t12;
t7 = [Ifges(9,2) * t114 ^ 2 + (-0.2e1 * Ifges(3,4) * t197 + Ifges(3,2) * t193) * t193 + Ifges(4,2) * t171 ^ 2 + 0.2e1 * (-mrSges(4,1) * t171 - mrSges(8,1) * t103 + mrSges(4,2) * t170 + mrSges(8,2) * t102) * t174 - 0.2e1 * (mrSges(3,1) * t193 - mrSges(9,1) * t114 + mrSges(3,2) * t197 + mrSges(9,2) * t113) * pkin(16) + (Ifges(9,1) * t113 + 0.2e1 * Ifges(9,4) * t114) * t113 + Ifges(3,1) * t197 ^ 2 + Ifges(8,2) * t103 ^ 2 + (Ifges(8,1) * t102 + 0.2e1 * Ifges(8,4) * t103) * t102 + (Ifges(4,1) * t170 + 0.2e1 * Ifges(4,4) * t171) * t170 + (t246 * t381 * t25 + 0.2e1 * t191 * t27 + 0.2e1 * t195 * t28) * t25 + (m(8) + m(4)) * t174 ^ 2 + (m(3) + m(9)) * pkin(16) ^ 2 + m(7) * pkin(15) ^ 2 + (0.2e1 * mrSges(5,2) * t154 + Ifges(5,1) * t39 + 0.2e1 * Ifges(5,4) * t38 + t195 * t22 + (-t21 + t337) * t191) * t39 + (-0.2e1 * mrSges(5,1) * t154 + Ifges(5,2) * t38 - t329) * t38 + m(5) * t154 ^ 2 + 0.2e1 * pkin(15) * (-mrSges(7,1) * t123 + mrSges(7,2) * t124) + t123 * (Ifges(7,4) * t124 + Ifges(7,2) * t123) + t124 * (Ifges(7,1) * t124 + t123 * Ifges(7,4)) + m(10) * t108 ^ 2 + 0.2e1 * t108 * (-mrSges(10,1) * t87 + mrSges(10,2) * t86) + t87 * (Ifges(10,4) * t86 + Ifges(10,2) * t87) + t86 * (Ifges(10,1) * t86 + Ifges(10,4) * t87) + m(11) * t78 ^ 2 + 0.2e1 * t78 * (-mrSges(11,1) * t33 + mrSges(11,2) * t32) + t32 * (Ifges(11,1) * t32 + Ifges(11,4) * t33) + t33 * (Ifges(11,4) * t32 + Ifges(11,2) * t33) + Ifges(2,3); t103 * Ifges(8,6) * t57 + t41 * t340 + t8 * t26 + t331 * t380 + Ifges(3,5) * t197 - Ifges(3,6) * t193 + t234 * t9 + (t30 * t33 - t31 * t32) * mrSges(11,3) + (-t40 * mrSges(5,3) + t2 * t358 + t3 * t356) * t39 + t228 * t19 + (t119 * t103 * mrSges(8,3) + (-t170 * t192 - t171 * t196) * mrSges(4,3)) * pkin(1) + t224 * t12 + t58 * (Ifges(7,5) * t124 + Ifges(7,6) * t123) + (-mrSges(8,3) * t351 + Ifges(8,5) * t57) * t102 + t287 + t387 + t396; 0.2e1 * t275 + (t246 * t9 ^ 2 + 0.2e1 * t8 ^ 2) * t381 + 0.2e1 * t274 + (t40 ^ 2 + t41 ^ 2) * m(5) + 0.2e1 * t8 * t1 + t2 * t322 + t3 * t323 + t58 ^ 2 * Ifges(7,3) + Ifges(3,3) + Ifges(4,3) + Ifges(10,3) + t9 * t339 * t246 + (0.2e1 * mrSges(8,1) * t351 - 0.2e1 * mrSges(8,2) * t350 + Ifges(8,3) * t57) * t57 + (-0.2e1 * t326 + 0.2e1 * t328 + t395) * t19 + (Ifges(5,3) * t12 - 0.2e1 * t333 + 0.2e1 * t334) * t12 + t389 + t394 + m(11) * (t30 ^ 2 + t31 ^ 2) + (m(4) * (t192 ^ 2 + t196 ^ 2) + m(8) * (t117 + t222)) * t218; t15 * t26 + t330 * t380 + t340 * t370 + t70 * t396 + t234 * t14 + (-mrSges(5,3) * t369 + t356 * t6 + t358 * t5) * t39 + t228 * t20 + t387 * t74 + t224 * t18 + t287; t275 + t274 + (t40 * t52 + t41 * t51) * pkin(5) * m(5) + (t384 * t8 + t9 * t402) * t381 + t5 * t322 / 0.2e1 + t6 * t323 / 0.2e1 + t2 * t320 / 0.2e1 + t3 * t321 / 0.2e1 + t15 * t1 + t8 * t4 + Ifges(4,3) + (-t333 + t334) * t18 + (Ifges(10,3) + t390) * t70 + (-t326 + t328 + t395) * t20 + (t277 - t278 + t397) * t12 + (t390 + t394) * t74 + t391 * (t14 * t339 + t9 * t338); t6 * t321 + t5 * t320 + t20 ^ 2 * Ifges(11,3) + (t14 ^ 2 * t246 + 0.2e1 * t15 ^ 2) * t381 + (t223 + t46) * t212 * m(5) + t4 * t384 + Ifges(4,3) + t394 * t74 ^ 2 + t338 * t402 + (Ifges(10,3) * t70 + t389 * t74) * t70 + (0.2e1 * t277 - 0.2e1 * t278 + t397) * t18; t25 * t327 + (-mrSges(6,2) * t25 - Ifges(6,6) * t39) * t191 + t329; -t237 * t9 + t331; -t14 * t237 + t330; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
