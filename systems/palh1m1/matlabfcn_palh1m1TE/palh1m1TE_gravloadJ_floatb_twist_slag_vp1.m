% Calculate Gravitation load on the joints for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1TE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:13:56
% EndTime: 2020-04-12 20:17:52
% DurationCPUTime: 54.17s
% Computational Cost: add. (948777->475), mult. (1438729->785), div. (59868->15), fcn. (915119->28), ass. (0->306)
t359 = sin(qJ(2));
t360 = sin(pkin(19));
t362 = cos(qJ(2));
t363 = cos(pkin(19));
t270 = t359 * t363 - t362 * t360;
t263 = pkin(7) * t270;
t253 = (-0.2e1 * t263 + pkin(1)) * pkin(1);
t384 = pkin(7) ^ 2;
t136 = t253 + t384;
t180 = pkin(3) ^ 2;
t325 = pkin(8) ^ 2 - t180;
t252 = t136 - t325;
t255 = -t263 + pkin(1);
t152 = t359 * t360 + t362 * t363;
t374 = pkin(7) + pkin(8);
t375 = pkin(7) - pkin(8);
t179 = sqrt(-((-pkin(3) + t374) * (pkin(3) + t375) + t253) * ((pkin(3) + t374) * (-pkin(3) + t375) + t253));
t334 = t152 * t179;
t314 = pkin(7) * t334;
t386 = 0.1e1 / pkin(3);
t246 = t386 * (t252 * t255 - t314);
t378 = 0.1e1 / t136;
t245 = t378 * t246;
t242 = t245 / 0.2e1;
t251 = pkin(7) * t252;
t247 = t386 * (t152 * t251 + t179 * t255);
t321 = t378 / 0.2e1;
t244 = t247 * t321;
t358 = sin(qJ(3));
t361 = cos(qJ(3));
t238 = t242 * t358 + t244 * t361;
t243 = -t245 / 0.2e1;
t239 = t243 * t361 + t244 * t358;
t315 = pkin(23) + pkin(22);
t295 = sin(t315);
t296 = cos(t315);
t231 = t238 * t296 - t295 * t239;
t230 = pkin(4) * t231;
t222 = (0.2e1 * t230 + pkin(5)) * pkin(5);
t372 = pkin(11) - pkin(9);
t373 = -pkin(9) - pkin(11);
t178 = sqrt(-((pkin(4) - t372) * (pkin(4) + t372) + t222) * ((pkin(4) - t373) * (pkin(4) + t373) + t222));
t388 = t238 * t295 + t296 * t239;
t409 = t178 * t388;
t413 = pkin(4) * t409;
t412 = pkin(5) * t409;
t182 = pkin(9) ^ 2;
t324 = pkin(11) ^ 2 - t182;
t229 = pkin(5) * t231;
t385 = pkin(5) ^ 2;
t72 = t385 + (0.2e1 * t229 + pkin(4)) * pkin(4);
t219 = t72 + t324;
t224 = t230 + pkin(5);
t382 = 0.1e1 / pkin(11);
t204 = t382 * (t219 * t224 - t413);
t217 = pkin(4) * t219;
t205 = t382 * (t178 * t224 + t217 * t388);
t340 = sin(pkin(21));
t342 = cos(pkin(21));
t411 = t340 * t204 + t342 * t205;
t171 = pkin(13) ^ 2;
t174 = pkin(6) ^ 2;
t175 = pkin(2) ^ 2;
t162 = sin(pkin(20));
t164 = cos(pkin(20));
t147 = t162 * t361 + t164 * t358;
t355 = pkin(6) * t147;
t367 = -t171 / 0.2e1 + t174 / 0.2e1 - t175 / 0.2e1 + (t355 + pkin(1) / 0.2e1) * pkin(1);
t380 = 0.2e1 * pkin(1);
t327 = t355 * t380 + t174;
t376 = -pkin(2) - pkin(13);
t119 = (pkin(1) - t376) * (pkin(1) + t376) + t327;
t371 = pkin(13) - pkin(2);
t120 = (pkin(1) - t371) * (pkin(1) + t371) + t327;
t148 = t162 * t358 - t164 * t361;
t354 = pkin(6) * t148;
t318 = pkin(1) * t354;
t177 = sqrt(-t120 * t119);
t365 = t177 / 0.2e1;
t338 = (t119 + t120) * t318 / t365;
t301 = -t338 / 0.4e1;
t126 = pkin(1) ^ 2 + t327;
t124 = 0.1e1 / t126;
t176 = 0.1e1 / pkin(2);
t298 = 0.1e1 / t126 ^ 2 * t318;
t122 = t126 - t171 + t175;
t141 = -pkin(1) - t355;
t111 = t122 * t354 - t141 * t177;
t303 = t359 * t111;
t335 = t148 * t177;
t110 = -pkin(6) * t335 - t122 * t141;
t306 = t362 * t110;
t311 = -t362 / 0.2e1;
t302 = -t338 / 0.2e1;
t86 = (-t147 * t177 + (t141 * t380 - t122 + t302) * t148) * pkin(6);
t381 = -0.2e1 * pkin(1);
t88 = t141 * t302 + t174 * t148 ^ 2 * t381 + (t122 * t147 - t335) * pkin(6);
t410 = t176 * (t124 * (t359 * t88 / 0.2e1 + t86 * t311) + (t303 - t306) * t298);
t165 = sin(qJ(4));
t168 = cos(qJ(4));
t169 = cos(qJ(1));
t151 = -t359 * t358 + t362 * t361;
t166 = sin(qJ(1));
t137 = t151 * t166;
t150 = -t358 * t362 - t359 * t361;
t138 = t150 * t166;
t379 = 0.1e1 / t72;
t330 = t379 / 0.2e1;
t45 = t411 * t330;
t199 = t342 * t204;
t200 = t340 * t205;
t331 = -t379 / 0.2e1;
t46 = t199 * t331 + t200 * t330;
t28 = -t137 * t46 - t138 * t45;
t408 = t169 * t165 + t168 * t28;
t407 = t165 * t28 - t168 * t169;
t406 = -0.2e1 * t388;
t258 = t270 * t179;
t357 = pkin(1) * t152;
t317 = pkin(7) * t357;
t337 = 0.4e1 / t179 * (t374 * t375 - t180 + t253) * t317;
t300 = -t337 / 0.2e1;
t87 = (t258 + (t300 + t325 - t384) * t152) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t263) * t317;
t405 = -t87 / 0.2e1;
t319 = t152 ^ 2 * t381;
t93 = -t314 + t255 * t337 / 0.2e1 - t270 * t251 + t384 * t319;
t404 = t93 / 0.2e1;
t161 = sin(pkin(22));
t403 = t161 / 0.2e1;
t163 = cos(pkin(22));
t402 = -t163 / 0.2e1;
t401 = -t388 / 0.2e1;
t400 = t340 / 0.2e1;
t399 = -t342 / 0.2e1;
t364 = -pkin(12) - rSges(6,3);
t304 = t359 * t110;
t305 = t362 * t111;
t336 = t124 * t176;
t95 = (t304 / 0.2e1 + t305 / 0.2e1) * t336;
t377 = pkin(4) * pkin(5);
t320 = 0.1e1 / t72 ^ 2 * t377;
t398 = t411 * t320;
t397 = (-t199 + t200) * t320;
t389 = 0.4e1 / t178 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t222 - t182) * t377;
t213 = t388 * t389;
t396 = -t178 * t231 + t213 * t401;
t131 = 0.1e1 / t136 ^ 2;
t240 = t246 * t317;
t241 = t247 * t317;
t322 = -t378 / 0.2e1;
t279 = t386 * t361 * t322;
t356 = t386 * t378;
t294 = t358 * t356;
t234 = t294 * t405 + t93 * t279 + (-t240 * t358 - t241 * t361) * t131;
t235 = t87 * t279 + t294 * t404 + (-t240 * t361 + t241 * t358) * t131;
t210 = -t234 * t295 + t235 * t296;
t60 = t234 * t296 + t235 * t295;
t214 = t60 * t389;
t395 = -t178 * t210 + t214 * t401;
t218 = t72 - t324;
t216 = pkin(5) * t218;
t223 = -t229 - pkin(4);
t394 = 0.2e1 * t223 * t377 - t216;
t393 = -0.2e1 * t224 * t377 - t217;
t383 = 0.1e1 / pkin(9);
t206 = t383 * (-t178 * t223 + t216 * t388);
t202 = t206 * t320;
t53 = -t218 * t223 - t412;
t368 = t53 * t383;
t291 = t320 * t368;
t392 = -t161 * t202 - t163 * t291;
t391 = t161 * t291 - t163 * t202;
t390 = g(1) * t169 + g(2) * t166;
t370 = t383 * t379;
t369 = t382 * t379;
t366 = -t177 / 0.2e1;
t212 = t214 / 0.2e1;
t225 = pkin(5) * pkin(4) ^ 2 * t406;
t312 = t382 * t330;
t344 = t178 * t60;
t184 = (-pkin(4) * t344 + t210 * t217 + t212 * t224 + t225 * t60) * t312;
t186 = (t395 * pkin(4) + t393 * t60) * t369;
t15 = t342 * t184 + t186 * t400 + t398 * t60;
t351 = t15 + t46;
t16 = t340 * t184 + t186 * t399 + t397 * t60;
t350 = t16 - t45;
t211 = t213 / 0.2e1;
t188 = (t211 * t224 + t217 * t231 + t225 * t388 - t413) * t312;
t190 = (t396 * pkin(4) + t388 * t393) * t369;
t19 = t342 * t188 + t190 * t400 + t388 * t398;
t349 = t19 + t46;
t20 = t340 * t188 + t190 * t399 + t388 * t397;
t348 = t20 - t45;
t347 = -t137 * t45 + t138 * t46;
t139 = t151 * t169;
t140 = t150 * t169;
t30 = -t139 * t45 + t140 * t46;
t346 = t150 * t45 + t151 * t46;
t345 = t166 * pkin(16);
t341 = cos(pkin(23));
t339 = sin(pkin(23));
t332 = 0.1e1 / pkin(13) * t176;
t329 = t138 * rSges(4,1) - t137 * rSges(4,2);
t328 = t140 * rSges(4,1) - t139 * rSges(4,2);
t326 = t151 * rSges(4,1) + t150 * rSges(4,2);
t313 = t383 * t331;
t310 = t362 * pkin(1);
t309 = t359 * pkin(1);
t308 = rSges(10,1) * t332;
t307 = rSges(10,2) * t332;
t167 = sin(pkin(18));
t299 = t167 * t322;
t297 = t166 * t309 - t345;
t293 = t166 * t310;
t292 = t169 * t310;
t290 = t341 * t356;
t289 = t332 * t367;
t288 = t332 * t365;
t173 = 0.1e1 / pkin(8);
t285 = t131 * t173 * t317;
t146 = t151 * pkin(5);
t284 = t146 - t309;
t283 = -t137 * pkin(5) + t297;
t278 = t386 * t339 * t321;
t277 = t167 * t285;
t170 = cos(pkin(18));
t276 = t170 * t285;
t133 = t138 * pkin(5);
t275 = t133 - t293;
t135 = t140 * pkin(5);
t274 = t135 - t292;
t160 = t169 * pkin(16);
t273 = -t169 * t309 + t160;
t123 = t136 + t325;
t143 = pkin(1) * t270 - pkin(7);
t112 = -pkin(1) * t334 - t143 * t123;
t113 = t123 * t357 - t143 * t179;
t272 = -pkin(15) + (rSges(7,1) * (t112 * t170 * t321 + t113 * t299) + rSges(7,2) * (t113 * t170 * t322 + t112 * t299)) * t173;
t271 = rSges(6,1) * t168 - rSges(6,2) * t165 + pkin(10);
t89 = t166 * t95;
t94 = (t303 / 0.2e1 - t306 / 0.2e1) * t336;
t90 = t166 * t94;
t269 = t366 * t89 + t367 * t90;
t91 = t169 * t95;
t92 = t169 * t94;
t268 = t92 * t365 + t367 * t91;
t267 = t139 * pkin(5) + t273;
t266 = g(1) * t271;
t265 = g(2) * t271;
t264 = g(3) * t271;
t262 = -rSges(3,1) * t359 - rSges(3,2) * t362;
t96 = t243 * t341 + t244 * t339;
t97 = t242 * t339 + t244 * t341;
t260 = t359 * t96 + t362 * t97;
t259 = t359 * t97 - t362 * t96;
t67 = t290 * t405 + t93 * t278 + (-t240 * t341 + t241 * t339) * t131;
t68 = t290 * t404 + t87 * t278 + (t240 * t339 + t241 * t341) * t131;
t58 = -t359 * t68 + t362 * t67 - t260;
t59 = -t359 * t67 - t362 * t68 + t259;
t250 = pkin(1) * t173 * (t258 + (0.2e1 * t143 * pkin(7) - t123 + t300) * t152) * t321;
t249 = t173 * t378 * (t143 * t300 + (pkin(7) * t319 - t123 * t270 - t334) * pkin(1));
t66 = ((t88 * t311 - t359 * t86 / 0.2e1) * t124 + (-t304 - t305) * t298) * t176;
t226 = pkin(4) * t385 * t406;
t203 = t206 * t331;
t191 = (t396 * pkin(5) + t388 * t394) * t370;
t189 = (-t211 * t223 + t216 * t231 + t226 * t388 - t412) * t313;
t187 = (t395 * pkin(5) + t394 * t60) * t370;
t185 = (-pkin(5) * t344 + t210 * t216 - t212 * t223 + t226 * t60) * t313;
t83 = t92 * t289;
t82 = t89 * t289;
t76 = t259 * t169;
t75 = t260 * t169;
t74 = t259 * t166;
t73 = t260 * t166;
t70 = t170 * t249 / 0.2e1 + t113 * t276 + t167 * t250 + t112 * t277;
t69 = t170 * t250 + t112 * t276 - t167 * t249 / 0.2e1 - t113 * t277;
t64 = t169 * t410;
t63 = t169 * t66;
t62 = t166 * t410;
t61 = t166 * t66;
t57 = t58 * t169;
t56 = t59 * t169;
t55 = t58 * t166;
t54 = t59 * t166;
t44 = t163 * t313 * t53 + t161 * t203;
t43 = t161 * t330 * t368 + t163 * t203;
t31 = t139 * t46 + t140 * t45;
t26 = t165 * t166 + t168 * t31;
t25 = -t165 * t31 + t166 * t168;
t18 = t161 * t189 + t191 * t402 + t388 * t392;
t17 = t163 * t189 + t191 * t403 + t388 * t391;
t14 = t161 * t185 + t187 * t402 + t392 * t60;
t13 = t163 * t185 + t187 * t403 + t391 * t60;
t12 = -t150 * t20 + t151 * t19 + t346;
t11 = t150 * t349 + t151 * t348;
t10 = t139 * t20 + t140 * t19 + t30;
t9 = -t139 * t349 + t140 * t348;
t8 = t137 * t20 + t138 * t19 + t347;
t7 = -t137 * t349 + t138 * t348;
t6 = t15 * t151 - t150 * t16 + t346;
t5 = t150 * t351 + t151 * t350;
t4 = t139 * t16 + t140 * t15 + t30;
t3 = -t139 * t351 + t140 * t350;
t2 = t137 * t16 + t138 * t15 + t347;
t1 = -t137 * t351 + t138 * t350;
t21 = [-m(2) * (g(1) * (-rSges(2,1) * t166 - rSges(2,2) * t169) + g(2) * (rSges(2,1) * t169 - rSges(2,2) * t166)) - m(3) * (g(2) * t160 + (g(1) * rSges(3,3) + g(2) * t262) * t169 + (g(1) * (-pkin(16) - t262) + g(2) * rSges(3,3)) * t166) - m(4) * (g(1) * (-rSges(4,1) * t137 - rSges(4,2) * t138 + rSges(4,3) * t169 + t297) + g(2) * (t139 * rSges(4,1) + t140 * rSges(4,2) + t166 * rSges(4,3) + t273)) - m(5) * (g(1) * (rSges(5,1) * t28 - rSges(5,2) * t347 + rSges(5,3) * t169 + t283) + g(2) * (t31 * rSges(5,1) + t30 * rSges(5,2) + t166 * rSges(5,3) + t267)) - m(6) * (g(1) * (t28 * pkin(10) + rSges(6,1) * t408 - rSges(6,2) * t407 - t364 * t347 + t283) + g(2) * (t31 * pkin(10) + t26 * rSges(6,1) + t25 * rSges(6,2) + t30 * t364 + t267)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t272) * t169 + (g(2) * rSges(7,3) - g(1) * t272) * t166) - m(8) * (g(1) * (rSges(8,1) * t73 - rSges(8,2) * t74 + rSges(8,3) * t169 + t297) + g(2) * (-t75 * rSges(8,1) + t76 * rSges(8,2) + t166 * rSges(8,3) + t273)) - m(9) * (g(1) * (rSges(9,1) * t89 - rSges(9,2) * t90 + rSges(9,3) * t169 - t345) + g(2) * (-rSges(9,1) * t91 + rSges(9,2) * t92 + rSges(9,3) * t166 + t160)) - m(10) * (g(1) * (t89 * pkin(2) + t82 * rSges(10,1) + t169 * rSges(10,3) - t345) + g(2) * (-t91 * pkin(2) + t83 * rSges(10,2) + t166 * rSges(10,3) + t160) + (g(1) * (rSges(10,1) * t365 * t90 - rSges(10,2) * t269) + g(2) * (-rSges(10,2) * t365 * t91 - rSges(10,1) * t268)) * t332) - m(11) * (g(1) * ((-t43 * t74 + t44 * t73) * rSges(11,1) + (-t43 * t73 - t44 * t74) * rSges(11,2) + t169 * rSges(11,3) + t297) + g(2) * ((t43 * t76 - t44 * t75) * rSges(11,1) + (t43 * t75 + t44 * t76) * rSges(11,2) + t166 * rSges(11,3) + t273) + (g(1) * (t161 * t74 + t163 * t73) + g(2) * (-t161 * t76 - t163 * t75)) * pkin(4)), -m(3) * (g(3) * t262 + t390 * (-rSges(3,1) * t362 + rSges(3,2) * t359)) - m(4) * (g(1) * (-t292 + t328) + g(2) * (-t293 + t329) + g(3) * (-t309 + t326)) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t274) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t275) + g(3) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t284)) - m(6) * (g(1) * (t364 * t3 + t274) + g(2) * (t364 * t1 + t275) + g(3) * (t364 * t5 + t284) + t6 * t264 + t4 * t266 + t2 * t265) - m(7) * (g(3) * (rSges(7,1) * t70 + rSges(7,2) * t69) + t390 * (rSges(7,1) * t69 - rSges(7,2) * t70)) - m(8) * (g(1) * (t56 * rSges(8,1) - t57 * rSges(8,2) - t292) + g(2) * (t54 * rSges(8,1) - t55 * rSges(8,2) - t293) + g(3) * (t58 * rSges(8,1) + t59 * rSges(8,2) - t309)) - m(9) * (g(1) * (rSges(9,1) * t92 + rSges(9,2) * t91) + g(2) * (rSges(9,1) * t90 + rSges(9,2) * t89) + g(3) * (-rSges(9,1) * t95 + rSges(9,2) * t94)) - m(10) * (g(1) * (t92 * pkin(2) + (-t288 * t91 + t83) * rSges(10,1) + t268 * t307) + g(2) * (t90 * pkin(2) + t269 * t308 + (t288 * t90 + t82) * rSges(10,2)) + g(3) * (-t95 * pkin(2) + (t366 * t94 - t367 * t95) * t308 + (-t365 * t95 + t367 * t94) * t307)) - m(11) * (g(1) * (-t292 + (t13 * t76 - t14 * t75 - t43 * t57 + t44 * t56) * rSges(11,1) + (t13 * t75 + t14 * t76 - t43 * t56 - t44 * t57) * rSges(11,2)) + g(2) * (-t293 + (t13 * t74 - t14 * t73 - t43 * t55 + t44 * t54) * rSges(11,1) + (t13 * t73 + t14 * t74 - t43 * t54 - t44 * t55) * rSges(11,2)) + g(3) * (-t309 + (-t13 * t260 - t14 * t259 + t43 * t59 + t44 * t58) * rSges(11,1) + (t13 * t259 - t14 * t260 - t43 * t58 + t44 * t59) * rSges(11,2)) + (g(1) * (t161 * t57 + t163 * t56) + g(2) * (t161 * t55 + t163 * t54) + g(3) * (-t161 * t59 + t163 * t58)) * pkin(4)), -m(4) * (g(1) * t328 + g(2) * t329 + g(3) * t326) - m(5) * (g(1) * (rSges(5,1) * t10 + rSges(5,2) * t9 + t135) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t7 + t133) + g(3) * (rSges(5,1) * t12 + rSges(5,2) * t11 + t146)) - m(6) * (g(1) * (t364 * t9 + t135) + g(2) * (t364 * t7 + t133) + g(3) * (t364 * t11 + t146) + t8 * t265 + t12 * t264 + t10 * t266) - m(9) * (g(1) * (rSges(9,1) * t63 + rSges(9,2) * t64) + g(2) * (rSges(9,1) * t61 + rSges(9,2) * t62) + g(3) * (-rSges(9,1) * t410 + rSges(9,2) * t66)) - m(10) * ((g(1) * t63 + g(2) * t61 - g(3) * t410) * pkin(2) + (g(1) * ((t301 * t92 + t318 * t91 + t366 * t64 + t367 * t63) * rSges(10,1) + (t301 * t91 - t92 * t318 + t63 * t365 + t64 * t367) * rSges(10,2)) + g(2) * ((t301 * t90 + t318 * t89 + t366 * t62 + t367 * t61) * rSges(10,1) + (t301 * t89 - t318 * t90 + t365 * t61 + t367 * t62) * rSges(10,2)) + g(3) * ((-t301 * t95 + t318 * t94 + t366 * t66 - t367 * t410) * rSges(10,1) + (t301 * t94 + t318 * t95 - t365 * t410 + t367 * t66) * rSges(10,2))) * t332) - m(11) * (g(1) * ((t17 * t76 - t18 * t75) * rSges(11,1) + (t17 * t75 + t18 * t76) * rSges(11,2)) + g(2) * ((t17 * t74 - t18 * t73) * rSges(11,1) + (t17 * t73 + t18 * t74) * rSges(11,2)) + g(3) * ((-t17 * t260 - t18 * t259) * rSges(11,1) + (t17 * t259 - t18 * t260) * rSges(11,2))), -m(6) * (g(1) * (rSges(6,1) * t25 - rSges(6,2) * t26) + g(2) * (rSges(6,1) * t407 + rSges(6,2) * t408) + g(3) * (-rSges(6,1) * t165 - rSges(6,2) * t168) * (-t150 * t46 + t151 * t45))];
taug = t21(:);
