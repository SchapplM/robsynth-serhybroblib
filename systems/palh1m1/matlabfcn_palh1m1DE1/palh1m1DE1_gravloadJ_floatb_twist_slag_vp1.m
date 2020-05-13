% Calculate Gravitation load on the joints for
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1DE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:52:20
% EndTime: 2020-04-13 15:00:06
% DurationCPUTime: 122.02s
% Computational Cost: add. (3266664->453), mult. (4931557->781), div. (217896->33), fcn. (3127181->46), ass. (0->344)
t211 = cos(pkin(21));
t416 = pkin(5) * pkin(4);
t232 = pkin(4) ^ 2;
t231 = pkin(5) ^ 2;
t229 = pkin(7) ^ 2;
t237 = pkin(1) ^ 2;
t215 = sin(qJ(2));
t220 = cos(qJ(2));
t222 = cos(pkin(19));
t398 = sin(pkin(19));
t192 = t215 * t222 - t220 * t398;
t392 = pkin(7) * t192;
t419 = -2 * pkin(1);
t345 = t392 * t419 + t237;
t175 = t229 + t345;
t342 = pkin(3) ^ 2 - pkin(8) ^ 2;
t161 = t175 + t342;
t182 = pkin(1) - t392;
t413 = pkin(7) + pkin(8);
t414 = pkin(7) - pkin(8);
t155 = (pkin(3) + t413) * (-pkin(3) + t414) + t345;
t156 = (-pkin(3) + t413) * (pkin(3) + t414) + t345;
t239 = sqrt(-t156 * t155);
t194 = t215 * t398 + t220 * t222;
t391 = pkin(7) * t194;
t141 = t161 * t391 + t182 * t239;
t219 = cos(qJ(3));
t353 = t219 * t141;
t361 = t194 * t239;
t140 = -pkin(7) * t361 + t161 * t182;
t214 = sin(qJ(3));
t360 = t214 * t140;
t169 = 0.1e1 / t175;
t234 = 0.1e1 / pkin(3);
t364 = t169 * t234;
t131 = (t360 / 0.2e1 + t353 / 0.2e1) * t364;
t354 = t219 * t140;
t359 = t214 * t141;
t132 = (-t354 / 0.2e1 + t359 / 0.2e1) * t364;
t204 = pkin(23) + pkin(22);
t201 = sin(t204);
t202 = cos(t204);
t101 = t131 * t202 - t132 * t201;
t272 = pkin(5) * t101;
t349 = 0.2e1 * pkin(4) * t272 + t231;
t98 = t232 + t349;
t340 = 0.1e1 / t98 ^ 2 * t416;
t316 = t211 * t340;
t339 = pkin(1) * t391;
t369 = 0.2e1 / t239 * (t155 + t156) * t339;
t417 = -0.2e1 * t194 ^ 2;
t122 = t182 * t369 / 0.2e1 + t229 * pkin(1) * t417 + (-t161 * t192 - t361) * pkin(7);
t312 = 0.1e1 / t175 ^ 2 * t339;
t402 = -t219 / 0.2e1;
t322 = -t369 / 0.2e1;
t362 = t192 * t239;
t116 = (t362 + (t182 * t419 - t161 + t322) * t194) * pkin(7);
t407 = -t116 / 0.2e1;
t85 = ((t122 * t402 + t214 * t407) * t169 + (-t353 - t360) * t312) * t234;
t406 = t122 / 0.2e1;
t86 = ((t116 * t402 + t214 * t406) * t169 + (-t354 + t359) * t312) * t234;
t74 = t201 * t86 + t202 * t85;
t305 = t74 * t316;
t207 = sin(pkin(21));
t317 = t207 * t340;
t306 = t74 * t317;
t412 = -pkin(9) - pkin(11);
t91 = (pkin(4) - t412) * (pkin(4) + t412) + t349;
t411 = pkin(11) - pkin(9);
t92 = (pkin(4) - t411) * (pkin(4) + t411) + t349;
t240 = sqrt(-t92 * t91);
t310 = t131 * t201 + t202 * t132;
t78 = 0.1e1 / t240;
t428 = -t78 / 0.2e1;
t329 = t310 * t428;
t429 = 0.2e1 * t416;
t279 = (t91 + t92) * t429;
t62 = t74 * t279;
t75 = -t201 * t85 + t202 * t86;
t273 = -t75 * t240 + t329 * t62;
t418 = -0.2e1 * pkin(5);
t341 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t93 = t98 + t341;
t99 = pkin(4) * t101 + pkin(5);
t324 = t418 * t99 - t93;
t31 = (t324 * t74 + t273) * pkin(4);
t315 = t310 * t232 * t418;
t332 = t78 * t99 / 0.2e1;
t32 = t62 * t332 + t74 * t315 + (-t240 * t74 + t75 * t93) * pkin(4);
t96 = 0.1e1 / t98;
t379 = t211 * t96;
t326 = t379 / 0.2e1;
t327 = -t379 / 0.2e1;
t403 = t207 / 0.2e1;
t328 = t96 * t403;
t226 = 0.1e1 / pkin(11);
t376 = t226 * t96;
t375 = t310 * t240;
t71 = -pkin(4) * t375 + t93 * t99;
t396 = pkin(4) * t310;
t72 = t240 * t99 + t396 * t93;
t60 = (-t71 * t211 / 0.2e1 + t72 * t403) * t376;
t58 = 0.1e1 / t60 ^ 2;
t59 = (t71 * t403 + t72 * t211 / 0.2e1) * t376;
t377 = t226 / (t58 * t59 ^ 2 + 0.1e1);
t387 = t58 * t59;
t57 = 0.1e1 / t60;
t427 = ((t305 * t72 + t306 * t71 + t31 * t328 + t32 * t326) * t57 - (-t305 * t71 + t306 * t72 + t31 * t327 + t32 * t328) * t387) * t377 + 0.1e1;
t303 = t310 * t316;
t304 = t310 * t317;
t73 = t310 * t279;
t271 = -t101 * t240 + t329 * t73;
t54 = (t310 * t324 + t271) * pkin(4);
t55 = t73 * t332 + t310 * t315 + (t101 * t93 - t375) * pkin(4);
t18 = ((t303 * t72 + t304 * t71 + t326 * t55 + t328 * t54) * t57 - (-t303 * t71 + t304 * t72 + t327 * t54 + t328 * t55) * t387) * t377;
t409 = t18 + 0.1e1;
t213 = sin(qJ(4));
t218 = cos(qJ(4));
t221 = cos(qJ(1));
t216 = sin(qJ(1));
t352 = t219 * t220;
t358 = t215 * t216;
t176 = t214 * t358 - t216 * t352;
t47 = atan2(t59, t60);
t45 = cos(t47);
t33 = t176 * t45;
t280 = t220 * t214 + t215 * t219;
t177 = t280 * t216;
t44 = sin(t47);
t384 = t177 * t44;
t27 = -t384 - t33;
t431 = t221 * t213 - t218 * t27;
t430 = t213 * t27 + t218 * t221;
t399 = pkin(12) + rSges(6,3);
t386 = t96 * t98;
t208 = sin(pkin(20));
t212 = cos(pkin(20));
t188 = t208 * t219 + t212 * t214;
t394 = pkin(6) * t188;
t338 = pkin(1) * t394;
t181 = 0.2e1 * t338;
t235 = pkin(2) ^ 2;
t230 = pkin(6) ^ 2;
t343 = t230 + t237;
t325 = -pkin(13) ^ 2 + t343;
t159 = t181 + t235 + t325;
t180 = -pkin(1) - t394;
t189 = t208 * t214 - t212 * t219;
t346 = t181 + t230;
t415 = -pkin(2) - pkin(13);
t153 = (pkin(1) - t415) * (pkin(1) + t415) + t346;
t410 = pkin(13) - pkin(2);
t154 = (pkin(1) - t410) * (pkin(1) + t410) + t346;
t368 = t154 * t153;
t238 = sqrt(-t368);
t363 = t189 * t238;
t337 = pkin(6) * t363;
t137 = -t159 * t180 - t337;
t393 = pkin(6) * t189;
t138 = t159 * t393 - t180 * t238;
t164 = t181 + t343;
t162 = 0.1e1 / t164;
t236 = 0.1e1 / pkin(2);
t400 = t236 / 0.2e1;
t321 = t162 * t400;
t126 = atan2(t138 * t321, t137 * t321);
t123 = sin(t126);
t124 = cos(t126);
t281 = t123 * t220 + t124 * t215;
t136 = 0.1e1 / t137 ^ 2;
t163 = 0.1e1 / t164 ^ 2;
t370 = 0.2e1 / t238 * (t153 + t154) * pkin(1) * t393;
t323 = -t370 / 0.2e1;
t395 = pkin(6) * t138;
t397 = pkin(1) * t189;
t405 = t162 / 0.2e1;
t76 = 0.2e1 * (((t180 * t323 + (t159 * t188 - t363) * pkin(6)) * t405 + (-t162 * t189 * t230 + t163 * t395) * t397) / t137 - ((-t188 * t238 + (t323 - t159) * t189) * t405 + (t137 * t163 + t162 * t180) * t397) * t136 * t395) * pkin(2) / (t136 * t138 ^ 2 + 0.1e1) * t164 * t236;
t67 = t281 * t76;
t35 = t177 * t45;
t25 = t176 * t44 - t35;
t179 = t280 * t221;
t193 = -t214 * t215 + t352;
t178 = t193 * t221;
t38 = t178 * t45;
t29 = -t179 * t44 + t38;
t426 = t193 * t44 + t280 * t45;
t267 = -t272 - pkin(4);
t283 = t98 - t341;
t259 = -pkin(5) * t375 - t267 * t283;
t258 = 0.1e1 / t259 ^ 2;
t274 = pkin(5) * t283;
t260 = -t240 * t267 + t274 * t310;
t257 = 0.1e1 / (t258 * t260 ^ 2 + 0.1e1);
t251 = t257 * t258 * t260;
t252 = t257 / t259;
t261 = t267 * t429 - t274;
t265 = t267 * t428;
t277 = -pkin(5) * t240 - 0.2e1 * t231 * t396;
t425 = (t101 * t274 + t265 * t73 + t277 * t310) * t252 - (pkin(5) * t271 + t261 * t310) * t251;
t424 = (t265 * t62 + t274 * t75 + t277 * t74) * t252 - (pkin(5) * t273 + t261 * t74) * t251;
t206 = sin(pkin(22));
t210 = cos(pkin(22));
t334 = t96 / pkin(9) / 0.2e1;
t256 = atan2(t260 * t334, t259 * t334);
t255 = sin(t256);
t61 = cos(t256);
t423 = -t206 * t61 + t210 * t255;
t43 = -t206 * t255 - t210 * t61;
t422 = g(1) * t221 + g(2) * t216;
t421 = t386 * t43;
t420 = t423 * t386;
t209 = cos(pkin(23));
t205 = sin(pkin(23));
t404 = t205 / 0.2e1;
t128 = (-t209 * t140 / 0.2e1 + t141 * t404) * t364;
t127 = 0.1e1 / t128 ^ 2;
t129 = (t209 * t141 / 0.2e1 + t140 * t404) * t364;
t296 = t140 * t312;
t299 = t141 * t312;
t320 = t169 * t404;
t367 = t169 * t209;
t69 = ((t116 * t320 + t205 * t296 + t209 * t299 + t367 * t406) / t128 - (t122 * t320 + t205 * t299 - t209 * t296 + t367 * t407) * t129 * t127) / (t127 * t129 ^ 2 + 0.1e1) * t234;
t408 = t69 + 0.1e1;
t223 = cos(pkin(18));
t401 = t223 / 0.2e1;
t388 = t215 * pkin(1);
t286 = -t178 * t44 - t179 * t45;
t284 = t193 * t45 - t280 * t44;
t378 = t216 * pkin(16);
t111 = atan2(t129, t128);
t106 = sin(t111);
t374 = t106 * t215;
t373 = t106 * t220;
t107 = cos(t111);
t372 = t107 * t215;
t371 = t107 * t220;
t217 = sin(pkin(18));
t366 = t169 * t217;
t228 = 0.1e1 / pkin(8);
t365 = t169 * t228;
t357 = t215 * t221;
t356 = t216 * t220;
t351 = t220 * t124;
t350 = t220 * t221;
t348 = -t177 * rSges(4,1) + t176 * rSges(4,2);
t347 = -t179 * rSges(4,1) - t178 * rSges(4,2);
t344 = t193 * rSges(4,1) - rSges(4,2) * t280;
t336 = pkin(1) * t356;
t335 = pkin(1) * t350;
t333 = t409 * t44;
t331 = t408 * t220;
t330 = t408 * t215;
t319 = t169 * t401;
t318 = 0.1e1 / pkin(13) * t400;
t187 = t193 * pkin(5);
t314 = t187 - t388;
t311 = pkin(1) * t358 - t378;
t308 = t107 * t331;
t307 = t106 * t330;
t172 = t177 * pkin(5);
t302 = -t172 - t336;
t174 = t179 * pkin(5);
t301 = -t174 - t335;
t203 = t221 * pkin(16);
t300 = -pkin(1) * t357 + t203;
t298 = t217 * t312;
t297 = t223 * t312;
t295 = t176 * pkin(5) + t311;
t293 = -rSges(3,1) * t215 - rSges(3,2) * t220;
t160 = t175 - t342;
t183 = pkin(1) * t192 - pkin(7);
t139 = -pkin(1) * t361 - t160 * t183;
t142 = pkin(1) * t194 * t160 - t183 * t239;
t133 = (t139 * t401 - t217 * t142 / 0.2e1) * t365;
t134 = (t142 * t401 + t139 * t217 / 0.2e1) * t365;
t113 = atan2(t134, t133);
t109 = sin(t113);
t110 = cos(t113);
t292 = rSges(7,1) * t110 - rSges(7,2) * t109;
t158 = t235 - t325 - 0.2e1 * t338;
t145 = atan2(t238 * t318, t158 * t318);
t143 = sin(t145);
t144 = cos(t145);
t87 = t123 * t358 - t216 * t351;
t88 = t281 * t216;
t290 = t143 * t88 + t144 * t87;
t289 = t143 * t87 - t144 * t88;
t94 = t123 * t215 - t351;
t89 = t94 * t221;
t90 = -t123 * t350 - t124 * t357;
t288 = t143 * t90 - t144 * t89;
t287 = t143 * t89 + t144 * t90;
t282 = t372 + t373;
t278 = t178 * pkin(5) + t300;
t276 = rSges(6,1) * t218 - rSges(6,2) * t213 + pkin(10);
t275 = -pkin(15) + t292;
t68 = t94 * t76;
t270 = g(1) * t276;
t269 = g(2) * t276;
t268 = g(3) * t276;
t266 = t107 * t330 + t373 * t408;
t264 = -rSges(10,1) * t290 + rSges(10,2) * t289;
t263 = rSges(10,1) * t288 + rSges(10,2) * t287;
t262 = (-t143 * t94 + t144 * t281) * rSges(10,1) + (-t143 * t281 - t144 * t94) * rSges(10,2);
t157 = 0.1e1 / t158 ^ 2;
t130 = 0.1e1 / t133 ^ 2;
t121 = t183 * t322 + t237 * pkin(7) * t417 + (-t160 * t192 - t361) * pkin(1);
t115 = (t362 + (0.2e1 * t183 * pkin(7) - t160 + t322) * t194) * pkin(1);
t105 = t106 * t357;
t104 = t107 * t356;
t84 = t371 - t374;
t82 = t282 * t221;
t81 = -t107 * t350 + t105;
t80 = t282 * t216;
t79 = t106 * t358 - t104;
t66 = t221 * t67;
t65 = t76 * t89;
t64 = t76 * t88;
t63 = t216 * t68;
t53 = -t106 * t331 - t372 * t408;
t52 = -t308 + t307;
t51 = t105 + (t374 * t69 - t308) * t221;
t50 = t266 * t221;
t49 = -t104 + (-t371 * t69 + t307) * t216;
t48 = t266 * t216;
t24 = t213 * t216 + t218 * t29;
t23 = -t213 * t29 + t216 * t218;
t17 = t425 * t420;
t16 = t425 * t421;
t14 = t424 * t420;
t13 = t424 * t421;
t12 = t409 * t284;
t11 = t409 * t426;
t10 = t409 * t286;
t9 = -t179 * t333 + t18 * t38 + t38;
t8 = t176 * t333 - t18 * t35 - t35;
t7 = -t18 * t33 - t384 * t409 - t33;
t6 = t427 * t284;
t5 = t427 * t426;
t4 = t427 * t286;
t3 = t427 * t29;
t2 = t427 * t25;
t1 = t427 * t27;
t15 = [-m(2) * (g(1) * (-rSges(2,1) * t216 - rSges(2,2) * t221) + g(2) * (rSges(2,1) * t221 - rSges(2,2) * t216)) - m(3) * (g(2) * t203 + (g(1) * rSges(3,3) + g(2) * t293) * t221 + (g(1) * (-pkin(16) - t293) + g(2) * rSges(3,3)) * t216) - m(4) * (g(1) * (rSges(4,1) * t176 + rSges(4,2) * t177 + rSges(4,3) * t221 + t311) + g(2) * (rSges(4,1) * t178 - rSges(4,2) * t179 + rSges(4,3) * t216 + t300)) - m(5) * (g(1) * (-rSges(5,1) * t27 - rSges(5,2) * t25 + rSges(5,3) * t221 + t295) + g(2) * (rSges(5,1) * t29 + rSges(5,2) * t286 + rSges(5,3) * t216 + t278)) - m(6) * (g(1) * (-t27 * pkin(10) + rSges(6,1) * t431 + rSges(6,2) * t430 + t399 * t25 + t295) + g(2) * (pkin(10) * t29 + rSges(6,1) * t24 + rSges(6,2) * t23 - t286 * t399 + t278)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t275) * t221 + (g(2) * rSges(7,3) - g(1) * t275) * t216) - m(8) * (g(1) * (rSges(8,1) * t80 - rSges(8,2) * t79 + rSges(8,3) * t221 + t311) + g(2) * (-rSges(8,1) * t82 + rSges(8,2) * t81 + rSges(8,3) * t216 + t300)) - m(9) * (g(1) * (rSges(9,1) * t88 - rSges(9,2) * t87 + rSges(9,3) * t221 - t378) + g(2) * (rSges(9,1) * t90 + rSges(9,2) * t89 + rSges(9,3) * t216 + t203)) - m(10) * (g(1) * (t88 * pkin(2) + rSges(10,1) * t289 + rSges(10,2) * t290 + t221 * rSges(10,3) - t378) + g(2) * (t90 * pkin(2) - rSges(10,1) * t287 + rSges(10,2) * t288 + t216 * rSges(10,3) + t203)) - m(11) * (g(1) * ((t423 * t79 + t43 * t80) * rSges(11,1) + (t423 * t80 - t43 * t79) * rSges(11,2) + t221 * rSges(11,3) + t311) + g(2) * ((-t423 * t81 - t43 * t82) * rSges(11,1) + (-t423 * t82 + t43 * t81) * rSges(11,2) + t216 * rSges(11,3) + t300) + (g(1) * (t79 * t206 + t80 * t210) + g(2) * (-t81 * t206 - t82 * t210)) * pkin(4)), -m(3) * (g(3) * t293 + t422 * (-rSges(3,1) * t220 + rSges(3,2) * t215)) - m(4) * (g(1) * (-t335 + t347) + g(2) * (-t336 + t348) + g(3) * (t344 - t388)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t301) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 + t302) + g(3) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t314)) - m(6) * (g(1) * (t3 * t399 + t301) + g(2) * (t1 * t399 + t302) + g(3) * (t399 * t5 + t314) + t6 * t268 + t4 * t270 + t2 * t269) - m(7) * (g(3) * t292 + t422 * (-rSges(7,1) * t109 - rSges(7,2) * t110)) * ((t121 * t319 + t142 * t297 + t115 * t366 / 0.2e1 + t139 * t298) / t133 - (t115 * t319 + t139 * t297 - t121 * t366 / 0.2e1 - t142 * t298) * t134 * t130) / (t130 * t134 ^ 2 + 0.1e1) * t228 - m(8) * (g(1) * (rSges(8,1) * t51 + rSges(8,2) * t50 - t335) + g(2) * (rSges(8,1) * t49 + rSges(8,2) * t48 - t336) + g(3) * (rSges(8,1) * t53 + rSges(8,2) * t52 - t388)) - m(9) * (g(1) * (rSges(9,1) * t89 - rSges(9,2) * t90) + g(2) * (rSges(9,1) * t87 + rSges(9,2) * t88) + g(3) * (-rSges(9,1) * t281 + rSges(9,2) * t94)) - m(10) * (g(1) * (t89 * pkin(2) + t263) + g(2) * (t87 * pkin(2) + t264) + g(3) * (-pkin(2) * t281 + t262)) - m(11) * (g(1) * (-t335 + (t13 * t81 - t14 * t82 - t423 * t50 + t43 * t51) * rSges(11,1) + (t13 * t82 + t14 * t81 + t423 * t51 + t43 * t50) * rSges(11,2)) + g(2) * (-t336 + (t13 * t79 - t14 * t80 - t423 * t48 + t43 * t49) * rSges(11,1) + (t13 * t80 + t14 * t79 + t423 * t49 + t43 * t48) * rSges(11,2)) + g(3) * (-t388 + (-t13 * t282 + t14 * t84 - t423 * t52 + t43 * t53) * rSges(11,1) + (-t13 * t84 - t14 * t282 + t423 * t53 + t43 * t52) * rSges(11,2)) + (g(1) * (-t206 * t50 + t210 * t51) + g(2) * (-t206 * t48 + t210 * t49) + g(3) * (-t206 * t52 + t210 * t53)) * pkin(4)), -m(4) * (g(1) * t347 + g(2) * t348 + g(3) * t344) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 - t174) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t7 - t172) + g(3) * (rSges(5,1) * t12 - rSges(5,2) * t11 + t187)) - m(6) * (g(1) * (t399 * t9 - t174) + g(2) * (t399 * t7 - t172) + g(3) * (t11 * t399 + t187) + t8 * t269 + t12 * t268 + t10 * t270) - m(9) * (g(1) * (rSges(9,1) * t65 + rSges(9,2) * t66) + g(2) * (rSges(9,1) * t63 + rSges(9,2) * t64) + g(3) * (-rSges(9,1) * t67 + rSges(9,2) * t68)) - m(10) * (g(1) * (t65 * pkin(2) + (-t143 * t66 - t144 * t65) * rSges(10,1) + (t143 * t65 - t144 * t66) * rSges(10,2)) + g(2) * (t63 * pkin(2) + (-t143 * t64 - t144 * t63) * rSges(10,1) + (t143 * t63 - t144 * t64) * rSges(10,2)) + g(3) * (-t67 * pkin(2) + (-t143 * t68 + t144 * t67) * rSges(10,1) + (-t143 * t67 - t144 * t68) * rSges(10,2)) + (g(1) * t263 + g(2) * t264 + g(3) * t262) * (0.1e1 / t158 * t370 / 0.2e1 + t157 * t337 * t419) / (-t157 * t368 + 0.1e1)) - m(11) * (g(1) * ((t16 * t81 - t17 * t82) * rSges(11,1) + (t16 * t82 + t17 * t81) * rSges(11,2)) + g(2) * ((t16 * t79 - t17 * t80) * rSges(11,1) + (t16 * t80 + t17 * t79) * rSges(11,2)) + g(3) * ((-t16 * t282 + t17 * t84) * rSges(11,1) + (-t16 * t84 - t17 * t282) * rSges(11,2))), -m(6) * (g(1) * (rSges(6,1) * t23 - rSges(6,2) * t24) + g(2) * (-rSges(6,1) * t430 + rSges(6,2) * t431) + g(3) * (-rSges(6,1) * t213 - t218 * rSges(6,2)) * t426)];
taug = t15(:);
