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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1DE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:52:20
% EndTime: 2020-04-13 15:00:10
% DurationCPUTime: 121.97s
% Computational Cost: add. (3266663->398), mult. (4931576->656), div. (217896->33), fcn. (3127181->46), ass. (0->336)
t304 = -pkin(12) * m(6) + mrSges(5,2) - mrSges(6,3);
t213 = cos(pkin(21));
t420 = pkin(4) * pkin(5);
t234 = pkin(4) ^ 2;
t233 = pkin(5) ^ 2;
t231 = pkin(7) ^ 2;
t239 = pkin(1) ^ 2;
t217 = sin(qJ(2));
t222 = cos(qJ(2));
t224 = cos(pkin(19));
t400 = sin(pkin(19));
t192 = t217 * t224 - t222 * t400;
t391 = pkin(7) * t192;
t423 = -2 * pkin(1);
t342 = t391 * t423 + t239;
t175 = t231 + t342;
t339 = pkin(3) ^ 2 - pkin(8) ^ 2;
t161 = t175 + t339;
t182 = pkin(1) - t391;
t416 = pkin(7) + pkin(8);
t417 = pkin(7) - pkin(8);
t155 = (pkin(3) + t416) * (-pkin(3) + t417) + t342;
t156 = (-pkin(3) + t416) * (pkin(3) + t417) + t342;
t241 = sqrt(-t156 * t155);
t194 = t217 * t400 + t222 * t224;
t390 = pkin(7) * t194;
t141 = t161 * t390 + t182 * t241;
t221 = cos(qJ(3));
t352 = t221 * t141;
t359 = t194 * t241;
t140 = -pkin(7) * t359 + t161 * t182;
t216 = sin(qJ(3));
t358 = t216 * t140;
t169 = 0.1e1 / t175;
t236 = 0.1e1 / pkin(3);
t362 = t169 * t236;
t131 = (t358 / 0.2e1 + t352 / 0.2e1) * t362;
t353 = t221 * t140;
t357 = t216 * t141;
t132 = (-t353 / 0.2e1 + t357 / 0.2e1) * t362;
t206 = pkin(23) + pkin(22);
t203 = sin(t206);
t204 = cos(t206);
t101 = t131 * t204 - t132 * t203;
t274 = pkin(5) * t101;
t348 = 0.2e1 * pkin(4) * t274 + t233;
t98 = t234 + t348;
t337 = 0.1e1 / t98 ^ 2 * t420;
t312 = t213 * t337;
t335 = pkin(1) * t390;
t367 = 0.2e1 / t241 * (t155 + t156) * t335;
t421 = -0.2e1 * t194 ^ 2;
t122 = t182 * t367 / 0.2e1 + t231 * pkin(1) * t421 + (-t161 * t192 - t359) * pkin(7);
t307 = 0.1e1 / t175 ^ 2 * t335;
t403 = -t221 / 0.2e1;
t318 = -t367 / 0.2e1;
t360 = t192 * t241;
t116 = (t360 + (t182 * t423 - t161 + t318) * t194) * pkin(7);
t408 = -t116 / 0.2e1;
t85 = ((t122 * t403 + t216 * t408) * t169 + (-t352 - t358) * t307) * t236;
t407 = t122 / 0.2e1;
t86 = ((t116 * t403 + t216 * t407) * t169 + (-t353 + t357) * t307) * t236;
t74 = t203 * t86 + t204 * t85;
t300 = t74 * t312;
t209 = sin(pkin(21));
t313 = t209 * t337;
t301 = t74 * t313;
t415 = -pkin(9) - pkin(11);
t91 = (pkin(4) - t415) * (pkin(4) + t415) + t348;
t414 = pkin(11) - pkin(9);
t92 = (pkin(4) - t414) * (pkin(4) + t414) + t348;
t242 = sqrt(-t92 * t91);
t306 = t131 * t203 + t132 * t204;
t78 = 0.1e1 / t242;
t442 = -t78 / 0.2e1;
t325 = t306 * t442;
t443 = 0.2e1 * t420;
t279 = (t91 + t92) * t443;
t62 = t74 * t279;
t75 = -t203 * t85 + t204 * t86;
t275 = -t242 * t75 + t325 * t62;
t422 = -0.2e1 * pkin(5);
t338 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t93 = t98 + t338;
t99 = pkin(4) * t101 + pkin(5);
t320 = t422 * t99 - t93;
t31 = (t320 * t74 + t275) * pkin(4);
t311 = t306 * t234 * t422;
t328 = t78 * t99 / 0.2e1;
t32 = t62 * t328 + t74 * t311 + (-t242 * t74 + t75 * t93) * pkin(4);
t96 = 0.1e1 / t98;
t376 = t213 * t96;
t322 = t376 / 0.2e1;
t323 = -t376 / 0.2e1;
t404 = t209 / 0.2e1;
t324 = t96 * t404;
t228 = 0.1e1 / pkin(11);
t374 = t228 * t96;
t373 = t306 * t242;
t71 = -pkin(4) * t373 + t93 * t99;
t397 = pkin(4) * t306;
t72 = t242 * t99 + t397 * t93;
t60 = (-t71 * t213 / 0.2e1 + t72 * t404) * t374;
t58 = 0.1e1 / t60 ^ 2;
t59 = (t71 * t404 + t72 * t213 / 0.2e1) * t374;
t375 = t228 / (t58 * t59 ^ 2 + 0.1e1);
t387 = t58 * t59;
t57 = 0.1e1 / t60;
t444 = ((t300 * t72 + t301 * t71 + t31 * t324 + t32 * t322) * t57 - (-t300 * t71 + t301 * t72 + t31 * t323 + t32 * t324) * t387) * t375 + 0.1e1;
t446 = t304 * t444;
t298 = t306 * t312;
t299 = t306 * t313;
t73 = t306 * t279;
t273 = -t101 * t242 + t325 * t73;
t54 = (t306 * t320 + t273) * pkin(4);
t55 = t73 * t328 + t306 * t311 + (t101 * t93 - t373) * pkin(4);
t18 = ((t298 * t72 + t299 * t71 + t322 * t55 + t324 * t54) * t57 - (-t298 * t71 + t299 * t72 + t323 * t54 + t324 * t55) * t387) * t375;
t410 = t18 + 0.1e1;
t351 = t221 * t222;
t193 = -t216 * t217 + t351;
t280 = t216 * t222 + t217 * t221;
t47 = atan2(t59, t60);
t44 = sin(t47);
t45 = cos(t47);
t438 = t193 * t44 + t280 * t45;
t445 = t304 * t438;
t215 = sin(qJ(4));
t220 = cos(qJ(4));
t439 = -pkin(10) * m(6) - mrSges(6,1) * t220 + mrSges(6,2) * t215 - mrSges(5,1);
t429 = t439 * t444;
t441 = m(6) + m(5);
t412 = -m(8) - m(4);
t440 = pkin(4) * m(11);
t386 = t96 * t98;
t218 = sin(qJ(1));
t356 = t217 * t218;
t176 = t216 * t356 - t218 * t351;
t177 = t280 * t218;
t35 = t177 * t45;
t287 = t176 * t44 - t35;
t223 = cos(qJ(1));
t179 = t280 * t223;
t178 = t193 * t223;
t38 = t178 * t45;
t29 = -t179 * t44 + t38;
t33 = t176 * t45;
t381 = t177 * t44;
t27 = -t381 - t33;
t437 = -pkin(2) * m(10) - mrSges(9,1);
t271 = -t274 - pkin(4);
t283 = t98 - t338;
t261 = -pkin(5) * t373 - t271 * t283;
t260 = 0.1e1 / t261 ^ 2;
t277 = pkin(5) * t283;
t262 = -t242 * t271 + t277 * t306;
t259 = 0.1e1 / (t260 * t262 ^ 2 + 0.1e1);
t253 = t259 * t260 * t262;
t254 = t259 / t261;
t263 = t271 * t443 - t277;
t268 = t271 * t442;
t278 = -pkin(5) * t242 - 0.2e1 * t233 * t397;
t436 = (t101 * t277 + t268 * t73 + t278 * t306) * t254 - (pkin(5) * t273 + t263 * t306) * t253;
t435 = (t268 * t62 + t277 * t75 + t278 * t74) * t254 - (pkin(5) * t275 + t263 * t74) * t253;
t208 = sin(pkin(22));
t212 = cos(pkin(22));
t330 = t96 / pkin(9) / 0.2e1;
t258 = atan2(t262 * t330, t261 * t330);
t257 = sin(t258);
t61 = cos(t258);
t434 = -t208 * t61 + t212 * t257;
t43 = -t208 * t257 - t212 * t61;
t432 = m(11) - t412;
t431 = t386 * t43;
t430 = t434 * t386;
t428 = t439 * t410;
t237 = pkin(2) ^ 2;
t232 = pkin(6) ^ 2;
t340 = t232 + t239;
t321 = -pkin(13) ^ 2 + t340;
t210 = sin(pkin(20));
t214 = cos(pkin(20));
t188 = t210 * t221 + t214 * t216;
t393 = pkin(6) * t188;
t334 = pkin(1) * t393;
t158 = t237 - t321 - 0.2e1 * t334;
t181 = 0.2e1 * t334;
t343 = t181 + t232;
t418 = -pkin(2) - pkin(13);
t153 = (pkin(1) - t418) * (pkin(1) + t418) + t343;
t411 = pkin(13) - pkin(2);
t154 = (pkin(1) - t411) * (pkin(1) + t411) + t343;
t366 = t154 * t153;
t240 = sqrt(-t366);
t238 = 0.1e1 / pkin(2);
t401 = t238 / 0.2e1;
t314 = 0.1e1 / pkin(13) * t401;
t145 = atan2(t240 * t314, t158 * t314);
t143 = sin(t145);
t144 = cos(t145);
t427 = mrSges(10,1) * t143 + mrSges(10,2) * t144 - mrSges(9,2);
t426 = mrSges(10,1) * t144 - mrSges(10,2) * t143 + t437;
t159 = t181 + t237 + t321;
t180 = -pkin(1) - t393;
t189 = t210 * t216 - t214 * t221;
t361 = t189 * t240;
t331 = pkin(6) * t361;
t137 = -t159 * t180 - t331;
t136 = 0.1e1 / t137 ^ 2;
t392 = pkin(6) * t189;
t138 = t159 * t392 - t180 * t240;
t164 = t181 + t340;
t162 = 0.1e1 / t164;
t163 = 0.1e1 / t164 ^ 2;
t368 = 0.2e1 / t240 * (t153 + t154) * pkin(1) * t392;
t319 = -t368 / 0.2e1;
t394 = pkin(6) * t138;
t399 = pkin(1) * t189;
t406 = t162 / 0.2e1;
t76 = 0.2e1 * (((t180 * t319 + (t159 * t188 - t361) * pkin(6)) * t406 + (-t162 * t189 * t232 + t163 * t394) * t399) / t137 - ((-t188 * t240 + (t319 - t159) * t189) * t406 + (t137 * t163 + t162 * t180) * t399) * t136 * t394) * pkin(2) / (t136 * t138 ^ 2 + 0.1e1) * t164 * t238;
t425 = t76 * t427;
t424 = t76 * t426;
t413 = -m(3) - m(9);
t211 = cos(pkin(23));
t207 = sin(pkin(23));
t405 = t207 / 0.2e1;
t128 = (-t211 * t140 / 0.2e1 + t141 * t405) * t362;
t127 = 0.1e1 / t128 ^ 2;
t129 = (t211 * t141 / 0.2e1 + t140 * t405) * t362;
t294 = t211 * t307;
t295 = t207 * t307;
t316 = t169 * t405;
t365 = t169 * t211;
t69 = ((t116 * t316 + t140 * t295 + t141 * t294 + t365 * t407) / t128 - (t122 * t316 - t140 * t294 + t141 * t295 + t365 * t408) * t129 * t127) / (t127 * t129 ^ 2 + 0.1e1) * t236;
t409 = t69 + 0.1e1;
t225 = cos(pkin(18));
t402 = t225 / 0.2e1;
t398 = pkin(1) * t217;
t396 = pkin(4) * t208;
t395 = pkin(4) * t212;
t172 = t177 * pkin(5);
t174 = t179 * pkin(5);
t187 = t193 * pkin(5);
t286 = -t178 * t44 - t179 * t45;
t284 = t193 * t45 - t280 * t44;
t385 = mrSges(3,1) * t217;
t384 = mrSges(3,2) * t222;
t383 = mrSges(6,2) * t220;
t111 = atan2(t129, t128);
t106 = sin(t111);
t372 = t106 * t217;
t371 = t106 * t222;
t107 = cos(t111);
t370 = t107 * t217;
t369 = t107 * t222;
t219 = sin(pkin(18));
t364 = t169 * t219;
t230 = 0.1e1 / pkin(8);
t363 = t169 * t230;
t355 = t217 * t223;
t354 = t218 * t222;
t317 = t162 * t401;
t126 = atan2(t138 * t317, t137 * t317);
t124 = cos(t126);
t350 = t222 * t124;
t349 = t222 * t223;
t347 = -mrSges(4,1) * t177 + mrSges(4,2) * t176;
t346 = -mrSges(4,1) * t179 - mrSges(4,2) * t178;
t205 = t223 * pkin(16);
t344 = pkin(5) * t178 + t205;
t341 = mrSges(4,1) * t193 - mrSges(4,2) * t280;
t336 = m(11) + t441;
t333 = pkin(1) * t354;
t332 = pkin(1) * t349;
t329 = t410 * t44;
t327 = t409 * t222;
t326 = t409 * t217;
t315 = t169 * t402;
t303 = t107 * t327;
t302 = t106 * t326;
t293 = t219 * t307;
t292 = t225 * t307;
t160 = t175 - t339;
t183 = pkin(1) * t192 - pkin(7);
t139 = -pkin(1) * t359 - t160 * t183;
t142 = pkin(1) * t160 * t194 - t183 * t241;
t133 = (t139 * t402 - t219 * t142 / 0.2e1) * t363;
t134 = (t142 * t402 + t139 * t219 / 0.2e1) * t363;
t113 = atan2(t134, t133);
t109 = sin(t113);
t110 = cos(t113);
t290 = mrSges(7,1) * t110 - mrSges(7,2) * t109;
t123 = sin(t126);
t94 = t123 * t217 - t350;
t89 = t94 * t223;
t90 = -t123 * t349 - t124 * t355;
t289 = t143 * t90 - t144 * t89;
t288 = t143 * t89 + t144 * t90;
t282 = t370 + t371;
t281 = t123 * t222 + t124 * t217;
t68 = t94 * t76;
t88 = t281 * t218;
t270 = -mrSges(8,3) - mrSges(7,3) - mrSges(3,3) + mrSges(2,2) - mrSges(5,3) - mrSges(11,3) - mrSges(10,3) - mrSges(9,3) - mrSges(4,3);
t269 = t107 * t326 + t371 * t409;
t267 = -m(7) * pkin(15) + mrSges(2,1) + t290 - t384;
t87 = t123 * t356 - t218 * t350;
t266 = (-t143 * t88 - t144 * t87) * mrSges(10,1) + (t143 * t87 - t144 * t88) * mrSges(10,2);
t265 = mrSges(10,1) * t289 + mrSges(10,2) * t288;
t264 = (-t143 * t94 + t144 * t281) * mrSges(10,1) + (-t143 * t281 - t144 * t94) * mrSges(10,2);
t202 = t223 * t215;
t198 = pkin(1) * t356;
t157 = 0.1e1 / t158 ^ 2;
t130 = 0.1e1 / t133 ^ 2;
t121 = t183 * t318 + t239 * pkin(7) * t421 + (-t160 * t192 - t359) * pkin(1);
t115 = (t360 + (0.2e1 * pkin(7) * t183 - t160 + t318) * t194) * pkin(1);
t114 = (0.1e1 / t158 * t368 / 0.2e1 + t157 * t331 * t423) / (-t157 * t366 + 0.1e1);
t105 = t106 * t355;
t104 = t107 * t354;
t84 = t369 - t372;
t82 = t282 * t223;
t81 = -t107 * t349 + t105;
t80 = t282 * t218;
t79 = t106 * t356 - t104;
t70 = ((t121 * t315 + t142 * t292 + t115 * t364 / 0.2e1 + t139 * t293) / t133 - (t115 * t315 + t139 * t292 - t121 * t364 / 0.2e1 - t142 * t293) * t134 * t130) / (t130 * t134 ^ 2 + 0.1e1) * t230;
t53 = -t106 * t327 - t370 * t409;
t52 = -t303 + t302;
t51 = t105 + (t372 * t69 - t303) * t223;
t50 = t269 * t223;
t49 = -t104 + (-t369 * t69 + t302) * t218;
t48 = t269 * t218;
t24 = t215 * t218 + t220 * t29;
t23 = -t215 * t29 + t218 * t220;
t17 = t436 * t430;
t16 = t436 * t431;
t14 = t435 * t430;
t13 = t435 * t431;
t1 = [(mrSges(8,1) * t82 - mrSges(8,2) * t81 - m(6) * (pkin(10) * t29 + t344) - m(5) * t344 - m(11) * (-t395 * t82 - t396 * t81 + t205) - t29 * mrSges(5,1) - m(10) * (pkin(2) * t90 + t205) - t24 * mrSges(6,1) - t23 * mrSges(6,2) + t288 * mrSges(10,1) - t289 * mrSges(10,2) - mrSges(9,1) * t90 - mrSges(9,2) * t89 - mrSges(4,1) * t178 + mrSges(4,2) * t179 - (-t43 * t82 - t434 * t81) * mrSges(11,1) - (t43 * t81 - t434 * t82) * mrSges(11,2) - t304 * t286 + t412 * (-pkin(1) * t355 + t205) + (t413 * pkin(16) + (pkin(1) * t336 + mrSges(3,1)) * t217 - t267) * t223 + t270 * t218) * g(2) + (-t202 * mrSges(6,1) - m(11) * (t395 * t80 + t396 * t79 + t198) - mrSges(4,1) * t176 - mrSges(4,2) * t177 - mrSges(8,1) * t80 + mrSges(8,2) * t79 - (t43 * t80 + t434 * t79) * mrSges(11,1) - (-t43 * t79 + t434 * t80) * mrSges(11,2) - t427 * t87 + t426 * t88 - t439 * t27 + t304 * t287 + t412 * (-pkin(16) * t218 + t198) + (-t385 + (m(10) + t336 - t413) * pkin(16) + t267) * t218 + (t270 - t383) * t223 - t441 * (pkin(5) * t176 + t198)) * g(1), (-t266 - (-t208 * t48 + t212 * t49) * t440 - (t13 * t79 - t14 * t80 + t43 * t49 - t434 * t48) * mrSges(11,1) - (t13 * t80 + t14 * t79 + t43 * t48 + t434 * t49) * mrSges(11,2) - mrSges(8,1) * t49 - mrSges(8,2) * t48 - t347 - mrSges(9,2) * t88 + t437 * t87 + t432 * t333 - t441 * (-t172 - t333) + t287 * t429 + t27 * t446) * g(2) + (-t265 - mrSges(8,1) * t51 - mrSges(8,2) * t50 - (-t208 * t50 + t212 * t51) * t440 - (t13 * t81 - t14 * t82 + t43 * t51 - t434 * t50) * mrSges(11,1) - (t13 * t82 + t14 * t81 + t43 * t50 + t434 * t51) * mrSges(11,2) - t346 + mrSges(9,2) * t90 + t437 * t89 + t286 * t429 + t432 * t332 - t441 * (-t174 - t332) + t29 * t446) * g(1) + (-t264 - t290 * t70 - mrSges(9,2) * t94 - mrSges(8,1) * t53 - mrSges(8,2) * t52 - t341 - (-t208 * t52 + t212 * t53) * t440 - (-t13 * t282 + t14 * t84 + t43 * t53 - t434 * t52) * mrSges(11,1) - (-t13 * t84 - t14 * t282 + t43 * t52 + t434 * t53) * mrSges(11,2) + t384 + t385 - t437 * t281 + t284 * t429 - t441 * (t187 - t398) + t432 * t398 + t444 * t445) * g(3) + (mrSges(3,1) * t222 - mrSges(3,2) * t217 - (-mrSges(7,1) * t109 - mrSges(7,2) * t110) * t70) * (g(1) * t223 + g(2) * t218), (-(-t16 * t282 + t17 * t84) * mrSges(11,1) - (-t16 * t84 - t17 * t282) * mrSges(11,2) - t264 * t114 - t341 + t427 * t68 - t281 * t424 - t441 * t187 + t284 * t428 + t410 * t445) * g(3) + (-t347 - (t16 * t79 - t17 * t80) * mrSges(11,1) - (t16 * t80 + t17 * t79) * mrSges(11,2) - t266 * t114 + t439 * (t176 * t329 - t18 * t35 - t35) + t441 * t172 + t304 * (-t18 * t33 - t381 * t410 - t33) + t88 * t425 + t426 * t218 * t68) * g(2) + (-(t16 * t81 - t17 * t82) * mrSges(11,1) - (t16 * t82 + t17 * t81) * mrSges(11,2) - t346 - t265 * t114 + t441 * t174 + t304 * (-t179 * t329 + t18 * t38 + t38) + t281 * t223 * t425 + t89 * t424 + t286 * t428) * g(1), -g(1) * (mrSges(6,1) * t23 - mrSges(6,2) * t24) - g(2) * ((-t215 * t27 - t220 * t223) * mrSges(6,1) + (-t220 * t27 + t202) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t215 - t383) * t438];
taug = t1(:);
