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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1TE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:13:56
% EndTime: 2020-04-12 20:17:35
% DurationCPUTime: 53.18s
% Computational Cost: add. (948776->426), mult. (1438748->679), div. (59868->15), fcn. (915119->28), ass. (0->293)
t358 = sin(qJ(2));
t359 = sin(pkin(19));
t361 = cos(qJ(2));
t362 = cos(pkin(19));
t273 = t358 * t362 - t359 * t361;
t265 = pkin(7) * t273;
t255 = (-0.2e1 * t265 + pkin(1)) * pkin(1);
t383 = pkin(7) ^ 2;
t136 = t255 + t383;
t182 = pkin(3) ^ 2;
t324 = pkin(8) ^ 2 - t182;
t254 = t136 - t324;
t258 = -t265 + pkin(1);
t152 = t358 * t359 + t361 * t362;
t372 = pkin(7) + pkin(8);
t373 = pkin(7) - pkin(8);
t181 = sqrt(-((-pkin(3) + t372) * (pkin(3) + t373) + t255) * ((pkin(3) + t372) * (-pkin(3) + t373) + t255));
t333 = t152 * t181;
t313 = pkin(7) * t333;
t385 = 0.1e1 / pkin(3);
t248 = t385 * (t254 * t258 - t313);
t377 = 0.1e1 / t136;
t247 = t377 * t248;
t244 = t247 / 0.2e1;
t253 = pkin(7) * t254;
t249 = t385 * (t152 * t253 + t181 * t258);
t320 = t377 / 0.2e1;
t246 = t249 * t320;
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t240 = t244 * t357 + t246 * t360;
t245 = -t247 / 0.2e1;
t241 = t245 * t360 + t246 * t357;
t314 = pkin(23) + pkin(22);
t295 = sin(t314);
t296 = cos(t314);
t233 = t240 * t296 - t241 * t295;
t232 = pkin(4) * t233;
t224 = (0.2e1 * t232 + pkin(5)) * pkin(5);
t370 = pkin(11) - pkin(9);
t371 = -pkin(9) - pkin(11);
t180 = sqrt(-((pkin(4) - t370) * (pkin(4) + t370) + t224) * ((pkin(4) - t371) * (pkin(4) + t371) + t224));
t389 = t240 * t295 + t241 * t296;
t416 = t180 * t389;
t420 = pkin(4) * t416;
t419 = pkin(5) * t416;
t184 = pkin(9) ^ 2;
t323 = pkin(11) ^ 2 - t184;
t231 = pkin(5) * t233;
t384 = pkin(5) ^ 2;
t72 = t384 + (0.2e1 * t231 + pkin(4)) * pkin(4);
t221 = t72 + t323;
t226 = t232 + pkin(5);
t381 = 0.1e1 / pkin(11);
t206 = t381 * (t221 * t226 - t420);
t219 = pkin(4) * t221;
t207 = t381 * (t180 * t226 + t219 * t389);
t339 = sin(pkin(21));
t341 = cos(pkin(21));
t418 = t339 * t206 + t341 * t207;
t173 = pkin(13) ^ 2;
t176 = pkin(6) ^ 2;
t177 = pkin(2) ^ 2;
t164 = sin(pkin(20));
t166 = cos(pkin(20));
t147 = t164 * t360 + t166 * t357;
t354 = pkin(6) * t147;
t365 = -t173 / 0.2e1 + t176 / 0.2e1 - t177 / 0.2e1 + (t354 + pkin(1) / 0.2e1) * pkin(1);
t379 = 0.2e1 * pkin(1);
t326 = t354 * t379 + t176;
t374 = -pkin(2) - pkin(13);
t119 = (pkin(1) - t374) * (pkin(1) + t374) + t326;
t369 = pkin(13) - pkin(2);
t120 = (pkin(1) - t369) * (pkin(1) + t369) + t326;
t148 = t164 * t357 - t166 * t360;
t353 = pkin(6) * t148;
t317 = pkin(1) * t353;
t179 = sqrt(-t120 * t119);
t363 = t179 / 0.2e1;
t337 = (t119 + t120) * t317 / t363;
t300 = -t337 / 0.4e1;
t126 = pkin(1) ^ 2 + t326;
t124 = 0.1e1 / t126;
t178 = 0.1e1 / pkin(2);
t297 = 0.1e1 / t126 ^ 2 * t317;
t122 = t126 - t173 + t177;
t141 = -pkin(1) - t354;
t111 = t122 * t353 - t141 * t179;
t302 = t358 * t111;
t334 = t148 * t179;
t110 = -pkin(6) * t334 - t122 * t141;
t305 = t361 * t110;
t310 = -t361 / 0.2e1;
t301 = -t337 / 0.2e1;
t86 = (-t147 * t179 + (t141 * t379 - t122 + t301) * t148) * pkin(6);
t380 = -0.2e1 * pkin(1);
t88 = t141 * t301 + t176 * t148 ^ 2 * t380 + (t122 * t147 - t334) * pkin(6);
t417 = t178 * (t124 * (t358 * t88 / 0.2e1 + t86 * t310) - (t305 - t302) * t297);
t167 = sin(qJ(4));
t170 = cos(qJ(4));
t171 = cos(qJ(1));
t151 = -t357 * t358 + t360 * t361;
t168 = sin(qJ(1));
t137 = t151 * t168;
t150 = -t357 * t361 - t358 * t360;
t138 = t150 * t168;
t378 = 0.1e1 / t72;
t329 = t378 / 0.2e1;
t45 = t418 * t329;
t201 = t341 * t206;
t202 = t339 * t207;
t330 = -t378 / 0.2e1;
t46 = t201 * t330 + t202 * t329;
t28 = -t137 * t46 - t138 * t45;
t415 = t171 * t167 + t170 * t28;
t414 = t167 * t28 - t170 * t171;
t413 = -0.2e1 * t389;
t261 = t273 * t181;
t356 = pkin(1) * t152;
t316 = pkin(7) * t356;
t336 = 0.4e1 / t181 * (t372 * t373 - t182 + t255) * t316;
t299 = -t336 / 0.2e1;
t87 = (t261 + (t299 + t324 - t383) * t152) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t265) * t316;
t412 = -t87 / 0.2e1;
t318 = t152 ^ 2 * t380;
t93 = -t313 + t258 * t336 / 0.2e1 - t273 * t253 + t383 * t318;
t411 = t93 / 0.2e1;
t410 = m(6) + m(5);
t409 = pkin(4) * m(11);
t163 = sin(pkin(22));
t408 = t163 / 0.2e1;
t165 = cos(pkin(22));
t407 = -t165 / 0.2e1;
t406 = -t389 / 0.2e1;
t405 = t339 / 0.2e1;
t404 = -t341 / 0.2e1;
t375 = m(10) * pkin(2);
t403 = -mrSges(9,1) - t375;
t402 = -pkin(10) * m(6) - mrSges(6,1) * t170 + mrSges(6,2) * t167 - mrSges(5,1);
t303 = t358 * t110;
t304 = t361 * t111;
t335 = t124 * t178;
t95 = (t303 / 0.2e1 + t304 / 0.2e1) * t335;
t376 = pkin(4) * pkin(5);
t319 = 0.1e1 / t72 ^ 2 * t376;
t401 = t418 * t319;
t400 = (-t201 + t202) * t319;
t390 = 0.4e1 / t180 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t224 - t184) * t376;
t215 = t389 * t390;
t399 = -t180 * t233 + t215 * t406;
t131 = 0.1e1 / t136 ^ 2;
t242 = t248 * t316;
t243 = t249 * t316;
t321 = -t377 / 0.2e1;
t279 = t385 * t360 * t321;
t355 = t385 * t377;
t294 = t357 * t355;
t236 = t294 * t412 + t93 * t279 + (-t242 * t357 - t243 * t360) * t131;
t237 = t87 * t279 + t294 * t411 + (-t242 * t360 + t243 * t357) * t131;
t212 = -t236 * t295 + t237 * t296;
t60 = t236 * t296 + t237 * t295;
t216 = t60 * t390;
t398 = -t180 * t212 + t216 * t406;
t220 = t72 - t323;
t218 = pkin(5) * t220;
t225 = -t231 - pkin(4);
t397 = 0.2e1 * t225 * t376 - t218;
t396 = -0.2e1 * t226 * t376 - t219;
t382 = 0.1e1 / pkin(9);
t208 = t382 * (-t180 * t225 + t218 * t389);
t204 = t208 * t319;
t53 = -t220 * t225 - t419;
t366 = t53 * t382;
t291 = t319 * t366;
t395 = -t163 * t204 - t165 * t291;
t394 = t163 * t291 - t165 * t204;
t392 = m(8) + m(11) + m(4);
t123 = t136 + t324;
t143 = pkin(1) * t273 - pkin(7);
t112 = -pkin(1) * t333 - t123 * t143;
t113 = t123 * t356 - t143 * t181;
t172 = cos(pkin(18));
t175 = 0.1e1 / pkin(8);
t262 = -mrSges(3,1) * t358 - mrSges(3,2) * t361;
t169 = sin(pkin(18));
t298 = t169 * t321;
t391 = -m(3) * pkin(16) - mrSges(2,1) - t262 + m(7) * pkin(15) + (-mrSges(7,1) * (t112 * t172 * t320 + t113 * t298) - mrSges(7,2) * (t113 * t172 * t321 + t112 * t298)) * t175;
t388 = pkin(12) * m(6) - mrSges(5,2) + mrSges(6,3);
t386 = -mrSges(4,3) - mrSges(10,3) - mrSges(5,3) - mrSges(11,3) - mrSges(3,3) - mrSges(8,3) - mrSges(9,3) + mrSges(2,2) - mrSges(7,3);
t368 = t382 * t378;
t367 = t381 * t378;
t364 = -t179 / 0.2e1;
t133 = t138 * pkin(5);
t140 = t150 * t171;
t135 = t140 * pkin(5);
t146 = t151 * pkin(5);
t214 = t216 / 0.2e1;
t227 = pkin(5) * pkin(4) ^ 2 * t413;
t311 = t381 * t329;
t343 = t180 * t60;
t186 = (-pkin(4) * t343 + t212 * t219 + t214 * t226 + t227 * t60) * t311;
t188 = (pkin(4) * t398 + t396 * t60) * t367;
t15 = t341 * t186 + t188 * t405 + t401 * t60;
t350 = t15 + t46;
t16 = t339 * t186 + t188 * t404 + t400 * t60;
t349 = t16 - t45;
t213 = t215 / 0.2e1;
t190 = (t213 * t226 + t219 * t233 + t227 * t389 - t420) * t311;
t192 = (pkin(4) * t399 + t389 * t396) * t367;
t19 = t341 * t190 + t192 * t405 + t389 * t401;
t348 = t19 + t46;
t20 = t339 * t190 + t192 * t404 + t389 * t400;
t347 = t20 - t45;
t346 = -t137 * t45 + t138 * t46;
t139 = t151 * t171;
t30 = -t139 * t45 + t140 * t46;
t345 = t150 * t45 + t151 * t46;
t344 = pkin(16) * t168;
t162 = t171 * pkin(16);
t340 = cos(pkin(23));
t338 = sin(pkin(23));
t331 = 0.1e1 / pkin(13) * t178;
t328 = mrSges(4,1) * t138 - mrSges(4,2) * t137;
t327 = mrSges(4,1) * t140 - mrSges(4,2) * t139;
t325 = mrSges(4,1) * t151 + mrSges(4,2) * t150;
t312 = t382 * t330;
t309 = pkin(1) * t361;
t308 = t358 * pkin(1);
t307 = mrSges(10,1) * t331;
t306 = mrSges(10,2) * t331;
t153 = t168 * t308 - t344;
t293 = t168 * t309;
t292 = t171 * t309;
t290 = t340 * t355;
t289 = t331 * t365;
t288 = t331 * t363;
t285 = t131 * t175 * t316;
t283 = -pkin(5) * t137 + t153;
t278 = t385 * t338 * t320;
t277 = t169 * t285;
t276 = t172 * t285;
t154 = -t171 * t308 + t162;
t89 = t168 * t95;
t94 = (t302 / 0.2e1 - t305 / 0.2e1) * t335;
t90 = t168 * t94;
t272 = t288 * t90 + t289 * t89;
t91 = t171 * t95;
t92 = t171 * t94;
t271 = -t288 * t91 + t289 * t92;
t270 = t364 * t89 + t365 * t90;
t269 = t92 * t363 + t365 * t91;
t268 = pkin(5) * t139 + t154;
t96 = t245 * t340 + t246 * t338;
t97 = t244 * t338 + t246 * t340;
t78 = -t358 * t97 + t361 * t96;
t263 = t358 * t96 + t361 * t97;
t67 = t290 * t412 + t93 * t278 + (-t242 * t340 + t243 * t338) * t131;
t68 = t290 * t411 + t87 * t278 + (t242 * t338 + t243 * t340) * t131;
t58 = -t358 * t68 + t361 * t67 - t263;
t59 = -t358 * t67 - t361 * t68 - t78;
t252 = pkin(1) * t175 * (t261 + (0.2e1 * pkin(7) * t143 - t123 + t299) * t152) * t320;
t251 = t175 * t377 * (t143 * t299 + (pkin(7) * t318 - t123 * t273 - t333) * pkin(1));
t66 = ((t88 * t310 - t358 * t86 / 0.2e1) * t124 + (-t304 - t303) * t297) * t178;
t228 = pkin(4) * t384 * t413;
t205 = t208 * t330;
t193 = (pkin(5) * t399 + t389 * t397) * t368;
t191 = (-t213 * t225 + t218 * t233 + t228 * t389 - t419) * t312;
t189 = (pkin(5) * t398 + t397 * t60) * t368;
t187 = (-pkin(5) * t343 + t212 * t218 - t214 * t225 + t228 * t60) * t312;
t76 = t78 * t171;
t75 = t263 * t171;
t74 = t78 * t168;
t73 = t263 * t168;
t70 = t172 * t251 / 0.2e1 + t113 * t276 + t169 * t252 + t112 * t277;
t69 = t172 * t252 + t112 * t276 - t169 * t251 / 0.2e1 - t113 * t277;
t64 = t171 * t417;
t63 = t171 * t66;
t62 = t168 * t417;
t61 = t168 * t66;
t57 = t58 * t171;
t56 = t59 * t171;
t55 = t58 * t168;
t54 = t59 * t168;
t44 = t165 * t312 * t53 + t163 * t205;
t43 = t163 * t329 * t366 + t165 * t205;
t31 = t139 * t46 + t140 * t45;
t26 = t167 * t168 + t170 * t31;
t25 = -t167 * t31 + t168 * t170;
t18 = t163 * t191 + t193 * t407 + t389 * t395;
t17 = t165 * t191 + t193 * t408 + t389 * t394;
t14 = t163 * t187 + t189 * t407 + t395 * t60;
t13 = t165 * t187 + t189 * t408 + t394 * t60;
t1 = [(-t26 * mrSges(6,1) + t269 * t307 - (t163 * t76 - t165 * t75) * t409 - mrSges(4,1) * t139 - mrSges(4,2) * t140 - m(10) * (-pkin(2) * t91 + t162) - t31 * mrSges(5,1) + mrSges(8,1) * t75 + mrSges(8,2) * t76 - t271 * mrSges(10,2) + mrSges(9,1) * t91 - mrSges(9,2) * t92 - m(9) * t162 - m(6) * (pkin(10) * t31 + t268) - m(5) * t268 - t25 * mrSges(6,2) - (-t43 * t76 - t44 * t75) * mrSges(11,1) - (t43 * t75 - t44 * t76) * mrSges(11,2) + t388 * t30 - t392 * t154 + t391 * t171 + t386 * t168) * g(2) + (-m(6) * (pkin(10) * t28 + t283) - m(5) * t283 + t270 * t306 - t272 * mrSges(10,1) + mrSges(4,1) * t137 + mrSges(4,2) * t138 - mrSges(9,1) * t89 + mrSges(9,2) * t90 + t414 * mrSges(6,2) - mrSges(8,1) * t73 - mrSges(8,2) * t74 - t415 * mrSges(6,1) - (t43 * t74 + t44 * t73) * mrSges(11,1) - (-t43 * t73 + t44 * t74) * mrSges(11,2) + m(9) * t344 - m(10) * (pkin(2) * t89 - t344) - (-t163 * t74 + t165 * t73) * t409 - t28 * mrSges(5,1) - t388 * t346 - t392 * t153 - t391 * t168 + t386 * t171) * g(1), (-mrSges(9,2) * t89 - t272 * mrSges(10,2) - t270 * t307 - t54 * mrSges(8,1) + t55 * mrSges(8,2) - (t163 * t55 + t165 * t54) * t409 - (-t13 * t74 - t14 * t73 - t43 * t55 + t44 * t54) * mrSges(11,1) - (t13 * t73 - t14 * t74 - t43 * t54 - t44 * t55) * mrSges(11,2) - t328 + t403 * t90 + t392 * t293 - t410 * (t133 - t293) + t402 * (t137 * t16 + t138 * t15 + t346) + t388 * (-t137 * t350 + t138 * t349)) * g(2) + (-t56 * mrSges(8,1) + t57 * mrSges(8,2) - (t163 * t57 + t165 * t56) * t409 - (-t13 * t76 - t14 * t75 - t43 * t57 + t44 * t56) * mrSges(11,1) - (t13 * t75 - t14 * t76 - t43 * t56 - t44 * t57) * mrSges(11,2) - t271 * mrSges(10,1) - t269 * t306 - mrSges(9,2) * t91 - t327 + t403 * t92 + t402 * (t139 * t16 + t140 * t15 + t30) - t410 * (t135 - t292) + t388 * (-t139 * t350 + t140 * t349) + t392 * t292) * g(1) + (-mrSges(7,1) * t70 - mrSges(7,2) * t69 - t58 * mrSges(8,1) - t59 * mrSges(8,2) - (-t163 * t59 + t165 * t58) * t409 - (-t13 * t263 + t14 * t78 + t43 * t59 + t44 * t58) * mrSges(11,1) - (-t13 * t78 - t14 * t263 - t43 * t58 + t44 * t59) * mrSges(11,2) - t262 - t325 + (-t306 * t365 - t307 * t364 - mrSges(9,2)) * t94 - (-t306 * t363 - t307 * t365 + t403) * t95 + t402 * (t15 * t151 - t150 * t16 + t345) - t410 * (t146 - t308) + t388 * (t150 * t350 + t151 * t349) + t392 * t308) * g(3) + (mrSges(3,1) * t361 - mrSges(7,1) * t69 - mrSges(3,2) * t358 + mrSges(7,2) * t70) * (g(1) * t171 + g(2) * t168), (-(-t17 * t263 + t18 * t78) * mrSges(11,1) - (-t17 * t78 - t18 * t263) * mrSges(11,2) - t325 + t417 * t375 - ((-t300 * t95 + t317 * t94 + t364 * t66 - t365 * t417) * mrSges(10,1) + (t300 * t94 + t317 * t95 - t363 * t417 + t365 * t66) * mrSges(10,2)) * t331 + mrSges(9,1) * t417 - mrSges(9,2) * t66 - t410 * t146 + t402 * (-t150 * t20 + t151 * t19 + t345) + t388 * (t150 * t348 + t151 * t347)) * g(3) + (-t61 * t375 - ((t300 * t90 + t317 * t89 + t364 * t62 + t365 * t61) * mrSges(10,1) + (t300 * t89 - t317 * t90 + t363 * t61 + t365 * t62) * mrSges(10,2)) * t331 - mrSges(9,1) * t61 - mrSges(9,2) * t62 - (-t17 * t74 - t18 * t73) * mrSges(11,1) - (t17 * t73 - t18 * t74) * mrSges(11,2) - t328 + t402 * (t137 * t20 + t138 * t19 + t346) - t410 * t133 + t388 * (-t137 * t348 + t138 * t347)) * g(2) + (-(-t17 * t76 - t18 * t75) * mrSges(11,1) - (t17 * t75 - t18 * t76) * mrSges(11,2) - t327 - mrSges(9,1) * t63 - mrSges(9,2) * t64 - t63 * t375 - ((t300 * t92 + t317 * t91 + t364 * t64 + t365 * t63) * mrSges(10,1) + (t300 * t91 - t92 * t317 + t63 * t363 + t64 * t365) * mrSges(10,2)) * t331 - t410 * t135 + t388 * (-t139 * t348 + t140 * t347) + t402 * (t139 * t20 + t140 * t19 + t30)) * g(1), -g(1) * (mrSges(6,1) * t25 - mrSges(6,2) * t26) - g(2) * (mrSges(6,1) * t414 + mrSges(6,2) * t415) - g(3) * (-mrSges(6,1) * t167 - mrSges(6,2) * t170) * (-t150 * t46 + t151 * t45)];
taug = t1(:);
