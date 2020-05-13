% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:11:24
% EndTime: 2020-05-03 00:12:06
% DurationCPUTime: 9.87s
% Computational Cost: add. (3491->545), mult. (6542->843), div. (0->0), fcn. (1564->8), ass. (0->335)
t230 = pkin(4) * m(6);
t192 = mrSges(5,1) + t230;
t219 = cos(qJ(5));
t409 = t219 * mrSges(6,1);
t215 = sin(qJ(5));
t413 = t215 * mrSges(6,2);
t266 = t409 - t413;
t120 = t266 + t192;
t218 = sin(qJ(2));
t239 = qJD(1) ^ 2;
t382 = t218 * t239;
t106 = pkin(2) * t120 * t382;
t190 = mrSges(6,2) * pkin(6) - Ifges(6,6);
t163 = t190 * t215;
t191 = pkin(6) * mrSges(6,1) - Ifges(6,5);
t164 = t191 * t219;
t364 = qJD(5) * qJD(1);
t181 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t391 = t181 * t239;
t309 = -pkin(1) * t391 + (0.2e1 * t163 - 0.2e1 * t164) * t364;
t494 = t106 + t309;
t350 = 2 * qJD(4);
t217 = sin(qJ(3));
t238 = qJD(2) ^ 2;
t435 = pkin(2) * t238;
t65 = t309 * t218;
t402 = t120 * t435 - t65;
t493 = t217 * t402;
t205 = t219 ^ 2;
t302 = t205 * t364;
t274 = Ifges(6,4) * t302;
t212 = Ifges(6,1) - Ifges(6,2);
t368 = qJD(5) * t215;
t303 = t212 * t368;
t373 = qJD(1) * t215;
t201 = qJD(2) + qJD(3);
t284 = t201 * t190;
t417 = pkin(6) * qJD(4);
t292 = mrSges(6,2) * t417 - Ifges(6,6) * qJD(4);
t400 = mrSges(6,1) * qJD(5);
t349 = pkin(4) * t400;
t85 = t349 / 0.2e1 + t284 + t292;
t310 = t85 * t373;
t371 = qJD(1) * t219;
t268 = t191 * qJD(4);
t399 = mrSges(6,2) * qJD(5);
t86 = -pkin(4) * t399 / 0.2e1 + t201 * t191 + t268;
t337 = 0.8e1 * t274 + (0.8e1 * t86 + 0.4e1 * t303) * t371 - 0.8e1 * t310;
t216 = sin(qJ(4));
t431 = pkin(3) * t216;
t213 = qJD(4) / 0.2e1;
t178 = t213 + t201;
t253 = qJD(1) * t120;
t60 = t178 * t253;
t357 = t60 * t431;
t161 = pkin(4) * mrSges(6,3) + pkin(6) * t230 + Ifges(5,4);
t145 = -Ifges(4,4) + t161;
t136 = t145 * qJD(3);
t307 = t161 * qJD(4);
t396 = qJD(5) * Ifges(6,4);
t276 = -t396 / 0.2e1 + t307;
t264 = t136 + t276;
t398 = qJD(1) * (qJD(2) * t145 + t264);
t492 = 0.8e1 * t357 + 0.8e1 * t398 + t337;
t59 = t201 * t161 + t276;
t397 = qJD(1) * t59;
t491 = t337 + 0.8e1 * t397;
t419 = pkin(1) * qJD(1);
t124 = mrSges(6,1) * t419 - t292;
t174 = t212 * t205;
t200 = qJD(3) + qJD(4);
t185 = qJD(2) + t200;
t427 = -qJD(2) / 0.2e1;
t299 = t427 - qJD(3) / 0.2e1;
t420 = Ifges(6,4) * t215;
t328 = t219 * t420;
t348 = mrSges(6,2) * t419;
t351 = 0.2e1 * qJD(2);
t235 = qJD(5) ^ 2;
t236 = qJD(4) ^ 2;
t377 = -t235 + t236;
t452 = 0.2e1 * qJD(3);
t237 = qJD(3) ^ 2;
t473 = -t237 - t238;
t477 = pkin(6) * t377;
t487 = t192 * t419;
t489 = t419 * t350;
t1 = (-0.2e1 * t487 - 0.2e1 * (-qJD(4) / 0.2e1 + t299) * Ifges(5,6) + (0.4e1 * t328 - 0.2e1 * t174 - Ifges(6,3) + t212) * qJD(5)) * t185 + (mrSges(6,1) * t477 + mrSges(6,2) * t489 - t377 * Ifges(6,5) - t473 * t191 + (mrSges(6,1) * t417 - Ifges(6,5) * qJD(4) + t348) * t351) * t215 - (mrSges(6,1) * t489 - mrSges(6,2) * t477 + t377 * Ifges(6,6) + t124 * t351 + t473 * t190) * t219 + ((t191 * qJD(2) + t268 + t348) * t215 + (t190 * qJD(2) - t124) * t219) * t452;
t450 = -mrSges(5,3) / 0.2e1;
t421 = mrSges(6,2) * t219;
t422 = mrSges(6,1) * t215;
t472 = t421 + t422;
t490 = -t472 * pkin(3) * (qJD(5) + t201) * (-qJD(5) + t201) + 0.2e1 * t201 * (mrSges(4,2) * t419 + t201 * (pkin(3) * t450 + Ifges(4,5) / 0.2e1)) - t1 * t216;
t312 = t181 * t382;
t135 = pkin(2) * t312;
t170 = 0.2e1 * qJD(5) * pkin(4) + t419;
t330 = t170 * t413;
t426 = Ifges(6,3) - Ifges(6,2);
t476 = -(t487 + qJD(5) * (Ifges(6,1) + t426)) * qJD(1) - (mrSges(6,1) * t170 + 0.4e1 * Ifges(6,4) * t368) * t371;
t18 = -qJD(1) * t330 - 0.2e1 * t212 * t302 - t476;
t251 = t135 + t18;
t486 = t217 * t251;
t220 = cos(qJ(4));
t485 = t494 * t220;
t375 = qJD(1) * t185;
t275 = t375 * t174;
t293 = -Ifges(6,5) * qJD(5) + pkin(6) * t400;
t429 = t185 * pkin(4);
t358 = mrSges(6,2) * t429;
t98 = t293 - 0.2e1 * t358;
t315 = t98 * t373;
t390 = t185 * t215;
t329 = Ifges(6,4) * t390;
t196 = pkin(6) * t399;
t291 = -Ifges(6,6) * qJD(5) + t196;
t359 = mrSges(6,1) * t429;
t97 = t291 + 0.2e1 * t359;
t339 = 0.2e1 * t275 + (-0.2e1 * t97 - 0.4e1 * t329) * t371 - 0.2e1 * t315;
t376 = qJD(1) * t181;
t379 = -t472 * t364 / 0.2e1;
t76 = t178 * t376 + t379;
t356 = t76 * t431;
t445 = pkin(6) * mrSges(6,3);
t202 = 0.2e1 * t445;
t240 = pkin(6) ^ 2;
t229 = t240 * m(6);
t298 = Ifges(5,2) - Ifges(5,1) + t426;
t241 = pkin(4) ^ 2;
t428 = t241 * m(6);
t255 = -t229 + t298 + t428;
t123 = t202 - t255;
t290 = t123 + t174;
t411 = t215 * t98;
t11 = (-t411 - t97 * t219 + (t290 - 0.2e1 * t328) * t185) * qJD(1);
t206 = t220 ^ 2;
t415 = t11 * t206;
t484 = 0.4e1 * t356 + t339 - 0.4e1 * t415;
t134 = (mrSges(6,1) * pkin(4) + t420) * t219;
t355 = pkin(4) * t413;
t483 = t174 - 0.2e1 * t134 + 0.2e1 * t355;
t372 = qJD(1) * t216;
t381 = t219 * t216;
t447 = pkin(3) * t76;
t482 = 0.8e1 * (t185 * t290 - t411) * t372 - (0.8e1 * t97 + 0.16e2 * t329) * qJD(1) * t381 + 0.8e1 * t447;
t481 = -0.2e1 * t201;
t480 = 0.2e1 * t205;
t479 = 0.2e1 * qJD(1);
t221 = cos(qJ(3));
t222 = cos(qJ(2));
t383 = t217 * t218;
t112 = (-t221 * t222 + t383) * qJD(1);
t113 = (-t217 * t222 - t218 * t221) * qJD(1);
t51 = t112 * t216 + t113 * t220;
t25 = t185 * t219 - t215 * t51;
t449 = -t25 / 0.2e1;
t478 = pkin(2) * t239;
t369 = qJD(5) * t205;
t318 = Ifges(6,4) * t369;
t412 = t215 * t85;
t475 = 0.16e2 * (t318 + t59 - t412) * qJD(1);
t169 = mrSges(4,1) + (m(5) + m(6)) * pkin(3);
t392 = t181 * t216;
t289 = t169 + t392;
t242 = pkin(3) ^ 2;
t232 = t242 * m(5);
t471 = Ifges(4,2) + t232 - Ifges(4,1);
t465 = -0.2e1 * t445 - t483;
t287 = 0.2e1 * t303;
t386 = t206 * t239;
t90 = -t163 + t164 + t161;
t321 = t90 * t386;
t370 = qJD(5) * t185;
t463 = -0.2e1 * t321 + (-t190 * t239 + qJD(5) * (t291 - 0.2e1 * t359)) * t215 + (-qJD(5) * (t293 + 0.2e1 * t358) + t185 * t287 + t191 * t239) * t219 + (0.2e1 * t480 - 0.2e1) * Ifges(6,4) * t370;
t462 = 0.2e1 * pkin(2);
t461 = 0.4e1 * pkin(2);
t393 = t120 * t216;
t341 = pkin(3) * t393;
t89 = -Ifges(4,4) + t90;
t460 = -0.4e1 * (t89 + t341) * t239;
t56 = t255 + t465;
t459 = -0.2e1 * t56;
t458 = 0.2e1 * pkin(1) * t364 + 0.2e1 * (-t238 / 0.2e1 - t200 * qJD(2) - t237 / 0.2e1 - qJD(3) * qJD(4) - t236 / 0.2e1 + t235 / 0.2e1) * pkin(4);
t457 = 0.2e1 * t206;
t456 = 0.8e1 * t206;
t455 = 0.4e1 * t220;
t454 = 0.8e1 * t220;
t453 = -0.2e1 * qJD(2);
t451 = -mrSges(4,2) / 0.2e1;
t448 = pkin(3) * t60;
t444 = (-t393 / 0.2e1 + t451) * t478;
t443 = t491 * t216 + 0.4e1 * t448;
t26 = t219 * t51 + t390;
t440 = Ifges(6,4) * t26;
t439 = pkin(1) * t235;
t438 = pkin(1) * t239;
t436 = pkin(2) * t217;
t434 = pkin(3) * t120;
t433 = pkin(3) * t181;
t432 = pkin(3) * t201;
t430 = t289 * t478;
t338 = -0.4e1 * t275 + (0.4e1 * t97 + 0.8e1 * t329) * t371 + 0.4e1 * t315;
t425 = (-0.4e1 * t123 * t375 + t338) * t216 - 0.4e1 * t447;
t16 = t18 * t218;
t108 = t239 * t434;
t384 = t216 * t239;
t325 = t90 * t384;
t37 = 0.2e1 * t108 + 0.4e1 * t325;
t424 = -t37 * t217 + t16;
t407 = t239 * t90;
t78 = 0.4e1 * t321;
t423 = t78 - 0.2e1 * t407;
t418 = pkin(2) * qJD(2);
t416 = qJD(3) * pkin(3);
t410 = t216 * t90;
t408 = t239 * t89;
t288 = t108 * t216;
t405 = -0.2e1 * t288 + t78;
t154 = pkin(3) * t391;
t143 = 0.2e1 * t154;
t246 = t134 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t205 + Ifges(6,3) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - t355 - t445;
t46 = t428 / 0.2e1 - t229 / 0.2e1 + t246;
t322 = t46 * t384;
t404 = -0.4e1 * t322 + t143;
t320 = t56 * t384;
t403 = -0.2e1 * t320 + t154;
t101 = t242 * m(6) + t123 + t471;
t401 = t101 * qJD(3) + qJD(4) * t123;
t193 = qJD(2) + qJD(3) / 0.2e1;
t173 = t213 + t193;
t116 = 0.2e1 * t173 * t200 + t239;
t395 = t116 * t181;
t394 = t116 * t192;
t389 = (qJD(2) * t101 + t401) * qJD(1);
t388 = t200 * t216;
t387 = t200 * t220;
t207 = t221 ^ 2;
t385 = t207 * t218;
t137 = -0.2e1 * qJD(4) * t178 + t239;
t378 = t241 / 0.2e1 - t240 / 0.2e1;
t374 = qJD(1) * t193;
t367 = qJD(5) * t216;
t366 = qJD(5) * t220;
t53 = t173 * t253;
t365 = t53 * t461;
t363 = -0.4e1 * t436;
t362 = 0.2e1 * t434;
t354 = t218 * t444;
t353 = t217 * t456;
t352 = 0.2e1 * t385;
t346 = t181 * t435;
t345 = pkin(2) * t391;
t343 = t217 * t430;
t342 = t120 * t436;
t340 = pkin(3) * t392;
t336 = 0.4e1 * t274 + (0.4e1 * t86 + t287) * t371 - 0.4e1 * t310;
t199 = mrSges(4,2) * t438;
t296 = t218 * t353;
t62 = t296 * t407;
t334 = t18 * t216 + t199 + t62;
t333 = -0.2e1 * t408 + t405;
t80 = -0.8e1 * t321;
t332 = t80 + 0.4e1 * t288 + 0.4e1 * t408;
t327 = qJD(1) * t432;
t326 = pkin(6) * t370;
t42 = t378 * m(6) + t246;
t324 = t42 * t386;
t323 = t42 * t384;
t316 = pkin(2) * t364;
t314 = qJD(1) * t418;
t311 = t206 * t382;
t243 = pkin(2) ^ 2;
t306 = -t241 + t242 - t243;
t305 = mrSges(6,1) * t364;
t304 = mrSges(6,2) * t364;
t301 = -0.2e1 * t340;
t297 = 0.2e1 * t314;
t295 = -0.4e1 * t311;
t283 = qJD(4) * t220 * t432;
t282 = t367 * t432;
t281 = pkin(3) * t304;
t280 = pkin(3) * t305;
t279 = t169 * t438 - 0.2e1 * t215 * t281 - t216 * t309 + 0.2e1 * t219 * t280;
t278 = t218 * t314;
t277 = t222 * t316;
t273 = (-t242 / 0.2e1 + t378) * m(6) - t232 / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + t246 - t340;
t270 = -t491 * t206 + t336 + 0.4e1 * t357;
t34 = t239 * (t101 + t483);
t22 = t239 * t301 - t34;
t269 = -0.4e1 * t46 * t386 + t22;
t267 = -0.4e1 * t410 - t434;
t263 = t112 * t220 - t113 * t216;
t262 = -t181 * t220 + t393;
t261 = (t137 * t192 + t266 * (-t235 + t137)) * t431 + t463;
t258 = -t409 / 0.2e1 + t413 / 0.2e1;
t252 = t472 * t370;
t92 = -pkin(4) * t304 + t191 * t375;
t99 = -Ifges(6,6) * t364 + (t196 + t359) * qJD(1);
t249 = t216 * t99 - t220 * t92 + t281;
t100 = -Ifges(6,5) * t364 + pkin(6) * t305 - qJD(1) * t358;
t91 = (t190 * t185 + t349) * qJD(1);
t248 = t100 * t216 + t220 * t91 + t280;
t224 = m(4) + m(5);
t247 = -t224 * t243 + Ifges(3,1) - Ifges(3,2) - t298 + t471;
t179 = t185 ^ 2;
t5 = Ifges(6,4) * t179 * t480 + (t179 * t212 * t215 + mrSges(6,2) * t458) * t219 + t422 * t458 + (-0.2e1 * pkin(1) * t376 + (-Ifges(6,4) + Ifges(5,5)) * t185) * t185;
t208 = t222 ^ 2;
t176 = pkin(2) * t221 + pkin(3);
t168 = t201 ^ 2 + t239;
t153 = mrSges(3,1) + (m(6) + t224) * pkin(2);
t151 = qJD(2) * t452 + t237 + t239;
t142 = -0.4e1 * t154;
t140 = t154 / 0.2e1;
t138 = Ifges(3,4) + t145;
t126 = t138 * qJD(2);
t111 = t235 + t116;
t105 = t239 * t342;
t74 = t173 * t376 + t379;
t68 = (t176 * t216 + t220 * t436) * qJD(2) + t216 * t416 + t185 * pkin(6);
t54 = (0.16e2 * t86 + 0.8e1 * t303) * qJD(1);
t40 = -0.4e1 * pkin(3) * (t137 * t181 / 0.4e1 + (t421 / 0.2e1 + t422 / 0.2e1) * t370);
t38 = t108 + 0.2e1 * t325;
t32 = t169 * t374 + t216 * t74;
t30 = (0.4e1 * t108 + 0.8e1 * t325) * t217;
t29 = -0.2e1 * t38 * t218;
t28 = -mrSges(4,2) * t374 - t216 * t53;
t21 = -mrSges(6,1) * t327 - t216 * t92 - t220 * t99;
t20 = mrSges(6,2) * t327 - t100 * t220 + t216 * t91;
t17 = (0.2e1 * t212 * t369 + t330) * qJD(1) + t476;
t14 = t279 * t218;
t13 = t273 * t239;
t12 = (t17 * t216 - t199) * t218;
t9 = ((pkin(3) * t221 + pkin(2)) * t222 - pkin(3) * t383 + pkin(1)) * qJD(1) - t51 * pkin(6) - t263 * pkin(4);
t8 = (-0.2e1 * t318 + (-0.2e1 * t86 - t303) * t219 + 0.2e1 * t412 - 0.2e1 * t59) * qJD(1);
t2 = t5 * t216 + (t299 * Ifges(4,6) + t169 * t419) * t481;
t3 = [(((t219 * t54 + t475) * t206 + t482 * t220 - t492) * t207 + (t11 * t353 + ((-t216 * t475 - t54 * t381 - 0.8e1 * t448) * t217 + t74 * t461) * t220 + (t338 - 0.8e1 * t356 - 0.4e1 * t389) * t217 + t28 * t461) * t221 + (-t217 * t365 + t425) * t220 + t32 * t363 + 0.4e1 * qJD(1) * (t126 + t264) + t270) * t208 + (0.8e1 * (t415 + (t216 * t8 - t448) * t220 - t356 - t275 / 0.2e1 + (t329 + t97 / 0.2e1) * t371 + t315 / 0.2e1 - t389 / 0.2e1) * t385 + (-t5 * t220 + t8 * t296 - t490) * t221 + (mrSges(3,2) * t419 + ((-mrSges(4,3) / 0.2e1 + t450) * pkin(2) + Ifges(3,5) / 0.2e1) * qJD(2)) * t453 + t472 * pkin(2) * (-t235 + t238) + (t1 * t220 + t2) * t217 + (((-t482 * t217 - t365) * t220 + t492 * t217 - 0.4e1 * t32 * pkin(2)) * t221 + (t74 * t363 + t443) * t220 + t28 * t363 + ((t306 * m(6) + t202 + t229 + t247) * qJD(2) + t401) * t479 + t484) * t218) * t222 + (t425 * t220 + t270 + 0.4e1 * t398) * t207 + ((t1 * t218 - 0.2e1 * t181 * t314) * t220 + t2 * t218 + (mrSges(4,2) + t393) * t297 + (t443 * t220 + 0.2e1 * t389 + t484) * t217) * t221 + (t336 + 0.4e1 * t397) * t206 + (t5 * t383 + t339 * t216 + (qJD(2) * t342 + t123 * t185 * t216 + pkin(3) * (-mrSges(6,1) * t368 + t181 * t185 - t219 * t399)) * t479) * t220 + (t490 * t217 + (Ifges(3,6) * t427 + t153 * t419) * t453) * t218 + t289 * t217 * t297 + 0.4e1 * pkin(3) * (-t192 / 0.2e1 + t258) * t185 * t372 - 0.4e1 * t274 + (-mrSges(6,2) * t439 + (-0.2e1 * t303 + (t481 - t350) * t191) * qJD(1)) * t219 + (-mrSges(6,1) * t439 + (t190 * t350 + 0.2e1 * t284) * qJD(1)) * t215 - 0.2e1 * (t126 + t136 + t307 - t396) * qJD(1); (((t142 + 0.8e1 * t323) * t220 + t332) * t207 + ((t30 - 0.2e1 * t345) * t220 - 0.4e1 * t444 + (-0.4e1 * t13 + 0.8e1 * t324) * t217) * t221 + 0.2e1 * (t105 + (-0.2e1 * t216 * t42 + t433) * t239) * t220 + 0.2e1 * t343 - 0.2e1 * t239 * (Ifges(3,4) + t89) + t405) * t208 + (0.4e1 * (t220 * t38 + (t42 * t457 - t273) * t239) * t385 + ((0.2e1 * t106 + t309) * t220 + (0.2e1 * t430 + ((t140 - t323) * t454 + t460) * t217) * t218 + t334) * t221 + t42 * t295 + ((t18 + 0.2e1 * t135) * t217 + t29) * t220 + (t279 + 0.4e1 * t354) * t217 + ((-t240 - t306) * m(6) + t301 - t247 + t465) * t382 + mrSges(3,2) * t438) * t222 + ((t143 - 0.4e1 * t323) * t220 + t333) * t207 + ((pkin(2) * t395 - t252 * t462 + t424) * t220 + t14 + ((-t394 / 0.2e1 + t258 * t111) * t216 + t151 * t451) * t462 + (-0.4e1 * t324 + 0.2e1 * t13) * t217) * t221 + ((-t65 - pkin(2) * (t266 * t111 + t394)) * t217 + 0.2e1 * t323 + t40) * t220 + (t12 - pkin(2) * ((-0.2e1 * t252 + t395) * t216 + t151 * t169)) * t217 + (t153 * t438 + 0.2e1 * t266 * t316) * t218 + t138 * t239 + t261; (((t142 + 0.8e1 * t322) * t220 + t332) * t207 + ((t30 - t345) * t220 - 0.2e1 * t444 + (0.2e1 * t34 + (t46 * t456 + 0.4e1 * t340) * t239) * t217) * t221 + (t105 + t404) * t220 + t343 + t333) * t208 + ((t220 * t37 - t269) * t352 + (t485 + (t430 + ((t140 - t322) * t454 + t460) * t217) * t218 + t334) * t221 + t46 * t295 + (t29 + t486) * t220 + (t279 + 0.2e1 * t354) * t217 + t22 * t218) * t222 + (t404 * t220 + t333) * t207 + (t424 * t220 + t14 + (mrSges(4,2) + t262) * t435 + t269 * t217) * t221 + (0.2e1 * t322 + t40 + t493) * t220 + (t289 * t435 + t12) * t217 + t145 * t239 + t261; ((-0.2e1 * t154 * t220 + t80 + (0.4e1 * t90 + (t56 * t455 + t362) * t216) * t239) * t207 + (t262 * pkin(2) + ((t362 + 0.8e1 * t410) * t220 + 0.2e1 * t340 + (0.4e1 * t206 - 0.2e1) * t56) * t217) * t239 * t221 + (t105 + t403) * t220 + (t217 * t345 - t108) * t216 + t423) * t208 + ((t108 * t220 + ((t457 - 0.1e1) * t56 + (t90 * t455 + t433) * t216) * t239) * t352 + (t62 + t485 + t251 * t216 + ((-0.4e1 * t216 * t56 + 0.2e1 * t433) * t220 - 0.2e1 * t341 - 0.4e1 * t90) * t217 * t382) * t221 + t311 * t459 + (t267 * t382 + t486) * t220 + (-pkin(3) * t312 - t494 * t217) * t216 + t56 * t382) * t222 + (t220 * t403 - t288 + t423) * t207 + ((t16 - t346) * t220 + t402 * t216 + (t206 * t459 + t220 * t267 - t340 + t56) * t239 * t217) * t221 + (-t168 * t433 + t320 + t493) * t220 + ((t17 * t218 + t346) * t217 + t168 * t434) * t216 + t161 * t239 + t463; ((t21 * t218 - t249 * t222 + (-mrSges(6,1) * t367 - mrSges(6,2) * t387) * t418) * t221 + (t249 * t218 + t21 * t222 + (-mrSges(6,1) * t366 + mrSges(6,2) * t388) * t418) * t217 - mrSges(6,1) * t278 - mrSges(6,2) * t277 - mrSges(6,1) * t282 - mrSges(6,2) * t283 - mrSges(6,1) * t326 - pkin(1) * t304 + Ifges(6,5) * t370) * t219 + ((t20 * t218 - t248 * t222 + (-mrSges(6,1) * t387 + mrSges(6,2) * t367) * t418) * t221 + (t248 * t218 + t20 * t222 + (mrSges(6,1) * t388 + mrSges(6,2) * t366) * t418) * t217 + mrSges(6,2) * t278 - mrSges(6,1) * t277 + mrSges(6,2) * t282 - mrSges(6,1) * t283 - pkin(1) * t305 + mrSges(6,2) * t326 - Ifges(6,6) * t370) * t215 - ((t216 * t222 + t220 * t218) * t221 + (-t216 * t218 + t220 * t222) * t217) * Ifges(6,3) * t375 - (-(t176 * t220 - t216 * t436) * qJD(2) - t220 * t416 - t429) * (mrSges(6,1) * t26 + mrSges(6,2) * t25) - t26 * (Ifges(6,1) * t25 - t440) / 0.2e1 + t26 * (Ifges(6,2) * t25 + t440) / 0.2e1 + (0.2e1 * Ifges(6,4) * t25 + t212 * t26) * t449 + ((t215 * t9 + t219 * t68) * mrSges(6,1) + (-t215 * t68 + t219 * t9) * mrSges(6,2) + Ifges(6,6) * t26 + 0.2e1 * Ifges(6,5) * t449) * (qJD(5) - t263);];
tauc = t3(:);
