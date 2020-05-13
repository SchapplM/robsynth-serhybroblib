% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% tauc [13x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m2OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:43
% EndTime: 2020-05-02 21:19:36
% DurationCPUTime: 26.16s
% Computational Cost: add. (8800->718), mult. (22821->1080), div. (0->0), fcn. (18570->20), ass. (0->355)
t396 = qJD(3) + qJD(4);
t295 = qJD(2) + t396;
t314 = cos(qJ(4));
t315 = cos(qJ(3));
t283 = t315 * t314 * pkin(1);
t307 = sin(qJ(3));
t479 = pkin(1) * t307;
t288 = pkin(5) + t479;
t306 = sin(qJ(4));
t238 = t288 * t306 - t283;
t468 = pkin(5) * qJD(3);
t215 = t238 * qJD(2) + t306 * t468;
t205 = pkin(11) * t295 + t215;
t305 = sin(qJ(5));
t313 = cos(qJ(5));
t308 = sin(qJ(2));
t411 = t308 * t315;
t316 = cos(qJ(2));
t412 = t307 * t316;
t263 = t411 + t412;
t248 = t263 * qJD(1);
t410 = t315 * t316;
t264 = -t307 * t308 + t410;
t251 = t264 * qJD(1);
t348 = t314 * t248 + t251 * t306;
t369 = -t248 * t306 + t314 * t251;
t287 = pkin(5) * t307 + pkin(1);
t229 = -pkin(5) * t410 + t287 * t308 - pkin(15);
t408 = qJD(1) * t229;
t79 = -pkin(9) * t369 - pkin(11) * t348 + t408;
t40 = -t205 * t305 + t313 * t79;
t41 = t205 * t313 + t305 * t79;
t166 = qJD(5) - t369;
t494 = -t166 / 0.2e1;
t136 = t295 * t313 - t305 * t348;
t497 = -t136 / 0.2e1;
t522 = Ifges(5,2) * t369;
t461 = Ifges(5,4) * t348;
t137 = t295 * t305 + t313 * t348;
t496 = -t137 / 0.2e1;
t553 = Ifges(6,5) * t496;
t548 = t553 + t461 / 0.2e1;
t554 = 0.2e1 * Ifges(6,6) * t497 + 0.2e1 * Ifges(6,3) * t494 + t295 * Ifges(5,6) + t522 / 0.2e1 + t548 - t40 * mrSges(6,1) + t41 * mrSges(6,2);
t550 = -mrSges(5,1) * t408 + t215 * mrSges(5,3) + t554;
t304 = sin(qJ(6));
t547 = -t304 / 0.2e1;
t312 = cos(qJ(6));
t546 = t312 / 0.2e1;
t303 = sin(qJ(7));
t311 = cos(qJ(7));
t261 = -t303 * t316 - t308 * t311;
t247 = t261 * qJD(1);
t406 = qJD(1) * t316;
t407 = qJD(1) * t308;
t250 = -t303 * t407 + t311 * t406;
t300 = sin(pkin(19));
t434 = cos(pkin(19));
t470 = sin(qJ(10));
t471 = cos(qJ(10));
t258 = t300 * t470 + t434 * t471;
t325 = -t300 * t471 + t434 * t470;
t114 = -t247 * t325 - t250 * t258;
t545 = Ifges(11,4) * t114;
t543 = Ifges(5,2) / 0.2e1;
t393 = t306 * t479;
t403 = qJD(4) * t306;
t189 = -qJD(3) * t393 + t283 * t396 - t288 * t403;
t384 = qJD(4) * t468;
t149 = qJD(2) * t189 - t306 * t384;
t299 = qJD(2) + qJD(3);
t211 = t299 * t264;
t195 = t211 * qJD(1);
t212 = t299 * t263;
t323 = qJD(1) * t212;
t67 = qJD(4) * t369 + t314 * t195 - t306 * t323;
t26 = qJD(5) * t136 + t313 * t67;
t27 = -qJD(5) * t137 - t305 * t67;
t8 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t541 = -m(6) * t149 + t8;
t442 = t295 * Ifges(5,5);
t540 = t305 * t41;
t370 = -t258 * t247 + t250 * t325;
t109 = Ifges(11,4) * t370;
t298 = qJD(2) + qJD(7);
t293 = qJD(10) + t298;
t44 = Ifges(11,1) * t114 + t293 * Ifges(11,5) + t109;
t539 = t44 + t109;
t282 = pkin(5) * t411;
t254 = pkin(5) * t412 + t282;
t537 = mrSges(5,2) * t408 + t442 / 0.2e1;
t535 = (Ifges(7,5) * t546 + Ifges(7,6) * t547) * qJD(6);
t165 = Ifges(5,4) * t369;
t414 = t306 * t315;
t284 = pkin(1) * t414;
t416 = t288 * t314;
t239 = t284 + t416;
t217 = qJD(2) * t239 + t314 * t468;
t206 = -pkin(9) * t295 - t217;
t356 = mrSges(6,1) * t305 + mrSges(6,2) * t313;
t334 = t206 * t356;
t449 = t137 * Ifges(6,4);
t38 = t136 * Ifges(6,2) + t166 * Ifges(6,6) + t449;
t135 = Ifges(6,4) * t136;
t39 = t137 * Ifges(6,1) + t166 * Ifges(6,5) + t135;
t438 = t313 * t39;
t481 = t305 / 0.2e1;
t517 = t348 * Ifges(5,1);
t93 = t165 + t442 + t517;
t534 = -t438 / 0.2e1 + t38 * t481 - t93 / 0.2e1 - t165 / 0.2e1 - t334;
t302 = sin(qJ(8));
t310 = cos(qJ(8));
t259 = -t302 * t316 - t308 * t310;
t246 = t259 * qJD(1);
t249 = -t302 * t407 + t310 * t406;
t301 = sin(qJ(9));
t309 = cos(qJ(9));
t350 = -t246 * t309 + t249 * t301;
t164 = Ifges(10,4) * t350;
t430 = qJD(1) * pkin(15);
t225 = -pkin(2) * t246 - t430;
t297 = qJD(2) + qJD(8);
t294 = qJD(9) + t297;
t349 = t246 * t301 + t249 * t309;
t415 = t297 * t309;
t392 = pkin(2) * t415;
t361 = qJD(9) * t392;
t477 = pkin(2) * t301;
t387 = qJD(9) * t477;
t362 = t297 * t387;
t454 = Ifges(10,4) * t349;
t492 = -t349 / 0.2e1;
t207 = t297 * t259;
t191 = t207 * qJD(1);
t346 = t302 * t308 - t310 * t316;
t208 = t297 * t346;
t192 = t208 * qJD(1);
t65 = qJD(9) * t350 - t191 * t309 - t192 * t301;
t66 = qJD(9) * t349 + t191 * t301 - t192 * t309;
t90 = Ifges(10,2) * t350 + t294 * Ifges(10,6) - t454;
t92 = -Ifges(10,1) * t349 + t294 * Ifges(10,5) + t164;
t532 = mrSges(10,1) * t362 + mrSges(10,2) * t361 + Ifges(10,5) * t65 + Ifges(10,6) * t66 - (Ifges(10,5) * t350 + Ifges(10,6) * t349) * t294 / 0.2e1 - t225 * (-mrSges(10,1) * t349 + mrSges(10,2) * t350) + t90 * t492 + (Ifges(10,1) * t350 + t454) * t349 / 0.2e1 - (Ifges(10,2) * t349 + t164 + t92) * t350 / 0.2e1;
t531 = m(8) + m(4);
t530 = t114 / 0.2e1;
t529 = t247 / 0.2e1;
t528 = t251 / 0.2e1;
t527 = -t348 / 0.2e1;
t526 = t350 / 0.2e1;
t499 = -t370 / 0.2e1;
t525 = t370 / 0.2e1;
t520 = t246 * Ifges(9,2);
t519 = t297 * Ifges(9,5);
t518 = t297 * Ifges(9,6);
t235 = t258 * qJD(10);
t236 = t325 * qJD(10);
t339 = pkin(1) * (t258 * t303 + t311 * t325);
t516 = -qJD(2) * t339 + (t235 * t300 + t236 * t434) * pkin(4);
t338 = pkin(1) * (-t258 * t311 + t303 * t325);
t515 = -qJD(2) * t338 + (-t235 * t434 + t236 * t300) * pkin(4);
t514 = -mrSges(4,1) * t251 - mrSges(8,1) * t247 + mrSges(4,2) * t248 + mrSges(8,2) * t250;
t204 = t263 * t314 + t264 * t306;
t401 = qJD(5) * t313;
t203 = -t306 * t263 + t264 * t314;
t78 = qJD(4) * t203 + t211 * t314 - t212 * t306;
t343 = t204 * t401 + t305 * t78;
t513 = -t305 * t40 + t313 * t41;
t404 = qJD(2) * t316;
t190 = qJD(2) * t282 + qJD(3) * t254 + t287 * t404;
t409 = qJD(1) * t190;
t68 = qJD(4) * t348 + t195 * t306 + t314 * t323;
t11 = pkin(9) * t68 - pkin(11) * t67 + t409;
t413 = t307 * t314;
t188 = qJD(4) * t416 + (t315 * t403 + (t413 + t414) * qJD(3)) * pkin(1);
t148 = qJD(2) * t188 + t314 * t384;
t6 = qJD(5) * t40 + t11 * t305 + t148 * t313;
t429 = qJD(5) * t41;
t7 = t11 * t313 - t148 * t305 - t429;
t512 = -t7 * mrSges(6,1) + t6 * mrSges(6,2);
t511 = pkin(9) * t348 - pkin(11) * t369;
t508 = t26 / 0.2e1;
t507 = t27 / 0.2e1;
t43 = Ifges(11,2) * t370 + t293 * Ifges(11,6) + t545;
t504 = -t43 / 0.2e1;
t503 = t68 / 0.2e1;
t500 = pkin(2) * (-mrSges(10,1) * t350 - mrSges(10,2) * t349);
t498 = -t114 / 0.2e1;
t495 = t137 / 0.2e1;
t488 = t248 / 0.2e1;
t487 = t249 / 0.2e1;
t486 = t250 / 0.2e1;
t484 = -t293 / 0.2e1;
t482 = t304 / 0.2e1;
t478 = pkin(2) * t192;
t476 = pkin(4) * t300;
t474 = t305 * t7;
t473 = t313 * t6;
t469 = pkin(1) * qJD(2);
t467 = mrSges(4,3) * t307;
t466 = mrSges(6,3) * t313;
t465 = mrSges(8,3) * t250;
t464 = Ifges(7,1) * t304;
t463 = Ifges(3,4) * t308;
t462 = Ifges(3,4) * t316;
t460 = Ifges(6,4) * t305;
t459 = Ifges(6,4) * t313;
t458 = Ifges(7,4) * t304;
t457 = Ifges(7,4) * t312;
t456 = Ifges(8,4) * t250;
t455 = Ifges(9,4) * t249;
t452 = Ifges(7,2) * t312;
t446 = t369 * mrSges(5,3);
t444 = t217 * mrSges(5,3);
t443 = t248 * Ifges(4,4);
t437 = t313 * t40;
t286 = t308 * pkin(1) - pkin(15);
t435 = mrSges(5,1) * t295 + mrSges(6,1) * t136 - mrSges(6,2) * t137 - mrSges(5,3) * t348;
t433 = mrSges(4,2) * qJD(3);
t390 = t303 * t469;
t255 = t298 * t476 + t390;
t380 = t434 * pkin(4);
t388 = t311 * t469;
t256 = t298 * t380 + t388;
t115 = t255 * t325 - t256 * t258;
t432 = mrSges(11,3) * t115;
t116 = t255 * t258 + t256 * t325;
t431 = mrSges(11,3) * t116;
t428 = t114 * t116;
t427 = t370 * mrSges(11,3);
t426 = (-mrSges(10,2) * t294 + mrSges(10,3) * t350) * t309;
t424 = t348 * t215;
t423 = t204 * t305;
t422 = t204 * t313;
t421 = (mrSges(4,1) * t299 - mrSges(4,3) * t248) * t315;
t420 = t229 * (mrSges(5,1) * t348 + mrSges(5,2) * t369);
t317 = qJD(1) ^ 2;
t419 = t229 * t317;
t253 = pkin(1) * t413 + t284;
t405 = qJD(2) * t253;
t402 = qJD(5) * t305;
t400 = t286 * qJD(1);
t395 = pkin(2) * m(10) * t225;
t394 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t68;
t389 = t315 * t469;
t386 = pkin(1) * t404;
t385 = pkin(1) * t406;
t381 = t438 / 0.2e1;
t377 = -qJD(2) * t308 / 0.2e1;
t376 = -t404 / 0.2e1;
t375 = -t402 / 0.2e1;
t372 = mrSges(10,3) * t297 * t477;
t368 = mrSges(4,3) * t299;
t367 = t467 * t469;
t366 = mrSges(4,3) * t389;
t365 = mrSges(8,3) * t388;
t360 = qJD(3) * t389;
t359 = qJD(2) * t385;
t355 = Ifges(6,1) * t313 - t460;
t354 = -Ifges(6,2) * t305 + t459;
t353 = Ifges(6,5) * t313 - Ifges(6,6) * t305;
t352 = t437 + t540;
t209 = t298 * t261;
t193 = t209 * qJD(1);
t345 = t303 * t308 - t311 * t316;
t210 = t298 * t345;
t194 = t210 * qJD(1);
t31 = -t193 * t258 - t194 * t325 - t235 * t247 + t236 * t250;
t32 = t193 * t325 - t194 * t258 + t235 * t250 + t236 * t247;
t326 = qJD(7) * t338;
t72 = qJD(2) * t326 - t235 * t256 + t236 * t255;
t327 = qJD(7) * t339;
t73 = qJD(2) * t327 + t235 * t255 + t236 * t256;
t351 = t73 * mrSges(11,1) - t72 * mrSges(11,2) + Ifges(11,5) * t31 + Ifges(11,6) * t32;
t201 = -t259 * t309 - t301 * t346;
t347 = t259 * t301 - t309 * t346;
t342 = t204 * t402 - t313 * t78;
t341 = pkin(14) * (mrSges(7,1) * t304 + mrSges(7,2) * t312);
t340 = pkin(15) * (mrSges(3,1) * t316 - mrSges(3,2) * t308);
t337 = t136 * t354;
t336 = t137 * t355;
t335 = t166 * t353;
t333 = t304 * (Ifges(7,1) * t312 - t458);
t332 = t308 * (-Ifges(3,2) * t316 - t463);
t331 = t316 * (-Ifges(3,1) * t308 - t462);
t330 = (Ifges(3,1) * t316 - t463) * qJD(1);
t329 = (-Ifges(3,2) * t308 + t462) * qJD(1);
t328 = (t452 + t458) * qJD(1);
t153 = (-t247 * t300 + t250 * t434) * pkin(4);
t324 = -qJD(5) * t352 - t474;
t4 = Ifges(6,4) * t26 + Ifges(6,2) * t27 + t68 * Ifges(6,6);
t5 = Ifges(6,1) * t26 + Ifges(6,4) * t27 + t68 * Ifges(6,5);
t322 = -t148 * mrSges(5,2) + t313 * t4 / 0.2e1 + t5 * t481 + (Ifges(6,1) * t305 + t459) * t508 + (Ifges(6,2) * t313 + t460) * t507 + t6 * t466 + t38 * t375 + (Ifges(6,5) * t305 + Ifges(6,6) * t313) * t503 - Ifges(5,6) * t68 + Ifges(5,5) * t67 + (mrSges(6,1) * t313 - mrSges(6,2) * t305 + mrSges(5,1)) * t149 + (t334 + t381) * qJD(5) + (t337 + t336 + t335) * qJD(5) / 0.2e1;
t159 = t247 * Ifges(8,2) + t298 * Ifges(8,6) + t456;
t233 = Ifges(8,4) * t247;
t162 = Ifges(8,1) * t250 + t298 * Ifges(8,5) + t233;
t321 = -t114 * t504 + t159 * t486 - t250 * (Ifges(8,1) * t247 - t456) / 0.2e1 - t298 * (Ifges(8,5) * t247 - Ifges(8,6) * t250) / 0.2e1 + (Ifges(11,1) * t370 - t545) * t498 + Ifges(8,6) * t194 + Ifges(8,5) * t193 + t247 * t365 + t351 + (Ifges(11,5) * t370 - Ifges(11,6) * t114) * t484 + t370 * t432 + (-Ifges(11,2) * t114 + t539) * t499 - (-Ifges(8,2) * t250 + t162 + t233) * t247 / 0.2e1;
t10 = -mrSges(6,2) * t68 + t27 * mrSges(6,3);
t83 = -mrSges(6,2) * t166 + mrSges(6,3) * t136;
t84 = mrSges(6,1) * t166 - mrSges(6,3) * t137;
t9 = mrSges(6,1) * t68 - t26 * mrSges(6,3);
t320 = -t84 * t401 - t83 * t402 + m(6) * (-t40 * t401 - t402 * t41 + t473 - t474) + t313 * t10 - t305 * t9;
t144 = mrSges(10,1) * t294 + mrSges(10,3) * t349;
t158 = t455 + t518 + t520;
t232 = Ifges(9,4) * t246;
t161 = t249 * Ifges(9,1) + t232 + t519;
t319 = t349 * t372 + t158 * t487 - t249 * t500 - t249 * (Ifges(9,1) * t246 - t455) / 0.2e1 - t249 * t395 + t144 * t387 - t297 * (Ifges(9,5) * t246 - Ifges(9,6) * t249) / 0.2e1 + Ifges(9,5) * t191 + Ifges(9,6) * t192 + (mrSges(9,1) * t249 + mrSges(9,2) * t246) * t430 - (-Ifges(9,2) * t249 + t161 + t232) * t246 / 0.2e1 + (-qJD(9) * t426 + (-t301 * t66 + t309 * t65 - t350 * t415) * mrSges(10,3)) * pkin(2) + t532;
t160 = t251 * Ifges(4,2) + t299 * Ifges(4,6) + t443;
t234 = Ifges(4,4) * t251;
t163 = t248 * Ifges(4,1) + t299 * Ifges(4,5) + t234;
t318 = t322 + t160 * t488 - t461 * t527 - t248 * (Ifges(4,1) * t251 - t443) / 0.2e1 - Ifges(4,6) * t323 - t299 * (Ifges(4,5) * t251 - Ifges(4,6) * t248) / 0.2e1 + Ifges(4,5) * t195 + mrSges(4,1) * t360 - t248 * t366 + t251 * t367 - (-Ifges(4,2) * t248 + t163 + t234) * t251 / 0.2e1 + (t553 + t554) * t348 + (mrSges(6,3) * t540 + t40 * t466 + t353 * t494 + t355 * t496 + t354 * t497 - t442 / 0.2e1 + t444 + t348 * t543 + Ifges(5,1) * t527 + t534) * t369;
t291 = qJD(1) * t457;
t270 = pkin(1) * t311 + t380;
t269 = pkin(1) * t303 + t476;
t252 = t283 - t393;
t245 = Ifges(3,5) * qJD(2) + t330;
t244 = Ifges(7,5) * qJD(6) + qJD(1) * t464 + t291;
t243 = Ifges(3,6) * qJD(2) + t329;
t242 = Ifges(7,6) * qJD(6) + t328;
t237 = t287 * t316 + t282;
t231 = -pkin(9) - t239;
t230 = pkin(11) + t238;
t221 = -mrSges(4,2) * t299 + mrSges(4,3) * t251;
t220 = mrSges(8,1) * t298 - t465;
t218 = -mrSges(8,2) * t298 + mrSges(8,3) * t247;
t199 = mrSges(4,1) * t248 + mrSges(4,2) * t251;
t196 = mrSges(8,1) * t250 + mrSges(8,2) * t247;
t187 = (-t258 * t300 - t325 * t434) * pkin(4);
t186 = (-t258 * t434 + t300 * t325) * pkin(4);
t157 = (-t261 * t434 + t300 * t345) * pkin(4) + t286;
t147 = t153 + t385;
t143 = -mrSges(5,2) * t295 + t446;
t140 = -t258 * t269 - t270 * t325;
t139 = -t258 * t270 + t269 * t325;
t138 = t400 + (-t247 * t434 - t250 * t300) * pkin(4);
t126 = t258 * t345 - t261 * t325;
t125 = -t258 * t261 - t325 * t345;
t105 = mrSges(11,1) * t293 - mrSges(11,3) * t114;
t104 = -mrSges(11,2) * t293 + t427;
t103 = t386 + (-t209 * t300 - t210 * t434) * pkin(4);
t100 = -mrSges(5,1) * t369 + mrSges(5,2) * t348;
t87 = t359 + (-t193 * t300 - t194 * t434) * pkin(4);
t86 = t235 * t269 + t236 * t270 + t327;
t85 = -t235 * t270 + t236 * t269 + t326;
t82 = qJD(1) * t254 + t511;
t80 = qJD(1) * t237 + t511;
t77 = qJD(4) * t204 + t211 * t306 + t314 * t212;
t76 = qJD(9) * t347 + t207 * t301 - t208 * t309;
t75 = qJD(9) * t201 - t207 * t309 - t208 * t301;
t70 = t217 * t313 + t305 * t511;
t69 = -t217 * t305 + t313 * t511;
t64 = t305 * t82 + t313 * t405;
t63 = -t305 * t405 + t313 * t82;
t51 = -mrSges(11,1) * t370 + mrSges(11,2) * t114;
t50 = mrSges(11,1) * t114 + mrSges(11,2) * t370;
t36 = t209 * t325 - t210 * t258 - t235 * t345 + t236 * t261;
t35 = -t209 * t258 - t210 * t325 - t235 * t261 - t236 * t345;
t1 = [t166 * (-Ifges(6,5) * t342 - Ifges(6,6) * t343) / 0.2e1 + t136 * (-Ifges(6,4) * t342 - Ifges(6,2) * t343) / 0.2e1 + (mrSges(5,2) * t409 + Ifges(5,1) * t67 + t353 * t503 + t354 * t507 + t355 * t508 + t375 * t39) * t204 + (-Ifges(6,1) * t342 - Ifges(6,4) * t343) * t495 + (t195 * t263 + t211 * t488) * Ifges(4,1) + (t342 * t40 - t343 * t41 - t422 * t7 - t423 * t6) * mrSges(6,3) + (t125 * t72 - t126 * t73) * mrSges(11,3) - t191 * Ifges(9,1) * t346 + (m(6) * t352 + t305 * t83 + t313 * t84) * (pkin(9) * t77 - pkin(11) * t78 + t190) + (-m(10) * t478 - mrSges(10,1) * t66 + mrSges(10,2) * t65) * (-pkin(2) * t259 - pkin(15)) + (t330 + t245) * t377 + (Ifges(6,5) * t495 - t522 / 0.2e1 - t550) * t77 + (t329 + t243) * t376 + (-t193 * t345 + t209 * t486) * Ifges(8,1) + (-mrSges(4,1) * t264 - mrSges(8,1) * t261 + mrSges(4,2) * t263 - mrSges(8,2) * t345) * t359 + ((t311 * t261 - t303 * t345) * qJD(7) * t469 + t210 * t390) * mrSges(8,3) + (t193 * t261 - t194 * t345 + t209 * t529 + t210 * t486) * Ifges(8,4) + (t208 * t487 + t259 * t191 - t192 * t346 + t246 * t207 / 0.2e1) * Ifges(9,4) + (mrSges(4,1) * t212 - mrSges(8,1) * t210 + mrSges(4,2) * t211 + mrSges(8,2) * t209) * t400 + t299 * (Ifges(4,5) * t211 - Ifges(4,6) * t212) / 0.2e1 + (mrSges(4,1) * t323 - mrSges(8,1) * t194 + t195 * mrSges(4,2) + mrSges(8,2) * t193 + 0.2e1 * t531 * t359) * t286 + (-t212 * t528 - t264 * t323) * Ifges(4,2) + (t264 * t195 + t211 * t528 - t212 * t488 - t263 * t323) * Ifges(4,4) + (-t347 * t65 + t492 * t75) * Ifges(10,1) + (-t201 * t361 + t347 * t362 + t392 * t75) * mrSges(10,3) - (-mrSges(10,1) * t201 - mrSges(10,2) * t347) * t478 + (t201 * t65 - t347 * t66 + t76 * t492 + t526 * t75) * Ifges(10,4) + (0.2e1 * m(5) * t408 + t100) * t190 + ((t244 + qJD(1) * (t457 + t464)) * t546 + (t328 + t242) * t547 + (t333 + t312 * (-Ifges(7,2) * t304 + t457) + 0.2e1 * t341) * qJD(1) + t535) * qJD(6) + (t381 - t444 + t93 / 0.2e1 + t517 / 0.2e1 + t537) * t78 + (qJD(3) * t264 - t211) * t367 + t206 * (mrSges(6,1) * t343 - mrSges(6,2) * t342) + (-t356 - mrSges(5,3)) * t149 * t204 + (t331 - t332 - 0.2e1 * t340) * qJD(1) * qJD(2) + (t125 * t31 + t32 * t126 + t35 * t525 + t36 * t530) * Ifges(11,4) + (t31 * t126 + t35 * t530) * Ifges(11,1) + (t125 * t32 + t36 * t525) * Ifges(11,2) + (t201 * t66 + t526 * t76) * Ifges(10,2) + (-t68 * t204 + t67 * t203 + t77 * t527 + t369 * t78 / 0.2e1) * Ifges(5,4) + (t83 * t401 + t313 * t9 + t305 * t10 - t84 * t402 + m(6) * (qJD(5) * t513 + t305 * t6 + t313 * t7)) * (-pkin(9) * t203 - pkin(11) * t204 + t229) - t343 * t38 / 0.2e1 + t514 * t386 + t259 * Ifges(9,2) * t192 + (t520 / 0.2e1 - t500 + mrSges(9,1) * t430 - t395 + t518 / 0.2e1 + t158 / 0.2e1) * t208 + (Ifges(9,1) * t487 - mrSges(9,2) * t430 + t519 / 0.2e1 + t161 / 0.2e1) * t207 + (-mrSges(5,1) * t409 - Ifges(6,3) * t503 - Ifges(6,6) * t507 - Ifges(6,5) * t508 + t148 * mrSges(5,3) - t394 / 0.2e1 - Ifges(5,2) * t68 + t512) * t203 + (t103 * t138 + t157 * t87) * m(11) + (t194 * t261 + t210 * t529) * Ifges(8,2) - t35 * t432 - t4 * t423 / 0.2e1 + t5 * t422 / 0.2e1 + t298 * (Ifges(8,5) * t209 + Ifges(8,6) * t210) / 0.2e1 + t293 * (Ifges(11,5) * t35 + Ifges(11,6) * t36) / 0.2e1 + t294 * (Ifges(10,5) * t75 + Ifges(10,6) * t76) / 0.2e1 + t229 * (mrSges(5,1) * t68 + mrSges(5,2) * t67) + t225 * (-mrSges(10,1) * t76 + mrSges(10,2) * t75) - t76 * t372 + t209 * t162 / 0.2e1 + t210 * t159 / 0.2e1 + t211 * t163 / 0.2e1 - t212 * t160 / 0.2e1 + t212 * t366 - pkin(15) * (-mrSges(9,1) * t192 + mrSges(9,2) * t191) - t209 * t365 + t157 * (-mrSges(11,1) * t32 + mrSges(11,2) * t31) + t138 * (-mrSges(11,1) * t36 + mrSges(11,2) * t35) + t87 * (-mrSges(11,1) * t125 + mrSges(11,2) * t126) + t103 * t51 + qJD(2) ^ 2 * (-Ifges(3,5) * t308 - Ifges(3,6) * t316) / 0.2e1 + t76 * t90 / 0.2e1 + t75 * t92 / 0.2e1 + t36 * t43 / 0.2e1 + t35 * t44 / 0.2e1 - t36 * t431 - t263 * mrSges(4,3) * t360; t321 + t318 + t319 + (-t188 * t84 - t80 * t83 + (-qJD(5) * t83 - t9) * t230 + (-t7 - t429) * mrSges(6,3)) * t305 + (-t195 * t467 - t531 * t317 * t316 * t286 + (t218 * t311 - t220 * t303) * qJD(7) + (t221 * t307 + t421) * qJD(3) + (-t193 * t311 + t194 * t303) * mrSges(8,3) + (-t307 * t433 + t303 * t465 + (-mrSges(8,1) * t303 - mrSges(8,2) * t311) * qJD(7)) * qJD(2)) * pkin(1) + (t230 * t10 + t188 * t83 - t80 * t84 + (-t40 * mrSges(6,3) - t230 * t84) * qJD(5)) * t313 + m(11) * (t115 * t86 - t116 * t85 + t139 * t73 + t140 * t72) + (-m(11) * t147 - t50) * t138 + t435 * t189 + (-t238 * t68 - t239 * t67 + t424) * mrSges(5,3) + (-t139 * t31 + t140 * t32 - t428) * mrSges(11,3) + t231 * t8 + t188 * t143 - t147 * t51 + t85 * t104 + t86 * t105 + (t316 * t243 / 0.2e1 + t308 * t245 / 0.2e1 - t237 * t100 - t420 + Ifges(3,6) * t376 + (-t331 / 0.2e1 + t332 / 0.2e1 + t340) * qJD(1) + (-t196 - t199) * t286 + (t315 ^ 2 * t308 * t368 + (t307 * t315 * t368 - t514) * t316) * pkin(1) + Ifges(3,5) * t377) * qJD(1) + (t148 * t238 + t149 * t239 + t188 * t215 + t189 * t217 - t237 * t419) * m(5) + (-t352 * t80 - t149 * t231 - t189 * t206 + (t324 + t473) * t230 + t513 * t188) * m(6); (-t100 * t254 - t286 * t199 - t420) * qJD(1) + t318 - m(5) * (t254 * t419 + (t215 * t253 + t217 * t252) * qJD(2)) + t324 * mrSges(6,3) + (-t253 * t143 - t435 * t252 + (-t421 + (-t221 - t433) * t307) * pkin(1)) * qJD(2) + mrSges(5,3) * t424 - t64 * t83 - t63 * t84 - m(6) * (-qJD(2) * t206 * t252 + t40 * t63 + t41 * t64) + (m(5) * (t148 * t306 + t149 * t314) + (-t306 * t68 - t314 * t67) * mrSges(5,3) + ((-m(5) * t217 + m(6) * t206 - t435) * t306 + (m(5) * t215 + m(6) * t513 - t305 * t84 + t313 * t83 + t143) * t314) * qJD(4)) * pkin(5) + t320 * (pkin(5) * t306 + pkin(11)) + t541 * (-pkin(5) * t314 - pkin(9)); t322 + t320 * pkin(11) + (-t166 * t437 + (-t166 * t41 - t7) * t305) * mrSges(6,3) + (-t335 / 0.2e1 - t337 / 0.2e1 - t336 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + t543) * t348 + t534 - t537) * t369 + (t548 + t550) * t348 + t435 * t215 - m(6) * (t206 * t215 + t40 * t69 + t41 * t70) - t70 * t83 - t69 * t84 + (t446 - t143) * t217 - t541 * pkin(9); -t206 * (mrSges(6,1) * t137 + mrSges(6,2) * t136) + (Ifges(6,1) * t136 - t449) * t496 + t38 * t495 + (Ifges(6,5) * t136 - Ifges(6,6) * t137) * t494 - t40 * t83 + t41 * t84 + (t136 * t40 + t137 * t41) * mrSges(6,3) + t394 + (-Ifges(6,2) * t137 + t135 + t39) * t497 - t512; (t242 * t482 + (-t341 - t333 / 0.2e1 + t452 * t482) * qJD(1) + t535 - (t244 + t291) * t312 / 0.2e1) * qJD(1); t321 + t516 * t105 + t515 * t104 + (-t186 * t31 + t187 * t32 - t428) * mrSges(11,3) - t196 * t400 + ((-mrSges(8,2) * qJD(7) - t218) * t311 + (-mrSges(8,1) * qJD(7) + t220 + t465) * t303) * t469 - t153 * t51 - t138 * t50 + (t115 * t516 - t116 * t515 - t138 * t153 + t186 * t73 + t187 * t72) * m(11); t319; (t426 - t144 * t301 + (t301 * t349 - t309 * t350) * mrSges(10,3)) * t297 * pkin(2) + t532; -t116 * t105 + (-t104 + t427) * t115 + t351 + t539 * t499 + (-t138 * mrSges(11,2) + Ifges(11,1) * t498 + Ifges(11,5) * t484) * t370 - (t138 * mrSges(11,1) + Ifges(11,4) * t498 + Ifges(11,2) * t499 + Ifges(11,6) * t484 + t431 + t504) * t114; 0; 0; 0;];
tauc = t1(:);
