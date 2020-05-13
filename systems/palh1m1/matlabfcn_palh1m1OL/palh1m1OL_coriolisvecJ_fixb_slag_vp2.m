% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m1OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:36
% EndTime: 2020-04-15 19:29:48
% DurationCPUTime: 15.50s
% Computational Cost: add. (8812->695), mult. (23180->1042), div. (0->0), fcn. (18560->20), ass. (0->325)
t304 = sin(qJ(4));
t313 = cos(qJ(3));
t430 = pkin(1) * qJD(2);
t374 = t313 * t430;
t275 = t304 * t374;
t312 = cos(qJ(4));
t377 = qJD(2) + qJD(3);
t365 = t377 * pkin(5);
t305 = sin(qJ(3));
t372 = t305 * t430;
t326 = t365 + t372;
t218 = t312 * t326 + t275;
t376 = qJD(3) + qJD(4);
t294 = qJD(2) + t376;
t213 = -t294 * pkin(9) - t218;
t303 = sin(qJ(5));
t311 = cos(qJ(5));
t355 = mrSges(6,1) * t303 + mrSges(6,2) * t311;
t330 = t213 * t355;
t306 = sin(qJ(2));
t314 = cos(qJ(2));
t263 = t305 * t314 + t313 * t306;
t252 = t263 * qJD(1);
t266 = -t305 * t306 + t313 * t314;
t255 = t266 * qJD(1);
t362 = -t252 * t304 + t312 * t255;
t344 = t312 * t252 + t255 * t304;
t136 = t294 * t311 - t303 * t344;
t167 = qJD(5) - t362;
t137 = t294 * t303 + t311 * t344;
t414 = t137 * Ifges(6,4);
t38 = t136 * Ifges(6,2) + t167 * Ifges(6,6) + t414;
t135 = Ifges(6,4) * t136;
t39 = t137 * Ifges(6,1) + t167 * Ifges(6,5) + t135;
t410 = t344 * Ifges(5,1);
t445 = -t311 / 0.2e1;
t447 = t303 / 0.2e1;
t507 = t218 * mrSges(5,3) + t38 * t447 + t39 * t445 - t330 - t410 / 0.2e1 - t294 * Ifges(5,5) - Ifges(5,4) * t362;
t256 = t304 * t326;
t219 = t312 * t374 - t256;
t214 = pkin(11) * t294 - t219;
t386 = qJD(1) * t306;
t398 = qJD(1) * pkin(15);
t273 = pkin(1) * t386 - t398;
t217 = -pkin(5) * t255 + t273;
t72 = -pkin(9) * t362 - pkin(11) * t344 + t217;
t40 = -t214 * t303 + t311 * t72;
t41 = t214 * t311 + t303 * t72;
t411 = t362 * Ifges(5,2);
t466 = -t167 / 0.2e1;
t467 = -t137 / 0.2e1;
t468 = -t136 / 0.2e1;
t506 = 0.2e1 * Ifges(6,5) * t467 + 0.2e1 * Ifges(6,6) * t468 + 0.2e1 * Ifges(6,3) * t466 + t294 * Ifges(5,6) + Ifges(5,4) * t344 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t411 / 0.2e1;
t301 = sin(qJ(7));
t309 = cos(qJ(7));
t262 = -t301 * t314 - t309 * t306;
t251 = t262 * qJD(1);
t385 = qJD(1) * t314;
t254 = -t301 * t386 + t309 * t385;
t298 = sin(pkin(19));
t401 = cos(pkin(19));
t431 = sin(qJ(10));
t432 = cos(qJ(10));
t260 = t298 * t431 + t401 * t432;
t324 = -t298 * t432 + t401 * t431;
t114 = -t251 * t324 - t254 * t260;
t504 = Ifges(11,4) * t114;
t388 = t312 * t313;
t375 = pkin(1) * t388;
t383 = qJD(3) * t305;
t487 = -t304 * pkin(1) * t383 + t376 * t375;
t149 = qJD(2) * t487 - qJD(4) * t256;
t211 = t377 * t263;
t196 = t211 * qJD(1);
t212 = t377 * t266;
t197 = t212 * qJD(1);
t67 = qJD(4) * t362 - t196 * t304 + t197 * t312;
t26 = qJD(5) * t136 + t311 * t67;
t27 = -qJD(5) * t137 - t303 * t67;
t8 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t502 = m(6) * t149 - t8;
t501 = t41 * t303;
t363 = -t260 * t251 + t254 * t324;
t109 = Ifges(11,4) * t363;
t297 = qJD(2) + qJD(7);
t292 = qJD(10) + t297;
t44 = t114 * Ifges(11,1) + t292 * Ifges(11,5) + t109;
t500 = t109 + t44;
t498 = (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t344;
t300 = sin(qJ(8));
t308 = cos(qJ(8));
t261 = -t300 * t314 - t308 * t306;
t250 = t261 * qJD(1);
t253 = -t300 * t386 + t308 * t385;
t299 = sin(qJ(9));
t307 = cos(qJ(9));
t346 = -t250 * t307 + t253 * t299;
t165 = Ifges(10,4) * t346;
t228 = -pkin(2) * t250 - t398;
t296 = qJD(2) + qJD(8);
t293 = qJD(9) + t296;
t345 = t250 * t299 + t253 * t307;
t441 = pkin(2) * t299;
t373 = qJD(9) * t441;
t378 = qJD(9) * t307;
t416 = Ifges(10,4) * t345;
t442 = pkin(2) * t296;
t463 = -t345 / 0.2e1;
t207 = t296 * t261;
t192 = t207 * qJD(1);
t342 = t300 * t306 - t308 * t314;
t208 = t296 * t342;
t193 = t208 * qJD(1);
t65 = qJD(9) * t346 - t192 * t307 - t193 * t299;
t66 = qJD(9) * t345 + t192 * t299 - t193 * t307;
t91 = Ifges(10,2) * t346 + Ifges(10,6) * t293 - t416;
t93 = -Ifges(10,1) * t345 + Ifges(10,5) * t293 + t165;
t495 = t296 * mrSges(10,1) * t373 + mrSges(10,2) * t378 * t442 + Ifges(10,5) * t65 + Ifges(10,6) * t66 - (Ifges(10,5) * t346 + Ifges(10,6) * t345) * t293 / 0.2e1 - t228 * (-mrSges(10,1) * t345 + mrSges(10,2) * t346) + t91 * t463 + (Ifges(10,1) * t346 + t416) * t345 / 0.2e1 - (Ifges(10,2) * t345 + t165 + t93) * t346 / 0.2e1;
t472 = -t363 / 0.2e1;
t490 = t307 * t346;
t395 = t363 * mrSges(11,3);
t239 = t260 * qJD(10);
t240 = t324 * qJD(10);
t335 = pkin(1) * (t260 * t301 + t309 * t324);
t489 = -qJD(2) * t335 + (t239 * t298 + t240 * t401) * pkin(4);
t334 = pkin(1) * (-t260 * t309 + t301 * t324);
t488 = -qJD(2) * t334 + (-t239 * t401 + t240 * t298) * pkin(4);
t486 = -t303 * t40 + t311 * t41;
t290 = pkin(1) * t385;
t282 = qJD(2) * t290;
t150 = pkin(5) * t196 + t282;
t68 = qJD(4) * t344 + t312 * t196 + t197 * t304;
t11 = pkin(9) * t68 - pkin(11) * t67 + t150;
t389 = t304 * t313;
t340 = t305 * t312 + t389;
t329 = t340 * qJD(3);
t381 = qJD(4) * t312;
t148 = t365 * t381 + (qJD(4) * t340 + t329) * t430;
t2 = qJD(5) * t40 + t11 * t303 + t148 * t311;
t397 = qJD(5) * t41;
t3 = t11 * t311 - t148 * t303 - t397;
t485 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t26 + Ifges(6,6) * t27;
t473 = m(4) + m(8);
t483 = -mrSges(4,1) * t255 - mrSges(8,1) * t251 + mrSges(4,2) * t252 + mrSges(8,2) * t254 + t273 * t473;
t102 = pkin(9) * t344 - pkin(11) * t362;
t350 = Ifges(6,5) * t311 - Ifges(6,6) * t303;
t331 = t167 * t350;
t421 = Ifges(6,4) * t303;
t354 = Ifges(6,1) * t311 - t421;
t332 = t137 * t354;
t420 = Ifges(6,4) * t311;
t352 = -Ifges(6,2) * t303 + t420;
t333 = t136 * t352;
t348 = t40 * t311 + t501;
t482 = t348 * mrSges(6,3) - t333 / 0.2e1 - t332 / 0.2e1 - t331 / 0.2e1 - t217 * mrSges(5,2) + t507;
t480 = t26 / 0.2e1;
t479 = t27 / 0.2e1;
t43 = Ifges(11,2) * t363 + t292 * Ifges(11,6) + t504;
t477 = -t43 / 0.2e1;
t476 = t68 / 0.2e1;
t471 = t363 / 0.2e1;
t470 = -t114 / 0.2e1;
t469 = t114 / 0.2e1;
t464 = t346 / 0.2e1;
t460 = t250 / 0.2e1;
t458 = t251 / 0.2e1;
t457 = -t252 / 0.2e1;
t456 = t252 / 0.2e1;
t455 = t253 / 0.2e1;
t454 = t254 / 0.2e1;
t452 = t255 / 0.2e1;
t451 = -t292 / 0.2e1;
t302 = sin(qJ(6));
t449 = -t302 / 0.2e1;
t448 = -t303 / 0.2e1;
t310 = cos(qJ(6));
t446 = t310 / 0.2e1;
t444 = t311 / 0.2e1;
t440 = pkin(4) * t298;
t439 = pkin(5) * t252;
t437 = pkin(15) * mrSges(3,1);
t436 = pkin(15) * mrSges(3,2);
t435 = t2 * t311;
t434 = t3 * t303;
t429 = mrSges(6,3) * t311;
t428 = mrSges(8,3) * t251;
t427 = mrSges(8,3) * t254;
t426 = mrSges(10,3) * t296;
t425 = Ifges(3,4) * t306;
t424 = Ifges(3,4) * t314;
t423 = Ifges(4,4) * t252;
t419 = Ifges(7,4) * t302;
t418 = Ifges(8,4) * t254;
t417 = Ifges(9,4) * t253;
t408 = t252 * mrSges(4,3);
t407 = t255 * mrSges(4,3);
t284 = t306 * pkin(1) - pkin(15);
t402 = mrSges(5,1) * t294 + mrSges(6,1) * t136 - mrSges(6,2) * t137 - mrSges(5,3) * t344;
t400 = Ifges(3,5) * qJD(2);
t399 = Ifges(3,6) * qJD(2);
t257 = t297 * t440 + t301 * t430;
t367 = t401 * pkin(4);
t371 = t309 * t430;
t258 = t297 * t367 + t371;
t116 = t257 * t260 + t258 * t324;
t396 = t114 * t116;
t394 = t344 * t219;
t393 = (mrSges(4,1) * t377 - t408) * t313;
t390 = t302 * t310;
t384 = qJD(2) * t314;
t382 = qJD(4) * t304;
t380 = qJD(5) * t303;
t379 = qJD(5) * t311;
t370 = Ifges(7,2) * t390;
t291 = pkin(1) * t384;
t181 = pkin(5) * t211 + t291;
t229 = -pkin(5) * t266 + t284;
t286 = pkin(1) * t305 + pkin(5);
t242 = t304 * t286 - t375;
t358 = -t2 * t303 - t3 * t311;
t356 = mrSges(6,1) * t311 - mrSges(6,2) * t303;
t353 = Ifges(6,1) * t303 + t420;
t351 = Ifges(6,2) * t311 + t421;
t349 = Ifges(6,5) * t303 + Ifges(6,6) * t311;
t209 = t297 * t262;
t194 = t209 * qJD(1);
t341 = t301 * t306 - t309 * t314;
t210 = t297 * t341;
t195 = t210 * qJD(1);
t31 = -t194 * t260 - t195 * t324 - t239 * t251 + t240 * t254;
t32 = t194 * t324 - t195 * t260 + t239 * t254 + t240 * t251;
t327 = qJD(7) * t334;
t73 = qJD(2) * t327 - t239 * t258 + t240 * t257;
t328 = qJD(7) * t335;
t74 = qJD(2) * t328 + t239 * t257 + t240 * t258;
t347 = t74 * mrSges(11,1) - t73 * mrSges(11,2) + Ifges(11,5) * t31 + Ifges(11,6) * t32;
t203 = -t261 * t307 - t299 * t342;
t343 = t261 * t299 - t307 * t342;
t206 = t263 * t312 + t266 * t304;
t204 = t263 * t304 - t312 * t266;
t241 = pkin(1) * t389 + t286 * t312;
t338 = m(6) * t348;
t337 = pkin(14) * (mrSges(7,1) * t302 + mrSges(7,2) * t310);
t288 = Ifges(7,4) * t310 * qJD(1);
t336 = (Ifges(7,6) * qJD(6) + (Ifges(7,2) * t310 + t419) * qJD(1)) * t449 + (Ifges(7,1) * qJD(1) * t302 + Ifges(7,5) * qJD(6) + t288) * t446;
t82 = t102 + t439;
t325 = (Ifges(7,5) * t446 + Ifges(7,6) * t449) * qJD(6);
t154 = (-t251 * t298 + t254 * t401) * pkin(4);
t323 = -qJD(5) * t348 - t434;
t322 = t323 * mrSges(6,3);
t6 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t68 * Ifges(6,6);
t7 = t26 * Ifges(6,1) + t27 * Ifges(6,4) + t68 * Ifges(6,5);
t321 = -t148 * mrSges(5,2) + t2 * t429 + t351 * t479 + t353 * t480 - t38 * t380 / 0.2e1 + t6 * t444 + t7 * t447 + qJD(5) * t330 + t39 * t379 / 0.2e1 + t349 * t476 - Ifges(5,6) * t68 + Ifges(5,5) * t67 + (t356 + mrSges(5,1)) * t149 + (t333 + t332 + t331) * qJD(5) / 0.2e1;
t115 = t257 * t324 - t258 * t260;
t160 = t251 * Ifges(8,2) + t297 * Ifges(8,6) + t418;
t237 = Ifges(8,4) * t251;
t163 = t254 * Ifges(8,1) + t297 * Ifges(8,5) + t237;
t320 = -t114 * t477 + t160 * t454 - t254 * (Ifges(8,1) * t251 - t418) / 0.2e1 - t297 * (Ifges(8,5) * t251 - Ifges(8,6) * t254) / 0.2e1 + (Ifges(11,1) * t363 - t504) * t470 + Ifges(8,6) * t195 + Ifges(8,5) * t194 + t371 * t428 + t347 + (Ifges(11,5) * t363 - Ifges(11,6) * t114) * t451 + t115 * t395 + (-Ifges(11,2) * t114 + t500) * t472 - (-Ifges(8,2) * t254 + t163 + t237) * t251 / 0.2e1;
t10 = -mrSges(6,2) * t68 + mrSges(6,3) * t27;
t83 = -mrSges(6,2) * t167 + mrSges(6,3) * t136;
t84 = mrSges(6,1) * t167 - mrSges(6,3) * t137;
t9 = mrSges(6,1) * t68 - mrSges(6,3) * t26;
t319 = m(6) * (-t379 * t40 - t380 * t41 - t434 + t435) + t311 * t10 - t303 * t9 - t84 * t379 - t83 * t380;
t318 = -t217 * mrSges(5,1) - t219 * mrSges(5,3) + t506;
t142 = -mrSges(10,2) * t293 + mrSges(10,3) * t346;
t144 = mrSges(10,1) * t293 + mrSges(10,3) * t345;
t159 = t250 * Ifges(9,2) + t296 * Ifges(9,6) + t417;
t236 = Ifges(9,4) * t250;
t162 = t253 * Ifges(9,1) + t296 * Ifges(9,5) + t236;
t99 = -mrSges(10,1) * t346 - mrSges(10,2) * t345;
t317 = t345 * t426 * t441 - t253 * (Ifges(9,1) * t250 - t417) / 0.2e1 - t296 * (Ifges(9,5) * t250 - Ifges(9,6) * t253) / 0.2e1 + Ifges(9,5) * t192 + Ifges(9,6) * t193 + (mrSges(9,1) * t253 + mrSges(9,2) * t250) * t398 + t144 * t373 + t159 * t455 - (-Ifges(9,2) * t253 + t162 + t236) * t250 / 0.2e1 + t495 + ((-m(10) * t228 - t99) * t253 - t142 * t378 + (-t296 * t490 - t299 * t66 + t307 * t65) * mrSges(10,3)) * pkin(2);
t161 = Ifges(4,2) * t255 + Ifges(4,6) * t377 + t423;
t238 = Ifges(4,4) * t255;
t164 = Ifges(4,1) * t252 + Ifges(4,5) * t377 + t238;
t316 = t321 - t377 * (Ifges(4,5) * t255 - Ifges(4,6) * t252) / 0.2e1 - Ifges(4,6) * t196 + Ifges(4,5) * t197 + t372 * t407 + t161 * t456 + (Ifges(4,1) * t255 - t423) * t457 - (-Ifges(4,2) * t252 + t164 + t238) * t255 / 0.2e1 + (mrSges(4,1) * qJD(3) - t408) * t374 + t506 * t344 + (mrSges(6,3) * t501 + t350 * t466 + t352 * t468 + t354 * t467 + t40 * t429 + t498 + t507) * t362;
t272 = t309 * pkin(1) + t367;
t271 = pkin(1) * t301 + t440;
t249 = t400 + (Ifges(3,1) * t314 - t425) * qJD(1);
t247 = t399 + (-Ifges(3,2) * t306 + t424) * qJD(1);
t244 = -pkin(2) * t261 - pkin(15);
t235 = t312 * t372 + t275;
t234 = (t304 * t305 - t388) * t430;
t233 = pkin(11) + t242;
t232 = -pkin(9) - t241;
t224 = -mrSges(4,2) * t377 + t407;
t223 = mrSges(8,1) * t297 - t427;
t221 = -mrSges(8,2) * t297 + t428;
t220 = t290 + t439;
t201 = mrSges(4,1) * t252 + mrSges(4,2) * t255;
t198 = mrSges(8,1) * t254 + mrSges(8,2) * t251;
t191 = -t286 * t382 + t487;
t190 = t286 * t381 + (t313 * t382 + t329) * pkin(1);
t189 = (-t260 * t298 - t324 * t401) * pkin(4);
t188 = (-t260 * t401 + t298 * t324) * pkin(4);
t158 = (-t262 * t401 + t298 * t341) * pkin(4) + t284;
t147 = t290 + t154;
t143 = -mrSges(5,2) * t294 + mrSges(5,3) * t362;
t140 = -t260 * t271 - t272 * t324;
t139 = -t260 * t272 + t271 * t324;
t138 = (-t251 * t401 - t254 * t298) * pkin(4) + t273;
t126 = t260 * t341 - t262 * t324;
t125 = -t260 * t262 - t324 * t341;
t105 = mrSges(11,1) * t292 - mrSges(11,3) * t114;
t104 = -mrSges(11,2) * t292 + t395;
t103 = t291 + (-t209 * t298 - t210 * t401) * pkin(4);
t101 = mrSges(5,1) * t344 + mrSges(5,2) * t362;
t100 = -mrSges(5,1) * t362 + mrSges(5,2) * t344;
t88 = t282 + (-t194 * t298 - t195 * t401) * pkin(4);
t86 = t239 * t271 + t240 * t272 + t328;
t85 = -t239 * t272 + t240 * t271 + t327;
t80 = qJD(4) * t206 + t312 * t211 + t212 * t304;
t79 = -qJD(4) * t204 - t211 * t304 + t212 * t312;
t78 = qJD(9) * t343 + t207 * t299 - t208 * t307;
t77 = t203 * qJD(9) - t207 * t307 - t208 * t299;
t76 = t290 + t82;
t70 = t102 * t303 + t218 * t311;
t69 = t102 * t311 - t218 * t303;
t59 = Ifges(6,3) * t68;
t54 = t235 * t311 + t303 * t82;
t53 = -t235 * t303 + t311 * t82;
t51 = -mrSges(11,1) * t363 + mrSges(11,2) * t114;
t50 = mrSges(11,1) * t114 + mrSges(11,2) * t363;
t36 = t209 * t324 - t210 * t260 - t239 * t341 + t240 * t262;
t35 = -t209 * t260 - t210 * t324 - t239 * t262 - t240 * t341;
t1 = [(t410 / 0.2e1 - t482) * t79 + (-t196 * t266 - t211 * t452) * Ifges(4,2) + (-t263 * t196 + t197 * t266 - t211 * t456 + t212 * t452) * Ifges(4,4) + (mrSges(4,1) * t196 - mrSges(8,1) * t195 + mrSges(4,2) * t197 + mrSges(8,2) * t194) * t284 + (t7 * t444 + t6 * t448 + t350 * t476 + t354 * t480 + t352 * t479 + t150 * mrSges(5,2) + Ifges(5,1) * t67 - Ifges(5,4) * t68 + (-mrSges(5,3) - t355) * t149 + t358 * mrSges(6,3) + (-mrSges(6,3) * t486 + t213 * t356 + t349 * t466 + t351 * t468 + t353 * t467 + t38 * t445 + t39 * t448) * qJD(5)) * t206 + (-t148 * mrSges(5,3) - Ifges(5,4) * t67 + t150 * mrSges(5,1) + t59 / 0.2e1 + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t68 + t485) * t204 + (t325 + t336) * qJD(6) + (t103 * t138 + t158 * t88) * m(11) + (-t411 / 0.2e1 - t318) * t80 + (m(6) * (t379 * t41 - t380 * t40 - t358) + t311 * t9 + t303 * t10 + t83 * t379 - t84 * t380) * (pkin(9) * t204 - pkin(11) * t206 + t229) - t247 * t384 / 0.2e1 + (t303 * t83 + t311 * t84 + t338) * (pkin(9) * t80 - pkin(11) * t79 + t181) + qJD(2) ^ 2 * (-Ifges(3,5) * t306 - Ifges(3,6) * t314) / 0.2e1 + (t195 * t262 + t210 * t458) * Ifges(8,2) + (t193 * t261 + t208 * t460) * Ifges(9,2) + (t66 * t203 + t464 * t78) * Ifges(10,2) + (t31 * t126 + t35 * t469) * Ifges(11,1) + (t125 * t32 + t36 * t471) * Ifges(11,2) + (t125 * t31 + t32 * t126 + t35 * t471 + t36 * t469) * Ifges(11,4) + (t150 * t229 + t181 * t217) * m(5) - qJD(2) * t306 * t249 / 0.2e1 + t293 * (Ifges(10,5) * t77 + Ifges(10,6) * t78) / 0.2e1 + t296 * (Ifges(9,5) * t207 + Ifges(9,6) * t208) / 0.2e1 + t297 * (Ifges(8,5) * t209 + Ifges(8,6) * t210) / 0.2e1 + t292 * (Ifges(11,5) * t35 + Ifges(11,6) * t36) / 0.2e1 + t244 * (-mrSges(10,1) * t66 + mrSges(10,2) * t65) + t228 * (-mrSges(10,1) * t78 + mrSges(10,2) * t77) + t229 * (mrSges(5,1) * t68 + mrSges(5,2) * t67) + t210 * t160 / 0.2e1 - t211 * t161 / 0.2e1 + t212 * t164 / 0.2e1 + t207 * t162 / 0.2e1 + t208 * t159 / 0.2e1 + t209 * t163 / 0.2e1 + (t194 * t262 - t195 * t341 + t209 * t458 + t210 * t454) * Ifges(8,4) + (-t194 * t341 + t209 * t454) * Ifges(8,1) + (mrSges(4,1) * t211 - mrSges(8,1) * t210 + mrSges(4,2) * t212 + mrSges(8,2) * t209) * t273 + t377 * (Ifges(4,5) * t212 - Ifges(4,6) * t211) / 0.2e1 + (t483 * t314 + (-t209 * t309 + t210 * t301 + (t262 * t309 - t301 * t341) * qJD(7)) * mrSges(8,3) + (t211 * t313 - t212 * t305 + (-t263 * t313 + t266 * t305) * qJD(3)) * mrSges(4,3)) * t430 - pkin(15) * (-mrSges(9,1) * t193 + mrSges(9,2) * t192) + (-t115 * t35 - t116 * t36 + t125 * t73 - t126 * t74) * mrSges(11,3) + t181 * t100 + (-t343 * t65 + t463 * t77) * Ifges(10,1) + (t65 * t203 - t343 * t66 + t463 * t78 + t464 * t77) * Ifges(10,4) + (-t193 * (-mrSges(10,1) * t203 - mrSges(10,2) * t343) - t208 * t99 + (-t193 * t244 - t228 * t208) * m(10) + (-t299 * t78 + t307 * t77 + (-t203 * t307 + t299 * t343) * qJD(9)) * t426) * pkin(2) + (t192 * t261 - t193 * t342 + t207 * t460 + t208 * t455) * Ifges(9,4) + (-t192 * t342 + t207 * t455) * Ifges(9,1) + (-pkin(15) * (-mrSges(9,1) * t208 + mrSges(9,2) * t207) + ((0.2e1 * t436 + 0.3e1 / 0.2e1 * t425) * t306 + (-0.2e1 * t437 - 0.3e1 / 0.2e1 * t424 + (0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(3,1)) * t306 + (-mrSges(4,1) * t266 - mrSges(8,1) * t262 + mrSges(4,2) * t263 - mrSges(8,2) * t341 + t284 * t473) * pkin(1)) * t314) * qJD(2) + (0.2e1 * t337 + 0.3e1 / 0.2e1 * Ifges(7,1) * t390 - 0.3e1 / 0.2e1 * t370 + (-0.3e1 / 0.2e1 * t302 ^ 2 + 0.3e1 / 0.2e1 * t310 ^ 2) * Ifges(7,4)) * qJD(6)) * qJD(1) + t158 * (-mrSges(11,1) * t32 + mrSges(11,2) * t31) + t138 * (-mrSges(11,1) * t36 + mrSges(11,2) * t35) + t88 * (-mrSges(11,1) * t125 + mrSges(11,2) * t126) + t103 * t51 + t78 * t91 / 0.2e1 + t77 * t93 / 0.2e1 + t36 * t43 / 0.2e1 + t35 * t44 / 0.2e1 + (t263 * t197 + t212 * t456) * Ifges(4,1); t320 + (-m(11) * t147 - t50) * t138 + m(11) * (t115 * t86 - t116 * t85 + t139 * t74 + t140 * t73) + m(5) * (t148 * t242 + t149 * t241 - t190 * t219 + t191 * t218) + t402 * t191 + (-t241 * t67 - t242 * t68 - t394) * mrSges(5,3) + t316 - t76 * t338 + ((t249 / 0.2e1 - t400 / 0.2e1 + (-t436 - t425 / 0.2e1) * qJD(1)) * t306 + (t247 / 0.2e1 - t399 / 0.2e1 + (t437 + t424 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t306) * qJD(1) - t483 * pkin(1)) * t314) * qJD(1) + t317 + (-t198 - t201) * t273 + (-t139 * t31 + t140 * t32 - t396) * mrSges(11,3) + ((t221 * t309 - t223 * t301) * qJD(7) + (t224 * t305 + t393) * qJD(3) + (-t194 * t309 + t195 * t301) * mrSges(8,3) + (t196 * t313 - t197 * t305) * mrSges(4,3) + (-mrSges(4,2) * t383 + t301 * t427 + (-mrSges(8,1) * t301 - mrSges(8,2) * t309) * qJD(7)) * qJD(2)) * pkin(1) + (-t190 * t84 - t76 * t83 + (-qJD(5) * t83 - t9) * t233 + (-t3 - t397) * mrSges(6,3)) * t303 + (t233 * t10 + t190 * t83 - t76 * t84 + (-t40 * mrSges(6,3) - t233 * t84) * qJD(5)) * t311 + t232 * t8 - t220 * t100 + t190 * t143 - t147 * t51 + t85 * t104 + t86 * t105 + (-m(5) * t220 - t101) * t217 + ((t323 + t435) * t233 - t149 * t232 - t191 * t213 + t486 * t190) * m(6); t402 * t234 - m(6) * (t213 * t234 + t40 * t53 + t41 * t54) + t316 + (-t393 + (-mrSges(4,2) * qJD(3) - t224) * t305) * t430 - m(5) * (-t218 * t234 - t219 * t235) + t322 + (-t252 * t100 + (-t304 * t68 - t312 * t67) * mrSges(5,3) + (-t402 * t304 + (-t303 * t84 + t311 * t83 + t143) * t312 + m(6) * (t213 * t304 + t312 * t486)) * qJD(4) + (t148 * t304 + t149 * t312 + 0.2e1 * t217 * t457 + (-t218 * t304 - t219 * t312) * qJD(4)) * m(5)) * pkin(5) + t319 * (pkin(5) * t304 + pkin(11)) - mrSges(5,3) * t394 - t273 * t201 - t235 * t143 - t217 * t101 - t54 * t83 - t53 * t84 - t502 * (-pkin(5) * t312 - pkin(9)); (t498 + t482) * t362 + t318 * t344 + t319 * pkin(11) + t322 - t402 * t219 + t321 - m(6) * (-t213 * t219 + t40 * t69 + t41 * t70) - t218 * t143 - t70 * t83 - t69 * t84 + t502 * pkin(9); t59 - t213 * (mrSges(6,1) * t137 + mrSges(6,2) * t136) + (Ifges(6,1) * t136 - t414) * t467 + t137 * t38 / 0.2e1 + (Ifges(6,5) * t136 - Ifges(6,6) * t137) * t466 - t40 * t83 + t41 * t84 + (t136 * t40 + t137 * t41) * mrSges(6,3) + (-Ifges(6,2) * t137 + t135 + t39) * t468 + t485; (-t310 * t288 / 0.2e1 + (-t337 + (Ifges(7,1) * t310 - t419) * t449 + t370 / 0.2e1) * qJD(1) + t325 - t336) * qJD(1); t489 * t105 + t488 * t104 + (-t188 * t31 + t189 * t32 - t396) * mrSges(11,3) + t320 - t273 * t198 + ((-mrSges(8,2) * qJD(7) - t221) * t309 + (-mrSges(8,1) * qJD(7) + t223 + t427) * t301) * t430 - t154 * t51 - t138 * t50 + (t115 * t489 - t116 * t488 - t138 * t154 + t188 * t74 + t189 * t73) * m(11); t317; (t142 * t307 - t144 * t299 + (t299 * t345 - t490) * mrSges(10,3)) * t442 + t495; -t116 * t105 + (-t104 + t395) * t115 + t347 + t500 * t472 + (-t138 * mrSges(11,2) + Ifges(11,1) * t470 + Ifges(11,5) * t451) * t363 - (t138 * mrSges(11,1) + mrSges(11,3) * t116 + Ifges(11,4) * t470 + Ifges(11,2) * t472 + Ifges(11,6) * t451 + t477) * t114; 0; 0; 0;];
tauc = t1(:);
