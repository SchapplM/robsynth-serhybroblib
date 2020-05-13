% Calculate vector of cutting forces with Newton-Euler
% palh2m2OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = palh2m2OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 02:21:15
% EndTime: 2020-05-03 02:22:23
% DurationCPUTime: 31.11s
% Computational Cost: add. (12193->776), mult. (15135->1033), div. (0->0), fcn. (10586->12), ass. (0->422)
t322 = sin(qJ(2));
t328 = cos(qJ(2));
t330 = m(6) + m(7);
t306 = m(5) + t330;
t302 = m(4) + t306;
t588 = pkin(4) * t302;
t321 = sin(qJ(3));
t318 = mrSges(7,3) + mrSges(6,2);
t319 = sin(qJ(5));
t527 = t318 * t319;
t589 = sin(qJ(6));
t301 = mrSges(7,2) * t589;
t324 = cos(qJ(6));
t564 = t324 * mrSges(7,1);
t268 = -t301 + t564;
t574 = pkin(3) * m(7);
t243 = mrSges(6,1) + t268 + t574;
t325 = cos(qJ(5));
t541 = t243 * t325;
t213 = t541 - t527;
t334 = (pkin(5) * m(7));
t335 = (pkin(5) * m(6));
t427 = mrSges(5,1) + t334 + t335;
t184 = t427 + t213;
t320 = sin(qJ(4));
t637 = t320 * t184;
t525 = t319 * t243;
t526 = t318 * t325;
t216 = t526 + t525;
t209 = mrSges(5,2) + t216;
t326 = cos(qJ(4));
t638 = t209 * t326;
t128 = t637 + t638;
t617 = mrSges(4,2) + t128;
t114 = t617 * t321;
t327 = cos(qJ(3));
t483 = pkin(2) * t306;
t258 = mrSges(4,1) + t483;
t191 = t209 * t320;
t641 = t184 * t326 - t191;
t649 = t258 + t641;
t658 = t649 * t327 - t114;
t669 = t658 + t588;
t48 = mrSges(3,1) + t669;
t45 = t48 * t328;
t655 = t617 * t327;
t93 = t321 * t649;
t57 = t655 + t93;
t53 = mrSges(3,2) + t57;
t12 = -t53 * t322 + t45;
t618 = m(3) + t302;
t622 = pkin(1) * t618;
t673 = -mrSges(2,1) - t622 - t12;
t570 = t322 * t48;
t13 = t328 * t53 + t570;
t51 = t57 * t322;
t19 = t328 * t658 - t51;
t569 = t322 * t658;
t20 = t328 * t57 + t569;
t657 = t128 * t321;
t67 = t641 * t327 - t657;
t567 = t322 * t67;
t553 = t128 * t327;
t656 = t321 * t641;
t65 = t656 + t553;
t28 = t328 * t65 + t567;
t568 = t322 * t65;
t29 = t328 * t67 - t568;
t341 = qJD(5) ^ 2;
t342 = qJD(4) ^ 2;
t343 = qJD(3) ^ 2;
t344 = qJD(2) ^ 2;
t345 = qJD(1) ^ 2;
t203 = t216 * t320;
t142 = t213 * t326 - t203;
t204 = t216 * t326;
t628 = t320 * t213;
t140 = t204 + t628;
t551 = t140 * t321;
t73 = t142 * t327 - t551;
t565 = t322 * t73;
t550 = t140 * t327;
t642 = t321 * t142;
t72 = t642 + t550;
t36 = t328 * t72 + t565;
t566 = t322 * t72;
t37 = t328 * t73 - t566;
t672 = -t13 * qJDD(2) - t20 * qJDD(3) - t28 * qJDD(4) - t36 * qJDD(5) - t12 * t344 - t19 * t343 - t29 * t342 - t37 * t341 + t673 * t345;
t665 = mrSges(3,2) + t655;
t455 = 0.2e1 * t527;
t459 = 0.2e1 * t541;
t206 = -t459 + t455;
t163 = (2 * mrSges(5,1)) + (2 * t334) + (2 * t335) - t206;
t159 = t163 * t320;
t453 = 0.2e1 * t525;
t611 = t453 + 0.2e1 * t526;
t189 = 0.2e1 * mrSges(5,2) + t611;
t630 = t189 * t326;
t107 = -t159 - t630 - 0.2e1 * mrSges(4,2);
t644 = t107 * t327;
t46 = (t644 - 0.2e1 * t93) * t328;
t661 = (-t46 + 0.2e1 * t569) * qJD(3);
t109 = (t630 + 0.2e1 * t637) * t327;
t660 = (t109 + 0.2e1 * t656) * t328;
t586 = pkin(4) * t321;
t105 = t302 * t586 - t617;
t446 = t327 * t588;
t653 = -t446 - t649;
t652 = t446 - t649;
t160 = t163 * t326;
t177 = t189 * t320;
t632 = -t177 + 0.2e1 * t258 + t160;
t90 = t632 * t321;
t651 = (t644 - t90 - 0.2e1 * mrSges(3,2)) * t328 - 0.2e1 * t570;
t648 = (t660 + 0.2e1 * t567) * qJD(4);
t598 = -0.2e1 * t322;
t646 = t73 * t598;
t187 = t611 * t326;
t606 = (t187 + 0.2e1 * t628) * t327;
t645 = (0.2e1 * t642 + t606) * t328;
t640 = (t645 + 0.2e1 * t565) * qJD(5);
t283 = t320 * t483;
t169 = t283 - t209;
t481 = pkin(2) * t326;
t284 = t306 * t481;
t634 = t284 + t184;
t323 = sin(qJ(1));
t329 = cos(qJ(1));
t610 = g(1) * t329 + g(2) * t323;
t443 = t306 * t586;
t262 = t320 * t443;
t529 = t306 * t326;
t445 = pkin(4) * t529;
t626 = t327 * t445 - t262 + t634;
t625 = qJD(1) * t268;
t269 = mrSges(7,1) * t589 + mrSges(7,2) * t324;
t259 = -mrSges(6,3) + t269;
t255 = -mrSges(5,3) + t259;
t624 = -t255 * qJDD(1) - t610 * t306;
t188 = 0.4e1 * t209;
t395 = qJD(6) * t625;
t385 = -0.2e1 * t395;
t251 = -t319 * t320 + t325 * t326;
t250 = t326 * t319 + t320 * t325;
t539 = t250 * t327;
t172 = t251 * t321 + t539;
t236 = t251 * t327;
t386 = -t250 * t321 + t236;
t623 = t322 * t172 - t328 * t386;
t621 = -t259 * qJDD(1) - t610 * t330;
t308 = t326 ^ 2;
t620 = -0.4e1 * t308;
t192 = t209 * t308;
t615 = t209 - 0.2e1 * t192;
t307 = t325 ^ 2;
t542 = t243 * t307;
t225 = 0.2e1 * t542;
t614 = -t225 + t243;
t340 = qJD(6) ^ 2;
t487 = -t268 * qJDD(6) + t269 * t340;
t528 = t318 * t307;
t279 = 0.2e1 * t528;
t613 = t279 - t318;
t480 = pkin(2) * t330;
t402 = t320 * t480;
t612 = t453 - t402;
t464 = t269 * qJDD(6);
t535 = t268 * t340;
t366 = -t464 - t535;
t475 = qJD(6) * t269;
t400 = t623 * t475;
t609 = qJD(1) * (qJD(4) * t28 + qJD(5) * t36 - t400);
t277 = 0.4e1 * t528;
t523 = t319 * t330;
t442 = pkin(5) * t523;
t377 = -0.2e1 * t442 + (2 * mrSges(7,3)) + (2 * mrSges(6,2)) - t277;
t524 = t319 * t325;
t424 = t243 * t524;
t608 = (t377 - 0.4e1 * t424) * t308 + t612 * t325 + t279;
t245 = -mrSges(4,3) + t255;
t530 = t302 * t329;
t531 = t302 * t323;
t607 = -g(1) * t530 - g(2) * t531 - t245 * qJDD(1);
t425 = t318 * t524;
t591 = 0.8e1 * t307 - 0.4e1;
t604 = t591 * t243 - 0.8e1 * t425;
t101 = t172 * t328 + t322 * t386;
t559 = qJD(5) * t37;
t369 = -qJD(4) * t29 - t559;
t394 = qJD(5) * t475;
t382 = 0.2e1 * t394;
t594 = 0.2e1 * qJD(4);
t595 = 0.2e1 * qJD(3);
t596 = 2 * qJD(2);
t603 = t385 + (-qJD(3) * t19 + t369) * t596 + t369 * t595 - t559 * t594 + t101 * t382 + t672;
t602 = -0.2e1 * t184;
t460 = 0.2e1 * t191;
t601 = 0.4e1 * t319;
t600 = -0.2e1 * t320;
t599 = -0.2e1 * t321;
t597 = 0.2e1 * t326;
t593 = 0.2e1 * qJD(5);
t592 = 0.4e1 * t307 - 0.2e1;
t590 = t620 + 0.2e1;
t587 = pkin(4) * t306;
t585 = pkin(5) * t330;
t582 = g(3) * t302;
t581 = g(3) * t306;
t580 = g(3) * t328;
t579 = g(3) * t330;
t399 = t101 * t475;
t139 = t206 * t326 + 0.2e1 * t203;
t71 = t139 * t327 + 0.2e1 * t551;
t578 = (t328 * t71 + 0.2e1 * t566) * qJD(5) + 0.2e1 * t399;
t577 = 0.2e1 * t400 - t640;
t305 = qJD(2) + qJD(3);
t300 = qJD(4) + t305;
t290 = qJD(5) + t300;
t479 = qJD(1) * t290;
t575 = -qJDD(1) * t623 - t101 * t479;
t433 = -0.2e1 * t475;
t561 = (t139 * t321 - 0.2e1 * t550) * qJD(5) + t386 * t433;
t432 = 0.2e1 * t475;
t560 = t71 * qJD(5) + t172 * t432;
t309 = t327 ^ 2;
t558 = t649 * t309;
t310 = t328 ^ 2;
t557 = t649 * t310;
t112 = -t483 - t641;
t556 = t112 * t321;
t197 = -t213 - t585;
t133 = t197 * t326 + t203;
t124 = -t480 + t133;
t555 = t124 * t321;
t303 = -mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t157 = t528 + (t541 + t585 / 0.2e1) * t319 + t303;
t548 = t157 * t308;
t547 = (t191 - t483 / 0.2e1) * t326;
t546 = t184 * t308;
t545 = t197 * t320;
t520 = t320 * t243;
t452 = -0.4e1 * t520;
t196 = (t319 * t452 + t480) * t325;
t242 = -mrSges(3,3) + t245;
t241 = mrSges(2,2) + t242;
t544 = t241 * t323;
t543 = t241 * t329;
t538 = (t527 - t585 / 0.2e1) * t325;
t537 = t268 * t323;
t536 = t268 * t329;
t534 = t269 * t323;
t533 = t269 * t329;
t532 = t618 * t322;
t519 = t320 * t306;
t518 = t320 * t318;
t517 = t320 * t326;
t516 = t320 * t330;
t513 = t321 * t243;
t512 = t321 * t308;
t511 = t321 * t318;
t510 = t321 * t320;
t509 = t321 * t325;
t508 = t321 * t327;
t506 = t322 * t243;
t505 = t322 * t308;
t504 = t322 * t309;
t503 = t322 * t325;
t502 = t322 * t345;
t501 = t325 * t330;
t500 = t330 * (pkin(4) * t320 + pkin(5) * t599);
t499 = t330 * (pkin(5) * t600 + t586);
t451 = -0.4e1 * t510;
t498 = t330 * (pkin(5) * t451 + pkin(4));
t275 = pkin(4) * t510 - pkin(5);
t497 = t330 * t275;
t496 = t112 * t327 + t657;
t134 = t204 - t545;
t127 = t134 * t321;
t495 = t124 * t327 + t127;
t494 = (t206 * t320 - 0.2e1 * t204) * qJD(5) + t251 * t433;
t493 = t139 * qJD(5) + t250 * t432;
t456 = -0.4e1 * t528;
t491 = t320 * t456 + t196;
t441 = pkin(5) * t501;
t218 = t243 + t441;
t490 = -t218 * t594 - t243 * t593;
t489 = t206 * qJD(5) + t319 * t432;
t265 = t442 - t318;
t488 = t265 * t594 - t318 * t593;
t485 = pkin(1) * t302;
t484 = pkin(1) * t322;
t482 = pkin(2) * t322;
t478 = qJD(1) * t323;
t477 = qJD(1) * t329;
t470 = qJD(1) * qJD(2);
t469 = qJD(4) * qJD(5);
t468 = t242 * qJDD(1);
t463 = 0.4e1 * t308;
t462 = -0.4e1 * t548;
t458 = -0.2e1 * t538;
t457 = 0.4e1 * t538;
t450 = 0.2e1 * t506;
t449 = qJDD(2) + qJDD(3);
t118 = -t160 + t460;
t62 = t118 * t327 + 0.2e1 * t657;
t448 = (t328 * t62 + 0.2e1 * t568) * qJD(4) + t578;
t447 = t577 - t648;
t444 = pkin(4) * t519;
t401 = t325 * t480;
t431 = t320 * t401 - t318;
t430 = (t118 * t321 - 0.2e1 * t553) * qJD(4) + t561;
t429 = t62 * qJD(4) + t560;
t428 = -0.2e1 * t469;
t426 = t184 * t517;
t233 = t574 / 0.2e1 + t564 / 0.2e1 - t301 / 0.2e1 + mrSges(6,1) / 0.2e1;
t148 = t233 + t538 - t542;
t423 = t148 * t512;
t422 = t307 * t513;
t421 = t319 * t511;
t420 = t322 * t528;
t419 = t321 * t506;
t260 = pkin(1) * t321 - t482;
t418 = t260 * t523;
t266 = t482 * t599 + pkin(1);
t417 = t266 * t523;
t404 = t321 * t484;
t267 = -pkin(2) + t404;
t416 = t267 * t523;
t415 = t319 * t500;
t414 = t319 * t499;
t413 = t319 * t498;
t412 = t325 * t498;
t411 = (-t159 - 0.2e1 * t638) * qJD(4) + t494;
t410 = t118 * qJD(4) + t493;
t409 = -qJD(4) * t163 + t489;
t398 = t325 * t475;
t408 = -qJD(4) * t189 - qJD(5) * t611 - 0.2e1 * t398;
t407 = t591 * t318;
t406 = t592 * t243;
t405 = t327 * t485;
t403 = t319 * t480;
t396 = t322 * t470;
t393 = t480 / 0.4e1;
t392 = 0.2e1 * t508;
t391 = -0.2e1 * t426;
t390 = 0.4e1 * t424;
t389 = -0.4e1 * t422;
t388 = -0.8e1 * t321 * t505;
t85 = -0.2e1 * t400;
t387 = t487 + (t640 + t648 + t85 + t661) * qJD(1);
t252 = -t321 * t322 + t327 * t328;
t237 = t252 * qJD(1);
t249 = t321 * t328 + t322 * t327;
t238 = t249 * qJD(1);
t155 = t237 * t326 - t238 * t320;
t156 = t237 * t320 + t238 * t326;
t88 = t325 * t155 - t156 * t319;
t384 = 0.2e1 * t395;
t383 = -0.2e1 * t394;
t381 = qJDD(4) + t449;
t378 = 0.2e1 * t283 - t189;
t295 = m(2) + t618;
t89 = t155 * t319 + t156 * t325;
t78 = -t589 * t290 + t324 * t89;
t375 = t323 * g(1) - t329 * g(2);
t374 = t192 + t426;
t371 = -t189 * t308 - t169;
t370 = (pkin(5) + t481) * t330;
t247 = -t327 * t326 + t510;
t248 = t320 * t327 + t321 * t326;
t368 = t247 * t322 - t248 * t328;
t367 = t319 * t497 + t318;
t365 = t307 * t452 - t403 + 0.2e1 * t520;
t364 = -0.2e1 * mrSges(6,1) + 0.4e1 * t542 + (-0.4e1 * t527 + 0.2e1 * t585) * t325 - 0.2e1 * t564 + 0.2e1 * t301 - 0.2e1 * t574;
t77 = -t324 * t290 - t589 * t89;
t361 = (-t177 + t483) * t326 + t163 * t308 - t184;
t360 = t265 + 0.2e1 * t424 + t279;
t246 = -pkin(1) * t345 - t610;
t358 = -qJDD(1) * pkin(1) - t375;
t263 = t319 * t402;
t356 = t243 - t263;
t353 = t364 * t308 - t225 + t356;
t2 = t447 - t661;
t352 = t2 * qJD(2) + t447 * qJD(3) + t577 * qJD(4) + t12 * qJDD(2) + t19 * qJDD(3) + t29 * qJDD(4) + t37 * qJDD(5) + t366 * t101 - t13 * t344 - t20 * t343 - t28 * t342 - t341 * t36 - t383 * t623;
t288 = m(1) + t295;
t287 = qJDD(5) + t381;
t285 = t303 * t319;
t261 = -pkin(2) * t321 + t484;
t253 = qJDD(1) * t328 - t396;
t235 = (t455 - t585) * t325;
t220 = t233 * t319;
t186 = t319 * t370 + t431;
t185 = -g(3) * t322 + t246 * t328 + (-t310 * t345 - t344) * pkin(4);
t175 = -t580 - t246 * t322 + (t328 * t502 + qJDD(2)) * pkin(4);
t170 = -t247 * t328 - t248 * t322;
t154 = t325 * t370 + t356;
t146 = t248 * t587 + t169;
t123 = (t187 - 0.2e1 * t545) * t327;
t122 = (pkin(4) * t539 + (pkin(2) * t319 + pkin(4) * t509) * t326 - t275 * t319) * t330 + t431;
t97 = (pkin(4) * t236 + (pkin(2) * t325 - t319 * t586) * t326 - t275 * t325) * t330 + t356;
t84 = qJD(6) + t88;
t82 = t329 * t399;
t81 = t323 * t399;
t79 = t175 * t321 + t185 * t327 + (-t237 ^ 2 - t305 ^ 2) * pkin(2);
t76 = t327 * t175 - t321 * t185 + (t237 * t238 + t449) * pkin(2);
t64 = -pkin(4) * t330 + t495;
t63 = (t123 - 0.2e1 * t555) * t328;
t59 = t496 - t587;
t52 = (t109 - 0.2e1 * t556) * t328;
t49 = -t327 * t632 + 0.2e1 * t114;
t43 = t101 * qJDD(1) - t623 * t479;
t42 = qJDD(6) + t575;
t41 = mrSges(7,1) * t84 - mrSges(7,3) * t78;
t40 = -mrSges(7,2) * t84 + mrSges(7,3) * t77;
t39 = -mrSges(7,1) * t77 + mrSges(7,2) * t78;
t31 = t320 * t76 + t326 * t79 + (-t155 ^ 2 - t300 ^ 2) * pkin(5);
t27 = t77 * qJD(6) - t589 * t287 + t324 * t43;
t26 = -t78 * qJD(6) - t324 * t287 - t589 * t43;
t18 = -t320 * t79 + t326 * t76 + (t155 * t156 + t381) * pkin(5);
t6 = (-t170 * qJDD(1) + (-t368 * qJD(1) + t156) * t300) * pkin(5) + (-t253 + t396) * pkin(4) + (t290 * t89 - t575) * pkin(3) + (t321 * (qJDD(1) * t322 + t328 * t470) - t327 * t253 + (qJD(3) + t305) * t238) * pkin(2) + t358;
t3 = t18 * t319 + t31 * t325 + (-t290 ^ 2 - t88 ^ 2) * pkin(3);
t1 = [-t288 * g(1) + t82 * t594 + (t20 * t478 + t82) * t595 + (t13 * t478 + t82) * t596 + t544 * t345 + (t323 * t673 - t543) * qJDD(1) + (t533 * t623 - t537) * qJDD(6) + 0.2e1 * t323 * t609 + (t536 * t623 + t534) * t340 + t603 * t329, -t241 * qJDD(1) + t578 * qJD(4) + ((t328 * t49 + 0.2e1 * t51) * qJD(3) + t448) * qJD(2) + t448 * qJD(3) - t610 * t295 - t366 * t623 + 0.2e1 * (qJD(5) * t101 * t269 - t625) * qJD(6) + t672, -t67 * t342 - g(3) * t532 - t53 * qJDD(2) - t48 * t344 - t65 * qJDD(4) - t658 * t343 - t57 * qJDD(3) - t73 * t341 - t72 * qJDD(5) + t560 * qJD(4) + (qJD(3) * t49 + t429) * qJD(2) + t429 * qJD(3) + t172 * t382 + t366 * t386 + (t246 * t618 - t45 * t345 + t53 * t502 + t385 - t468) * t328, -t249 * t582 + t250 * t382 + t105 * qJDD(2) + t653 * t344 - t649 * t343 - t617 * qJDD(3) - t142 * t341 - t140 * qJDD(5) - t641 * t342 - t128 * qJDD(4) + ((-t309 * t632 + t392 * t617 - t652) * t310 + (-t405 + (0.2e1 * t309 * t617 + t392 * t649 + t105) * t322) * t328 + t558 - t617 * t508 + t302 * t404 - t649) * t345 + (-qJD(3) * t632 + t410) * qJD(2) + t410 * qJD(3) + t493 * qJD(4) + t366 * t251 + (t385 + t607) * t252, t368 * t581 - t209 * qJDD(4) - t213 * t341 - t325 * t464 - t325 * t535 + t319 * t382 + t409 * qJD(3) + ((-0.2e1 * t284 - t163) * qJD(3) + t409) * qJD(2) + t489 * qJD(4) + (((-t427 * t463 + (t188 * t320 - 0.2e1 * t483) * t326 + t213 * t620 + t163) * t309 + (t512 * t188 + (-t184 * t451 - t587) * t326 - t169 * t599) * t327 + t262 + t361) * t310 + t361 * t309 + (t261 * t519 + (t391 + t615) * t321) * t327 + t326 * (t306 * t267 - t641) + (0.4e1 * (t283 / 0.2e1 + t303 * t325 - t220 - mrSges(5,2) / 0.2e1 + t374) * t504 - t266 * t529 * t327 + t260 * t519 + ((t444 + (-t184 * t590 - 0.4e1 * t209 * t517) * t321) * t327 + (-t637 + t443 / 0.2e1) * t597 + t615) * t322) * t328) * t345 - t626 * t344 - t184 * t342 - t634 * t343 + t169 * qJDD(3) + t146 * qJDD(2) - t216 * qJDD(5) + (t385 + t624) * t170, -t154 * t343 - t97 * t344 - t218 * t342 - t243 * t341 + ((((t243 + 0.2e1 * t425 - t441 - 0.2e1 * t542) * t463 + (-0.2e1 * t401 + (0.4e1 * t265 + 0.8e1 * t424 + 0.8e1 * t528) * t320) * t326 + 0.2e1 * t263 + t364) * t309 + (0.8e1 * t157 * t512 + (-t412 + (t320 * t604 + 0.2e1 * t403) * t321) * t326 + t321 * t456 - 0.2e1 * t612 * t509 + t415 + 0.2e1 * t511) * t327 + (t414 + t491 + 0.2e1 * t518) * t326 + (t455 + t497) * t325 + t353) * t310 + (0.8e1 * (t548 + (t319 * t393 + (-t538 + (t307 - 0.1e1 / 0.2e1) * t243) * t320) * t326 + t303 * t307 + (t320 * t393 - t220) * t325 - t442 / 0.4e1 + mrSges(7,3) / 0.4e1 + mrSges(6,2) / 0.4e1) * t504 + (t148 * t388 - t266 * t501 * t326 + 0.4e1 * (t421 + t500 / 0.4e1) * t503 + t320 * t417 + 0.2e1 * t419) * t327 + t418 * t326 + 0.2e1 * t420 + (t260 * t516 + t319 * t450) * t325 + (((t413 + (-t407 - 0.8e1 * t424) * t510) * t326 + t389) * t327 + t462 + (t325 * t499 + (-t406 + 0.4e1 * t425) * t320) * t326 - t367) * t322) * t328 + ((t265 * t600 + t491) * t326 + t235 + t353) * t309 + (t250 * t330 * t261 + (t462 + (-t406 + t457) * t517 + t360) * t321) * t327 + (t235 + t614) * t308 + (t267 * t501 + t360 * t320) * t326 + t542 - t425 - t320 * t416 - t243) * t345 - t101 * t579 + t122 * qJDD(2) + t186 * qJDD(3) + t490 * qJD(3) + (-0.2e1 * qJD(3) * t154 + t490) * qJD(2) + t243 * t428 + t265 * qJDD(4) - t318 * qJDD(5) + t366 - (t385 + t621) * t623, m(7) * (t324 * t3 - t589 * t6) + t26 * mrSges(7,3) - t42 * mrSges(7,2) + t77 * t39 - t84 * t41; -t288 * g(2) + t81 * t594 + (-t13 * t477 + t81) * t596 + (-t20 * t477 + t81) * t595 - t543 * t345 + (-t329 * t673 - t544) * qJDD(1) + (t534 * t623 + t536) * qJDD(6) - 0.2e1 * t329 * t609 + (t537 * t623 - t533) * t340 + t603 * t323, -t673 * qJDD(1) - t241 * t345 + t375 * t295 + (t651 * qJD(2) + t2) * qJD(1) - t487, t658 * qJDD(3) + ((-t649 + t557) * t321 + t665 * t310 - t665) * t345 + t73 * qJDD(5) - t618 * t580 + ((-t90 - 0.2e1 * t655) * qJD(3) + t430) * qJD(2) + t561 * qJD(4) + t430 * qJD(3) + t48 * qJDD(2) - t53 * t344 - t65 * t342 + t67 * qJDD(4) - t57 * t343 - t72 * t341 + t386 * t383 + t610 * t532 + ((t622 + t45) * t345 + t468 + t384) * t322 + t366 * t172, -t252 * t582 - t128 * t342 - t140 * t341 + t142 * qJDD(5) + t411 * qJD(3) + t494 * qJD(4) + (qJD(3) * t107 + t411) * qJD(2) + t251 * t383 - t617 * t343 + t649 * qJDD(3) + t641 * qJDD(4) + ((-t107 * t309 + t105) * t310 - t655 * t327 + ((0.2e1 * t558 + t652) * t328 + t405) * t322 + (t328 * t485 + (t328 * t598 * t617 + 0.2e1 * t557 - t649) * t327) * t321) * t345 - t653 * qJDD(2) + t105 * t344 + t366 * t250 + (t384 - t607) * t249, -t611 * t469 + t146 * t344 + t169 * t343 + t634 * qJDD(3) - t170 * t581 + 0.2e1 * (-qJD(4) - qJD(5)) * t398 + t408 * qJD(3) + (t378 * qJD(3) + t408) * qJD(2) + (((t188 * t308 + t378 + 0.4e1 * t426) * t309 + t327 * t444 - t159 * t326 + ((t602 + 0.4e1 * t546 - 0.4e1 * t547) * t327 + t445) * t321 + t371) * t310 + (-0.4e1 * (-t546 + t547 + t285 + t233 * t325 + t335 / 0.2e1 + t334 / 0.2e1 + mrSges(5,1) / 0.2e1) * t504 + (t266 * t519 + (t445 + (t590 * t209 - 0.4e1 * t426) * t321) * t322) * t327 + t505 * t602 + (t306 * t260 + t322 * t460) * t326 - t322 * (t262 - t184)) * t328 + (t371 + t391) * t309 + ((t261 * t306 + t321 * t460) * t326 - (0.2e1 * t308 - 0.1e1) * t321 * t184) * t327 - t267 * t519 + t374 - t209) * t345 + t626 * qJDD(2) + t184 * qJDD(4) - t216 * t341 - t209 * t342 + t213 * qJDD(5) + t366 * t319 - (t384 - t624) * t368, t243 * qJDD(5) + t218 * qJDD(4) + t154 * qJDD(3) + t623 * t579 + t97 * qJDD(2) + t265 * t342 + t122 * t344 - t318 * t341 + t488 * qJD(3) + (t186 * t595 + t488) * qJD(2) + t318 * t428 + t186 * t343 + (((t613 * t463 + ((0.8e1 * t541 + 0.4e1 * t585) * t308 + t480 * t597 - 0.4e1 * t541) * t319 + ((-0.4e1 * t243 - 0.8e1 * t538 + 0.8e1 * t542) * t326 + 0.2e1 * t401) * t320 + t377) * t309 + (-0.8e1 * t423 + (t413 + (-t320 * t407 + 0.2e1 * t196) * t321) * t326 + t389 + (0.4e1 * t421 + t500) * t325 + 0.2e1 * t321 * t356) * t327 + ((t518 * t601 + t499) * t325 + t365) * t326 - t367 + t608) * t310 + (-0.8e1 * (t148 * t308 + (-t401 / 0.4e1 + (t528 + t424 + t265 / 0.2e1) * t320) * t326 + t233 * t307 + (t285 + t585 / 0.4e1) * t325 + t263 / 0.4e1 - t574 / 0.4e1 - t564 / 0.4e1 + t301 / 0.4e1 - mrSges(6,1) / 0.4e1) * t504 + (t157 * t388 + (t417 + (-t510 * t604 + t412) * t322) * t326 + 0.4e1 * t321 * t420 + (t266 * t516 + t419 * t601) * t325 + (t415 / 0.2e1 + t511) * t598) * t327 + 0.4e1 * t148 * t505 + (t260 * t501 + (-t414 + (t592 * t318 + t390) * t320) * t322) * t326 + t307 * t450 - 0.2e1 * (t527 + t497 / 0.2e1) * t503 - t320 * t418 - t506) * t328 + ((t320 * t457 + t365) * t326 + t265 + t608) * t309 + (0.4e1 * t423 + (t261 * t501 + (0.2e1 * t265 + t390 + t277) * t510) * t326 + 0.2e1 * t422 + t321 * t458 - t261 * t319 * t516 - t513) * t327 + ((t459 + t585) * t319 + t613) * t308 + (-t416 + (t458 - t614) * t320) * t326 - t325 * (t267 * t516 + t216)) * t345 + t290 * t433 + (t384 - t621) * t101, m(7) * (-t589 * t3 - t324 * t6) - t27 * mrSges(7,3) + t42 * mrSges(7,1) - t78 * t39 + t84 * t40; -g(3) * t288 + t352, -g(3) * t295 + t352, -qJDD(1) * t12 + t242 * t345 + t358 * t618 - t651 * t470 + t387, (-t669 * t328 - t485 + t51) * qJDD(1) - g(1) * t531 + g(2) * t530 + (0.2e1 * t322 * t669 - t46) * t470 + t245 * t345 + t387, (t59 * t328 + (-t556 + t553) * t322) * qJDD(1) + t255 * t345 + t358 * t306 + ((t59 * t598 + t52) * qJD(2) + (t496 * t598 + t52) * qJD(3) + (-t67 * t598 + t660) * qJD(4) + ((-t142 * t599 + t606) * t328 - t646) * qJD(5) + t85) * qJD(1) + t487, (t64 * t328 + (t134 * t327 - t555) * t322) * qJDD(1) + t259 * t345 + ((t64 * t598 + t63) * qJD(2) + (t495 * t598 + t63) * qJD(3) + ((t133 * t599 + t123) * t328 + (t133 * t327 + t127) * t598) * qJD(4) + (t645 - t646) * qJD(5) + t85) * qJD(1) + t358 * t330 + t487, m(7) * (t18 * t325 - t31 * t319 + (t88 * t89 + t287) * pkin(3)) + t27 * mrSges(7,2) - t26 * mrSges(7,1) + t78 * t41 - t77 * t40;];
f_new = t1;
