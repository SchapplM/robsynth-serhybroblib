% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% palh2m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1IC_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:57
% EndTime: 2020-05-03 01:04:19
% DurationCPUTime: 22.35s
% Computational Cost: add. (5832->809), mult. (9340->936), div. (0->0), fcn. (3788->10), ass. (0->435)
t341 = sin(qJ(3));
t346 = cos(qJ(3));
t339 = sin(qJ(5));
t602 = mrSges(6,1) * t339;
t283 = pkin(4) * t602;
t344 = cos(qJ(5));
t326 = t344 ^ 2;
t291 = Ifges(6,4) * t326;
t338 = Ifges(6,1) - Ifges(6,2);
t627 = mrSges(6,2) * pkin(4);
t667 = (-t338 * t339 + t627) * t344;
t133 = t667 + Ifges(6,4) + t283 - 0.2e1 * t291;
t129 = -Ifges(5,5) + t133;
t340 = sin(qJ(4));
t345 = cos(qJ(4));
t290 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t253 = t290 * t339;
t626 = mrSges(6,2) * pkin(6);
t289 = -Ifges(6,6) + t626;
t173 = -t289 * t344 - t253;
t701 = Ifges(5,6) - t173;
t51 = -t129 * t340 + t701 * t345;
t47 = Ifges(4,6) + t51;
t600 = mrSges(6,2) * t344;
t284 = pkin(3) * t600;
t285 = pkin(3) * t602;
t350 = mrSges(5,3) * pkin(3);
t139 = t701 * t340;
t50 = t129 * t345 + t139;
t715 = t284 + t285 + t350 - Ifges(4,5) + t50;
t719 = -t341 * t715 + t47 * t346;
t718 = Ifges(3,6) + t719;
t716 = -t341 * t50 + t346 * t51;
t292 = mrSges(6,2) * t339;
t228 = mrSges(6,1) * t344 - t292;
t629 = pkin(4) * m(6);
t193 = mrSges(5,1) + t228 + t629;
t545 = t340 * t193;
t153 = t545 + mrSges(4,2);
t293 = m(6) * pkin(6) + mrSges(6,3);
t282 = -mrSges(5,2) + t293;
t555 = t282 * t345;
t463 = -t153 + t555;
t561 = t193 * t345;
t557 = t282 * t340;
t351 = m(5) + m(6);
t672 = -pkin(3) * t351 - mrSges(4,1);
t695 = t557 - t672;
t705 = t695 + t561;
t711 = t346 * t705;
t714 = t463 * t341 + t711;
t552 = t290 * t344;
t554 = t289 * t339;
t174 = t552 - t554;
t681 = -Ifges(5,4) - t174;
t481 = 0.8e1 * t681;
t687 = -mrSges(6,3) * pkin(4) - pkin(6) * t629;
t120 = t481 + 0.8e1 * t687;
t106 = t120 * t340;
t686 = -0.2e1 * t552 + 0.2e1 * t554;
t477 = -0.2e1 * Ifges(5,4) + t686;
t122 = t477 + 0.2e1 * t687;
t110 = t122 * t340;
t712 = t341 * t705;
t327 = t345 ^ 2;
t134 = -t687 - t681;
t165 = pkin(3) * t545;
t671 = -Ifges(4,4) + t165;
t73 = t134 + t671;
t709 = t122 * t327 + t73;
t708 = pkin(2) * t695;
t608 = t695 * pkin(1);
t541 = t341 * t345;
t202 = t340 * t346 + t541;
t615 = pkin(2) * t341;
t166 = t193 * t615;
t265 = pkin(3) * t282;
t132 = t166 - t265;
t188 = pkin(3) * t193;
t230 = t282 * t615;
t531 = t230 + t188;
t651 = t555 - t545;
t572 = t651 * t346;
t707 = -pkin(2) * t572 + t132 * t345 + t340 * t531;
t517 = pkin(4) * t292;
t274 = 0.2e1 * t517;
t625 = mrSges(6,3) * pkin(6);
t315 = 0.2e1 * t625;
t604 = Ifges(5,1) + Ifges(6,2);
t673 = Ifges(5,2) + Ifges(6,3);
t468 = -t604 + t673;
t599 = Ifges(6,4) * t339;
t628 = mrSges(6,1) * pkin(4);
t206 = (t599 + t628) * t344;
t547 = t338 * t326;
t652 = t547 - 0.2e1 * t206;
t706 = t274 + t315 + t652 - t468;
t532 = -t561 - t557;
t700 = t346 * t532;
t704 = t651 * t341 - t700;
t111 = t463 * t346;
t703 = t111 - t712;
t399 = pkin(6) ^ 2;
t366 = m(6) * t399;
t400 = pkin(4) ^ 2;
t367 = m(6) * t400;
t423 = 0.2e1 * t706;
t81 = 0.2e1 * t366 - 0.2e1 * t367 + t423;
t71 = t81 * t340;
t568 = t134 * t327;
t125 = t134 * t340;
t566 = t153 * t341;
t699 = t468 - 0.2e1 * t625;
t279 = t293 * pkin(4);
t698 = -t477 + 0.2e1 * t279;
t479 = 0.4e1 * t681;
t697 = t479 - 0.4e1 * t279;
t443 = t479 + 0.4e1 * t687;
t696 = t481 - 0.8e1 * t279;
t401 = pkin(3) ^ 2;
t662 = -t366 / 0.2e1 + t367 / 0.2e1;
t694 = t662 - (m(5) / 0.2e1 + m(6) / 0.2e1) * t401;
t692 = 0.2e1 * t345;
t691 = (Ifges(4,2) - Ifges(4,1));
t590 = t47 * t341;
t18 = t346 * t715 + t590;
t259 = -0.4e1 * t265;
t677 = -0.2e1 * t517;
t424 = 0.4e1 * t677 - 0.4e1 * t652 + 0.4e1 * t699;
t80 = -0.4e1 * t366 + 0.4e1 * t367 + t424;
t70 = t80 * t340;
t578 = t70 + t259;
t690 = t345 * t578;
t229 = pkin(3) * t557;
t689 = 0.4e1 * t229 + (2 * t691);
t688 = -0.2e1 * t229 - t691;
t343 = sin(qJ(1));
t348 = cos(qJ(1));
t649 = g(1) * t348 + g(2) * t343;
t685 = Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 - t229;
t658 = t367 - t366;
t405 = t340 * (-t658 + t706);
t678 = 0.8e1 * t405;
t683 = t678 + 0.8e1 * t265;
t312 = qJD(2) + qJD(3);
t682 = (qJD(4) * t312);
t680 = 0.2e1 * pkin(2);
t614 = pkin(3) * t651;
t676 = 0.2e1 * t614;
t634 = 0.2e1 * t265;
t577 = t71 + t634;
t570 = t531 * t345;
t156 = 0.4e1 * t165;
t665 = t156 - 0.4e1 * Ifges(4,4);
t664 = t81 + t689;
t515 = 0.2e1 * t206;
t482 = t515 + t677 - t547;
t427 = t482 + t699;
t87 = t427 + t658;
t663 = t87 + t688;
t352 = m(4) + m(5);
t490 = mrSges(3,1) + (m(6) + t352) * pkin(2);
t369 = m(5) * t401;
t635 = 0.2e1 * t229;
t657 = t369 + t635;
t353 = Ifges(6,3) / 0.2e1;
t356 = -Ifges(6,2) / 0.2e1;
t295 = Ifges(6,1) / 0.2e1 + t356;
t639 = t206 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - t295 * t326 - t517 - t625;
t413 = t353 + t356 + t639;
t655 = t413 + t685;
t460 = t547 + t315 + t604;
t428 = t274 - t515 + t460 - t673;
t620 = m(6) * (-t399 + t400);
t84 = t428 - t620;
t588 = t84 * t327;
t654 = -t369 + t427 + 0.2e1 * t588 + t688;
t227 = t600 + t602;
t560 = t227 * t340;
t112 = pkin(3) * t560 - t173;
t288 = qJD(4) + t312;
t442 = qJD(5) * t288;
t402 = pkin(2) ^ 2;
t303 = t352 * t402;
t648 = Ifges(3,2) + t303 - Ifges(3,1);
t646 = Ifges(4,2) + t657;
t644 = Ifges(6,4) / 0.4e1 - t291 / 0.2e1;
t258 = -0.2e1 * t265;
t643 = (t70 + t258) * t345;
t342 = sin(qJ(2));
t542 = t341 * t342;
t494 = pkin(1) * t542;
t241 = -pkin(3) + t494;
t642 = -t698 * t327 + (t241 * t282 + t87 * t340) * t345;
t563 = t695 * t341;
t142 = mrSges(3,2) + t563;
t113 = t490 - t566;
t26 = t282 * t541 + t113 + t711;
t313 = m(3) + t352;
t347 = cos(qJ(2));
t641 = t342 * (t193 * t541 - t111 + t142) - t26 * t347 - mrSges(2,1) - (m(6) + t313) * pkin(1);
t302 = t351 * t401;
t636 = 0.2e1 * t188;
t440 = Ifges(6,1) + Ifges(5,3) + t315 + t482;
t88 = t366 + t367 + t440;
t640 = t636 * t345 + Ifges(4,3) + t302 + t635 + t88;
t638 = -0.2e1 * pkin(1);
t385 = -0.2e1 * Ifges(6,4);
t158 = -0.2e1 * t165;
t633 = 0.2e1 * t340;
t632 = 0.2e1 * t341;
t631 = 0.8e1 * t341;
t630 = 0.2e1 * t342;
t105 = t129 * t692;
t546 = t339 * t344;
t246 = t546 * t385;
t355 = Ifges(6,2) / 0.2e1;
t544 = t340 * t341;
t201 = -t345 * t346 + t544;
t90 = t201 * t342 - t202 * t347;
t42 = 0.2e1 * t90 * (t547 + t246 - Ifges(6,1) / 0.2e1 + t355 + t353) * qJD(5);
t594 = t342 * t716;
t595 = t341 * t51;
t623 = t42 + (((t105 + 0.2e1 * t139) * t346 + 0.2e1 * t595) * t347 + 0.2e1 * t594) * qJD(4);
t406 = -t369 / 0.2e1 + (t110 - t188) * t345 + t655;
t328 = t346 ^ 2;
t548 = t328 * t342;
t505 = 0.8e1 * t548;
t605 = t401 * m(6);
t549 = t327 * t342;
t507 = -0.4e1 * t549;
t86 = t428 - t658;
t64 = t86 * t507;
t75 = t86 * t327;
t622 = t64 + (-t605 / 0.2e1 + t75 + t406 + t662) * t505;
t528 = t399 + t401;
t280 = -t400 + t528;
t621 = m(6) * t280;
t619 = pkin(1) * t227;
t267 = pkin(1) * t282;
t618 = pkin(1) * t342;
t616 = pkin(2) * t193;
t266 = pkin(2) * t282;
t613 = pkin(3) * t227;
t612 = pkin(3) * t340;
t609 = t153 * pkin(2);
t606 = (-t401 + t402) * m(6);
t550 = t327 * t341;
t492 = t86 * t550;
t65 = 0.8e1 * t492;
t603 = (0.4e1 * Ifges(4,1) - 0.4e1 * Ifges(4,2) - 0.8e1 * t229 - 0.4e1 * t369 + t80 - 0.4e1 * t605) * t341 + t65;
t598 = t719 * t342;
t593 = t345 * (t340 * t86 + t265);
t182 = -0.2e1 * t188;
t62 = t443 * t340 + t182;
t592 = t345 * t62;
t74 = t413 + t662;
t589 = t74 * t340;
t586 = 0.8e1 * t671;
t149 = t153 * pkin(1);
t483 = t134 * t549;
t91 = -0.16e2 * t341 * t483;
t583 = -0.2e1 * t149 + t91;
t495 = t342 * t608;
t66 = -0.4e1 * t492;
t582 = -0.2e1 * t495 + t66;
t180 = 0.4e1 * t188;
t61 = -t106 + t180;
t580 = t61 * t341 - 0.2e1 * t266;
t579 = (t423 - 0.2e1 * t620) * t340 + t634;
t118 = -t547 + (0.2e1 * t599 + t628) * t344 - t517 + t353 + t295;
t575 = t118 * t340;
t124 = mrSges(4,2) - t651;
t574 = t124 * t341;
t569 = t133 * t345;
t135 = t279 - t681;
t567 = t135 * t340;
t565 = t153 * t342;
t564 = t173 * t327;
t562 = t193 * t342;
t559 = t227 * t341;
t558 = t227 * t345;
t556 = t282 * t342;
t540 = t342 * t345;
t539 = t342 * t346;
t538 = -t106 - 0.2e1 * t230;
t537 = pkin(2) * t228 + t174 * t541;
t151 = t174 * t340;
t535 = pkin(3) * t228 + t151;
t161 = -0.2e1 * t166;
t534 = t161 + t258;
t167 = pkin(2) * t562;
t533 = t167 - t267;
t530 = (t626 / 0.4e1 - Ifges(6,6) / 0.4e1) * t344 + t253 / 0.4e1;
t529 = (-t626 / 0.2e1 + Ifges(6,6) / 0.2e1) * t344 - t253 / 0.2e1;
t527 = qJD(1) * t288;
t89 = t201 * t347 + t202 * t342;
t526 = qJDD(1) * t89;
t524 = qJD(1) * qJD(2);
t523 = qJD(1) * qJD(5);
t522 = qJD(2) * qJD(3);
t263 = 0.2e1 * t267;
t521 = t341 * t638;
t520 = -0.2e1 * t618;
t516 = -0.4e1 * t567;
t506 = 0.4e1 * t548;
t503 = qJDD(2) + qJDD(3);
t502 = pkin(1) * t556;
t501 = pkin(2) * t565;
t500 = pkin(2) * t563;
t499 = pkin(2) * t556;
t498 = pkin(2) * t545;
t497 = pkin(2) * t557;
t496 = pkin(2) * t566;
t493 = t74 * t550;
t491 = t612 / 0.4e1;
t470 = -0.2e1 * t502;
t487 = t341 * t470 + t577;
t486 = t193 * t544;
t485 = t327 * t542;
t484 = t340 * t569;
t475 = t400 + t528;
t474 = t564 / 0.2e1;
t473 = -t560 / 0.4e1;
t472 = t71 + t265;
t471 = t342 * t524;
t467 = -t73 + 0.2e1 * t568;
t464 = -0.2e1 * t282 * t544;
t462 = t153 * t494;
t461 = t193 * t494;
t294 = Ifges(6,1) / 0.4e1 - Ifges(6,2) / 0.4e1;
t455 = (-t294 * t339 + t627 / 0.4e1) * t344 + t283 / 0.4e1 + t644;
t126 = -0.2e1 * t283 + 0.4e1 * t291 + t385 - 0.2e1 * t667;
t454 = t291 + (t295 * t339 - t627 / 0.2e1) * t344 - Ifges(6,4) / 0.2e1 - t283 / 0.2e1;
t453 = t341 * t473;
t203 = pkin(4) * t228 + Ifges(6,3);
t451 = t203 * t345 + t535;
t449 = t327 * t696 - t697;
t448 = -t327 * t697 - t698;
t447 = 0.8e1 * t135 * t485 + t345 * t533;
t446 = -0.2e1 * pkin(3) * t558 + t126;
t177 = (-t346 * t347 + t542) * qJD(1);
t178 = (-t341 * t347 - t539) * qJD(1);
t82 = t177 * t345 - t178 * t340;
t439 = t135 + t671;
t438 = t454 * t327 + (t529 * t340 - t613 / 0.4e1) * t345 + t455;
t437 = t120 + 0.16e2 * t568;
t436 = t120 * t327 - t443;
t433 = mrSges(4,3) + mrSges(5,3) + t227;
t392 = 0.2e1 * Ifges(4,4);
t432 = t158 + t392 + t448;
t157 = 0.2e1 * t165;
t431 = t157 + t436;
t430 = -t327 * t443 + t122 + t158;
t426 = t436 + t665;
t425 = t392 + t430;
t421 = Ifges(4,3) + t440 + t657;
t420 = (t345 * t683 + t437 - t586) * t328 + t426;
t419 = -pkin(2) * t202 - t612;
t418 = 0.2e1 * t500 - (2 * Ifges(3,4)) + t425;
t337 = 0.2e1 * t369;
t410 = t337 + t664;
t409 = t81 * t327 - t369 + t663;
t408 = t410 + 0.2e1 * t605;
t189 = pkin(1) * t193;
t119 = (t189 + t499) * t341;
t176 = t188 / 0.4e1;
t232 = -t341 * pkin(3) + t618;
t329 = t347 ^ 2;
t37 = -pkin(3) * t532 + t88;
t370 = qJD(5) ^ 2;
t372 = qJD(3) ^ 2;
t373 = qJD(2) ^ 2;
t374 = qJD(1) ^ 2;
t38 = -t341 * t532 - t572;
t40 = (t345 * t686 + 0.2e1 * t575) * t346;
t152 = t174 * t345;
t41 = t341 * (t152 - t575);
t98 = t118 * t345;
t46 = t98 + t151;
t58 = 0.8e1 * t493;
t59 = -0.4e1 * t493;
t60 = t74 * t507;
t21 = t346 * t50 + t595;
t8 = t342 * t21 - t716 * t347;
t404 = -t126 * t442 - t8 * qJDD(1) - (-pkin(2) * t700 - t132 * t340 + t570 + t88) * qJDD(2) - t37 * qJDD(3) - t88 * qJDD(4) - t173 * qJDD(5) - (((t643 + t157 + t449) * t328 + (t498 - t266 * t345 + t58 + ((t636 + 0.8e1 * t567) * t345 + t635 - 0.4e1 * t74) * t341) * t346 + (t166 + t472) * t345 + (-t188 + t230) * t340 + t448) * t329 + ((t74 * t327 + (t176 + t567) * t345 + t282 * t491 + t294 * t326 + (-t599 / 0.2e1 - t628 / 0.2e1) * t344 + t517 / 0.2e1 - Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1 - Ifges(6,3) / 0.4e1 + Ifges(5,1) / 0.4e1 + t366 / 0.4e1 + t625 / 0.2e1 - t367 / 0.4e1) * t505 + (t189 * t340 + (t497 + ((t634 - 0.8e1 * t589) * t345 + t158 - 0.4e1 * t135) * t341) * t342 + t447) * t346 + t60 + (t119 + (t516 - t188) * t342) * t345 + (-pkin(3) * t556 - t341 * t533) * t340 + t74 * t630) * t347 + (t345 * t472 - t165 + t448) * t328 + (t59 + (t193 * t232 + t341 * t516) * t345 + t232 * t557 + t74 * t632) * t346 - t241 * t545 + t135 + t642) * t374 + t174 * t370 - t707 * t373 + t372 * t614 - t523 * ((t46 * t632 + t40) * t347 + (t346 * t46 + t41) * t630) + t522 * t676 + (t342 * t38 - t347 * t704) * g(3) - t649 * (t342 * t704 + t38 * t347);
t403 = pkin(1) ^ 2;
t371 = qJD(4) ^ 2;
t281 = qJDD(4) + t503;
t261 = 0.2e1 * t266;
t256 = 0.4e1 * t265;
t212 = -0.4e1 * t229;
t205 = -qJDD(1) * t347 + t471;
t204 = -pkin(1) * t374 - t649;
t192 = mrSges(2,2) + mrSges(3,3) + t433;
t184 = -0.2e1 * t616;
t183 = -0.4e1 * t188;
t164 = pkin(1) * t562;
t160 = 0.2e1 * t166;
t159 = -0.4e1 * t165;
t155 = -0.2e1 * t164;
t154 = t193 * t521;
t150 = t173 * t340;
t148 = -0.2e1 * t500;
t144 = -0.2e1 * t609;
t143 = 0.2e1 * t609;
t136 = -0.2e1 * t496;
t114 = t342 * g(3) + t204 * t347 + (-t329 * t374 - t373) * pkin(2);
t101 = -0.16e2 * t125;
t97 = g(3) * t347 - t204 * t342 + (t342 * t347 * t374 + qJDD(2)) * pkin(2);
t96 = -t150 + t613 / 0.2e1;
t83 = t177 * t340 + t178 * t345;
t79 = qJD(5) - t82;
t77 = t87 * t327;
t69 = t125 + t188 / 0.2e1;
t57 = t340 * t698 + t188;
t55 = (-t120 + t586) * t341;
t52 = (t101 - 0.8e1 * t188) * t341;
t49 = t288 * t339 + t344 * t83;
t48 = t288 * t344 - t339 * t83;
t45 = t683 * t341;
t36 = t98 + t535;
t34 = t124 * t346 + t712;
t33 = t487 * t345;
t31 = t36 * t346;
t29 = -mrSges(3,2) + t703;
t27 = qJDD(1) * t90 + t527 * t89;
t15 = (t36 * t632 + t40) * t347;
t14 = (t426 + t690) * t328;
t11 = t433 * pkin(2) - Ifges(3,5) + t18;
t5 = t18 * t342 - t719 * t347;
t4 = t11 * t342 - t718 * t347;
t3 = t345 * (t114 * t346 + t341 * t97 + (-t177 ^ 2 - t312 ^ 2) * pkin(3)) + t340 * (-t341 * t114 + t346 * t97 + (t177 * t178 + t503) * pkin(3)) + t281 * pkin(6) + t82 * (-pkin(4) * t82 - pkin(6) * t83) - t288 ^ 2 * pkin(4);
t2 = qJDD(1) * pkin(1) + g(1) * t343 - t348 * g(2) + (-t288 * t82 - t27) * pkin(6) + (-t526 + (qJD(1) * t90 + t83) * t288) * pkin(4) + (t341 * (-qJDD(1) * t342 - t347 * t524) - t346 * t205 + (qJD(3) + t312) * t178) * pkin(3) + (-t205 - t471) * pkin(2);
t1 = [(t21 * t347 + t594) * t371 + (t18 * t347 + t598) * t372 + (t11 * t347 + t342 * t718) * t373 + t5 * qJDD(3) + t8 * qJDD(4) + ((t409 + t592 - t605) * t328 + (t401 + t403) * m(6) + t313 * t403 + (-0.2e1 * t461 + t636 - t110) * t345 + t77 + t366 + t246 + t460 + ((t705 * t680 + (0.8e1 * t568 + (t256 + 0.4e1 * t405) * t345 + 0.4e1 * Ifges(4,4) + t159 + t443) * t341) * t346 + (0.2e1 * t230 + t62) * t345 + (t327 * t80 + t345 * t61 + t408) * t328 + t136 + t409 + t606 + t648) * t329 + ((t593 - t709) * t506 + (0.2e1 * t608 + t189 * t692 + (t261 * t345 + t144 + ((t183 + t106) * t345 - 0.2e1 * t605 - t81 - 0.2e1 * t369 + t212 + 0.4e1 * t75 + 0.2e1 * Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t341) * t342) * t346 - 0.4e1 * t483 + ((t534 - t71) * t342 + t341 * t263) * t345 + ((2 * Ifges(3,4)) - 0.2e1 * Ifges(4,4) + t148 + t157 - t122) * t342 + 0.2e1 * pkin(1) * t113) * t347 + 0.2e1 * (t463 * t618 + (-t467 - t593) * t341) * t346 + t142 * t520 + Ifges(3,1) + Ifges(2,3) + t646) * qJDD(1) + t4 * qJDD(2) + ((-t203 * t544 + t346 * t451 + t537) * t347 + pkin(1) * t228 + ((-t203 * t340 + t152) * t346 - t451 * t341) * t342) * qJDD(5) - 0.8e1 * (((t133 * t327 + t96 * t345 + t454) * t328 + (pkin(2) * t558 / 0.4e1 + (-t564 - t484 - t112 / 0.2e1) * t341) * t346 + pkin(2) * t453 + t438) * t329 + ((-t564 + (-t569 - t613 / 0.2e1) * t340 + t529) * t548 + (-t133 * t485 + (-t96 * t542 + t619 / 0.4e1) * t345 - t342 * (pkin(2) * t560 + t126 * t341) / 0.4e1) * t346 + t342 * t474 - (pkin(2) * t559 - 0.2e1 * t133 * t340) * t540 / 0.4e1 + (t227 * t491 + t530) * t342 + pkin(1) * t453) * t347 + t438 * t328 + (t473 * t618 + (t474 + t484 / 0.2e1 + t112 / 0.4e1) * t341) * t346 + t455 * t327 + (t530 * t340 + (-t494 / 0.4e1 + pkin(3) / 0.4e1) * t227) * t345 + t338 * t546 / 0.4e1 - t644) * t523 + (((((t678 + t256) * t345 + t159 + t437) * t328 + (t65 + ((t101 + t183) * t341 + t261) * t345 + (t212 + t80) * t341 - 0.2e1 * t498) * t346 + (t70 + t534) * t345 + pkin(2) * t464 + t431) * t329 + ((t263 * t345 + t545 * t638 + t91) * t346 + t64 + t154 * t345 + pkin(1) * t464 + (t486 * t680 + (t636 + t538) * t345 + t81 + (((-t678 + t259) * t341 + t184) * t345 + (t156 - t120) * t341 - 0.2e1 * t497) * t346 + t635 + 0.8e1 * (t75 + (t110 - t188 / 0.2e1) * t345 - t229 / 0.2e1 + t74) * t328) * t342) * t347 + (t643 + t431) * t328 + (t66 + (t155 + (t125 + t176) * t631) * t345 + t340 * t470 + (t229 + t86) * t632) * t346 + t33 + t461 * t633 + t430) * qJD(1) + t42) * qJD(4) + (((((t52 + t261) * t345 + t144 + t603) * t346 + (t161 + t578) * t345 + t148 + t420) * t329 + ((((-t45 + t184) * t342 + t263) * t345 + (t55 - 0.2e1 * t708) * t342 + t583) * t346 + ((t180 + t538) * t342 + t154) * t345 + (t408 + 0.2e1 * t496) * t342 + t695 * t521 + t622) * t347 + t14 + ((t631 * t69 + t155) * t345 + (-Ifges(4,1) + t86 + t605 + t646) * t632 + t582) * t346 + t33 + 0.2e1 * t462 + t425) * qJD(1) + t623) * qJD(3) - (-t192 * t348 + t343 * t641) * g(1) + ((((t633 * t701 - 0.2e1 * Ifges(4,5) + t105 + 0.2e1 * t284 + 0.2e1 * t285 + 0.2e1 * t350) * t346 + 0.2e1 * t590) * t347 + 0.2e1 * t598) * qJD(3) + (((((-t45 - 0.4e1 * t616) * t342 + t263) * t345 + (t55 - 0.4e1 * t708) * t342 + t583) * t346 + ((-0.4e1 * t230 + t61) * t342 + t154) * t345 + (0.2e1 * Ifges(3,1) - 0.2e1 * Ifges(3,2) + t410 + 0.4e1 * t496) * t342 + t142 * t638 + (-t303 - t606) * t630 + t622) * t347 + t418 + (t160 + t487) * t345 + (((t52 + 0.4e1 * t266) * t345 - 0.4e1 * t609 + t603) * t346 + (-0.4e1 * t166 + t578) * t345 - 0.4e1 * t500 + (4 * Ifges(3,4)) + t420) * t329 + ((t155 + t580) * t345 + t408 * t341 + t143 + t582) * t346 + t14 + t113 * t520) * qJD(1) + t623) * qJD(2) + ((((-pkin(4) * t345 - pkin(3)) * t346 + pkin(4) * t544 - pkin(2)) * t227 + t202 * t173) * t347 - (-pkin(4) * t560 - t173 * t345) * t539 + pkin(4) * t540 * t559 + (-t150 + t613) * t542 - t619) * t370 - (-t192 * t343 - t348 * t641) * g(2); (t227 * t419 + t173) * qJDD(5) + (t228 * t419 - t174) * t370 + (t570 + t229 + (t399 + t400) * m(6) + (-t486 - t700) * pkin(2) + t440) * qJDD(4) + ((t230 + t636) * t345 + t475 * m(6) + (-t566 + t711) * pkin(2) + t421) * qJDD(3) - ((-t490 - t714) * t347 - t29 * t342) * g(3) + t404 + (t201 * t227 * t680 + t446) * t442 + t4 * qJDD(1) + ((t402 + t475) * m(6) + 0.2e1 * t570 + t303 + t136 + t421 + t711 * t680 + Ifges(3,3)) * qJDD(2) + (((((t424 + 0.4e1 * t620) * t340 + t259) * t345 + t426) * t328 + (t580 * t345 + t143 + (t337 + t423 - 0.4e1 * t588 + 0.2e1 * t621 + t689) * t341) * t346 + (t160 + t579) * t345 + t418) * t329 + (-0.4e1 * ((t400 / 0.2e1 - t401 / 0.2e1 - t399 / 0.2e1) * m(6) + t588 + t406) * t548 + ((0.2e1 * t167 - t267) * t345 + t149) * t346 + mrSges(3,2) * pkin(1) + ((t189 + 0.2e1 * t499) * t345 - 0.2e1 * t501 + t608) * t341 + ((0.2e1 * t708 + 0.4e1 * ((t340 * t84 + t265) * t345 + t467) * t341) * t346 - 0.4e1 * t69 * t345 + (t402 - t280) * m(6) + t648 + t654) * t342) * t347 + (t345 * t579 + t425) * t328 + ((t164 + t266) * t345 + t495 - t609 + (t592 - t621 + t654) * t341) * t346 + ((t502 - t616) * t341 + (t427 + t620) * t340 - t265) * t345 + (-pkin(1) * t565 - t708) * t341 + t490 * t618 + Ifges(3,4) + t709) * t374 + (t15 + (-t118 * t544 + t31 + t537) * t630) * t523 + (0.2e1 * t522 + t372) * pkin(2) * t703 - (t371 + (2 * t682)) * t707 - t649 * (-t26 * t342 + t29 * t347); t446 * t442 - t112 * qJDD(5) + (-t228 * t612 - t174) * t370 + t640 * qJDD(3) + (pkin(2) * t714 + t640) * qJDD(2) - (((t532 + t672) * t346 + t574) * t347 + t34 * t342) * g(3) + t404 + t371 * t614 + t34 * pkin(2) * t373 + t5 * qJDD(1) + (((t449 + t665 + t690) * t328 + (t58 + ((-t340 * t696 + t180) * t341 - t266) * t345 + (0.2e1 * t302 + t664) * t341 + t609) * t346 + (t166 + t577) * t345 + t341 * t708 + t432) * t329 + ((t77 + t57 * t345 + t355 - Ifges(6,3) / 0.2e1 - t639 - t685 - t694) * t506 + (t149 + (-0.8e1 * (t589 - t265 / 0.2e1) * t541 - 0.4e1 * t439 * t341 + t708) * t342 + t447) * t346 + t60 + (-0.2e1 * t342 * t57 + t119) * t345 + (-t501 + t608) * t341 + (t655 + t694) * t630) * t347 + (t345 * t577 + t432) * t328 + (t59 + ((t340 * t697 + t182) * t341 + t164) * t345 + (-t302 + t663) * t341 + t695 * t618) * t346 - t462 + t439 + t642) * t374 + (t15 + (t31 + t41) * t630) * t523 + t37 * qJDD(4) + t676 * t682 + t649 * (t34 * t347 + t342 * (-t574 + t711)); Ifges(6,5) * (qJD(5) * t48 + t27 * t344 + t281 * t339) + Ifges(6,6) * (-qJD(5) * t49 - t27 * t339 + t281 * t344) + Ifges(6,3) * (t527 * t90 + qJDD(5) - t526) + t49 * (Ifges(6,4) * t49 + Ifges(6,2) * t48 + Ifges(6,6) * t79) - t48 * (Ifges(6,1) * t49 + Ifges(6,4) * t48 + Ifges(6,5) * t79) + mrSges(6,1) * (t2 * t344 - t3 * t339) - mrSges(6,2) * (t2 * t339 + t3 * t344);];
tau = t1(:);
