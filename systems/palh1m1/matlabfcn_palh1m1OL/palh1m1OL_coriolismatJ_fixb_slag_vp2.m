% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [13x13]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m1OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:29:01
% EndTime: 2020-04-15 19:29:49
% DurationCPUTime: 11.43s
% Computational Cost: add. (13766->568), mult. (29349->777), div. (0->0), fcn. (34650->20), ass. (0->360)
t370 = sin(pkin(19));
t371 = cos(pkin(19));
t372 = cos(qJ(10));
t590 = sin(qJ(10));
t334 = t370 * t372 - t371 * t590;
t335 = t370 * t590 + t371 * t372;
t375 = sin(qJ(7));
t599 = sin(qJ(2));
t601 = cos(qJ(7));
t605 = cos(qJ(2));
t339 = t375 * t599 - t601 * t605;
t340 = -t375 * t605 - t599 * t601;
t188 = t334 * t340 + t335 * t339;
t468 = t334 * t339 - t335 * t340;
t735 = Ifges(11,5) * t468 - Ifges(11,6) * t188;
t733 = qJD(10) * t735;
t374 = sin(qJ(8));
t600 = cos(qJ(8));
t337 = t374 * t599 - t600 * t605;
t338 = -t374 * t605 - t599 * t600;
t373 = sin(qJ(9));
t378 = cos(qJ(9));
t435 = -t337 * t373 - t338 * t378;
t436 = -t337 * t378 + t338 * t373;
t728 = Ifges(10,5) * t435 + Ifges(10,6) * t436;
t734 = qJD(9) * t728;
t377 = sin(qJ(4));
t593 = pkin(5) * t377;
t363 = pkin(11) + t593;
t598 = sin(qJ(3));
t604 = cos(qJ(3));
t341 = t598 * t605 + t599 * t604;
t401 = -t598 * t599 + t604 * t605;
t603 = cos(qJ(4));
t390 = t341 * t603 + t377 * t401;
t273 = t341 * t377 - t401 * t603;
t376 = sin(qJ(5));
t697 = t273 * t376;
t705 = -mrSges(6,2) * t390 + mrSges(6,3) * t697;
t731 = t363 * t705;
t457 = t604 * t603;
t505 = t598 * pkin(1);
t459 = t505 + pkin(5);
t309 = pkin(1) * t457 - t377 * t459;
t304 = pkin(11) - t309;
t379 = cos(qJ(5));
t724 = t379 * t705;
t696 = t273 * t379;
t704 = mrSges(6,1) * t390 + mrSges(6,3) * t696;
t725 = t376 * t704;
t692 = t724 / 0.2e1 - t725 / 0.2e1;
t730 = t692 * t304;
t631 = t188 / 0.2e1;
t729 = -t631 - t188 / 0.2e1;
t655 = -t725 + t724;
t629 = t435 / 0.2e1;
t727 = 0.2e1 * t629;
t700 = t273 / 0.2e1;
t673 = Ifges(5,4) * t390;
t726 = -t673 / 0.2e1;
t576 = Ifges(6,2) * t376;
t581 = Ifges(6,4) * t379;
t446 = -t576 + t581;
t582 = Ifges(6,4) * t376;
t586 = Ifges(6,1) * t379;
t447 = -t582 + t586;
t608 = -t379 / 0.2e1;
t609 = -t376 / 0.2e1;
t350 = Ifges(6,1) * t376 + t581;
t515 = t379 * t350;
t349 = Ifges(6,2) * t379 + t582;
t519 = t376 * t349;
t706 = t519 / 0.2e1 - t515 / 0.2e1;
t385 = t446 * t608 + t447 * t609 + t706;
t684 = Ifges(6,6) * t390 - t273 * t446;
t709 = t379 * t684;
t683 = Ifges(6,5) * t390 - t273 * t447;
t710 = t376 * t683;
t723 = t709 / 0.4e1 + t710 / 0.4e1;
t348 = Ifges(6,5) * t376 + Ifges(6,6) * t379;
t662 = t390 * t348;
t669 = Ifges(5,6) * t390;
t722 = t662 / 0.2e1 - t669 + t709 / 0.2e1 + t710 / 0.2e1;
t721 = -t376 * t705 - t379 * t704;
t482 = t363 * t609;
t607 = t379 / 0.2e1;
t720 = t482 * t704 + t607 * t731;
t719 = Ifges(5,2) * t700 + t726;
t717 = -t683 / 0.4e1;
t716 = -t684 / 0.4e1;
t715 = -t696 / 0.2e1;
t714 = t697 / 0.2e1;
t699 = Ifges(5,5) * t273;
t713 = t699 / 0.2e1;
t368 = t376 ^ 2;
t369 = t379 ^ 2;
t512 = t368 + t369;
t712 = t512 * mrSges(6,3);
t562 = t369 * mrSges(6,3);
t563 = t368 * mrSges(6,3);
t708 = -t563 / 0.2e1 - t562 / 0.2e1;
t707 = -t662 / 0.4e1 + t669 / 0.2e1 + t713;
t702 = -pkin(9) / 0.2e1;
t701 = -t273 / 0.2e1;
t270 = Ifges(5,4) * t273;
t698 = -mrSges(5,2) + t712;
t365 = Ifges(6,5) * t379;
t573 = Ifges(6,6) * t376;
t695 = (t365 / 0.2e1 - t573 / 0.2e1) * t273;
t561 = t376 * mrSges(6,1);
t588 = mrSges(6,2) * t379;
t448 = t561 + t588;
t694 = t448 * t273;
t596 = pkin(1) * t375;
t345 = pkin(4) * t370 + t596;
t506 = t601 * pkin(1);
t347 = pkin(4) * t371 + t506;
t197 = -t334 * t345 - t335 * t347;
t199 = t334 * t347 - t335 * t345;
t33 = -mrSges(11,1) * t199 - t197 * mrSges(11,2);
t689 = qJD(10) * t33;
t249 = (-t334 * t370 - t335 * t371) * pkin(4);
t678 = (t334 * t371 - t335 * t370) * pkin(4);
t39 = -mrSges(11,1) * t678 - t249 * mrSges(11,2);
t688 = t39 * qJD(10);
t660 = t365 - t573;
t78 = Ifges(6,3) * t273 + t390 * t660;
t687 = t78 / 0.2e1 + t719 + t726;
t171 = pkin(9) * t390 + pkin(11) * t273;
t682 = t729 * t188 * Ifges(11,4) + (0.2e1 * Ifges(11,1) * t631 + Ifges(11,4) * t468 + Ifges(11,2) * t729) * t468;
t310 = -pkin(2) * t338 - pkin(15);
t26 = (t310 * mrSges(10,2) + Ifges(10,4) * t727) * t435 + (-t310 * mrSges(10,1) + Ifges(10,2) * t727 - Ifges(10,4) * t436 + (-t435 / 0.2e1 - t629) * Ifges(10,1)) * t436;
t624 = t390 / 0.2e1;
t680 = t401 / 0.2e1;
t367 = t605 * pkin(1);
t675 = mrSges(6,3) * (t369 / 0.2e1 + t368 / 0.2e1);
t170 = Ifges(5,1) * t390 - t270;
t420 = t447 * t390;
t578 = Ifges(6,5) * t273;
t83 = t420 + t578;
t554 = t379 * t83;
t419 = t446 * t390;
t574 = Ifges(6,6) * t273;
t80 = t419 + t574;
t559 = t376 * t80;
t659 = -t170 / 0.2e1 - t554 / 0.2e1 + t559 / 0.2e1 + t270 / 0.2e1;
t360 = pkin(1) * t599 - pkin(15);
t291 = -pkin(5) * t401 + t360;
t125 = pkin(9) * t273 - pkin(11) * t390 + t291;
t532 = t390 * t379;
t162 = mrSges(6,1) * t273 - mrSges(6,3) * t532;
t516 = t379 * t162;
t533 = t390 * t376;
t159 = -mrSges(6,2) * t273 - mrSges(6,3) * t533;
t523 = t376 * t159;
t438 = t516 + t523;
t475 = m(6) * t512;
t658 = t125 * t475 + t438;
t154 = t448 * t390;
t356 = t377 * t604 * pkin(1);
t308 = t459 * t603 + t356;
t303 = -pkin(9) - t308;
t622 = t303 / 0.2e1;
t657 = -t309 * t154 / 0.2e1 - t694 * t622;
t560 = t376 * mrSges(6,2);
t589 = mrSges(6,1) * t379;
t431 = -t560 / 0.2e1 + t589 / 0.2e1;
t481 = -t516 / 0.2e1;
t653 = t390 * t708 + t481;
t651 = (mrSges(4,1) * t604 - mrSges(4,2) * t598) * pkin(1);
t541 = t678 * t188;
t543 = t249 * t468;
t650 = (-t541 / 0.2e1 - t543 / 0.2e1) * mrSges(11,3);
t343 = (mrSges(10,1) * t373 + mrSges(10,2) * t378) * pkin(2);
t510 = qJD(2) + qJD(8);
t38 = t510 * t343;
t66 = mrSges(11,1) * t188 + mrSges(11,2) * t468;
t649 = t694 * t702 + t707;
t646 = -pkin(11) / 0.2e1;
t645 = pkin(5) * m(6);
t644 = -mrSges(6,1) / 0.2e1;
t643 = mrSges(6,1) / 0.2e1;
t642 = -mrSges(6,2) / 0.2e1;
t641 = mrSges(6,2) / 0.2e1;
t640 = -Ifges(6,3) / 0.2e1;
t638 = -t162 / 0.2e1;
t627 = -t273 / 0.4e1;
t626 = t273 / 0.4e1;
t621 = -t304 / 0.2e1;
t620 = t304 / 0.2e1;
t619 = t308 / 0.2e1;
t507 = t603 * pkin(5);
t364 = -t507 - pkin(9);
t313 = t364 * t448;
t618 = -t313 / 0.2e1;
t315 = (t377 * t598 - t457) * pkin(1);
t617 = t315 / 0.2e1;
t316 = t505 * t603 + t356;
t616 = -t316 / 0.2e1;
t615 = -t349 / 0.4e1;
t614 = t349 / 0.4e1;
t613 = -t350 / 0.4e1;
t612 = t350 / 0.4e1;
t611 = t363 / 0.2e1;
t610 = -t364 / 0.2e1;
t602 = cos(qJ(6));
t597 = sin(qJ(6));
t594 = pkin(5) * t341;
t46 = Ifges(9,5) * t338 + Ifges(9,6) * t337 + (-t373 * t436 + t378 * t435) * pkin(2) * mrSges(10,3) + t728;
t591 = qJD(8) * t46 + t734;
t587 = Ifges(5,1) * t273;
t585 = Ifges(4,4) * t341;
t571 = Ifges(6,3) * t390;
t128 = t171 + t594;
t126 = t128 + t367;
t19 = t26 + (-pkin(15) * mrSges(9,2) + Ifges(9,4) * t338) * t338 + (pkin(15) * mrSges(9,1) + (-Ifges(9,1) + Ifges(9,2)) * t338 + (-m(10) * t310 + mrSges(10,1) * t435 + mrSges(10,2) * t436) * pkin(2) - Ifges(9,4) * t337) * t337;
t231 = (t339 * t370 - t340 * t371) * pkin(4) + t360;
t264 = (-t339 * t371 - t340 * t370) * pkin(4);
t234 = t367 + t264;
t325 = Ifges(4,4) * t401;
t382 = t340 ^ 2 * Ifges(8,4) + t360 * (-mrSges(8,1) * t339 + mrSges(8,2) * t340) + t682;
t393 = Ifges(4,2) * t401 + t585;
t397 = Ifges(4,1) * t401;
t427 = -Ifges(8,4) * t339 + (-Ifges(8,1) + Ifges(8,2)) * t340;
t466 = -(mrSges(5,1) * t390 - mrSges(5,2) * t273) * t291 - t360 * (mrSges(4,1) * t341 + mrSges(4,2) * t401) + t721 * t125;
t474 = m(5) * t291 + mrSges(5,1) * t273 + mrSges(5,2) * t390;
t67 = -mrSges(11,1) * t468 + mrSges(11,2) * t188;
t1 = -t466 + (-t273 * t660 + t571) * t700 + t390 * t719 + (-t401 * mrSges(4,1) - mrSges(8,1) * t340 + (m(4) + m(8)) * t360) * t367 + 0.2e1 * t680 * t325 + t19 + t382 - pkin(15) * (mrSges(3,1) * t605 - mrSges(3,2) * t599) + (t599 ^ 2 - t605 ^ 2) * Ifges(3,4) + t474 * (t367 + t594) + (m(11) * t234 + t66) * t231 + t683 * t532 / 0.2e1 - t684 * t533 / 0.2e1 + t83 * t715 + t80 * t714 + (-mrSges(8,2) * t367 + t427) * t339 + t234 * t67 + t658 * t126 + (-t393 / 0.2e1 + mrSges(4,2) * t367 + t397 / 0.2e1 - t585 / 0.2e1 + (-Ifges(4,2) + Ifges(4,1)) * t680) * t341 + (-t673 - t587 + t78) * t624 + (Ifges(3,2) - Ifges(3,1)) * t605 * t599 + (-Ifges(5,2) * t390 + t170 - t270) * t701;
t570 = t1 * qJD(1);
t503 = t640 - Ifges(5,2) / 0.2e1;
t2 = -t401 * t325 - (t660 * t701 + t659) * t273 + (-pkin(5) * t474 + t393 - t397) * t341 + t466 - t658 * t128 + (t503 * t273 + t587 / 0.2e1 + t683 * t608 + t376 * t684 / 0.2e1 - t687) * t390;
t569 = t2 * qJD(1);
t3 = (t291 * mrSges(5,1) + t607 * t683 + t609 * t684 + t687) * t390 + (-t291 * mrSges(5,2) - t695 + (-Ifges(5,1) / 0.2e1 - t503) * t390 + t659) * t273 + t438 * t171 + (t171 * t475 - t721) * t125;
t568 = t3 * qJD(1);
t567 = t308 * mrSges(5,2);
t566 = t309 * mrSges(5,1);
t565 = t315 * mrSges(5,1);
t564 = t316 * mrSges(5,2);
t551 = qJD(1) * t19;
t11 = t264 * t67 + (m(11) * t264 + t66) * t231 + t427 * t339 + t382;
t549 = t11 * qJD(1);
t156 = t390 * t349;
t157 = t390 * t350;
t517 = t379 * t159;
t521 = t376 * t162;
t17 = (t517 - t521) * t125 + (-t157 * t607 + t348 * t701 + t80 * t608 + (-t156 + t83) * t609) * t390;
t546 = t17 * qJD(1);
t18 = t231 * t66 + t682;
t545 = t18 * qJD(1);
t252 = (-t334 * t375 - t335 * t601) * pkin(1);
t540 = t252 * mrSges(11,2);
t539 = t26 * qJD(1);
t535 = t390 * t377;
t530 = t303 * t694;
t529 = t308 * t273;
t528 = t309 * t390;
t449 = -t560 + t589;
t527 = t309 * t449;
t526 = t315 * t449;
t524 = t364 * t694;
t172 = -pkin(14) * (mrSges(7,1) * t597 + mrSges(7,2) * t602) + (t597 ^ 2 - t602 ^ 2) * Ifges(7,4) + (-Ifges(7,1) + Ifges(7,2)) * t602 * t597;
t511 = t172 * qJD(1);
t509 = t449 * t593;
t508 = -t603 / 0.2e1;
t504 = -t593 / 0.2e1;
t500 = t197 * t468 * mrSges(11,3);
t499 = t199 * t188 * mrSges(11,3);
t497 = t571 / 0.2e1;
t491 = t603 * t273;
t484 = t154 * t617;
t483 = -t523 / 0.2e1;
t479 = -t249 / 0.2e1 - t197 / 0.2e1;
t477 = t611 + t621;
t472 = t512 * t308;
t471 = t512 * t316;
t470 = t512 * t363;
t465 = -t507 / 0.2e1;
t464 = t507 / 0.2e1;
t461 = Ifges(8,5) * t340 + Ifges(8,6) * t339 + t735;
t454 = t376 * t508;
t451 = t614 - t586 / 0.4e1;
t429 = pkin(9) * t448;
t444 = -t429 / 0.2e1 - t385;
t395 = mrSges(6,3) * t472 + t527 + t566 - t567;
t34 = m(6) * (-t303 * t309 + t304 * t472) + t395;
t418 = t515 * t627 + t519 * t626 - t707 + t723;
t4 = (t159 * t619 - t273 * t613 + t716 + (t620 + t646) * t705) * t379 + (-t273 * t614 + t308 * t638 + t717 + (t621 + pkin(11) / 0.2e1) * t704) * t376 + t418 + t649 + t657;
t443 = t4 * qJD(1) + t34 * qJD(2);
t423 = t303 * t448;
t138 = t423 - t385;
t404 = (t700 + t626) * Ifges(6,6) + t157 / 0.4e1 + t80 / 0.4e1;
t408 = t612 + t581 / 0.2e1 - t576 / 0.4e1;
t428 = t365 * t627 + t497;
t430 = t156 / 0.4e1 - t83 / 0.4e1 + Ifges(6,5) * t701;
t9 = (t126 * t642 + t159 * t620 + t404) * t376 + (t126 * t643 + t162 * t620 + t430) * t379 + ((t303 * t644 + t451) * t379 + (mrSges(6,2) * t622 + t408) * t376 + t304 * t675) * t390 + t428;
t442 = -qJD(1) * t9 + qJD(2) * t138;
t394 = t390 * t617 + t273 * t616 + t529 / 0.2e1 + t528 / 0.2e1;
t16 = t484 - (t622 + t610) * t694 + (t316 * t159 / 0.2e1 - t477 * t705) * t379 + (t162 * t616 + t477 * t704) * t376 + ((t535 / 0.2e1 - t491 / 0.2e1) * pkin(5) + t394) * mrSges(5,3);
t35 = (-mrSges(5,1) - t449) * t315 + t651 + t698 * t316 + m(6) * (t303 * t315 + t304 * t471) + m(5) * (-t308 * t315 - t309 * t316);
t441 = qJD(1) * t16 + qJD(2) * t35;
t413 = pkin(1) * (-t334 * t601 + t335 * t375);
t384 = t500 / 0.2e1 + t499 / 0.2e1 + (-t252 * t468 / 0.2e1 + t413 * t631) * mrSges(11,3);
t21 = t384 + t650;
t247 = mrSges(11,1) * t413;
t398 = -mrSges(8,1) * t596 - mrSges(8,2) * t506 + t247 - t540;
t36 = -m(11) * (t197 * t413 + t199 * t252) - t398;
t440 = -qJD(1) * t21 - qJD(2) * t36;
t439 = t33 * qJD(2);
t432 = t561 / 0.2e1 + t588 / 0.2e1;
t425 = -t521 / 0.2e1 + t517 / 0.2e1;
t417 = t154 * t504 - t610 * t694;
t416 = t449 * t702;
t415 = t429 / 0.2e1;
t414 = t512 * t603;
t412 = t432 * t308;
t411 = t432 * t316;
t410 = t449 * t624;
t386 = -mrSges(5,1) * t593 + t507 * t698 - t509;
t139 = (t363 * t414 + t364 * t377) * t645 + t386;
t380 = m(6) * (-t364 * t309 + t308 * t470 + (t303 * t377 + t304 * t414) * pkin(5)) / 0.2e1 - t567 / 0.2e1 + t566 / 0.2e1 + t527 / 0.2e1 + mrSges(5,1) * t504 - t509 / 0.2e1 + mrSges(5,2) * t465 + (t464 + t619) * t712;
t383 = -m(6) * (-pkin(9) * t315 + pkin(11) * t471) / 0.2e1 + t565 / 0.2e1 + t526 / 0.2e1 + t564 / 0.2e1 + t708 * t316;
t25 = t380 + t383;
t400 = -t649 + (t519 / 0.4e1 - t515 / 0.4e1) * t273 + t692 * pkin(11) + t723;
t6 = t400 + t713 + (t273 * t612 - t731 / 0.2e1 + t159 * t465 + t716) * t379 + (t162 * t464 + t273 * t615 + t611 * t704 + t717) * t376 + (-t348 / 0.4e1 + Ifges(5,6) / 0.2e1) * t390 + t417;
t407 = -qJD(1) * t6 + qJD(2) * t25 + qJD(3) * t139;
t389 = (-t678 / 0.2e1 - t199 / 0.2e1) * mrSges(11,1);
t22 = -t247 / 0.2e1 + (t252 / 0.2e1 + t479) * mrSges(11,2) + t389;
t406 = qJD(2) * t22 + qJD(7) * t39;
t12 = (t128 * t642 + t159 * t611 + t404) * t376 + (t128 * t643 + t162 * t611 + t430) * t379 + ((mrSges(6,1) * t610 + t451) * t379 + (t364 * t641 + t408) * t376 + t363 * t675) * t390 + t428;
t164 = -t313 + t385;
t381 = -t423 / 0.2e1 + t385;
t62 = t618 - t411 + t381;
t405 = -t12 * qJD(1) - t62 * qJD(2) - t164 * qJD(3);
t403 = t533 * t613 + t532 * t615 + t660 * t626 - t559 / 0.4e1 + t554 / 0.4e1 - (t157 + t419) * t376 / 0.4e1 + (-t156 + t420) * t379 / 0.4e1;
t402 = Ifges(4,5) * t401 - Ifges(4,6) * t341 + t273 * t706 - t699 + t722;
t14 = (pkin(11) * t638 + t171 * t644 + t578 / 0.2e1) * t379 + (t159 * t646 + t171 * t641 - t574 / 0.2e1) * t376 + (-pkin(11) * t675 + t416 + t640) * t390 + t403;
t173 = t429 + t385;
t64 = t415 - t412 + t381;
t396 = (mrSges(6,1) * t454 + t508 * t588) * pkin(5);
t97 = t415 + t618 + t396 + t385;
t399 = qJD(1) * t14 - qJD(2) * t64 - qJD(3) * t97 - qJD(4) * t173;
t388 = t400 + t418;
t387 = t403 + t497 - t695;
t336 = t343 * qJD(9);
t307 = t313 / 0.2e1;
t282 = t423 / 0.2e1;
t98 = t307 + t396 + t444;
t65 = t282 - t412 + t444;
t63 = t282 + t307 - t411 - t385;
t37 = t336 + t38;
t24 = t380 - t383;
t23 = -t540 / 0.2e1 + t247 / 0.2e1 + t479 * mrSges(11,2) + t389;
t20 = t650 - t384 + t461;
t15 = (-t390 * t675 + t481 + t483) * pkin(11) + Ifges(6,5) * t715 + Ifges(6,6) * t714 + t403 + t390 * t416 + Ifges(6,3) * t624 + t431 * t171;
t13 = t431 * t128 + t159 * t482 + t363 * t653 + t364 * t410 + t387;
t10 = t431 * t126 + t303 * t410 + t387 + (t483 + t653) * t304;
t8 = t425 * t316 + t730 + t402 - t524 / 0.2e1 + t484 - t530 / 0.2e1 + (-t273 * t465 + t390 * t504 + t394) * mrSges(5,3) + t720;
t7 = pkin(5) * t162 * t454 + t464 * t517 + t388 - t417 + t720;
t5 = t425 * t308 + t388 + t657 + t730;
t27 = [qJD(10) * t18 + qJD(2) * t1 - qJD(3) * t2 + qJD(4) * t3 + qJD(5) * t17 - qJD(6) * t172 + qJD(7) * t11 + qJD(8) * t19 + qJD(9) * t26, t733 + t8 * qJD(3) + t5 * qJD(4) + t10 * qJD(5) + t20 * qJD(7) + t570 + t591 + (-Ifges(3,5) * t599 - Ifges(3,6) * t605 + t402 + t46 + t461 - t499 - t500 - t530 + ((t339 * t375 - t340 * t601) * mrSges(8,3) + (t341 * t604 - t401 * t598) * mrSges(4,3)) * pkin(1) + t655 * t304 + (t528 + t529) * mrSges(5,3)) * qJD(2), -t569 + t8 * qJD(2) + (-t524 + t655 * t363 + (t491 - t535) * mrSges(5,3) * pkin(5) + t402) * qJD(3) + t7 * qJD(4) + t13 * qJD(5), t5 * qJD(2) + t7 * qJD(3) + t15 * qJD(5) + t568 + (pkin(9) * t694 + (-Ifges(5,5) + t706) * t273 + t655 * pkin(11) + t722) * qJD(4), t546 + t10 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (-t125 * t448 - t662) * qJD(5), -t511 + (Ifges(7,5) * t602 - Ifges(7,6) * t597) * qJD(6), t549 + t20 * qJD(2) + ((-t541 - t543) * mrSges(11,3) + t461) * qJD(7) + t733, qJD(2) * t46 + t551 + t591, t510 * t728 + t539 + t734, t733 + t545 + (qJD(2) + qJD(7)) * t735, 0, 0, 0; qJD(3) * t16 + qJD(4) * t4 - qJD(5) * t9 - qJD(7) * t21 - t570, qJD(3) * t35 + qJD(4) * t34 + qJD(5) * t138 - qJD(7) * t36 + t336 + t689, (t316 * t562 + t316 * t563 + m(6) * (t315 * t364 + t316 * t470) - t526 + m(5) * (-t315 * t603 + t316 * t377) * pkin(5) - t564 - t565 + t651) * qJD(3) + t24 * qJD(4) + t63 * qJD(5) + t441, t24 * qJD(3) + (m(6) * (pkin(9) * t309 + pkin(11) * t472) + t395) * qJD(4) + t65 * qJD(5) + t443, t63 * qJD(3) + t65 * qJD(4) + (-t304 * t449 + t660) * qJD(5) + t442, 0, (m(11) * (t249 * t413 + t252 * t678) + t398) * qJD(7) + t23 * qJD(10) + t440, t336, t37, qJD(7) * t23 + t439 + t689, 0, 0, 0; -qJD(2) * t16 - qJD(4) * t6 - qJD(5) * t12 + t569, qJD(4) * t25 - qJD(5) * t62 - t441, qJD(4) * t139 - qJD(5) * t164, ((-pkin(9) * t377 + pkin(11) * t414) * t645 + t386) * qJD(4) + t98 * qJD(5) + t407, t98 * qJD(4) + (-t363 * t449 + t660) * qJD(5) + t405, 0, 0, 0, 0, 0, 0, 0, 0; -qJD(2) * t4 + qJD(3) * t6 + qJD(5) * t14 - t568, -qJD(3) * t25 - qJD(5) * t64 - t443, -qJD(5) * t97 - t407, -t173 * qJD(5), (-pkin(11) * t449 + t660) * qJD(5) + t399, 0, 0, 0, 0, 0, 0, 0, 0; qJD(2) * t9 + qJD(3) * t12 - qJD(4) * t14 - t546, qJD(3) * t62 + qJD(4) * t64 - t442, qJD(4) * t97 - t405, -t399, 0, 0, 0, 0, 0, 0, 0, 0, 0; t511, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; qJD(2) * t21 - t549, qJD(10) * t22 - t440, 0, 0, 0, 0, t688, 0, 0, t406 + t688, 0, 0, 0; -t551, t336, 0, 0, 0, 0, 0, t336, t37, 0, 0, 0, 0; -t539, -t38, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0; -t545, -qJD(7) * t22 - t439, 0, 0, 0, 0, -t406, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t27;
