% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [13x13]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:18:02
% EndTime: 2020-05-02 21:18:46
% DurationCPUTime: 10.05s
% Computational Cost: add. (13766->594), mult. (29379->785), div. (0->0), fcn. (34561->20), ass. (0->378)
t373 = sin(qJ(8));
t379 = sin(qJ(2));
t385 = cos(qJ(2));
t602 = cos(qJ(8));
t326 = t373 * t379 - t602 * t385;
t327 = -t373 * t385 - t602 * t379;
t372 = sin(qJ(9));
t380 = cos(qJ(9));
t442 = -t326 * t372 - t327 * t380;
t443 = -t326 * t380 + t327 * t372;
t711 = Ifges(10,5) * t442 + Ifges(10,6) * t443;
t712 = qJD(9) * t711;
t382 = cos(qJ(5));
t363 = Ifges(6,5) * t382;
t384 = cos(qJ(3));
t525 = t379 * t384;
t378 = sin(qJ(3));
t526 = t378 * t385;
t330 = t525 + t526;
t519 = t384 * t385;
t331 = -t379 * t378 + t519;
t377 = sin(qJ(4));
t383 = cos(qJ(4));
t441 = -t377 * t330 + t331 * t383;
t376 = sin(qJ(5));
t579 = Ifges(6,6) * t376;
t702 = (t363 / 0.2e1 - t579 / 0.2e1) * t441;
t369 = sin(pkin(19));
t370 = cos(pkin(19));
t371 = cos(qJ(10));
t595 = sin(qJ(10));
t323 = t369 * t371 - t370 * t595;
t324 = t369 * t595 + t370 * t371;
t374 = sin(qJ(7));
t340 = pkin(1) * t374 + pkin(4) * t369;
t603 = cos(qJ(7));
t341 = pkin(1) * t603 + t370 * pkin(4);
t193 = t323 * t341 - t324 * t340;
t553 = t193 * mrSges(11,1);
t191 = -t323 * t340 - t324 * t341;
t555 = t191 * mrSges(11,2);
t32 = -t553 - t555;
t710 = qJD(10) * t32;
t680 = (t323 * t370 - t324 * t369) * pkin(4);
t550 = t680 * mrSges(11,1);
t243 = (-t323 * t369 - t324 * t370) * pkin(4);
t551 = t243 * mrSges(11,2);
t38 = -t550 - t551;
t709 = t38 * qJD(10);
t600 = pkin(1) * t384;
t350 = t383 * t600;
t601 = pkin(1) * t378;
t360 = pkin(5) + t601;
t302 = t360 * t377 - t350;
t296 = pkin(11) + t302;
t351 = t377 * t600;
t303 = t360 * t383 + t351;
t297 = -pkin(9) - t303;
t367 = t376 ^ 2;
t368 = t382 ^ 2;
t516 = t367 + t368;
t704 = t516 * t383;
t707 = t296 * t704 + t297 * t377;
t328 = t374 * t379 - t603 * t385;
t329 = -t374 * t385 - t603 * t379;
t181 = t323 * t329 + t324 * t328;
t632 = -t181 / 0.2e1;
t706 = t376 / 0.2e1;
t604 = t382 / 0.2e1;
t705 = Ifges(11,4) * t181;
t272 = t330 * t383 + t331 * t377;
t581 = Ifges(6,6) * t272;
t365 = Ifges(6,4) * t382;
t664 = -Ifges(6,2) * t376 + t365;
t74 = t441 * t664 + t581;
t588 = Ifges(6,4) * t376;
t339 = Ifges(6,1) * t382 - t588;
t585 = Ifges(6,5) * t272;
t76 = t339 * t441 + t585;
t665 = t363 - t579;
t683 = -t665 / 0.2e1;
t686 = t441 / 0.2e1;
t703 = -(Ifges(5,2) + Ifges(6,3)) * t686 + (-Ifges(5,4) - t683) * t272;
t701 = t193 / 0.2e1;
t566 = t382 * mrSges(6,2);
t570 = t376 * mrSges(6,1);
t454 = t566 + t570;
t667 = t454 * t441;
t699 = pkin(9) * t667;
t698 = Ifges(9,1) - Ifges(9,2);
t697 = mrSges(6,1) * t272;
t696 = mrSges(6,2) * t272;
t260 = Ifges(5,6) * t272;
t646 = t516 * mrSges(6,3);
t695 = -mrSges(5,2) + t646;
t335 = t376 * Ifges(6,5) + t382 * Ifges(6,6);
t547 = t272 * t335;
t338 = Ifges(6,1) * t376 + t365;
t521 = t382 * t338;
t336 = t382 * Ifges(6,2) + t588;
t530 = t376 * t336;
t423 = -t530 / 0.2e1 + t521 / 0.2e1;
t531 = t376 * t441;
t152 = -mrSges(6,3) * t531 - t696;
t524 = t382 * t152;
t546 = t441 * t382;
t155 = -mrSges(6,3) * t546 + t697;
t533 = t376 * t155;
t694 = t524 / 0.2e1 - t533 / 0.2e1;
t598 = pkin(5) * t383;
t362 = -pkin(9) - t598;
t599 = pkin(5) * t377;
t361 = pkin(11) + t599;
t469 = t516 * t361;
t693 = t302 * t362 + t303 * t469;
t586 = Ifges(10,4) * t443;
t692 = Ifges(10,1) * t442 + t586;
t671 = Ifges(11,6) * t181;
t465 = t323 * t328 - t324 * t329;
t673 = Ifges(11,5) * t465;
t691 = -t671 + t673;
t609 = -t336 / 0.4e1;
t476 = t339 / 0.4e1 + t609;
t149 = t272 * t336;
t584 = Ifges(6,5) * t441;
t77 = t272 * t339 - t584;
t501 = -t149 / 0.4e1 + t77 / 0.4e1;
t400 = t476 * t272 + t501 - t584 / 0.2e1;
t652 = -t673 / 0.2e1 + t671 / 0.2e1;
t165 = pkin(9) * t272 - pkin(11) * t441;
t460 = t339 * t706 + t664 * t604 + t423;
t304 = -pkin(2) * t327 - pkin(15);
t689 = (-t304 * mrSges(10,2) - Ifges(10,4) * t442) * t442 + (t304 * mrSges(10,1) - (-Ifges(10,1) + Ifges(10,2)) * t442 / 0.2e1) * t443;
t633 = t465 / 0.2e1;
t688 = (Ifges(11,1) * t465 - t705) * t181 / 0.2e1 + (Ifges(11,2) * t465 + t705) * t632 + (0.2e1 * Ifges(11,4) * t465 + (Ifges(11,1) - Ifges(11,2)) * t181) * t633;
t619 = t272 / 0.2e1;
t610 = t327 / 0.2e1;
t687 = -t441 / 0.4e1;
t621 = t441 / 0.4e1;
t684 = -t443 / 0.2e1;
t679 = mrSges(6,3) * (t368 / 0.2e1 + t367 / 0.2e1);
t567 = t382 * mrSges(6,1);
t569 = t376 * mrSges(6,2);
t455 = t567 - t569;
t674 = -t455 - mrSges(5,1);
t591 = mrSges(6,3) * t382;
t154 = -t441 * t591 + t697;
t156 = -mrSges(6,1) * t441 - t272 * t591;
t512 = -t598 / 0.2e1;
t607 = -t361 / 0.2e1;
t663 = t154 * t607 + t156 * t512;
t592 = mrSges(6,3) * t376;
t151 = -t441 * t592 - t696;
t153 = mrSges(6,2) * t441 - t272 * t592;
t511 = t598 / 0.2e1;
t662 = t361 * t151 / 0.2e1 + t153 * t511;
t148 = t454 * t272;
t513 = t599 / 0.2e1;
t606 = t362 / 0.2e1;
t661 = t148 * t513 + t606 * t667;
t618 = -t296 / 0.2e1;
t636 = -t156 / 0.2e1;
t659 = t154 * t618 + t303 * t636;
t614 = t303 / 0.2e1;
t617 = t296 / 0.2e1;
t658 = t151 * t617 + t153 * t614;
t615 = t302 / 0.2e1;
t616 = t297 / 0.2e1;
t657 = t148 * t615 + t616 * t667;
t656 = -t569 / 0.2e1 + t567 / 0.2e1;
t655 = t524 - t533;
t605 = -t376 / 0.2e1;
t580 = Ifges(6,6) * t441;
t75 = t272 * t664 - t580;
t653 = Ifges(5,1) * t619 + Ifges(5,4) * t441 + t77 * t604 + t75 * t605;
t577 = Ifges(6,3) * t272;
t650 = t577 / 0.2e1 + t702;
t246 = (-t603 * t323 + t324 * t374) * pkin(1);
t241 = t246 * mrSges(11,1);
t247 = (-t323 * t374 - t603 * t324) * pkin(1);
t549 = t247 * mrSges(11,2);
t648 = (t374 * mrSges(8,1) + t603 * mrSges(8,2)) * pkin(1) - t241 + t549;
t647 = (m(10) * t304 - mrSges(10,1) * t442 - mrSges(10,2) * t443) * pkin(2) + Ifges(9,4) * t326;
t333 = (mrSges(10,1) * t372 + mrSges(10,2) * t380) * pkin(2);
t514 = qJD(2) + qJD(8);
t37 = t514 * t333;
t308 = t383 * t601 + t351;
t645 = t695 * t308 + (mrSges(4,1) * t384 - mrSges(4,2) * t378) * pkin(1);
t62 = mrSges(11,1) * t181 + mrSges(11,2) * t465;
t262 = Ifges(5,5) * t441;
t644 = t260 / 0.2e1 - t262 / 0.2e1 - t547 / 0.4e1 + t699 / 0.2e1;
t641 = -pkin(11) / 0.2e1;
t640 = m(6) / 0.2e1;
t639 = -mrSges(6,1) / 0.2e1;
t638 = mrSges(6,2) / 0.2e1;
t628 = t241 / 0.2e1;
t627 = t243 / 0.2e1;
t626 = t680 / 0.2e1;
t623 = t443 / 0.2e1;
t307 = -t377 * t601 + t350;
t613 = -t307 / 0.2e1;
t612 = -t308 / 0.2e1;
t611 = t308 / 0.2e1;
t608 = -t338 / 0.4e1;
t45 = Ifges(9,5) * t327 + Ifges(9,6) * t326 + (-t372 * t443 + t380 * t442) * pkin(2) * mrSges(10,3) + t711;
t596 = t45 * qJD(8) + t712;
t359 = pkin(5) * t378 + pkin(1);
t295 = -pkin(5) * t519 + t359 * t379 - pkin(15);
t117 = -pkin(9) * t441 - t272 * pkin(11) + t295;
t347 = t359 * t385;
t349 = pkin(5) * t525;
t436 = t165 + t349;
t119 = t347 + t436;
t159 = mrSges(5,1) * t272 + mrSges(5,2) * t441;
t160 = -mrSges(5,1) * t441 + mrSges(5,2) * t272;
t356 = t379 * pkin(1) - pkin(15);
t225 = (t328 * t369 - t329 * t370) * pkin(4) + t356;
t258 = (-t328 * t370 - t329 * t369) * pkin(4);
t228 = pkin(1) * t385 + t258;
t273 = -mrSges(9,1) * t326 + mrSges(9,2) * t327;
t274 = -mrSges(8,1) * t328 + mrSges(8,2) * t329;
t275 = mrSges(4,1) * t330 + mrSges(4,2) * t331;
t301 = t349 + t347;
t316 = Ifges(9,4) * t327;
t508 = Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1;
t386 = t331 ^ 2 * Ifges(4,4) + (-Ifges(4,4) * t330 + (Ifges(4,1) - Ifges(4,2)) * t331) * t330 + (t441 * t683 + t653) * t441 + (-t441 * t508 + Ifges(5,1) * t686 + (Ifges(6,1) * t546 + t585) * t604 + (-Ifges(6,2) * t531 + t581) * t605 - t531 * t365 + t703) * t272;
t387 = t329 ^ 2 * Ifges(8,4) + (-Ifges(8,4) * t328 + (-Ifges(8,1) + Ifges(8,2)) * t329) * t328 + t688;
t161 = Ifges(10,2) * t442 - t586;
t388 = t161 * t684 + t692 * t623 + t689;
t445 = t376 * t153 + t382 * t156;
t446 = t376 * t152 + t382 * t155;
t474 = m(6) * t516;
t61 = -mrSges(11,1) * t465 + mrSges(11,2) * t181;
t1 = -t388 + (-pkin(15) * mrSges(3,1) - Ifges(3,4) * t385 + (Ifges(3,2) - Ifges(3,1)) * t379 + (-mrSges(4,1) * t331 - mrSges(8,1) * t329 + mrSges(4,2) * t330 - mrSges(8,2) * t328 + (m(8) + m(4)) * t356) * pkin(1)) * t385 + t386 + (t274 + t275) * t356 + t445 * t119 + t387 + (m(11) * t228 + t62) * t225 + (t379 * mrSges(3,2) - t273) * pkin(15) + 0.2e1 * t316 * t610 + t301 * t160 + t228 * t61 + (t119 * t474 + t446) * t117 + t379 ^ 2 * Ifges(3,4) + (m(5) * t301 + t159) * t295 + (-0.2e1 * t698 * t610 - t647) * t326;
t576 = t1 * qJD(1);
t352 = pkin(5) * t526;
t120 = t352 + t436;
t309 = t352 + t349;
t2 = t309 * t160 + t356 * t275 + (m(5) * t309 + t159) * t295 + t445 * t120 + (t120 * t474 + t446) * t117 + t386;
t575 = t2 * qJD(1);
t3 = (t295 * mrSges(5,1) + t76 * t604 + t74 * t605 + t703) * t272 + (t295 * mrSges(5,2) - t702 + (Ifges(5,1) / 0.2e1 - t508) * t272 + t653) * t441 + t445 * t165 + (t376 * t151 + t382 * t154 + t165 * t474) * t117;
t574 = t3 * qJD(1);
t573 = t302 * mrSges(5,1);
t572 = t303 * mrSges(5,2);
t571 = t307 * mrSges(5,1);
t568 = t377 * mrSges(5,1);
t565 = t383 * mrSges(5,2);
t19 = pkin(15) * t273 - t327 * t316 + (t698 * t327 + t647) * t326 + t388;
t562 = qJD(1) * t19;
t11 = t258 * t61 + t356 * t274 + (m(11) * t258 + t62) * t225 + t387;
t560 = t11 * qJD(1);
t150 = t272 * t338;
t523 = t382 * t153;
t532 = t376 * t156;
t17 = (t523 - t532) * t117 + (t335 * t686 - t150 * t604 - t382 * t75 / 0.2e1 + (t77 - t149) * t605) * t272;
t559 = t17 * qJD(1);
t558 = t181 * t680;
t18 = t225 * t62 + t688;
t557 = t18 * qJD(1);
t556 = t465 * t243;
t554 = t191 * t465;
t552 = t193 * t181;
t26 = t161 * t623 + t692 * t684 - t689;
t548 = t26 * qJD(1);
t545 = t297 * t667;
t543 = t302 * t272;
t542 = t302 * t455;
t540 = t303 * t441;
t538 = t307 * t455;
t536 = t362 * t667;
t528 = t377 * t272;
t527 = t377 * t455;
t520 = t383 * t441;
t375 = sin(qJ(6));
t381 = cos(qJ(6));
t166 = (mrSges(7,2) * pkin(14) + Ifges(7,4) * t381) * t381 + (pkin(14) * mrSges(7,1) - Ifges(7,4) * t375 + (Ifges(7,1) - Ifges(7,2)) * t381) * t375;
t515 = t166 * qJD(1);
t510 = t153 * t641;
t509 = pkin(11) * t636;
t504 = -t577 / 0.2e1;
t500 = -t150 / 0.4e1 - t75 / 0.4e1;
t498 = t363 * t687;
t496 = t153 * t618;
t494 = t156 * t618;
t489 = t148 * t613;
t487 = t153 * t607;
t485 = t156 * t607;
t478 = t617 + t607;
t477 = t608 - t664 / 0.4e1;
t471 = t516 * t303;
t470 = t516 * t308;
t461 = Ifges(8,5) * t329 + Ifges(8,6) * t328 + t691;
t457 = t74 / 0.4e1 + t338 * t621;
t456 = t76 / 0.4e1 + t441 * t609;
t428 = pkin(9) * t454;
t451 = -t428 / 0.2e1 + t460;
t399 = mrSges(6,3) * t471 - t542 - t572 - t573;
t33 = m(6) * (t296 * t471 + t297 * t302) + t399;
t389 = (t336 / 0.4e1 + t588 / 0.4e1 + (Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1) * t382) * t441 + pkin(11) * t155 / 0.2e1 - t585 / 0.4e1 + t456;
t411 = (-Ifges(5,6) / 0.2e1 + t335 / 0.4e1) * t272 + Ifges(5,5) * t686;
t392 = t411 + t644;
t393 = (t608 - t365 / 0.4e1) * t441 + t152 * t641 - t581 / 0.4e1 + t457;
t5 = (t393 + t658) * t382 + (t389 + t659) * t376 + t392 + t657;
t450 = t5 * qJD(1) + t33 * qJD(2);
t397 = -t543 / 0.2e1 - t540 / 0.2e1 + t441 * t611 + t272 * t613;
t16 = t489 + (t616 - t362 / 0.2e1) * t667 + (t478 * t152 + t153 * t611) * t382 + (-t478 * t155 + t156 * t612) * t376 + ((t528 / 0.2e1 + t520 / 0.2e1) * pkin(5) + t397) * mrSges(5,3);
t34 = -t674 * t307 + m(6) * (t296 * t470 - t297 * t307) + m(5) * (t302 * t308 + t303 * t307) + t645;
t449 = t16 * qJD(1) + t34 * qJD(2);
t426 = t246 * t632 + t247 * t633;
t21 = ((-t191 / 0.2e1 + t627) * t465 - (t701 - t680 / 0.2e1) * t181 + t426) * mrSges(11,3);
t35 = -m(11) * (t191 * t246 + t193 * t247) + t648;
t448 = t21 * qJD(1) - t35 * qJD(2);
t447 = t32 * qJD(2);
t440 = pkin(11) * t679;
t409 = t477 * t272 + t500;
t396 = (t621 + t686) * Ifges(6,6) + t409;
t145 = t455 * t272;
t425 = t145 * t616 + t498;
t438 = t296 * t679;
t10 = t504 - t272 * t438 + (t119 * t638 + t396 + t496) * t376 + (t119 * t639 + t400 + t494) * t382 + t425;
t422 = t297 * t454;
t130 = t422 + t460;
t439 = t10 * qJD(1) + t130 * qJD(2);
t437 = t361 * t679;
t435 = -pkin(9) * t145 / 0.2e1 + t498;
t432 = -t566 / 0.2e1 - t570 / 0.2e1;
t431 = t580 / 0.4e1 + t500;
t424 = t145 * t606 + t498;
t421 = t362 * t454;
t419 = t428 / 0.2e1;
t418 = t432 * t303;
t417 = t432 * t308;
t416 = -t422 / 0.2e1;
t398 = mrSges(6,3) * t704 - t527 - t565 - t568;
t131 = (m(6) * (t361 * t704 + t362 * t377) + t398) * pkin(5);
t407 = m(6) * (pkin(9) * t307 + pkin(11) * t470);
t25 = (t707 * pkin(5) + t693) * t640 - t407 / 0.2e1 + t674 * (t615 + t513 + t307 / 0.2e1) + t695 * (t511 + t614 + t612);
t7 = (t393 + t662) * t382 + (t389 + t663) * t376 + t392 + t661;
t414 = t7 * qJD(1) + t25 * qJD(2) + t131 * qJD(3);
t22 = t628 + (t701 + t626) * mrSges(11,1) + (-t247 / 0.2e1 + t191 / 0.2e1 + t627) * mrSges(11,2);
t394 = (t181 * t626 + t632 * t680) * mrSges(11,3) - t652;
t30 = t394 + t652;
t413 = t30 * qJD(1) - t22 * qJD(2) + t38 * qJD(7);
t13 = t504 - t272 * t437 + (t120 * t638 + t396 + t487) * t376 + (t120 * t639 + t400 + t485) * t382 + t424;
t157 = t421 + t460;
t390 = -t421 / 0.2e1 - t460;
t57 = t416 + t417 + t390;
t412 = t13 * qJD(1) - t57 * qJD(2) + t157 * qJD(3);
t410 = t432 * t598;
t401 = t382 * t74;
t402 = t376 * t76;
t406 = t547 / 0.2e1 - t260 + t262 + Ifges(4,5) * t331 + t402 / 0.2e1 + t401 / 0.2e1 - Ifges(4,6) * t330 + t423 * t441;
t15 = (-Ifges(6,3) / 0.2e1 - t440) * t272 + (t165 * t639 + t400 + t509) * t382 + (0.3e1 / 0.4e1 * t580 + t510 + t165 * t638 + t409) * t376 + t435;
t167 = -t428 + t460;
t59 = t416 + t419 + t418 - t460;
t90 = t419 + t410 + t390;
t404 = t15 * qJD(1) - t59 * qJD(2) - t90 * qJD(3) + t167 * qJD(4);
t403 = t477 * t376 + t476 * t382;
t391 = t530 * t687 + t521 * t621 + t402 / 0.4e1 + t401 / 0.4e1 + t411 - t644 + t694 * pkin(11);
t325 = t333 * qJD(9);
t300 = t421 / 0.2e1;
t276 = t422 / 0.2e1;
t91 = t300 + t410 + t451;
t60 = t276 + t418 + t451;
t58 = t276 + t300 + t417 + t460;
t36 = t325 + t37;
t29 = t394 - t652;
t27 = -0.2e1 * t652;
t24 = t693 * t640 - t542 / 0.2e1 - t572 / 0.2e1 - t573 / 0.2e1 + t407 / 0.2e1 + t538 / 0.2e1 + mrSges(5,2) * t612 + t571 / 0.2e1 + t303 * t679 + (t707 * t640 - t527 / 0.2e1 - t568 / 0.2e1 - t565 / 0.2e1 + t383 * t679) * pkin(5) + t611 * t646;
t23 = -t553 / 0.2e1 - t555 / 0.2e1 - t551 / 0.2e1 - t550 / 0.2e1 + t628 - t549 / 0.2e1;
t20 = t461 + (-t556 / 0.2e1 - t558 / 0.2e1 - t554 / 0.2e1 - t552 / 0.2e1 + t426) * mrSges(11,3);
t14 = Ifges(6,3) * t619 + (t510 + t431) * t376 + (t509 + t501) * t382 + (-t440 + t403) * t272 + t435 + t702 + t656 * t165;
t12 = (t487 + t431) * t376 + (t485 + t501) * t382 + (-t437 + t403) * t272 + t424 + t656 * t120 + t650;
t9 = (t496 + t431) * t376 + (t494 + t501) * t382 + (-t438 + t403) * t272 + t425 + t656 * t119 + t650;
t8 = t536 / 0.2e1 + t489 + t545 / 0.2e1 + t406 + (t523 / 0.2e1 - t532 / 0.2e1) * t308 + (-t272 * t513 + t441 * t512 + t397) * mrSges(5,3) + (t296 + t361) * t694;
t6 = t391 + (t456 + t663) * t376 + (t457 + t662) * t382 + t661;
t4 = t391 + (t457 + t658) * t382 + (t456 + t659) * t376 + t657;
t28 = [qJD(10) * t18 + qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t3 + qJD(5) * t17 + qJD(6) * t166 + qJD(7) * t11 - qJD(8) * t19 + qJD(9) * t26, t27 * qJD(10) + t8 * qJD(3) + t4 * qJD(4) + t9 * qJD(5) + t20 * qJD(7) + t576 + t596 + (-Ifges(3,5) * t379 - Ifges(3,6) * t385 + t406 + t45 + t461 + t545 + ((t328 * t374 - t603 * t329) * mrSges(8,3) + (t330 * t384 - t331 * t378) * mrSges(4,3)) * pkin(1) + t655 * t296 + (-t540 - t543) * mrSges(5,3) + (-t552 - t554) * mrSges(11,3)) * qJD(2), t575 + t8 * qJD(2) + (t536 + t655 * t361 + (-t520 - t528) * mrSges(5,3) * pkin(5) + t406) * qJD(3) + t6 * qJD(4) + t12 * qJD(5), t4 * qJD(2) + t6 * qJD(3) + t14 * qJD(5) + t574 + (t335 * t619 + t76 * t706 + t74 * t604 - t699 - t260 + (Ifges(5,5) + t423) * t441 + (t382 * t151 - t376 * t154) * pkin(11)) * qJD(4), t559 + t9 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + (-t117 * t454 - t547) * qJD(5), t515 + (Ifges(7,5) * t381 - Ifges(7,6) * t375) * qJD(6), t560 + t20 * qJD(2) + ((-t556 - t558) * mrSges(11,3) + t461) * qJD(7) + t29 * qJD(10), qJD(2) * t45 - t562 + t596, t514 * t711 + t548 + t712, qJD(10) * t691 + t27 * qJD(2) + t29 * qJD(7) + t557, 0, 0, 0; qJD(3) * t16 + qJD(4) * t5 + qJD(5) * t10 + qJD(7) * t21 - t576, qJD(3) * t34 + qJD(4) * t33 + qJD(5) * t130 - qJD(7) * t35 + t325 + t710, t24 * qJD(4) + t58 * qJD(5) + t449 + (t538 + t571 + 0.2e1 * (-t362 * t307 + t308 * t469) * t640 + m(5) * (t307 * t383 + t308 * t377) * pkin(5) + t645) * qJD(3), t24 * qJD(3) + (m(6) * (-pkin(9) * t302 + pkin(11) * t471) + t399) * qJD(4) + t60 * qJD(5) + t450, t58 * qJD(3) + t60 * qJD(4) + (-t296 * t455 + t665) * qJD(5) + t439, 0, (m(11) * (t243 * t246 + t247 * t680) - t648) * qJD(7) + t23 * qJD(10) + t448, t325, t36, t23 * qJD(7) + t447 + t710, 0, 0, 0; -qJD(2) * t16 + qJD(4) * t7 + qJD(5) * t13 - t575, qJD(4) * t25 - qJD(5) * t57 - t449, qJD(4) * t131 + qJD(5) * t157, t91 * qJD(5) + (m(6) * (-pkin(9) * t377 + pkin(11) * t704) + t398) * qJD(4) * pkin(5) + t414, t91 * qJD(4) + (-t361 * t455 + t665) * qJD(5) + t412, 0, 0, 0, 0, 0, 0, 0, 0; -qJD(2) * t5 - qJD(3) * t7 + qJD(5) * t15 - t574, -qJD(3) * t25 - qJD(5) * t59 - t450, -qJD(5) * t90 - t414, t167 * qJD(5), (-pkin(11) * t455 + t665) * qJD(5) + t404, 0, 0, 0, 0, 0, 0, 0, 0; -qJD(2) * t10 - qJD(3) * t13 - qJD(4) * t15 - t559, qJD(3) * t57 + qJD(4) * t59 - t439, qJD(4) * t90 - t412, -t404, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t515, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; qJD(10) * t30 - qJD(2) * t21 - t560, -qJD(10) * t22 - t448, 0, 0, 0, 0, t709, 0, 0, t413 + t709, 0, 0, 0; t562, t325, 0, 0, 0, 0, 0, t325, t36, 0, 0, 0, 0; -t548, -t37, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0; -qJD(7) * t30 - t557, qJD(7) * t22 - t447, 0, 0, 0, 0, -t413, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t28;
