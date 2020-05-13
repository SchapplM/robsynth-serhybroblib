% Calculate vector of inverse dynamics joint torques with ic for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m1IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1IC_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:01:34
% EndTime: 2020-04-15 20:02:25
% DurationCPUTime: 48.11s
% Computational Cost: add. (12924->982), mult. (29334->1389), div. (20->8), fcn. (23331->48), ass. (0->483)
t726 = pkin(2) * m(10);
t786 = mrSges(9,1) + t726;
t431 = qJ(2) + qJ(8);
t412 = sin(t431);
t416 = cos(t431);
t430 = qJ(8) + qJ(9);
t419 = qJ(2) + t430;
t382 = sin(t419);
t384 = cos(t419);
t595 = t382 * mrSges(10,1) + t384 * mrSges(10,2);
t785 = mrSges(9,2) * t416 + t786 * t412 - t595;
t659 = mrSges(10,2) * t382;
t784 = -mrSges(9,2) * t412 + t786 * t416 + t659;
t443 = sin(qJ(2));
t452 = cos(qJ(2));
t341 = -mrSges(3,1) * t443 - mrSges(3,2) * t452;
t428 = -qJ(7) + pkin(19);
t410 = -qJ(2) + t428;
t353 = pkin(4) * sin(t410);
t421 = t443 * pkin(1);
t433 = qJ(2) + qJ(3);
t414 = sin(t433);
t418 = cos(t433);
t599 = qJ(3) + qJ(4);
t420 = qJ(2) + t599;
t383 = sin(t420);
t385 = cos(t420);
t547 = -t385 * mrSges(5,1) + t383 * mrSges(5,2);
t764 = -t418 * mrSges(4,1) + mrSges(4,2) * t414 + t547;
t432 = qJ(2) + qJ(7);
t413 = sin(t432);
t417 = cos(t432);
t407 = -qJ(10) + t428;
t379 = -qJ(2) + t407;
t354 = sin(t379);
t355 = cos(t379);
t544 = -t354 * mrSges(11,1) + t355 * mrSges(11,2);
t765 = -mrSges(8,1) * t413 - mrSges(8,2) * t417 + t544;
t783 = -m(11) * (t353 - t421) - t341 + t764 - t765 + t785;
t775 = -m(8) - m(4);
t782 = m(9) + m(3);
t781 = t452 / 0.2e1;
t578 = qJD(3) + qJD(4);
t406 = qJD(2) + t578;
t440 = sin(qJ(5));
t449 = cos(qJ(5));
t442 = sin(qJ(3));
t451 = cos(qJ(3));
t309 = t442 * t452 + t443 * t451;
t286 = t309 * qJD(1);
t312 = -t442 * t443 + t451 * t452;
t289 = t312 * qJD(1);
t441 = sin(qJ(4));
t450 = cos(qJ(4));
t511 = t450 * t286 + t289 * t441;
t151 = t406 * t449 - t440 * t511;
t152 = t406 * t440 + t449 * t511;
t772 = -mrSges(5,1) * t406 - mrSges(6,1) * t151 + mrSges(6,2) * t152 + mrSges(5,3) * t511;
t662 = mrSges(5,2) * t385;
t780 = mrSges(4,1) * t414 + mrSges(5,1) * t383 + mrSges(4,2) * t418 + t662;
t661 = mrSges(6,2) * t440;
t779 = t383 * t661 + t385 * (pkin(11) * m(6) + mrSges(6,3));
t541 = -t286 * t441 + t450 * t289;
t192 = Ifges(5,4) * t541;
t636 = t406 * Ifges(5,5);
t102 = Ifges(5,1) * t511 + t192 + t636;
t592 = qJD(1) * t443;
t619 = qJD(1) * pkin(15);
t343 = pkin(1) * t592 - t619;
t234 = -pkin(5) * t289 + t343;
t193 = qJD(5) - t541;
t646 = t152 * Ifges(6,4);
t48 = t151 * Ifges(6,2) + t193 * Ifges(6,6) + t646;
t670 = pkin(1) * qJD(2);
t570 = t451 * t670;
t348 = t441 * t570;
t427 = qJD(2) + qJD(3);
t572 = t442 * t670;
t500 = pkin(5) * t427 + t572;
t236 = t450 * t500 + t348;
t230 = -t406 * pkin(9) - t236;
t529 = mrSges(6,1) * t440 + mrSges(6,2) * t449;
t488 = t230 * t529;
t150 = Ifges(6,4) * t151;
t49 = t152 * Ifges(6,1) + t193 * Ifges(6,5) + t150;
t626 = t449 * t49;
t642 = t236 * mrSges(5,3);
t698 = t440 / 0.2e1;
t778 = -t626 / 0.2e1 + t48 * t698 + t642 - t102 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t511 - t192 / 0.2e1 - t488 - t234 * mrSges(5,2) - t636 / 0.2e1;
t635 = t406 * Ifges(5,6);
t652 = Ifges(5,4) * t511;
t100 = Ifges(5,2) * t541 + t635 + t652;
t644 = t193 * Ifges(6,3);
t645 = t152 * Ifges(6,5);
t647 = t151 * Ifges(6,6);
t47 = t644 + t645 + t647;
t505 = -t441 * t500 + t450 * t570;
t641 = t505 * mrSges(5,3);
t231 = pkin(11) * t406 - t505;
t83 = -pkin(9) * t541 - pkin(11) * t511 + t234;
t60 = t231 * t449 + t440 * t83;
t739 = t60 * mrSges(6,2);
t59 = -t231 * t440 + t449 * t83;
t740 = t59 * mrSges(6,1);
t776 = t739 + t100 / 0.2e1 - t47 / 0.2e1 - t641 + t652 / 0.2e1 - t234 * mrSges(5,1) + t635 / 0.2e1 - t740;
t438 = sin(qJ(7));
t447 = cos(qJ(7));
t308 = -t438 * t452 - t443 * t447;
t285 = t308 * qJD(1);
t591 = qJD(1) * t452;
t288 = -t438 * t592 + t447 * t591;
t434 = sin(pkin(19));
t435 = cos(pkin(19));
t160 = (-t285 * t435 - t288 * t434) * pkin(4) + t343;
t671 = sin(qJ(10));
t672 = cos(qJ(10));
t299 = t434 * t671 + t435 * t672;
t468 = -t434 * t672 + t435 * t671;
t129 = -t285 * t468 - t288 * t299;
t542 = -t299 * t285 + t288 * t468;
t72 = -mrSges(11,1) * t542 + mrSges(11,2) * t129;
t771 = m(11) * t160 + t72;
t770 = -mrSges(4,1) * t289 - mrSges(8,1) * t285 + mrSges(4,2) * t286 + mrSges(8,2) * t288;
t444 = sin(qJ(1));
t664 = mrSges(10,1) * t384;
t331 = t444 * t664;
t769 = t444 * t784 - t331;
t453 = cos(qJ(1));
t334 = t453 * t664;
t768 = t453 * t784 - t334;
t762 = mrSges(11,1) * t355 + mrSges(11,2) * t354;
t598 = t762 * t444;
t660 = mrSges(8,2) * t413;
t767 = -t444 * t660 - t598;
t597 = t762 * t453;
t766 = -t453 * t660 - t597;
t666 = mrSges(6,1) * t449;
t761 = t661 - t666;
t760 = g(1) * t453 + g(2) * t444;
t109 = -mrSges(5,1) * t541 + mrSges(5,2) * t511;
t759 = m(5) * t234 + t109;
t532 = mrSges(3,1) * t452 - mrSges(3,2) * t443;
t653 = Ifges(3,4) * t452;
t758 = (-Ifges(3,1) * t443 - t653) * t781 - pkin(15) * t532;
t170 = -mrSges(5,2) * t406 + mrSges(5,3) * t541;
t92 = -mrSges(6,2) * t193 + mrSges(6,3) * t151;
t93 = mrSges(6,1) * t193 - mrSges(6,3) * t152;
t757 = -t440 * t93 + t449 * t92 + t170;
t517 = t440 * t60 + t449 * t59;
t756 = -m(6) * t517 - t440 * t92 - t449 * t93;
t424 = qJDD(2) + qJDD(3);
t401 = qJDD(4) + t424;
t581 = qJD(1) * qJD(2);
t325 = qJDD(1) * t452 - t443 * t581;
t326 = -qJDD(1) * t443 - t452 * t581;
t480 = t309 * qJD(3);
t165 = -qJD(1) * t480 + t325 * t451 + t326 * t442;
t479 = t312 * qJD(3);
t166 = qJD(1) * t479 + t325 * t442 - t326 * t451;
t57 = qJD(4) * t541 + t165 * t441 + t166 * t450;
t19 = qJD(5) * t151 + t401 * t440 + t449 * t57;
t58 = -qJD(4) * t511 + t165 * t450 - t166 * t441;
t54 = qJDD(5) - t58;
t11 = mrSges(6,1) * t54 - mrSges(6,3) * t19;
t20 = -qJD(5) * t152 + t401 * t449 - t440 * t57;
t12 = -mrSges(6,2) * t54 + mrSges(6,3) * t20;
t585 = qJD(5) * t449;
t586 = qJD(5) * t440;
t754 = -t440 * t11 + t449 * t12 - t93 * t585 - t92 * t586;
t685 = pkin(5) * t414;
t692 = pkin(1) * t452;
t329 = -t685 - t692;
t564 = t383 * t666;
t665 = mrSges(8,1) * t417;
t680 = pkin(9) * t383;
t689 = pkin(4) * cos(t410);
t409 = pkin(20) + t432;
t374 = sin(t409);
t377 = cos(t409);
t439 = sin(qJ(6));
t448 = cos(qJ(6));
t235 = 0.1e1 / (-t374 * t439 - t377 * t448) / pkin(7) / pkin(3);
t316 = pkin(3) * t374 + t421;
t317 = pkin(3) * t377 + t692;
t736 = (t316 * t439 + t317 * t448) * pkin(7) * t235;
t753 = (m(11) * t689 + t665) * t736 - m(6) * (t329 - t680) + t564 + m(8) * t692 + t665 - m(11) * (-t689 - t692);
t693 = pkin(1) * t451;
t568 = qJD(3) * t693;
t623 = pkin(1) * qJDD(2);
t303 = qJD(2) * t568 + t442 * t623;
t253 = pkin(5) * t424 + t303;
t589 = qJD(3) * t442;
t302 = (qJD(2) * t589 - qJDD(2) * t451) * pkin(1);
t119 = qJD(4) * t236 + t441 * t253 + t450 * t302;
t113 = pkin(11) * t401 + t119;
t582 = qJDD(1) * pkin(15);
t266 = -pkin(1) * t326 - t582;
t117 = -pkin(5) * t165 + t266;
t13 = -pkin(9) * t58 - pkin(11) * t57 + t117;
t5 = -qJD(5) * t60 - t113 * t440 + t13 * t449;
t675 = t440 * t5;
t463 = -qJD(5) * t517 - t675;
t4 = qJD(5) * t59 + t113 * t449 + t13 * t440;
t676 = t4 * t449;
t751 = m(6) * (t463 + t676) + t754;
t749 = pkin(3) * t235 * (t316 * t377 - t317 * t374);
t426 = qJD(2) + qJD(7);
t402 = qJD(10) + t426;
t701 = -t402 / 0.2e1;
t715 = -t129 / 0.2e1;
t748 = -t160 * mrSges(11,2) + Ifges(11,1) * t715 + Ifges(11,5) * t701;
t437 = sin(qJ(8));
t446 = cos(qJ(8));
t307 = -t437 * t452 - t443 * t446;
t284 = t307 * qJD(1);
t287 = -t437 * t592 + t446 * t591;
t436 = sin(qJ(9));
t445 = cos(qJ(9));
t513 = -t284 * t445 + t287 * t436;
t191 = Ifges(10,4) * t513;
t425 = qJD(2) + qJD(8);
t405 = qJD(9) + t425;
t512 = t284 * t436 + t287 * t445;
t101 = -Ifges(10,1) * t512 + t405 * Ifges(10,5) + t191;
t250 = -pkin(2) * t284 - t619;
t422 = qJDD(2) + qJDD(8);
t584 = qJD(9) * t436;
t262 = (-t422 * t445 + t425 * t584) * pkin(2);
t583 = qJD(9) * t445;
t263 = (-t422 * t436 - t425 * t583) * pkin(2);
t400 = qJDD(9) + t422;
t484 = t307 * qJD(8);
t161 = qJD(1) * t484 + t325 * t446 + t326 * t437;
t508 = t437 * t443 - t446 * t452;
t483 = t508 * qJD(8);
t162 = qJD(1) * t483 - t325 * t437 + t326 * t446;
t55 = qJD(9) * t513 - t161 * t445 - t162 * t436;
t56 = qJD(9) * t512 + t161 * t436 - t162 * t445;
t643 = t512 * Ifges(10,4);
t709 = -t512 / 0.2e1;
t99 = Ifges(10,2) * t513 + t405 * Ifges(10,6) - t643;
t747 = t262 * mrSges(10,1) - t263 * mrSges(10,2) + Ifges(10,5) * t55 + Ifges(10,6) * t56 + Ifges(10,3) * t400 - (Ifges(10,5) * t513 + Ifges(10,6) * t512) * t405 / 0.2e1 - t250 * (-mrSges(10,1) * t512 + mrSges(10,2) * t513) + t99 * t709 + (Ifges(10,1) * t513 + t643) * t512 / 0.2e1 - (Ifges(10,2) * t512 + t101 + t191) * t513 / 0.2e1;
t124 = Ifges(11,4) * t542;
t482 = t308 * qJD(7);
t163 = qJD(1) * t482 + t325 * t447 + t326 * t438;
t507 = t438 * t443 - t447 * t452;
t481 = t507 * qJD(7);
t164 = qJD(1) * t481 - t325 * t438 + t326 * t447;
t272 = t299 * qJD(10);
t273 = t468 * qJD(10);
t27 = -t163 * t299 - t164 * t468 - t272 * t285 + t273 * t288;
t28 = t163 * t468 - t164 * t299 + t272 * t288 + t273 * t285;
t423 = qJDD(2) + qJDD(7);
t393 = qJDD(10) + t423;
t573 = t438 * t670;
t688 = pkin(4) * t434;
t295 = t426 * t688 + t573;
t571 = t447 * t670;
t687 = pkin(4) * t435;
t296 = t426 * t687 + t571;
t131 = t295 * t299 + t296 * t468;
t617 = t131 * mrSges(11,3);
t669 = pkin(1) * qJD(7);
t563 = qJD(2) * t669;
t300 = -t438 * t563 + t447 * t623;
t248 = t423 * t687 + t300;
t301 = t438 * t623 + t447 * t563;
t249 = t423 * t688 + t301;
t62 = -t248 * t468 - t249 * t299 - t272 * t296 + t273 * t295;
t63 = -t248 * t299 + t249 * t468 + t272 * t295 + t273 * t296;
t65 = t129 * Ifges(11,4) + Ifges(11,2) * t542 + t402 * Ifges(11,6);
t66 = t129 * Ifges(11,1) + t402 * Ifges(11,5) + t124;
t716 = -t542 / 0.2e1;
t745 = t63 * mrSges(11,1) - t62 * mrSges(11,2) + Ifges(11,5) * t27 + Ifges(11,6) * t28 + Ifges(11,3) * t393 + (t124 + t66) * t716 - (t160 * mrSges(11,1) + Ifges(11,4) * t715 + Ifges(11,2) * t716 + Ifges(11,6) * t701 + t617 - t65 / 0.2e1) * t129;
t10 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t120 = qJD(4) * t505 + t450 * t253 - t441 * t302;
t114 = -t401 * pkin(9) - t120;
t519 = Ifges(6,5) * t449 - Ifges(6,6) * t440;
t489 = t193 * t519;
t651 = Ifges(6,4) * t440;
t524 = Ifges(6,1) * t449 - t651;
t490 = t152 * t524;
t650 = Ifges(6,4) * t449;
t522 = -Ifges(6,2) * t440 + t650;
t491 = t151 * t522;
t551 = -t586 / 0.2e1;
t556 = t626 / 0.2e1;
t657 = mrSges(6,3) * t449;
t720 = t54 / 0.2e1;
t723 = t20 / 0.2e1;
t724 = t19 / 0.2e1;
t8 = t19 * Ifges(6,4) + t20 * Ifges(6,2) + t54 * Ifges(6,6);
t9 = t19 * Ifges(6,1) + t20 * Ifges(6,4) + t54 * Ifges(6,5);
t458 = -t119 * mrSges(5,2) + t4 * t657 + t449 * t8 / 0.2e1 + t9 * t698 + t114 * t761 + t120 * mrSges(5,1) + (Ifges(6,1) * t440 + t650) * t724 + (Ifges(6,2) * t449 + t651) * t723 + Ifges(5,3) * t401 + (Ifges(6,5) * t440 + Ifges(6,6) * t449) * t720 + t48 * t551 + Ifges(5,6) * t58 + Ifges(5,5) * t57 + (t488 + t556) * qJD(5) + (t491 + t490 + t489) * qJD(5) / 0.2e1;
t725 = m(6) * pkin(9);
t461 = t662 + (mrSges(5,1) + t666 + t725) * t383;
t357 = t383 * mrSges(6,3);
t503 = -t385 * t761 + t357;
t533 = t779 * t453;
t534 = t779 * t444;
t594 = t385 * pkin(9) + t383 * pkin(11);
t733 = -t59 * t585 - t60 * t586;
t111 = pkin(9) * t511 - pkin(11) * t541;
t81 = t111 * t449 - t236 * t440;
t82 = t111 * t440 + t236 * t449;
t744 = pkin(8) * (-pkin(9) * t10 + (t453 * t461 - t533) * g(1) + (t444 * t461 - t534) * g(2) + (-t647 / 0.2e1 - t645 / 0.2e1 - t644 / 0.2e1 + t776) * t511 + (-t491 / 0.2e1 - t490 / 0.2e1 - t489 / 0.2e1 + t517 * mrSges(6,3) + t778) * t541 + t463 * mrSges(6,3) + (m(6) * (-t675 + t676 + t733) + t754) * pkin(11) - t114 * t725 + (-m(6) * t594 - t503 + t547) * g(3) - m(6) * (-t230 * t505 + t59 * t81 + t60 * t82) + t458 - t82 * t92 - t81 * t93 - t236 * t170 + t772 * t505);
t618 = t542 * mrSges(11,3);
t115 = -mrSges(11,2) * t402 + t618;
t116 = mrSges(11,1) * t402 - mrSges(11,3) * t129;
t130 = t295 * t468 - t296 * t299;
t743 = pkin(10) * (-t131 * t116 - g(3) * t544 - g(1) * t597 - g(2) * t598 + (-t115 + t618) * t130 + t748 * t542 + t745);
t699 = t439 / 0.2e1;
t697 = -t443 / 0.2e1;
t492 = pkin(1) * (-t299 * t447 + t438 * t468);
t735 = (-t272 * t435 + t273 * t434) * pkin(4) - qJD(2) * t492;
t493 = pkin(1) * (t299 * t438 + t447 * t468);
t734 = (t272 * t434 + t273 * t435) * pkin(4) - qJD(2) * t493;
t220 = t309 * t450 + t312 * t441;
t218 = t309 * t441 - t450 * t312;
t228 = -qJD(2) * t309 - t480;
t229 = qJD(2) * t312 + t479;
t88 = -qJD(4) * t218 + t228 * t441 + t229 * t450;
t498 = t220 * t585 + t440 * t88;
t340 = -mrSges(7,1) * t448 + mrSges(7,2) * t439;
t477 = m(7) * pkin(14) + t340;
t731 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t730 = -mrSges(9,3) - mrSges(11,3) - mrSges(5,3) - mrSges(4,3) - mrSges(7,3) - mrSges(3,3) - mrSges(10,3) - mrSges(8,3) + mrSges(2,2);
t380 = pkin(5) * t418;
t593 = t380 - t421;
t318 = pkin(15) + t593;
t381 = -t421 + pkin(15);
t729 = m(5) * t318 - t477 + m(6) * (t318 + t594) + mrSges(2,1) + t357 - t775 * t381 + (m(10) + m(11) + t782) * pkin(15) - t783;
t714 = -t151 / 0.2e1;
t713 = -t152 / 0.2e1;
t712 = t152 / 0.2e1;
t711 = -t193 / 0.2e1;
t705 = t286 / 0.2e1;
t704 = t287 / 0.2e1;
t703 = t288 / 0.2e1;
t696 = pkin(1) * t438;
t695 = pkin(1) * t442;
t694 = pkin(1) * t447;
t691 = pkin(2) * (-mrSges(10,1) * t513 - mrSges(10,2) * t512);
t686 = pkin(5) * t286;
t684 = pkin(5) * t441;
t683 = pkin(5) * t450;
t673 = -qJD(6) / 0.2e1;
t658 = mrSges(6,3) * t440;
t656 = mrSges(8,3) * t288;
t655 = mrSges(10,3) * t425;
t654 = Ifges(3,4) * t443;
t649 = Ifges(7,4) * t439;
t648 = Ifges(7,4) * t448;
t640 = t286 * Ifges(4,4);
t639 = t287 * Ifges(9,4);
t638 = t288 * Ifges(8,4);
t634 = t436 * (-mrSges(10,2) * t400 + mrSges(10,3) * t56);
t633 = t439 * Ifges(7,1);
t628 = t448 * Ifges(7,2);
t620 = mrSges(11,3) * t130;
t616 = t220 * t440;
t615 = t220 * t449;
t614 = (-mrSges(8,2) * t426 + mrSges(8,3) * t285) * t447;
t521 = t628 + t649;
t607 = t439 * (Ifges(7,6) * qJD(6) + qJD(1) * t521);
t606 = t440 * t444;
t605 = t440 * t453;
t604 = t441 * t451;
t603 = t444 * t449;
t392 = qJD(1) * t648;
t602 = t448 * (Ifges(7,5) * qJD(6) + qJD(1) * t633 + t392);
t601 = t449 * t453;
t600 = t450 * t451;
t590 = qJD(2) * t452;
t588 = qJD(4) * t441;
t587 = qJD(4) * t450;
t580 = qJD(1) * qJD(6);
t577 = t250 * t726;
t576 = pkin(2) * t655;
t575 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t54;
t574 = pkin(1) * t600;
t569 = pkin(1) * t589;
t396 = pkin(1) * t590;
t555 = t380 + t594;
t206 = -pkin(5) * t228 + t396;
t550 = -t581 / 0.2e1;
t546 = t436 * t576;
t545 = t445 * t576;
t540 = mrSges(4,3) * t572;
t539 = mrSges(4,3) * t570;
t538 = mrSges(8,3) * t573;
t537 = mrSges(8,3) * t571;
t251 = -pkin(5) * t312 - t381;
t387 = pkin(5) + t695;
t276 = t441 * t387 - t574;
t528 = mrSges(7,1) * t439 + mrSges(7,2) * t448;
t525 = t452 * Ifges(3,1) - t654;
t523 = -t443 * Ifges(3,2) + t653;
t520 = -Ifges(3,5) * t443 - Ifges(3,6) * t452;
t518 = Ifges(7,5) * t448 - Ifges(7,6) * t439;
t516 = -t440 * t59 + t449 * t60;
t217 = -t307 * t445 - t436 * t508;
t510 = t307 * t436 - t445 * t508;
t275 = pkin(1) * t604 + t387 * t450;
t497 = t220 * t586 - t449 * t88;
t495 = pkin(14) * t528;
t487 = t439 * (Ifges(7,1) * t448 - t649);
t486 = t443 * (-Ifges(3,2) * t452 - t654);
t180 = (-t285 * t434 + t288 * t435) * pkin(4);
t91 = t111 + t686;
t462 = m(6) * (-t680 - t685) - t564;
t171 = mrSges(10,1) * t405 + mrSges(10,3) * t512;
t185 = t284 * Ifges(9,2) + t425 * Ifges(9,6) + t639;
t269 = Ifges(9,4) * t284;
t188 = t287 * Ifges(9,1) + t425 * Ifges(9,5) + t269;
t457 = -t287 * t691 + t185 * t704 + (mrSges(9,1) * t287 + mrSges(9,2) * t284) * t619 + t512 * t546 - t287 * (Ifges(9,1) * t284 - t639) / 0.2e1 - t287 * t577 + pkin(2) * t171 * t584 + (-t262 * t445 - t263 * t436) * t726 + Ifges(9,5) * t161 + Ifges(9,6) * t162 + Ifges(9,3) * t422 - t425 * (Ifges(9,5) * t284 - Ifges(9,6) * t287) / 0.2e1 - (-Ifges(9,2) * t287 + t188 + t269) * t284 / 0.2e1 + t747;
t186 = t285 * Ifges(8,2) + t426 * Ifges(8,6) + t638;
t270 = Ifges(8,4) * t285;
t189 = t288 * Ifges(8,1) + t426 * Ifges(8,5) + t270;
t456 = -t301 * mrSges(8,2) - t343 * (mrSges(8,1) * t288 + mrSges(8,2) * t285) + t186 * t703 - t288 * (Ifges(8,1) * t285 - t638) / 0.2e1 + Ifges(8,6) * t164 + Ifges(8,5) * t163 - t426 * (Ifges(8,5) * t285 - Ifges(8,6) * t288) / 0.2e1 + t285 * t537 + t300 * mrSges(8,1) + Ifges(8,3) * t423 - (-Ifges(8,2) * t288 + t189 + t270) * t285 / 0.2e1 + (t620 + t748) * t542 + t745;
t187 = t289 * Ifges(4,2) + t427 * Ifges(4,6) + t640;
t271 = Ifges(4,4) * t289;
t190 = t286 * Ifges(4,1) + t427 * Ifges(4,5) + t271;
t455 = t187 * t705 - t5 * t658 + t733 * mrSges(6,3) - t286 * (Ifges(4,1) * t289 - t640) / 0.2e1 - t286 * t539 + t289 * t540 + t458 - t343 * (mrSges(4,1) * t286 + mrSges(4,2) * t289) + Ifges(4,6) * t165 + Ifges(4,5) * t166 - t302 * mrSges(4,2) + t303 * mrSges(4,1) + Ifges(4,3) * t424 - t427 * (Ifges(4,5) * t289 - Ifges(4,6) * t286) / 0.2e1 - (-Ifges(4,2) * t286 + t190 + t271) * t289 / 0.2e1 + (Ifges(6,5) * t713 + Ifges(6,6) * t714 + Ifges(6,3) * t711 + t776) * t511 + (t519 * t711 + t522 * t714 + t524 * t713 + t59 * t657 + t60 * t658 + t778) * t541;
t429 = qJ(3) + pkin(17);
t415 = cos(t430);
t411 = sin(t430);
t408 = pkin(18) + t599;
t404 = cos(t429);
t403 = sin(t429);
t388 = -pkin(9) - t683;
t376 = cos(t408);
t373 = sin(t408);
t371 = cos(t407);
t370 = sin(t407);
t342 = t687 + t694;
t339 = t688 + t696;
t328 = -pkin(2) * t446 + pkin(12) * t415;
t327 = pkin(2) * t437 - pkin(12) * t411;
t324 = qJDD(1) * t439 + t448 * t580;
t323 = qJDD(1) * t448 - t439 * t580;
t314 = pkin(5) * t451 + pkin(10) * t376;
t313 = pkin(5) * t442 + pkin(10) * t373;
t293 = t371 * pkin(8) - pkin(4) * cos(t428);
t292 = t370 * pkin(8) - pkin(4) * sin(t428);
t283 = Ifges(3,5) * qJD(2) + qJD(1) * t525;
t281 = Ifges(3,6) * qJD(2) + qJD(1) * t523;
t268 = t450 * t572 + t348;
t267 = (t441 * t442 - t600) * t670;
t261 = t385 * t601 + t606;
t260 = -t385 * t605 + t603;
t259 = -t385 * t603 + t605;
t258 = t385 * t606 + t601;
t242 = -mrSges(4,2) * t427 + mrSges(4,3) * t289;
t241 = mrSges(8,1) * t426 - t656;
t240 = mrSges(4,1) * t427 - mrSges(4,3) * t286;
t227 = qJD(2) * t507 + t481;
t226 = qJD(2) * t308 + t482;
t225 = qJD(2) * t508 + t483;
t224 = qJD(2) * t307 + t484;
t221 = 0.1e1 / (-t370 * t373 + t371 * t376) / pkin(10) / pkin(8);
t208 = (-t299 * t434 - t435 * t468) * pkin(4);
t207 = (-t299 * t435 + t434 * t468) * pkin(4);
t169 = -mrSges(10,2) * t405 + mrSges(10,3) * t513;
t168 = -t299 * t339 - t342 * t468;
t167 = -t299 * t342 + t339 * t468;
t141 = -pkin(2) * t162 - t582;
t95 = qJD(7) * t493 + t272 * t339 + t273 * t342;
t94 = qJD(7) * t492 - t272 * t342 + t273 * t339;
t89 = qJD(4) * t220 - t450 * t228 + t229 * t441;
t87 = qJD(9) * t510 + t224 * t436 - t225 * t445;
t86 = qJD(9) * t217 - t224 * t445 - t225 * t436;
t75 = t268 * t449 + t440 * t91;
t74 = -t268 * t440 + t449 * t91;
t73 = (-t163 * t434 - t164 * t435) * pkin(4) + t266;
t45 = t226 * t468 - t227 * t299 - t272 * t507 + t273 * t308;
t44 = -t226 * t299 - t227 * t468 - t272 * t308 - t273 * t507;
t41 = -mrSges(5,2) * t401 + mrSges(5,3) * t58;
t40 = mrSges(5,1) * t401 - mrSges(5,3) * t57;
t38 = mrSges(10,1) * t400 - mrSges(10,3) * t55;
t22 = -mrSges(11,2) * t393 + mrSges(11,3) * t28;
t21 = mrSges(11,1) * t393 - mrSges(11,3) * t27;
t1 = [(t117 * mrSges(5,2) - t120 * mrSges(5,3) + Ifges(5,1) * t57 + Ifges(5,4) * t58 + Ifges(5,5) * t401 + t114 * t529 + t49 * t551 + t519 * t720 + t522 * t723 + t524 * t724) * t220 + (t575 / 0.2e1 - t119 * mrSges(5,3) + Ifges(6,3) * t720 + Ifges(6,6) * t723 + Ifges(6,5) * t724 - Ifges(5,4) * t57 - Ifges(5,2) * t58 + t117 * mrSges(5,1) - Ifges(5,6) * t401 + t731) * t218 + t511 * (Ifges(5,1) * t88 - Ifges(5,4) * t89) / 0.2e1 + (-t259 * mrSges(6,1) - t258 * mrSges(6,2) + t444 * t729 + t453 * t730) * g(1) + (-t261 * mrSges(6,1) - t260 * mrSges(6,2) + t444 * t730 - t453 * t729) * g(2) + (mrSges(4,2) * t266 - mrSges(4,3) * t303 + Ifges(4,1) * t166 + Ifges(4,4) * t165 + Ifges(4,5) * t424) * t309 + (Ifges(3,4) * t325 + Ifges(3,2) * t326) * t697 + (-mrSges(4,1) * t266 + mrSges(4,3) * t302 + Ifges(4,4) * t166 + Ifges(4,2) * t165 + Ifges(4,6) * t424) * t312 + (-mrSges(8,1) * t266 + mrSges(8,3) * t301 + Ifges(8,4) * t163 + Ifges(8,2) * t164 + Ifges(8,6) * t423) * t308 + (mrSges(4,1) * t165 + mrSges(8,1) * t164 - mrSges(4,2) * t166 - mrSges(8,2) * t163) * t381 + t448 * (Ifges(7,4) * t324 + Ifges(7,2) * t323) / 0.2e1 + (m(5) * t117 - mrSges(5,1) * t58 + mrSges(5,2) * t57) * t251 + t89 * t740 + (Ifges(4,1) * t229 + Ifges(4,4) * t228) * t705 + t230 * (mrSges(6,1) * t498 - mrSges(6,2) * t497) + t151 * (-Ifges(6,4) * t497 - Ifges(6,2) * t498 + Ifges(6,6) * t89) / 0.2e1 + t193 * (-Ifges(6,5) * t497 - Ifges(6,6) * t498 + Ifges(6,3) * t89) / 0.2e1 + (pkin(15) ^ 2 * t782 + t477 * pkin(14) + Ifges(2,3)) * qJDD(1) - t281 * t590 / 0.2e1 + (Ifges(7,1) * t324 + Ifges(7,4) * t323) * t699 + t486 * t550 - t226 * t537 + (Ifges(8,1) * t226 + Ifges(8,4) * t227) * t703 + (Ifges(9,1) * t224 + Ifges(9,4) * t225) * t704 + t607 * t673 + t89 * t641 - t225 * t691 + t324 * (t633 + t648) / 0.2e1 + (Ifges(3,1) * t325 + Ifges(3,4) * t326) * t781 + t323 * t521 / 0.2e1 + t326 * t523 / 0.2e1 - (-mrSges(9,1) * t225 + mrSges(9,2) * t224) * t619 + (t283 * t697 + t520 * qJD(2) / 0.2e1) * qJD(2) - (-mrSges(9,2) * t582 + Ifges(9,1) * t161 + Ifges(9,4) * t162 + Ifges(9,5) * t422) * t508 + (m(11) * t73 - mrSges(11,1) * t28 + mrSges(11,2) * t27) * ((-t308 * t435 + t434 * t507) * pkin(4) - t381) + (mrSges(11,2) * t73 - mrSges(11,3) * t63 + Ifges(11,1) * t27 + Ifges(11,4) * t28 + Ifges(11,5) * t393) * (t299 * t507 - t308 * t468) + (-mrSges(11,1) * t73 + mrSges(11,3) * t62 + Ifges(11,4) * t27 + Ifges(11,2) * t28 + Ifges(11,6) * t393) * (-t299 * t308 - t468 * t507) - (mrSges(8,2) * t266 - mrSges(8,3) * t300 + Ifges(8,1) * t163 + Ifges(8,4) * t164 + Ifges(8,5) * t423) * t507 - t498 * t48 / 0.2e1 + t325 * t525 / 0.2e1 - t87 * t546 + (t487 + t448 * (-Ifges(7,2) * t439 + t648)) * t580 / 0.2e1 + t542 * (Ifges(11,4) * t44 + Ifges(11,2) * t45) / 0.2e1 + t513 * (Ifges(10,4) * t86 + Ifges(10,2) * t87) / 0.2e1 + t541 * (Ifges(5,4) * t88 - Ifges(5,2) * t89) / 0.2e1 + t9 * t615 / 0.2e1 - t8 * t616 / 0.2e1 - t228 * t539 + t88 * t556 + (-t4 * t616 + t497 * t59 - t498 * t60 - t5 * t615) * mrSges(6,3) + (mrSges(3,1) * t326 + mrSges(9,1) * t162 - mrSges(3,2) * t325 - mrSges(9,2) * t161) * pkin(15) + (-mrSges(4,1) * t228 - mrSges(8,1) * t227 + mrSges(4,2) * t229 + mrSges(8,2) * t226) * t343 + (mrSges(9,1) * t582 + Ifges(9,4) * t161 + Ifges(9,2) * t162 + Ifges(9,6) * t422) * t307 - t225 * t577 - t88 * t642 - t44 * t620 + t227 * t538 + t86 * t545 + (-mrSges(10,1) * t141 + mrSges(10,3) * t263 + Ifges(10,4) * t55 + Ifges(10,2) * t56 + Ifges(10,6) * t400) * t217 + (t602 / 0.2e1 + t518 * qJD(6) / 0.2e1) * qJD(6) + (Ifges(10,1) * t86 + Ifges(10,4) * t87) * t709 + t495 * t580 + (-Ifges(6,1) * t497 - Ifges(6,4) * t498 + Ifges(6,5) * t89) * t712 + t341 * t582 - t756 * (pkin(9) * t89 - pkin(11) * t88 + t206) + t758 * t581 + t759 * t206 + t770 * t396 - t89 * t739 + t45 * t65 / 0.2e1 + t44 * t66 / 0.2e1 + t771 * (t396 + (-t226 * t434 - t227 * t435) * pkin(4)) - t45 * t617 + (-t93 * t586 + m(6) * (qJD(5) * t516 + t4 * t440 + t449 * t5) + t92 * t585 + t440 * t12 + t449 * t11) * (pkin(9) * t218 - pkin(11) * t220 + t251) - (mrSges(10,2) * t141 - mrSges(10,3) * t262 + Ifges(10,1) * t55 + Ifges(10,4) * t56 + Ifges(10,5) * t400) * t510 + (m(10) * t141 - mrSges(10,1) * t56 + mrSges(10,2) * t55) * (-pkin(2) * t307 - pkin(15)) + (Ifges(3,5) * t452 + 0.2e1 * Ifges(3,6) * t697) * qJDD(2) + (0.2e1 * Ifges(7,5) * t699 + Ifges(7,6) * t448) * qJDD(6) + t89 * t47 / 0.2e1 + t87 * t99 / 0.2e1 - t89 * t100 / 0.2e1 + t86 * t101 / 0.2e1 + t88 * t102 / 0.2e1 + t129 * (Ifges(11,1) * t44 + Ifges(11,4) * t45) / 0.2e1 - t775 * (-t266 * t381 + t343 * t396) + t160 * (-mrSges(11,1) * t45 + mrSges(11,2) * t44) + t224 * t188 / 0.2e1 + t225 * t185 / 0.2e1 + t226 * t189 / 0.2e1 + t227 * t186 / 0.2e1 + t228 * t187 / 0.2e1 + t229 * t190 / 0.2e1 + t234 * (mrSges(5,1) * t89 + mrSges(5,2) * t88) + t250 * (-mrSges(10,1) * t87 + mrSges(10,2) * t86) - t229 * t540 + t284 * (Ifges(9,4) * t224 + Ifges(9,2) * t225) / 0.2e1 + t285 * (Ifges(8,4) * t226 + Ifges(8,2) * t227) / 0.2e1 + t289 * (Ifges(4,4) * t229 + Ifges(4,2) * t228) / 0.2e1 + pkin(14) * (-mrSges(7,1) * t323 + mrSges(7,2) * t324) + t402 * (Ifges(11,5) * t44 + Ifges(11,6) * t45) / 0.2e1 + t405 * (Ifges(10,5) * t86 + Ifges(10,6) * t87) / 0.2e1 + t406 * (Ifges(5,5) * t88 - Ifges(5,6) * t89) / 0.2e1 + t425 * (Ifges(9,5) * t224 + Ifges(9,6) * t225) / 0.2e1 + t426 * (Ifges(8,5) * t226 + Ifges(8,6) * t227) / 0.2e1 + t427 * (Ifges(4,5) * t229 + Ifges(4,6) * t228) / 0.2e1; (Ifges(7,5) * t324 + Ifges(7,6) * t323 + Ifges(7,3) * qJDD(6) + (-t602 / 0.2e1 + t607 / 0.2e1 - t448 * t392 / 0.2e1 + t518 * t673 + (-t495 - t487 / 0.2e1 + t628 * t699) * qJD(1)) * qJD(1)) * t749 - t513 * t545 + t457 + t456 + t455 + (t453 * t753 - t533 + t766 + t768) * g(1) + (t444 * t753 - t534 + t767 + t769) * g(2) + t756 * t91 - t759 * t686 - t771 * t180 + t751 * (pkin(11) + t276) + m(5) * (t119 * t276 + t120 * t275) + t520 * t550 + (mrSges(8,1) * t423 - mrSges(8,3) * t163) * t694 + (mrSges(4,1) * t424 - mrSges(4,3) * t166) * t695 + (-mrSges(8,2) * t423 + mrSges(8,3) * t164) * t696 + ((-(-t292 * t373 + t293 * t376) * t743 - (t292 * t371 - t293 * t370) * t744) * t221 + t456 + (-t614 + (t241 + t656) * t438) * t670 + t734 * t116 + t735 * t115 - t180 * t72 + t207 * t21 + t208 * t22 + (t130 * t734 - t131 * t735 - t160 * t180 + t207 * t63 + t208 * t62) * m(11) + t766 * g(1) + t767 * g(2) + (-m(11) * t353 - t765) * g(3)) * t736 - (-mrSges(4,2) * t424 + mrSges(4,3) * t165) * t693 + (-m(6) * t555 - t503 - m(5) * t593 + t340 * t749 + (m(6) - t775) * t421 + t783) * g(3) + (-t241 * t438 + t614) * t669 + (m(8) * (t300 * t447 + t301 * t438) + m(4) * (-t302 * t451 + t303 * t442)) * pkin(1) + t240 * t568 + t242 * t569 + t283 * t592 / 0.2e1 + (m(4) * t692 - m(5) * t329 + t528 * t749 + t532 + t780) * t760 + (-t169 * t583 - t38 * t445 - t634) * pkin(2) + ((t775 * t343 + t756 - t759 - t770 - t771) * pkin(1) + t281 / 0.2e1) * t591 + t288 * t538 + (m(6) * t114 + t10) * (-pkin(9) - t275) + m(11) * (t130 * t95 - t131 * t94 + t167 * t63 + t168 * t62) + (-m(5) * t505 + m(6) * t516 + t757) * (t387 * t587 + (t451 * t588 + (t442 * t450 + t604) * qJD(3)) * pkin(1)) + (t486 / 0.2e1 - t758) * qJD(1) ^ 2 + (m(5) * t236 - m(6) * t230 - t772) * (-t387 * t588 - t441 * t569 + t574 * t578) + Ifges(3,3) * qJDD(2) + t94 * t115 + t95 * t116 + t167 * t21 + t168 * t22 + t275 * t40 + t276 * t41 + Ifges(3,5) * t325 + Ifges(3,6) * t326; t455 - g(1) * (t453 * t462 + t533) - g(2) * (t444 * t462 + t534) + t40 * t683 + t41 * t684 - t240 * t570 - t109 * t686 - t75 * t92 - t74 * t93 - t268 * t170 - t242 * t572 + t388 * t10 + (-t503 + t764) * g(3) + ((-t327 * t403 + t328 * t404) * (-g(3) * t595 - g(1) * (-t453 * t659 + t334) - g(2) * (-t444 * t659 + t331) + (t169 * t445 - t171 * t436 + (t436 * t512 - t445 * t513) * mrSges(10,3)) * t425 * pkin(2) + t747) + (-t403 * t411 - t404 * t415) * pkin(12) * (t457 + (-t634 + (-qJD(9) * t169 - t513 * t655 - t38) * t445) * pkin(2) + t785 * g(3) + t769 * g(2) + t768 * g(1))) * pkin(6) / (t327 * t415 + t328 * t411) / pkin(12) + ((-t313 * t376 + t314 * t373) * t743 + (t313 * t370 - t314 * t371) * t744) * t221 + (-g(3) * t555 + t114 * t388 + (t230 * t441 + t450 * t516) * qJD(4) * pkin(5) - t230 * t267 - t59 * t74 - t60 * t75) * m(6) + ((t119 * t441 + t120 * t450 + (-t236 * t441 - t450 * t505) * qJD(4)) * pkin(5) - t234 * t686 + t236 * t267 + t268 * t505 - g(3) * t380) * m(5) + t757 * pkin(5) * t587 + t751 * (pkin(11) + t684) + t772 * (pkin(5) * t588 - t267) + (m(5) * t685 + t780) * t760; -t230 * (mrSges(6,1) * t152 + mrSges(6,2) * t151) + (Ifges(6,1) * t151 - t646) * t713 + t48 * t712 + (Ifges(6,5) * t151 - Ifges(6,6) * t152) * t711 - t59 * t92 + t60 * t93 - g(1) * (mrSges(6,1) * t260 - mrSges(6,2) * t261) - g(2) * (-mrSges(6,1) * t258 + mrSges(6,2) * t259) + g(3) * t529 * t383 + (t151 * t59 + t152 * t60) * mrSges(6,3) + t575 + (-Ifges(6,2) * t152 + t150 + t49) * t714 + t731;];
tau = t1(:);
