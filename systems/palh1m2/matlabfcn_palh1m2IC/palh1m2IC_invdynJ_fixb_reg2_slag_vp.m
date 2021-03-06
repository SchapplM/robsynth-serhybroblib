% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% palh1m2IC
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
% 
% Output:
% tau_reg [4x(13*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh1m2IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:47:01
% EndTime: 2020-05-02 23:47:29
% DurationCPUTime: 29.14s
% Computational Cost: add. (14033->887), mult. (28763->1253), div. (267->9), fcn. (24273->50), ass. (0->491)
t395 = cos(qJ(1));
t670 = t395 * g(1);
t387 = sin(qJ(1));
t671 = t387 * g(2);
t283 = t670 + t671;
t384 = sin(qJ(4));
t682 = cos(qJ(4));
t558 = t682 * g(3);
t233 = -t384 * t283 + t558;
t385 = sin(qJ(3));
t393 = cos(qJ(3));
t672 = t384 * g(3);
t466 = -t283 * t682 - t672;
t139 = -t233 * t385 + t466 * t393;
t386 = sin(qJ(2));
t394 = cos(qJ(2));
t528 = t233 * t393 + t466 * t385;
t736 = t139 * t386 + t528 * t394;
t677 = pkin(5) * t384;
t347 = qJD(3) * t677;
t337 = pkin(1) * t385 + pkin(5);
t539 = t682 * qJD(4);
t540 = t682 * qJD(3);
t620 = t384 * t385;
t565 = pkin(1) * t620;
t591 = qJD(4) * t384;
t167 = -t337 * t591 - qJD(3) * t565 + (t540 + t539) * t393 * pkin(1);
t619 = t384 * t393;
t328 = pkin(1) * t619;
t243 = t337 * t682 + t328;
t560 = pkin(5) * t682;
t549 = t167 * qJD(2) + t243 * qJDD(2) + qJDD(3) * t560;
t437 = qJD(4) * t347 - t549;
t368 = qJDD(2) + qJDD(3);
t516 = qJDD(4) + t368;
t82 = -pkin(9) * t516 + t437;
t766 = t736 + t82;
t765 = t139 * t394 - t528 * t386;
t402 = 0.1e1 / pkin(10);
t570 = qJ(4) + pkin(18);
t501 = pkin(19) + qJ(3) + t570;
t610 = qJ(6) - qJ(2);
t452 = t501 + t610;
t430 = -(2 * qJ(7)) - pkin(20) + t452;
t453 = t501 - t610;
t436 = pkin(20) + t453;
t696 = cos(qJ(10) - t430) + cos(qJ(10) - t436);
t194 = 0.1e1 / t696;
t405 = 0.1e1 / pkin(3);
t569 = qJ(7) + pkin(20);
t505 = t569 - t610;
t685 = pkin(4) * (pkin(1) * (sin(qJ(10) - t610) + sin(qJ(10) + t610)) + (sin(qJ(10) - t505) + sin(qJ(10) + t505)) * pkin(3)) * t194 * t405;
t493 = t402 * t685;
t751 = t493 - 0.1e1;
t278 = cos(qJ(7) + qJ(10) - t501);
t275 = 0.1e1 / t278;
t372 = pkin(19) - qJ(7);
t534 = -qJ(10) + t372;
t702 = (-pkin(10) * t278 - pkin(5) * cos(qJ(3) + t534)) * t275;
t552 = t402 * t702;
t750 = t552 + 0.1e1;
t536 = -qJ(8) + pkin(17) + qJ(3);
t313 = cos(qJ(9) - t536);
t380 = sin(qJ(9));
t563 = pkin(6) / t380 / pkin(2);
t634 = (-pkin(12) * t313 + pkin(2) * cos(t536)) / pkin(12);
t740 = t563 * (t313 + t634);
t390 = cos(qJ(7));
t381 = sin(qJ(7));
t597 = qJD(1) * t386;
t543 = t381 * t597;
t596 = qJD(1) * t394;
t254 = t390 * t596 - t543;
t379 = sin(pkin(19));
t659 = cos(pkin(19));
t668 = sin(qJ(10));
t669 = cos(qJ(10));
t264 = t379 * t668 + t659 * t669;
t429 = -t379 * t669 + t659 * t668;
t613 = t390 * t386;
t267 = t381 * t394 + t613;
t599 = qJD(1) * t267;
t477 = t254 * t429 + t264 * t599;
t99 = -t254 * t264 + t429 * t599;
t749 = t477 * t99;
t403 = 0.1e1 / pkin(8);
t438 = -qJ(7) + t452;
t439 = -qJ(7) + t453;
t553 = t194 * t403 * ((-(cos(t439) + cos(t438)) * pkin(1) + (-cos(t430) - cos(t436)) * pkin(3)) * pkin(4) + ((cos(qJ(10) - t439) + cos(qJ(10) - t438)) * pkin(1) + t696 * pkin(3)) * pkin(8));
t314 = cos(t505);
t306 = 0.1e1 / t314;
t633 = (-pkin(3) * t314 - pkin(1) * cos(t610)) * t306;
t758 = (t553 + t633) * t405 + 0.1e1;
t763 = t758 * t749;
t253 = t385 * t596 + t393 * t597;
t255 = -t385 * t597 + t393 * t596;
t151 = t253 * t384 - t682 * t255;
t458 = t682 * t253 + t384 * t255;
t336 = pkin(5) * t385 + pkin(1);
t612 = t393 * t394;
t695 = pkin(5) * t612 - t336 * t386 + pkin(15);
t600 = qJD(1) * t695;
t70 = pkin(9) * t151 - pkin(11) * t458 - t600;
t547 = t682 * t385;
t457 = t547 + t619;
t166 = t337 * t539 + (qJD(3) * t457 + t393 * t591) * pkin(1);
t546 = t682 * t393;
t327 = pkin(1) * t546;
t242 = t337 * t384 - t327;
t513 = pkin(5) * t540;
t89 = t166 * qJD(2) + qJD(4) * t513 + t242 * qJDD(2) + qJDD(3) * t677;
t81 = pkin(11) * t516 + t89;
t762 = -qJD(5) * t70 - t765 - t81;
t392 = cos(qJ(5));
t587 = qJD(5) * t392;
t726 = t151 * t392;
t761 = t587 + t726;
t197 = t243 * qJD(2) + t513;
t371 = qJD(2) + qJD(3);
t533 = qJD(4) + t371;
t186 = -pkin(9) * t533 - t197;
t383 = sin(qJ(5));
t195 = t242 * qJD(2) + t347;
t185 = pkin(11) * t533 + t195;
t55 = t185 * t392 + t383 * t70;
t422 = t186 * t587 + t766 * t383 + t55 * t458;
t333 = t385 * t386;
t476 = -t612 + t333;
t416 = t371 * t476;
t414 = qJD(1) * t416;
t614 = t386 * t393;
t618 = t385 * t394;
t268 = t614 + t618;
t574 = t268 * qJDD(1);
t116 = -t414 + t574;
t417 = t371 * t268;
t409 = -t417 * qJD(1) - t476 * qJDD(1);
t52 = -t682 * t116 + t253 * t591 - t255 * t539 - t384 * t409;
t529 = t384 * t116 - t682 * t409;
t53 = qJD(4) * t458 + t529;
t576 = qJDD(1) * t695;
t326 = pkin(5) * t614;
t260 = pkin(5) * t618 + t326;
t592 = qJD(2) * t394;
t168 = qJD(2) * t326 + t260 * qJD(3) + t336 * t592;
t603 = qJD(1) * t168;
t16 = pkin(9) * t53 + pkin(11) * t52 - t576 + t603;
t54 = -t185 * t383 + t392 * t70;
t7 = qJD(5) * t54 + t383 * t16 + t392 * t81;
t6 = t7 * t392;
t760 = t6 + t765;
t759 = -t437 - t736;
t754 = -t89 - t765;
t389 = cos(qJ(8));
t681 = sin(qJ(8));
t545 = t681 * t386;
t508 = qJD(1) * t545;
t249 = -t389 * t596 + t508;
t266 = t389 * t386 + t681 * t394;
t250 = t266 * qJD(1);
t388 = cos(qJ(9));
t154 = t249 * t380 - t250 * t388;
t369 = qJD(2) + qJD(8);
t356 = qJD(9) + t369;
t535 = qJDD(1) * t681;
t571 = t394 * qJDD(1);
t189 = t369 * t266;
t602 = qJD(1) * t189;
t112 = -t386 * t535 + t389 * t571 - t602;
t446 = -t369 * t508 + t394 * t535;
t518 = t394 * t369;
t572 = t386 * qJDD(1);
t113 = (qJD(1) * t518 + t572) * t389 + t446;
t584 = qJD(9) * t380;
t50 = t113 * t380 - t249 * t584 + (qJD(9) * t250 - t112) * t388;
t34 = t154 * t356 + t50;
t478 = -t249 * t388 - t250 * t380;
t724 = t478 * t154;
t756 = t724 * t740;
t673 = g(3) * t393;
t694 = t283 * t385 - t673;
t625 = t283 * t393;
t718 = g(3) * t385 + t625;
t753 = t718 * t386 + t694 * t394;
t583 = qJD(9) * t388;
t51 = t380 * t112 + t113 * t388 - t249 * t583 - t250 * t584;
t35 = -t356 * t478 + t51;
t746 = t151 * t458;
t551 = t405 * t633;
t274 = t283 * t394;
t361 = t386 * pkin(1);
t738 = pkin(1) * t274 + g(3) * t361;
t488 = t392 * t533;
t117 = t383 * t458 - t488;
t588 = qJD(5) * t383;
t30 = -qJD(5) * t488 - t383 * t516 + t392 * t52 + t458 * t588;
t119 = t383 * t533 + t392 * t458;
t590 = qJD(5) * t119;
t31 = -t383 * t52 - t392 * t516 + t590;
t733 = -t117 * t761 - t30 * t392 - t383 * t31;
t27 = t30 * t383;
t732 = t119 * t761 - t27;
t145 = qJD(5) + t151;
t47 = qJDD(5) + t53;
t731 = -t119 * t458 + t145 * t761 + t383 * t47;
t524 = pkin(5) * t275 * sin(t570) * t403;
t727 = t524 * t749;
t646 = t151 * t186;
t578 = qJD(1) * qJD(2);
t721 = t386 * t578 - t571;
t644 = t151 * t383;
t719 = t588 + t644;
t717 = -t477 ^ 2 + t99 ^ 2;
t716 = -t151 ^ 2 + t458 ^ 2;
t43 = -t154 ^ 2 + t478 ^ 2;
t715 = pkin(9) * t458 + pkin(11) * t151;
t492 = t381 * t572 - t390 * t571;
t370 = qJD(2) + qJD(7);
t191 = t370 * t267;
t601 = qJD(1) * t191;
t114 = t492 + t601;
t467 = -qJD(7) * t543 - t721 * t381;
t517 = t370 * t394;
t115 = (qJD(1) * t517 + t572) * t390 + t467;
t239 = t264 * qJD(10);
t240 = t429 * qJD(10);
t33 = -t114 * t429 + t264 * t115 + t239 * t254 - t240 * t599;
t350 = qJD(10) + t370;
t714 = t350 * t99 + t33;
t331 = t361 - pkin(15);
t580 = t331 * qJD(1);
t120 = t580 + (-t254 * t379 + t599 * t659) * pkin(4);
t330 = -qJ(2) + t534;
t309 = sin(t330);
t310 = cos(t330);
t658 = pkin(1) * qJDD(2);
t344 = t390 * t658;
t367 = qJDD(2) + qJDD(7);
t548 = t659 * pkin(4);
t579 = t381 * qJD(2);
t557 = pkin(1) * t579;
t199 = -qJD(7) * t557 + t367 * t548 + t344;
t586 = qJD(7) * t390;
t449 = qJD(2) * t586 + t381 * qJDD(2);
t678 = pkin(4) * t379;
t200 = pkin(1) * t449 + t367 * t678;
t261 = t370 * t678 + t557;
t593 = qJD(2) * t390;
t262 = pkin(1) * t593 + t370 * t548;
t57 = -t264 * t199 + t200 * t429 + t239 * t261 + t240 * t262;
t418 = g(3) * t309 - t120 * t99 - t283 * t310 + t57;
t713 = t151 * t533 - t52;
t365 = pkin(11) * t384;
t559 = pkin(9) * t682;
t503 = -t365 - t559;
t279 = pkin(5) - t503;
t567 = pkin(9) * t671;
t297 = -pkin(11) * g(3) + t567;
t338 = t365 + pkin(5);
t140 = pkin(9) * t672 + t279 * t670 + t297 * t682 + t338 * t671;
t296 = -pkin(9) * g(3) - pkin(11) * t671;
t285 = -pkin(9) * t384 + pkin(11) * t682;
t562 = t285 * t670;
t142 = -g(3) * t338 + t296 * t682 + t384 * t567 - t562;
t711 = t140 * t385 + t142 * t393;
t710 = t140 * t393 - t142 * t385;
t657 = qJD(1) * pkin(15);
t202 = pkin(2) * t250 - t657;
t378 = qJ(2) + qJ(8);
t358 = qJ(9) + t378;
t334 = sin(t358);
t335 = cos(t358);
t555 = pkin(2) * qJD(9) * t369;
t366 = qJDD(2) + qJDD(8);
t705 = pkin(2) * t366;
t709 = -g(3) * t335 + t154 * t202 + t283 * t334 + t380 * t705 + t388 * t555;
t511 = 0.1e1 + t551;
t699 = t511 * t599 * t254;
t519 = t255 * t371;
t88 = t116 - t519;
t697 = -t54 * t726 - t55 * t644 + t760;
t489 = -t383 * t54 + t392 * t55;
t675 = g(2) * t395;
t676 = g(1) * t387;
t499 = -t675 + t676;
t667 = pkin(1) * qJD(2);
t556 = t393 * t667;
t693 = qJD(3) * t556 - t253 * t580 + t385 * t658 + t753;
t45 = t392 * t47;
t692 = t117 * t458 + t45;
t532 = t186 * t588 - t54 * t458;
t32 = t264 * t114 + t115 * t429 + t239 * t599 + t240 * t254;
t691 = -t350 * t477 + t32;
t56 = -t199 * t429 - t264 * t200 - t239 * t262 + t240 * t261;
t410 = -g(3) * t310 - t120 * t477 - t283 * t309 - t56;
t687 = -g(3) * t334 + t478 * t202 - t283 * t335 + t380 * t555;
t686 = t458 * t371 - t529;
t190 = -t369 * t545 + t389 * t518;
t679 = pkin(2) * t190;
t29 = t31 * t392;
t664 = t54 * t145;
t663 = t55 * t145;
t656 = t119 * t117;
t649 = t145 * t383;
t647 = t458 * t695;
t645 = t151 * t695;
t182 = t682 * t268 - t384 * t476;
t642 = t182 * t383;
t641 = t182 * t392;
t187 = t457 * t499;
t640 = t187 * t386;
t399 = qJD(1) ^ 2;
t639 = t695 * t399;
t632 = t250 * t249;
t630 = t599 * t370;
t629 = t253 * t255;
t628 = t253 * t371;
t627 = t254 * t370;
t624 = t499 * t394;
t623 = t356 * t369;
t622 = t381 * t386;
t382 = sin(qJ(6));
t391 = cos(qJ(6));
t621 = t382 * t391;
t615 = t386 * t499;
t611 = t394 * t399;
t454 = pkin(1) * (-t264 * t390 + t381 * t429);
t609 = (-t239 * t659 + t240 * t379) * pkin(4) - qJD(2) * t454;
t455 = pkin(1) * (t264 * t381 + t390 * t429);
t608 = (t239 * t379 + t240 * t659) * pkin(4) - qJD(2) * t455;
t606 = g(3) * t386 + t274;
t374 = t382 ^ 2;
t376 = t391 ^ 2;
t605 = t374 - t376;
t375 = t386 ^ 2;
t377 = t394 ^ 2;
t604 = t375 - t377;
t598 = qJD(1) * t268;
t258 = t327 - t565;
t595 = qJD(2) * t258;
t259 = pkin(1) * t547 + t328;
t594 = qJD(2) * t259;
t589 = qJD(5) * t145;
t585 = qJD(7) * t394;
t582 = qJDD(1) * pkin(14);
t581 = qJDD(1) * pkin(15);
t577 = qJD(1) * qJD(6);
t575 = qJDD(2) * pkin(1) ^ 2;
t573 = t331 * qJDD(1);
t568 = 0.2e1 * pkin(15);
t554 = pkin(1) * t592;
t550 = t386 * t611;
t544 = t145 * t588;
t542 = t385 ^ 2 + t393 ^ 2;
t537 = t394 * t578;
t225 = pkin(11) + t242;
t241 = t336 * t394 + t326;
t71 = qJD(1) * t241 + t715;
t530 = qJD(5) * t225 + t71;
t525 = pkin(1) * t306 * sin(t569) / pkin(7);
t523 = t313 * t563;
t522 = 0.2e1 * pkin(14) * t577;
t521 = -qJD(5) * t185 - t499;
t515 = t116 + t574;
t514 = pkin(5) * t539;
t512 = pkin(1) * t537;
t346 = qJDD(10) + t367;
t349 = pkin(15) * t676;
t509 = -pkin(15) * t675 + t349;
t507 = t386 * t537;
t506 = t577 * t621;
t354 = sin(t372);
t355 = cos(t372);
t363 = g(3) * t394;
t504 = ((t283 * t386 - t363) * t354 + t606 * t355) * pkin(4);
t227 = -t681 * g(3) - t283 * t389;
t231 = t389 * g(3) - t681 * t283;
t497 = t227 * t386 + t231 * t394 - t250 * t657;
t496 = -t227 * t394 + t231 * t386 - t249 * t657;
t494 = t399 * t525;
t100 = t261 * t429 - t262 * t264;
t101 = t261 * t264 + t262 * t429;
t491 = t100 * t477 - t101 * t99;
t490 = t383 * t55 + t392 * t54;
t487 = qJDD(1) * t525;
t485 = -t225 * t47 + t646;
t339 = pkin(11) + t677;
t484 = -t339 * t47 + t646;
t480 = -t151 * t197 + t195 * t458;
t456 = t546 - t620;
t188 = t456 * t499;
t479 = t187 * t394 + t188 * t386;
t269 = t389 * t394 - t545;
t179 = t266 * t388 + t269 * t380;
t181 = -t266 * t380 + t269 * t388;
t474 = t402 * t516;
t473 = t719 * t117 - t29;
t471 = t750 * t458;
t470 = t523 * t632;
t464 = pkin(5) * t591 + t595;
t459 = t682 * t476;
t67 = qJD(4) * t459 + t268 * t591 + t384 * t417 + t682 * t416;
t463 = t182 * t587 - t383 * t67;
t462 = -t182 * t588 - t392 * t67;
t451 = t494 * t621;
t450 = -t537 - t572;
t445 = -pkin(14) * t399 + t283;
t444 = t751 * t458;
t443 = -t499 + 0.2e1 * t582;
t228 = -t381 * g(3) - t283 * t390;
t232 = g(3) * t390 - t283 * t381;
t442 = t228 * t386 + t232 * t394 + t580 * t599;
t441 = -pkin(1) * t331 * t611 + t738;
t440 = -t492 + t630;
t143 = (t254 * t659 + t379 * t599) * pkin(4);
t435 = -t145 * t644 - t544 + t692;
t433 = -t255 * t580 - t386 * t694 + t393 * t658 + t394 * t718;
t432 = -t228 * t394 + t232 * t386 - t254 * t580 + t344;
t431 = -t467 + t627;
t180 = t268 * t384 + t459;
t77 = pkin(9) * t180 - pkin(11) * t182 - t695;
t428 = t182 * t82 - t186 * t67 - t589 * t77 - t283;
t427 = -t371 * t598 - t628;
t426 = pkin(9) * t558 - t297 * t384 + t562;
t424 = -t719 * t119 + t733;
t68 = t182 * qJD(4) - t384 * t416 + t682 * t417;
t24 = pkin(9) * t68 + pkin(11) * t67 + t168;
t423 = qJD(5) * t182 * t186 + t145 * t24 + t47 * t77 - t640;
t15 = t392 * t16;
t8 = -qJD(5) * t55 - t383 * t81 + t15;
t421 = -qJD(5) * t490 - t8 * t383 + t6;
t398 = qJD(2) ^ 2;
t397 = qJD(6) ^ 2;
t357 = sin(t378);
t348 = qJDD(9) + t366;
t340 = -t560 - pkin(9);
t282 = pkin(1) * t390 + t548;
t281 = pkin(1) * t381 + t678;
t270 = t390 * t394 - t622;
t256 = t349 + (t581 - t675) * pkin(15);
t245 = pkin(2) * t266 - pkin(15);
t226 = -pkin(9) - t243;
t193 = t499 * t285;
t192 = -qJD(7) * t622 - t386 * t579 + t390 * t517;
t178 = t279 * t676 + (-t559 - t338) * t675;
t170 = t456 * t624;
t161 = (-t264 * t379 - t429 * t659) * pkin(4);
t160 = (-t264 * t659 + t379 * t429) * pkin(4);
t148 = -pkin(11) * t558 - t296 * t384 - t503 * t670 + t559 * t671;
t144 = (t267 * t659 - t270 * t379) * pkin(4) + t331;
t132 = (-t499 + 0.2e1 * t512 + t573) * t331;
t125 = pkin(1) * t596 + t143;
t124 = t253 ^ 2 - t255 ^ 2;
t123 = t249 ^ 2 - t250 ^ 2;
t122 = -t264 * t281 - t282 * t429;
t121 = -t264 * t282 + t281 * t429;
t107 = -t264 * t270 + t267 * t429;
t106 = t264 * t267 + t270 * t429;
t103 = pkin(2) * t113 - t581;
t87 = t409 + t628;
t86 = -t249 * t369 + (-t369 * t596 - t572) * t389 - t446;
t85 = t250 * t369 + t112;
t84 = t554 + (t191 * t379 + t192 * t659) * pkin(4);
t83 = (cos(t378) * t283 + g(3) * t357 + t202 * t249 + (t380 ^ 2 + t388 ^ 2) * t705) * pkin(2);
t74 = qJD(7) * t455 + t239 * t281 + t240 * t282;
t73 = qJD(7) * t454 - t239 * t282 + t240 * t281;
t72 = qJD(1) * t260 + t715;
t66 = qJD(9) * t181 - t189 * t380 + t190 * t388;
t65 = qJD(9) * t179 + t189 * t388 + t190 * t380;
t64 = t197 * t392 + t383 * t715;
t63 = -t197 * t383 + t392 * t715;
t61 = t383 * t72 + t392 * t594;
t60 = -t383 * t594 + t392 * t72;
t58 = t512 + t573 + (t114 * t379 + t115 * t659) * pkin(4);
t49 = -t191 * t429 + t192 * t264 + t239 * t270 - t240 * t267;
t48 = t191 * t264 + t192 * t429 + t239 * t267 + t240 * t270;
t42 = (-t249 * t478 + t348 * t380 + t356 * t583) * pkin(2) + t709;
t41 = (t356 * t584 + t154 * t249 + (-t348 - t366) * t388) * pkin(2) + t687;
t23 = t195 * t533 + t458 * t600 + t759;
t22 = -t151 * t600 + t197 * t533 + t754;
t18 = -t101 * t350 + t418;
t17 = t100 * t350 + t410;
t14 = t117 * t649 - t29;
t11 = -t145 * t649 + t692;
t9 = (t34 * t388 - t35 * t380) * pkin(2);
t5 = -t119 * t649 + t733;
t4 = -pkin(9) * t31 - t195 * t117 - t63 * t145 + (-pkin(11) * t47 + t646) * t383 + (-pkin(11) * t589 - t766) * t392 + t532;
t3 = pkin(9) * t30 + t64 * t145 - t195 * t119 + t186 * t726 + (t544 - t45) * pkin(11) + t422;
t2 = -t82 * pkin(9) - t55 * t64 - t54 * t63 - t186 * t195 + (t148 * t393 + t385 * t426) * t386 + (t148 * t385 - t393 * t426) * t394 + ((t386 * t547 - t394 * t546) * t671 + t421) * pkin(11);
t1 = t117 * t64 + t119 * t63 + (-t664 + (-t31 + t590) * pkin(11)) * t392 + (-t8 - t663 + (qJD(5) * t117 - t30) * pkin(11)) * t383 + t760;
t10 = [0, 0, 0, 0, 0, qJDD(1), t499, t283, 0, 0, qJDD(1) * t377 - 0.2e1 * t507, -0.2e1 * t386 * t571 + 0.2e1 * t578 * t604, qJDD(2) * t394 - t386 * t398, qJDD(1) * t375 + 0.2e1 * t507, -qJDD(2) * t386 - t394 * t398, 0, t450 * t568 - t615, t721 * t568 - t624, -t283, t256, t116 * t268 - t253 * t416, -t116 * t333 + ((t519 + t515) * t394 + t427 * t386) * t393 + ((-t574 - t519) * t386 + t427 * t394) * t385, t268 * t368 - t371 * t416, -t255 * t417 - t409 * t476, -t368 * t476 - t371 * t417, 0, -t255 * t554 + (0.2e1 * t331 * t417 + t476 * t554) * qJD(1) + (-t499 + 0.2e1 * t573) * t476, -t499 * t614 + (-t385 * t499 + (t253 + t598) * t667) * t394 + (-t414 + t515) * t331, ((-t385 * t268 + t393 * t476) * qJDD(2) + ((-t393 * t268 - t385 * t476) * qJD(3) + t371 * t386 * t542) * qJD(2)) * pkin(1) - t283, t132, -t182 * t52 - t458 * t67, t151 * t67 + t180 * t52 - t182 * t53 - t458 * t68, t182 * t516 - t533 * t67, t151 * t68 + t180 * t53, -t180 * t516 - t533 * t68, 0, -t640 + t170 + (qJD(1) * t180 + t151) * t168 - (qJD(1) * t68 + qJDD(1) * t180 + t53) * t695, (qJD(1) * t182 + t458) * t168 - (-qJD(1) * t67 + qJDD(1) * t182 - t52) * t695 - t479, -t180 * t89 + t182 * t437 - t195 * t68 + t197 * t67 - t283, (t576 - 0.2e1 * t603 + t499) * t695, t119 * t462 - t30 * t641, (t117 * t392 + t119 * t383) * t67 + (t27 - t29 + (t117 * t383 - t119 * t392) * qJD(5)) * t182, t119 * t68 + t145 * t462 - t180 * t30 + t47 * t641, t117 * t463 + t31 * t642, -t117 * t68 - t145 * t463 - t180 * t31 - t47 * t642, t145 * t68 + t180 * t47, t8 * t180 + t54 * t68 + (t170 + t423) * t392 + t428 * t383, -t7 * t180 - t55 * t68 + t428 * t392 + (-t188 * t394 - t423) * t383, (-t119 * t24 - t182 * t8 + t30 * t77 + t54 * t67 + (-t117 * t77 - t182 * t55) * qJD(5)) * t392 + (-t117 * t24 - t182 * t7 - t31 * t77 + t55 * t67 + (t119 * t77 + t182 * t54) * qJD(5)) * t383 + t479, (-pkin(1) * t499 - t178 * t385 + t193 * t393) * t386 + (t178 * t393 + t193 * t385) * t394 + t490 * t24 + (qJD(5) * t489 + t7 * t383 + t8 * t392) * t77 + t509, qJDD(1) * t374 + 0.2e1 * t506, 0.2e1 * qJDD(1) * t621 - 0.2e1 * t577 * t605, qJDD(6) * t382 + t391 * t397, qJDD(1) * t376 - 0.2e1 * t506, qJDD(6) * t391 - t382 * t397, 0, t382 * t522 - t391 * t443, t382 * t443 + t391 * t522, -t283, (-t499 + t582) * pkin(14), -t114 * t270 - t191 * t254, t114 * t267 - t115 * t270 + t191 * t599 - t192 * t254, -t191 * t370 + t270 * t367, t115 * t267 + t192 * t599, -t192 * t370 - t267 * t367, 0, -t499 * t613 + (-t381 * t499 + 0.2e1 * t599 * t667) * t394 + (qJD(1) * t192 + qJDD(1) * t267 + t115) * t331, t381 * t615 + (-t499 * t390 + (qJD(1) * t270 + t254) * t667) * t394 + (qJDD(1) * t270 - t114 - t601) * t331, ((-t267 * t381 - t270 * t390) * qJDD(2) + (t191 * t390 - t192 * t381 + (-t267 * t390 + t270 * t381) * qJD(7)) * qJD(2)) * pkin(1) - t283, t132, t112 * t269 + t189 * t249, -t112 * t266 - t113 * t269 + t189 * t250 + t190 * t249, -t189 * t369 + t269 * t366, t113 * t266 + t190 * t250, -t190 * t369 - t266 * t366, 0, -t266 * t499 + (-qJD(1) * t190 - qJDD(1) * t266 - t113) * pkin(15), -t269 * t499 + (-qJDD(1) * t269 - t112 + t602) * pkin(15), -t283, t256, -t181 * t50 - t478 * t65, -t154 * t65 + t179 * t50 - t181 * t51 - t478 * t66, -t181 * t348 + t356 * t65, -t154 * t66 + t179 * t51, t179 * t348 + t356 * t66, 0, -t103 * t179 + t154 * t679 - t202 * t66 - t245 * t51 + t334 * t499, -t103 * t181 + t202 * t65 + t245 * t50 + t335 * t499 - t478 * t679, ((-t179 * t380 - t181 * t388) * t366 + (-t380 * t66 + t388 * t65 + (-t179 * t388 + t181 * t380) * qJD(9)) * t369) * pkin(2) - t283, t103 * t245 + (t190 * t202 - t357 * t499) * pkin(2) + t509, t107 * t32 + t48 * t99, t106 * t32 + t107 * t33 + t477 * t48 + t49 * t99, t107 * t346 + t350 * t48, t106 * t33 + t477 * t49, t106 * t346 + t350 * t49, 0, -t106 * t58 - t120 * t49 - t144 * t33 - t309 * t499 - t477 * t84, t107 * t58 + t120 * t48 + t144 * t32 + t310 * t499 + t84 * t99, -t100 * t48 - t101 * t49 + t106 * t56 - t107 * t57 - t283, t120 * t84 + t144 * t58 - (t361 + (-t354 * t394 + t355 * t386) * pkin(4)) * t499 + t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t550, -t604 * t399, t571, -t550, -t572, qJDD(2), pkin(15) * t611 + t606, t363 + (-pkin(15) * t399 - t283) * t386, 0, 0, -t629, t124, t88, t629, t87, t368, (qJD(3) * t371 * t393 + t255 * t596 + t368 * t385) * pkin(1) + t693, (-t253 * t596 + t368 * t393 + (-qJD(2) - t371) * t385 * qJD(3)) * pkin(1) + t433, (((t371 * t597 - t571) * t393 - t628) * t393 + (-t393 * (-qJD(3) * t596 + t450) - t88) * t385) * pkin(1), t542 * t575 + t441, -t751 * t746, -t751 * t716, -t751 * t713, t444 * t151, -t751 * t686, -t474 * t685 + t516, t167 * t533 + t243 * t516 - t23 * t493 + (-t151 * t241 + t647) * qJD(1) + t759, -t166 * t533 - t242 * t516 - t22 * t493 + (-t241 * t458 - t645) * qJD(1) + t754, -t151 * t166 - t167 * t458 - t242 * t53 + t243 * t52 + t480, t89 * t242 + t195 * t166 - t437 * t243 + t197 * t167 + t241 * t639 + (-pkin(5) * t673 + t283 * t336) * t394 + (pkin(5) * t625 + g(3) * t336) * t386, -t751 * t732, -t493 * t5 + t424, -t751 * t731, -t14 * t493 + t473, -t11 * t493 + t435, t444 * t145, -t4 * t493 - t167 * t117 + t226 * t31 + (-t145 * t166 + t485) * t383 + (-t145 * t530 - t766) * t392 + t532, -t167 * t119 - t226 * t30 - t3 * t493 + t485 * t392 + (-t392 * t166 + t383 * t530) * t145 + t422, -t1 * t493 + (-t117 * t166 + t119 * t71 - t225 * t31 + (t119 * t225 - t54) * qJD(5)) * t392 + (t117 * t71 + t119 * t166 - t225 * t30 - t8 + (t117 * t225 - t55) * qJD(5)) * t383 + t697, t82 * t226 - t186 * t167 + (pkin(1) * t283 + t711) * t394 - (-pkin(1) * g(3) - t710) * t386 - t2 * t493 - t490 * t71 + t489 * t166 + t421 * t225, -t451, t605 * t494, t382 * t487, t451, t391 * t487, qJDD(6) * t525, (-g(3) * t391 + t382 * t445) * t525, (g(3) * t382 + t391 * t445) * t525, 0, 0, t699, t511 * (t254 ^ 2 - t599 ^ 2), t440 * t551 + t440 + t511 * qJD(1) * (-t381 * t585 - t394 * t579 + (-t586 - t593) * t386), -t699, t431 * t551 + t431 + t511 * t390 * (-qJD(1) * t585 + t450), t511 * t367, t432 * t551 + (-t599 * t596 + t390 * t367 + (-qJD(7) * t370 + (-qJD(7) + (-qJD(7) + t370) * t551) * qJD(2)) * t381) * pkin(1) + t432, t442 * t551 + (-t370 * t586 - t381 * t367 - t254 * t596 + (t370 * t593 - t449) * t551 - t449) * pkin(1) + t442, ((t114 - t630) * t390 + (-t115 + t627) * t381) * pkin(1), (t381 ^ 2 + t390 ^ 2) * t575 + t441, -t632, t123, t85, t632, t86, t366, t496, t497, 0, 0, -t724, t43, t34, t724, t35, t348, t41, t42, t9, t83, -t763, t758 * t717, t758 * t691, t763, t758 * t714, t758 * t346, t121 * t346 + t125 * t477 + t74 * t350 + ((t143 * t477 + t160 * t346 + t350 * t608 + t418) * t633 + t18 * t553) * t405 + t418, -t122 * t346 - t125 * t99 - t73 * t350 + ((-t143 * t99 - t161 * t346 - t350 * t609 + t410) * t633 + t17 * t553) * t405 + t410, -t121 * t32 + t122 * t33 + t73 * t477 - t74 * t99 + (-t160 * t32 + t161 * t33 + t477 * t609 - t608 * t99 + t491) * t551 + t491, t100 * t74 - t101 * t73 - t120 * t125 + t57 * t121 + t56 * t122 + (t100 * t608 - t101 * t609 - t120 * t143 + t160 * t57 + t161 * t56 + t504) * t551 + t504 + t738; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t629, t124, t88, t629, t87, t368, -t371 * t556 + t693, (-qJD(3) + t371) * t385 * t667 + t433, 0, 0, t750 * t746, t750 * t716, t750 * t713, -t471 * t151, t750 * t686, t474 * t702 + t516, -t533 * t595 + t23 * t552 + (-t151 * t260 + t647) * qJD(1) + (t682 * t516 + (-qJD(2) - 0.2e1 * qJD(3) - qJD(4)) * t591) * pkin(5) + t549 - t736, t533 * t594 + t22 * t552 + (-t260 * t458 - t645) * qJD(1) + (-t384 * t516 - t533 * t539) * pkin(5) + t754, (t151 * t259 + t258 * t458) * qJD(2) + (t682 * t52 - t384 * t53 + (-t682 * t151 + t384 * t458) * qJD(4)) * pkin(5) + t480, t260 * t639 + (-t195 * t259 - t197 * t258) * qJD(2) + (t89 * t384 - t437 * t682 + (t195 * t682 - t197 * t384) * qJD(4) + t753) * pkin(5), t750 * t732, t5 * t552 + t424, t750 * t731, t14 * t552 + t473, t11 * t552 + t435, -t471 * t145, t4 * t552 + t340 * t31 - t766 * t392 + t484 * t383 + t464 * t117 + (-t339 * t587 - t383 * t514 - t60) * t145 + t532, -t340 * t30 + t3 * t552 + t484 * t392 + t464 * t119 + (t339 * t588 - t392 * t514 + t61) * t145 + t422, t1 * t552 + t61 * t117 + t60 * t119 + (-t117 * t514 - t31 * t339 + (t119 * t339 - t54) * qJD(5)) * t392 + (t119 * t514 - t30 * t339 - t8 + (t117 * t339 - t55) * qJD(5)) * t383 + t697, t82 * t340 - t55 * t61 - t54 * t60 + t186 * t595 + t710 * t386 + t711 * t394 + t2 * t552 + (t186 * t384 + t489 * t682) * qJD(4) * pkin(5) + t421 * t339, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t470, t123 * t523, t85 * t523, t470, t86 * t523, t366 * t523, t496 * t523, t497 * t523, 0, 0, -t756, t43 * t740, t34 * t740, t756, t35 * t740, t348 * t740, (t313 * t41 + ((-t366 * t388 - t380 * t623) * pkin(2) + t687) * t634) * t563, (t313 * t42 + (-pkin(2) * t388 * t623 + t709) * t634) * t563, t9 * t523, t83 * t523, -t727, t717 * t524, t691 * t524, t727, t714 * t524, t346 * t524, t18 * t524, t17 * t524, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t656, -t117 ^ 2 + t119 ^ 2, t117 * t145 - t30, -t656, t119 * t145 - t31, t47, -t186 * t119 + t383 * t762 + t521 * t392 + t15 + t663, t186 * t117 + t664 + t762 * t392 + (-t16 - t521) * t383, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t10;
