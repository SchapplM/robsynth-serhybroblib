% Calculate time derivative of joint inertia matrix for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1DE1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1DE1_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 15:07:48
% EndTime: 2020-04-13 16:20:55
% DurationCPUTime: 1619.29s
% Computational Cost: add. (45247566->875), mult. (70383472->1708), div. (2716985->42), fcn. (44004788->44), ass. (0->700)
t385 = cos(pkin(20));
t387 = sin(qJ(3));
t391 = cos(qJ(3));
t726 = sin(pkin(20));
t361 = t385 * t387 + t391 * t726;
t778 = pkin(6) * t361;
t626 = pkin(1) * t778;
t353 = 0.2e1 * t626;
t409 = pkin(2) ^ 2;
t404 = pkin(6) ^ 2;
t411 = pkin(1) ^ 2;
t641 = t404 + t411;
t604 = -pkin(13) ^ 2 + t641;
t340 = t353 + t409 + t604;
t352 = -pkin(1) - t778;
t362 = -t385 * t391 + t387 * t726;
t645 = t353 + t404;
t826 = -pkin(2) - pkin(13);
t332 = (pkin(1) - t826) * (pkin(1) + t826) + t645;
t821 = pkin(13) - pkin(2);
t333 = (pkin(1) - t821) * (pkin(1) + t821) + t645;
t685 = t333 * t332;
t413 = sqrt(-t685);
t674 = t362 * t413;
t304 = -pkin(6) * t674 - t340 * t352;
t301 = 0.1e1 / t304 ^ 2;
t777 = pkin(6) * t362;
t305 = t340 * t777 - t352 * t413;
t303 = t305 ^ 2;
t290 = t301 * t303 + 0.1e1;
t288 = 0.1e1 / t290;
t868 = -0.2e1 * t288;
t300 = 0.1e1 / t304;
t347 = t353 + t641;
t789 = pkin(2) * t347;
t620 = t300 * t789;
t867 = 0.2e1 * t288 * t620;
t410 = 0.1e1 / pkin(2);
t619 = t301 * t789;
t560 = t305 * t619;
t511 = t560 * t868;
t866 = t410 * t511;
t865 = t410 * t867;
t339 = t409 - t604 - 0.2e1 * t626;
t337 = 0.1e1 / t339 ^ 2;
t834 = pkin(1) * pkin(6);
t630 = t337 * t834;
t864 = -0.2e1 * t413 * t630;
t406 = pkin(4) ^ 2;
t405 = pkin(5) ^ 2;
t403 = pkin(7) ^ 2;
t388 = sin(qJ(2));
t392 = cos(qJ(2));
t393 = cos(pkin(19));
t794 = sin(pkin(19));
t489 = -t388 * t393 + t392 * t794;
t776 = pkin(7) * t489;
t842 = -0.2e1 * pkin(1);
t644 = -t776 * t842 + t411;
t351 = t403 + t644;
t640 = pkin(3) ^ 2 - pkin(8) ^ 2;
t342 = t351 + t640;
t356 = pkin(1) + t776;
t367 = t388 * t794 + t392 * t393;
t824 = pkin(7) + pkin(8);
t825 = pkin(7) - pkin(8);
t334 = (pkin(3) + t824) * (-pkin(3) + t825) + t644;
t335 = (-pkin(3) + t824) * (pkin(3) + t825) + t644;
t683 = t335 * t334;
t414 = sqrt(-t683);
t308 = pkin(7) * t342 * t367 + t356 * t414;
t654 = t391 * t308;
t672 = t367 * t414;
t307 = -pkin(7) * t672 + t342 * t356;
t663 = t387 * t307;
t492 = t663 / 0.2e1 + t654 / 0.2e1;
t348 = 0.1e1 / t351;
t408 = 0.1e1 / pkin(3);
t680 = t348 * t408;
t283 = t492 * t680;
t655 = t391 * t307;
t662 = t387 * t308;
t490 = t655 / 0.2e1 - t662 / 0.2e1;
t284 = t490 * t680;
t376 = pkin(23) + pkin(22);
t373 = sin(t376);
t374 = cos(t376);
t236 = t283 * t374 + t284 * t373;
t781 = pkin(5) * t236;
t862 = 2 * pkin(4);
t648 = t781 * t862 + t405;
t232 = t406 + t648;
t229 = 0.1e1 / t232;
t398 = 0.1e1 / pkin(11);
t698 = t229 * t398;
t642 = pkin(9) ^ 2 - pkin(11) ^ 2;
t226 = t232 - t642;
t234 = pkin(4) * t236 + pkin(5);
t823 = (-pkin(9) - pkin(11));
t222 = ((pkin(4) - t823) * (pkin(4) + t823)) + t648;
t822 = (pkin(11) - pkin(9));
t223 = ((pkin(4) - t822) * (pkin(4) + t822)) + t648;
t701 = t223 * t222;
t412 = sqrt(-t701);
t562 = t283 * t373 - t284 * t374;
t173 = pkin(4) * t226 * t562 + t234 * t412;
t384 = cos(pkin(21));
t711 = t173 * t384;
t852 = t562 * t412;
t171 = -pkin(4) * t852 + t226 * t234;
t381 = sin(pkin(21));
t714 = t171 * t381;
t140 = (t714 / 0.2e1 + t711 / 0.2e1) * t698;
t136 = t140 ^ 2;
t712 = t173 * t381;
t713 = t171 * t384;
t141 = (-t713 / 0.2e1 + t712 / 0.2e1) * t698;
t138 = 0.1e1 / t141 ^ 2;
t118 = t136 * t138 + 0.1e1;
t116 = 0.1e1 / t118;
t137 = 0.1e1 / t141;
t360 = t367 * qJD(2);
t515 = -t655 + t662;
t349 = 0.1e1 / t351 ^ 2;
t628 = pkin(1) * pkin(7) * t349;
t451 = t515 * t628;
t806 = -t367 / 0.2e1;
t636 = 0.2e1 * pkin(1);
t508 = pkin(7) * (t334 + t335) * t636;
t315 = t360 * t508;
t326 = 0.1e1 / t414;
t853 = t315 * t326;
t548 = t806 * t853;
t293 = pkin(7) * t548;
t633 = t356 * t842;
t581 = -t342 + t633;
t359 = t489 * qJD(2);
t678 = t359 * t414;
t260 = t293 + (t360 * t581 - t678) * pkin(7);
t659 = t391 * t260;
t677 = t360 * t367;
t790 = pkin(1) * t403;
t557 = t677 * t790;
t808 = t326 / 0.2e1;
t592 = t356 * t808;
t676 = t360 * t414;
t261 = t315 * t592 - 0.2e1 * t557 + (t342 * t359 - t676) * pkin(7);
t666 = t387 * t261;
t194 = (t360 * t451 + (-t659 / 0.2e1 + t666 / 0.2e1 + t492 * qJD(3)) * t348) * t408;
t514 = t654 + t663;
t658 = t391 * t261;
t667 = t387 * t260;
t195 = (t514 * t360 * t628 + (t667 / 0.2e1 + t658 / 0.2e1 + t490 * qJD(3)) * t348) * t408;
t177 = t194 * t373 - t195 * t374;
t634 = 0.2e1 * pkin(5);
t507 = pkin(4) * (t222 + t223) * t634;
t143 = t177 * t507;
t563 = t194 * t374 + t195 * t373;
t780 = pkin(5) * t562;
t582 = -0.2e1 * t406 * t780;
t199 = 0.1e1 / t412;
t814 = t199 / 0.2e1;
t601 = t234 * t814;
t708 = t177 * t412;
t100 = t143 * t601 + t177 * t582 + (t226 * t563 - t708) * pkin(4);
t230 = 0.1e1 / t232 ^ 2;
t787 = pkin(4) * t230;
t627 = pkin(5) * t787;
t449 = (t711 + t714) * t627;
t800 = t384 / 0.2e1;
t804 = t381 / 0.2e1;
t812 = -t562 / 0.2e1;
t600 = t199 * t812;
t851 = t563 * t412;
t480 = t143 * t600 - t851;
t579 = -0.2e1 * pkin(5) * t234 - t226;
t98 = (t177 * t579 + t480) * pkin(4);
t72 = ((t100 * t800 + t804 * t98) * t229 + t177 * t449) * t398;
t722 = t138 * t140;
t450 = (t712 - t713) * t627;
t801 = -t384 / 0.2e1;
t73 = ((t100 * t804 + t801 * t98) * t229 + t177 * t450) * t398;
t49 = (t137 * t72 - t722 * t73) * t116;
t119 = atan2(t140, t141);
t115 = cos(t119);
t371 = pkin(1) * t387 + pkin(5);
t114 = sin(t119);
t621 = pkin(1) * t114 * t391;
t89 = t115 * t371 + t621;
t863 = t49 * t89;
t386 = sin(qJ(4));
t390 = cos(qJ(4));
t643 = t386 ^ 2 + t390 ^ 2;
t180 = t562 * t507;
t134 = t180 * t601 + t562 * t582 + (t226 * t236 - t852) * pkin(4);
t460 = t100 * t562 + t134 * t177 + t173 * t563;
t478 = t180 * t600 - t236 * t412;
t132 = (t562 * t579 + t478) * pkin(4);
t463 = t132 * t177 + t171 * t563 + t562 * t98;
t231 = t229 * t230;
t709 = t177 * t406;
t610 = t405 * t709;
t555 = t231 * t610;
t504 = 0.4e1 * t381 * t555;
t483 = t562 * t504;
t526 = t384 * t555;
t485 = 0.4e1 * t173 * t526;
t486 = -0.4e1 * t171 * t526;
t92 = ((t132 * t804 + t134 * t800) * t229 + t562 * t449) * t398;
t93 = ((t132 * t801 + t134 * t804) * t229 + t562 * t450) * t398;
t502 = -t137 * t92 + t722 * t93;
t746 = t137 * t138 * t73;
t568 = 0.2e1 * t140 * t746;
t572 = t116 * t627;
t699 = t229 * t384;
t597 = t699 / 0.2e1;
t598 = -t699 / 0.2e1;
t599 = t229 * t804;
t816 = -t180 / 0.2e1;
t818 = -t143 / 0.2e1;
t444 = (t177 * t816 + t562 * t818) * t199 - t851;
t564 = -0.8e1 * t610;
t130 = t507 * t563 + t562 * t564;
t721 = t143 * t199 / t701;
t603 = t721 / 0.4e1;
t472 = t130 * t814 + t180 * t603;
t858 = -0.2e1 * t562;
t521 = t563 * t858;
t585 = t406 * t634;
t71 = t472 * t234 + (-t177 * t236 + t521) * t585 + (-t177 * t226 + t444) * pkin(4);
t549 = -t562 * t721 / 0.4e1;
t436 = (t130 * t812 + t236 * t818 + t563 * t816) * t199 + t708 + t180 * t549;
t583 = 0.4e1 * pkin(5) * t709;
t74 = t562 * t583 + (t563 * t579 + t436) * pkin(4);
t747 = 0.2e1 / t118 ^ 2 * (-t136 * t746 + t72 * t722);
t18 = t502 * t747 + (t93 * t568 + (-t72 * t93 - t73 * t92) * t138) * t116 + (((t171 * t483 + t485 * t562 + t597 * t71 + t599 * t74) * t137 - (t173 * t483 + t486 * t562 + t598 * t74 + t599 * t71) * t722) * t116 + ((t137 * t460 + t463 * t722) * t384 + (t137 * t463 - t460 * t722) * t381) * t572) * t398;
t56 = -t116 * t502 + 0.1e1;
t637 = qJD(4) * t390;
t498 = t18 * t386 + t56 * t637;
t481 = t367 * t514;
t318 = t367 * t508;
t837 = -0.2e1 * t367 ^ 2;
t268 = t318 * t592 + t790 * t837 + (t342 * t489 - t672) * pkin(7);
t656 = t391 * t268;
t809 = -t326 / 0.2e1;
t593 = t318 * t809;
t545 = t593 - t342;
t673 = t489 * t414;
t265 = (-t673 + (t545 + t633) * t367) * pkin(7);
t665 = t387 * t265;
t493 = t665 / 0.2e1 + t656 / 0.2e1;
t216 = (-t348 * t493 - t481 * t628) * t408;
t657 = t391 * t265;
t664 = t387 * t268;
t491 = -t657 / 0.2e1 + t664 / 0.2e1;
t217 = (t348 * t491 + t367 * t451) * t408;
t185 = t216 * t374 + t217 * t373;
t148 = t185 * t507;
t186 = -t216 * t373 + t217 * t374;
t706 = t185 * t412;
t105 = t148 * t601 + t185 * t582 + (t186 * t226 - t706) * pkin(4);
t605 = t411 * t677;
t553 = t403 * t605;
t299 = t359 * t508 - 0.8e1 * t553;
t607 = 0.1e1 / t683 * t318 * t853;
t470 = t607 / 0.4e1 + t299 * t808;
t512 = -0.4e1 * t359 * t367 - 0.2e1 * t360 * t489;
t211 = t293 + t470 * t356 + t512 * t790 + (t360 * t545 - t678) * pkin(7);
t437 = t676 + (-t489 * t315 / 0.2e1 - t359 * t318 / 0.2e1 + t299 * t806) * t326 - t367 * t607 / 0.4e1;
t213 = 0.4e1 * t557 + (t359 * t581 + t437) * pkin(7);
t350 = t348 * t349;
t525 = t350 * t553;
t487 = 0.4e1 * t525;
t797 = -t391 / 0.2e1;
t149 = (t515 * t487 + (t213 * t797 + t387 * t211 / 0.2e1 + t493 * qJD(3)) * t348 + ((-t657 + t664) * t360 + t515 * t359 + (qJD(3) * t514 - t659 + t666) * t367) * t628) * t408;
t859 = -0.4e1 * t411;
t150 = (t403 * t360 * t350 * t481 * t859 + (-t387 * t213 / 0.2e1 + t211 * t797 + t491 * qJD(3)) * t348 + ((-t656 - t665) * t360 - t514 * t359 + (qJD(3) * t515 - t658 - t667) * t367) * t628) * t408;
t129 = t149 * t373 + t150 * t374;
t461 = t100 * t185 + t105 * t177 + t129 * t173;
t479 = t148 * t600 - t186 * t412;
t103 = (t185 * t579 + t479) * pkin(4);
t464 = t103 * t177 + t129 * t171 + t185 * t98;
t484 = t185 * t504;
t77 = ((t103 * t804 + t105 * t800) * t229 + t185 * t449) * t398;
t78 = ((t103 * t801 + t105 * t804) * t229 + t185 * t450) * t398;
t503 = -t137 * t77 + t722 * t78;
t128 = t149 * t374 - t150 * t373;
t817 = -t148 / 0.2e1;
t445 = (t177 * t817 + t185 * t818) * t199 - t129 * t412;
t91 = t129 * t507 + t185 * t564;
t475 = t148 * t603 + t814 * t91;
t522 = -t129 * t562 - t185 * t563;
t57 = t475 * t234 + (-t177 * t186 + t522) * t585 + (t128 * t226 + t445) * pkin(4);
t438 = (t186 * t818 + t563 * t817 + t812 * t91) * t199 - t128 * t412 + t148 * t549;
t58 = t185 * t583 + (t129 * t579 + t438) * pkin(4);
t11 = t503 * t747 + (t78 * t568 + (-t72 * t78 - t73 * t77) * t138) * t116 + (((t171 * t484 + t185 * t485 + t57 * t597 + t58 * t599) * t137 - (t173 * t484 + t185 * t486 + t57 * t599 + t58 * t598) * t722) * t116 + ((t137 * t461 + t464 * t722) * t384 + (t137 * t464 - t461 * t722) * t381) * t572) * t398;
t51 = -t116 * t503 + 0.1e1;
t500 = t11 * t386 + t51 * t637;
t109 = t114 * t371;
t723 = t115 * t391;
t36 = -pkin(1) * (qJD(3) * (t114 * t387 - t723) - t49 * t723) - t49 * t109;
t635 = 0.2e1 * pkin(2);
t791 = pkin(1) * t392;
t372 = qJD(2) * t791;
t365 = t387 * t392 + t388 * t391;
t477 = t365 * qJD(2);
t440 = -qJD(3) * t365 - t477;
t328 = -pkin(5) * t440 + t372;
t861 = 0.2e1 * t328;
t369 = pkin(1) * t388 - pkin(16);
t860 = 0.2e1 * t369;
t225 = t232 + t642;
t233 = -pkin(4) - t781;
t170 = -pkin(5) * t852 - t225 * t233;
t167 = 0.1e1 / t170 ^ 2;
t172 = t225 * t780 - t233 * t412;
t169 = t172 ^ 2;
t146 = t167 * t169 + 0.1e1;
t715 = t167 * t172;
t166 = 0.1e1 / t170;
t580 = t233 * t862 - t225;
t97 = (t177 * t580 + t480) * pkin(5);
t745 = t166 * t167 * t97;
t784 = pkin(4) * t405;
t584 = t784 * t858;
t602 = -t199 * t233 / 0.2e1;
t854 = t225 * t563;
t99 = t143 * t602 + t177 * t584 + (-t708 + t854) * pkin(5);
t857 = 0.2e1 / t146 ^ 2 * (-t169 * t745 + t715 * t99);
t355 = t362 * qJD(3);
t509 = (t332 + t333) * pkin(6) * t636;
t310 = t355 * t509;
t354 = t361 * qJD(3);
t513 = t340 * t354 - t355 * t413;
t323 = 0.1e1 / t413;
t811 = -t323 / 0.2e1;
t594 = t352 * t811;
t675 = t362 * t404;
t606 = t355 * t675;
t258 = pkin(6) * t513 + t310 * t594 - t606 * t636;
t689 = t301 * t305;
t686 = t310 * t362;
t291 = pkin(6) * t686 * t811;
t679 = t354 * t413;
t763 = t352 * pkin(1);
t257 = t291 + (-t679 - (t340 - 0.2e1 * t763) * t355) * pkin(6);
t692 = t257 * t300 * t301;
t856 = 0.2e1 * (t258 * t689 - t303 * t692) / t290 ^ 2;
t292 = pkin(1) * t548;
t341 = t351 - t640;
t357 = -pkin(1) * t489 - pkin(7);
t632 = 0.2e1 * t357 * pkin(7);
t578 = -t341 + t632;
t259 = t292 + (t360 * t578 - t678) * pkin(1);
t556 = pkin(7) * t605;
t591 = t357 * t809;
t262 = t315 * t591 - 0.2e1 * t556 + (t341 * t359 - t676) * pkin(1);
t402 = 0.1e1 / pkin(8);
t306 = -pkin(1) * t672 - t341 * t357;
t394 = cos(pkin(18));
t652 = t394 * t306;
t309 = pkin(1) * t341 * t367 - t357 * t414;
t389 = sin(pkin(18));
t660 = t389 * t309;
t455 = (t652 - t660) * t628;
t796 = t394 / 0.2e1;
t799 = -t389 / 0.2e1;
t206 = ((t259 * t796 + t262 * t799) * t348 + t360 * t455) * t402;
t681 = t348 * t402;
t285 = (t652 / 0.2e1 - t660 / 0.2e1) * t681;
t280 = 0.1e1 / t285 ^ 2;
t855 = t206 * t280;
t400 = 0.1e1 / pkin(9);
t813 = t229 / 0.2e1;
t596 = t400 * t813;
t446 = atan2(t172 * t596, t170 * t596);
t142 = cos(t446);
t380 = sin(pkin(22));
t383 = cos(pkin(22));
t443 = sin(t446);
t113 = -t142 * t383 - t380 * t443;
t382 = cos(pkin(23));
t670 = t382 * t307;
t379 = sin(pkin(23));
t671 = t379 * t308;
t517 = -t670 + t671;
t453 = t517 * t628;
t803 = -t382 / 0.2e1;
t805 = t379 / 0.2e1;
t204 = ((t260 * t803 + t261 * t805) * t348 + t360 * t453) * t408;
t669 = t382 * t308;
t687 = t307 * t379;
t516 = t669 + t687;
t452 = t516 * t628;
t802 = t382 / 0.2e1;
t205 = ((t260 * t805 + t261 * t802) * t348 + t360 * t452) * t408;
t214 = ((t265 * t803 + t268 * t805) * t348 + t367 * t453) * t408;
t215 = ((t265 * t805 + t268 * t802) * t348 + t367 * t452) * t408;
t850 = -t204 * t215 - t205 * t214;
t277 = (-t670 / 0.2e1 + t671 / 0.2e1) * t680;
t274 = 0.1e1 / t277 ^ 2;
t278 = (t669 / 0.2e1 + t687 / 0.2e1) * t680;
t276 = t278 ^ 2;
t244 = t274 * t276 + 0.1e1;
t242 = 0.1e1 / t244;
t273 = 0.1e1 / t277;
t691 = t274 * t278;
t159 = (-t204 * t691 + t205 * t273) * t242;
t848 = qJD(2) + t159;
t344 = 0.1e1 / t347 ^ 2;
t629 = t355 * t834;
t574 = t344 * t629;
t343 = 0.1e1 / t347;
t807 = t343 / 0.2e1;
t847 = qJD(2) + (t258 * t807 + t305 * t574) * t865 + (t257 * t807 + t304 * t574) * t866;
t795 = t410 / 0.2e1;
t590 = t343 * t795;
t272 = atan2(t305 * t590, t304 * t590);
t269 = sin(t272);
t270 = cos(t272);
t518 = t269 * t388 - t270 * t392;
t153 = t847 * t518;
t227 = -t269 * t392 - t270 * t388;
t154 = t847 * t227;
t684 = t333 * t337;
t321 = -t332 * t684 + 0.1e1;
t319 = 0.1e1 / t321;
t336 = 0.1e1 / t339;
t810 = t323 / 0.2e1;
t595 = t336 * t810;
t254 = (t310 * t595 + t355 * t864) * t319;
t588 = 0.1e1 / pkin(13) * t795;
t316 = atan2(t413 * t588, t339 * t588);
t313 = sin(t316);
t314 = cos(t316);
t520 = t227 * t313 - t314 * t518;
t124 = -t153 * t314 + t154 * t313 + t254 * t520;
t190 = -t227 * t314 - t313 * t518;
t125 = -t153 * t313 - t154 * t314 + t190 * t254;
t820 = pkin(2) * mrSges(10,3);
t846 = Ifges(9,5) * t154 + Ifges(9,6) * t153 + (-t124 * t313 + t125 * t314 + (-t190 * t314 + t313 * t520) * t254) * t820;
t144 = 0.1e1 / t146;
t550 = pkin(4) * pkin(5) * pkin(9) * t144 * t177;
t774 = pkin(9) * t232;
t618 = t144 * t774;
t559 = t167 * t618;
t617 = t166 * t774;
t845 = -0.2e1 * t166 * t550 - t559 * t97 - t617 * t857;
t616 = t172 * t774;
t558 = t167 * t616;
t843 = 0.2e1 * t144 * t616 * t745 + 0.2e1 * t550 * t715 + t558 * t857 - t559 * t99;
t841 = -0.2e1 * pkin(16);
t366 = -t387 * t388 + t391 * t392;
t331 = (qJD(2) + qJD(3)) * t366;
t567 = t366 * t49 + t331;
t725 = t115 * t365;
t33 = t114 * t567 - t115 * t440 + t49 * t725;
t34 = t567 * t115 + (-t365 * t49 + t440) * t114;
t24 = pkin(10) * t33 - pkin(12) * t34 + t328;
t840 = 0.2e1 * t24;
t279 = 0.1e1 / t285;
t247 = atan2(t278, t277);
t240 = sin(t247);
t241 = cos(t247);
t519 = t240 * t388 - t241 * t392;
t838 = -0.2e1 * t519;
t836 = m(5) / 0.2e1;
t835 = m(6) / 0.2e1;
t833 = -t34 / 0.2e1;
t832 = t34 / 0.2e1;
t749 = Ifges(6,4) * t390;
t87 = t114 * t365 - t115 * t366;
t756 = t87 * Ifges(6,6);
t88 = t114 * t366 + t725;
t59 = t756 + (-Ifges(6,2) * t386 + t749) * t88;
t831 = t59 / 0.2e1;
t750 = Ifges(6,4) * t386;
t772 = Ifges(6,5) * t87;
t60 = t772 + (Ifges(6,1) * t390 - t750) * t88;
t830 = t60 / 0.2e1;
t829 = t87 / 0.2e1;
t828 = -t88 / 0.2e1;
t827 = t88 / 0.2e1;
t317 = t362 * t509;
t547 = t317 * t810 + t340;
t263 = (-t361 * t413 + (-t547 + 0.2e1 * t763) * t362) * pkin(6);
t792 = pkin(1) * t344;
t573 = t777 * t792;
t252 = (t263 * t807 + t304 * t573) * t410;
t266 = t317 * t594 + t404 * t362 ^ 2 * t842 + (t340 * t361 - t674) * pkin(6);
t253 = (t266 * t807 + t305 * t573) * t410;
t188 = t252 * t511 + t253 * t867;
t448 = t317 * t595 + t362 * t864;
t182 = t319 * t448 + t188;
t815 = t182 / 0.2e1;
t798 = t389 / 0.2e1;
t793 = pkin(1) * t343;
t788 = pkin(4) * t229;
t786 = pkin(4) * t380;
t785 = pkin(4) * t383;
t783 = pkin(5) * t114;
t782 = pkin(5) * t115;
t779 = pkin(5) * t400;
t775 = pkin(7) * t411;
t773 = mrSges(6,3) * t88;
t771 = Ifges(5,3) * t11;
t770 = Ifges(5,3) * t18;
t768 = t33 * Ifges(6,5);
t766 = t33 * Ifges(6,6);
t639 = qJD(3) * t387;
t35 = t115 * pkin(1) * t639 + qJD(3) * t621 + t863;
t764 = t35 * mrSges(5,2);
t762 = t36 * mrSges(5,1);
t638 = qJD(4) * t386;
t734 = t390 * t11;
t499 = t51 * t638 - t734;
t3 = Ifges(6,1) * t500 - Ifges(6,4) * t499;
t761 = t386 * t3;
t743 = t18 * t390;
t497 = t56 * t638 - t743;
t6 = Ifges(6,1) * t498 - Ifges(6,4) * t497;
t760 = t386 * t6;
t2 = Ifges(6,4) * t500 - Ifges(6,2) * t499;
t759 = t390 * t2;
t5 = Ifges(6,4) * t498 - Ifges(6,2) * t497;
t758 = t390 * t5;
t754 = t89 * mrSges(5,1);
t90 = -pkin(1) * t723 + t109;
t753 = t90 * mrSges(5,2);
t733 = t390 * t34;
t752 = Ifges(6,5) * t733 + Ifges(6,3) * t33;
t751 = mrSges(6,1) * t390;
t554 = t411 * t606;
t298 = t354 * t509 - 0.8e1 * t554;
t608 = t323 / t685 * t317 * t310;
t471 = t608 / 0.4e1 + t298 * t810;
t543 = t619 * t868;
t551 = pkin(2) * t288 * t629;
t623 = t343 * t344 * t859;
t106 = ((-0.4e1 * pkin(1) * t354 * t675 - t352 * t471 + t291) * t807 - (t305 * t362 * t623 + t361 * t793) * t404 * t355 + ((-t355 * t547 - t679) * t807 + (t258 * t362 + t266 * t355 + t305 * t354) * t792) * pkin(6)) * t865 + (-(t304 * t623 - 0.2e1 * t793) * t606 + ((-t362 * t608 / 0.4e1 + (-t361 * t310 / 0.2e1 - t354 * t317 / 0.2e1 - t362 * t298 / 0.2e1) * t323 - t513) * t807 + (t352 * t354 * t343 + (t257 * t362 + t263 * t355 + t304 * t354) * t344) * pkin(1)) * pkin(6)) * t866 + (t257 * t543 - 0.4e1 * t300 * t551 - 0.2e1 * t620 * t856) * t253 + (0.4e1 * t288 * t305 * t692 * t789 + t258 * t543 + 0.4e1 * t551 * t689 + 0.2e1 * t560 * t856) * t252;
t338 = t336 * t337;
t101 = t448 / t321 ^ 2 * (-0.2e1 * t684 + (-0.4e1 * t333 * t338 - 0.2e1 * t337) * t332) * t629 + (0.8e1 * t413 * t338 * t554 + t471 * t336 + (-0.2e1 * t679 + (-t317 * t355 - t686) * t323) * t630) * t319 + t106;
t748 = Ifges(10,3) * t101;
t739 = t386 * t34;
t528 = Ifges(6,2) * t390 + t750;
t38 = t528 * t51;
t738 = t386 * t38;
t529 = Ifges(6,1) * t386 + t749;
t39 = t529 * t51;
t737 = t386 * t39;
t41 = t528 * t56;
t736 = t386 * t41;
t42 = t529 * t56;
t735 = t386 * t42;
t732 = t390 * t38;
t731 = t390 * t39;
t730 = t390 * t41;
t729 = t390 * t42;
t720 = t159 * t240;
t719 = t159 * t241;
t651 = t394 * t309;
t688 = t306 * t389;
t454 = (t651 + t688) * t628;
t207 = ((t259 * t798 + t262 * t796) * t348 + t360 * t454) * t402;
t286 = (t651 / 0.2e1 + t688 / 0.2e1) * t681;
t282 = t286 ^ 2;
t250 = t280 * t282 + 0.1e1;
t248 = 0.1e1 / t250;
t690 = t280 * t286;
t160 = (-t206 * t690 + t207 * t279) * t248;
t546 = t593 - t341;
t264 = (-t673 + (t546 + t632) * t367) * pkin(1);
t267 = t318 * t591 + t775 * t837 + (t341 * t489 - t672) * pkin(1);
t218 = ((t264 * t796 + t267 * t799) * t348 + t367 * t455) * t402;
t219 = ((t264 * t798 + t267 * t796) * t348 + t367 * t454) * t402;
t162 = (-t218 * t690 + t219 * t279) * t248;
t718 = t160 * t162;
t704 = t204 * t273 * t274;
t717 = 0.2e1 * (t205 * t691 - t276 * t704) / t244 ^ 2;
t702 = t279 * t855;
t716 = 0.2e1 * (t207 * t690 - t282 * t702) / t250 ^ 2;
t710 = t177 * t405;
t707 = t182 * t254;
t700 = t229 * t233;
t695 = t242 * t273;
t694 = t242 * t274;
t682 = t348 * t389;
t650 = t101 + t106;
t649 = Ifges(10,5) * t125 + Ifges(10,6) * t124;
t647 = Ifges(4,5) * t331 + Ifges(4,6) * t440;
t346 = -pkin(5) * t366 + t369;
t66 = pkin(10) * t87 - pkin(12) * t88 + t346;
t631 = 0.2e1 * t66;
t625 = Ifges(6,5) * t500 + Ifges(6,6) * t734;
t624 = 0.4e1 * t231 * t406;
t622 = Ifges(6,5) * t498 + Ifges(6,6) * t743;
t609 = t242 * t691;
t589 = t348 * t796;
t587 = -t59 - t756;
t586 = 0.2e1 * t784;
t577 = mrSges(6,3) * t643;
t571 = t177 * t627;
t570 = t185 * t627;
t569 = t562 * t627;
t566 = t254 * (t182 + t188);
t565 = t172 * t624;
t552 = 0.4e1 * t643;
t544 = 0.2e1 * t618;
t495 = t638 * t88 - t733;
t496 = t637 * t88 + t739;
t14 = -Ifges(6,4) * t495 - Ifges(6,2) * t496 + t766;
t542 = t14 / 0.2e1 + t766 / 0.2e1;
t15 = -Ifges(6,1) * t495 - Ifges(6,4) * t496 + t768;
t541 = t15 / 0.2e1 + t768 / 0.2e1;
t540 = -t59 / 0.2e1 - t756 / 0.2e1;
t534 = t144 * t400 * t617;
t533 = t144 * t558;
t532 = mrSges(5,1) * t115 - mrSges(5,2) * t114;
t531 = -mrSges(6,2) * t386 + t751;
t530 = mrSges(6,1) * t386 + mrSges(6,2) * t390;
t527 = 0.2e1 * t577;
t201 = -t240 * t392 - t241 * t388;
t510 = m(10) * t409 * (t313 ^ 2 + t314 ^ 2) / 0.2e1;
t501 = t389 * t525;
t161 = -t214 * t609 + t215 * t695 + 0.1e1;
t488 = qJD(3) * (mrSges(4,1) * t391 - mrSges(4,2) * t387);
t482 = t394 * t487;
t476 = (-mrSges(5,1) * t114 - mrSges(5,2) * t115) * t49 * pkin(5);
t474 = Ifges(5,5) * t34 - Ifges(5,6) * t33;
t473 = Ifges(5,5) * t88 - Ifges(5,6) * t87;
t102 = (t185 * t580 + t479) * pkin(5);
t94 = (t102 * t813 + t170 * t570) * t400;
t104 = t148 * t602 + t185 * t584 + (t186 * t225 - t706) * pkin(5);
t95 = (t104 * t813 + t172 * t570) * t400;
t64 = (t166 * t95 - t715 * t94) * t544 + t161;
t156 = pkin(1) * t240 + t161 * t786;
t157 = pkin(1) * t241 + t161 * t785;
t439 = -t142 * t380 + t383 * t443;
t79 = t113 * t157 + t156 * t439;
t80 = t113 * t156 - t157 * t439;
t469 = mrSges(11,1) * t79 - mrSges(11,2) * t80 + Ifges(11,3) * t64;
t126 = t848 * t519;
t127 = t848 * t201;
t63 = 0.2e1 * (t172 * t571 + t813 * t99) * t534 - 0.2e1 * (t170 * t571 + t813 * t97) * t400 * t533;
t47 = t113 * t63;
t48 = t439 * t63;
t29 = t113 * t126 + t127 * t439 + t201 * t48 + t47 * t519;
t30 = t113 * t127 - t126 * t439 + t201 * t47 - t48 * t519;
t468 = Ifges(11,5) * t30 + Ifges(11,6) * t29;
t83 = t113 * t201 - t439 * t519;
t84 = -t113 * t519 - t201 * t439;
t467 = Ifges(11,5) * t84 + Ifges(11,6) * t83;
t466 = t366 * Ifges(4,2) * t365;
t465 = t369 * t365 * mrSges(4,1);
t457 = t261 * t367 + t268 * t360 + t308 * t359;
t458 = t260 * t367 + t265 * t360 + t307 * t359;
t462 = -t215 * t273 * t717 + (t516 * t487 + (t211 * t802 + t213 * t805) * t348 + (t379 * t458 + t382 * t457) * t628) * t408 * t695 + (0.2e1 * t242 * t278 * t704 + t691 * t717) * t214;
t459 = t259 * t367 + t264 * t360 + t306 * t359;
t456 = t262 * t367 + t267 * t360 + t309 * t359;
t447 = t400 * (t170 * t624 + 0.2e1 * t788) * t710;
t251 = atan2(t286, t285);
t246 = cos(t251);
t245 = sin(t251);
t224 = -pkin(2) * t227 - pkin(16);
t212 = 0.4e1 * t556 + (t359 * t578 + t437) * pkin(1);
t210 = t292 - t470 * t357 + t512 * t775 + (t360 * t546 - t678) * pkin(1);
t158 = (t517 * t487 + (t211 * t805 + t213 * t803) * t348 + (t379 * t457 - t382 * t458) * t628) * t408;
t133 = t180 * t602 + t562 * t584 + (t225 * t236 - t852) * pkin(5);
t131 = (t562 * t580 + t478) * pkin(5);
t121 = (t133 * t813 + t172 * t569) * t400;
t120 = (t131 * t813 + t170 * t569) * t400;
t96 = t372 + (-t126 * t383 - t127 * t380) * pkin(4);
t86 = (-t248 * t855 - t279 * t716) * t219 + (t690 * t716 + (-t207 * t280 + 0.2e1 * t286 * t702) * t248) * t218 + ((t210 * t589 + t309 * t482 + t212 * t682 / 0.2e1 + 0.4e1 * t306 * t501) * t279 - (t212 * t589 + t306 * t482 - t210 * t682 / 0.2e1 - 0.4e1 * t309 * t501) * t690 + ((t279 * t456 - t459 * t690) * t394 + (t279 * t459 + t456 * t690) * t389) * t628) * t248 * t402;
t85 = (-t158 * t278 + t850) * t694 + t462;
t82 = pkin(1) * t719 + t786 * t85;
t81 = -pkin(1) * t720 + t785 * t85;
t76 = (-t120 * t715 + t121 * t166) * t544;
t69 = mrSges(6,1) * t87 - t390 * t773;
t68 = -mrSges(6,2) * t87 - t386 * t773;
t67 = t530 * t88;
t54 = -pkin(10) * t56 - t782;
t53 = pkin(12) * t56 + t783;
t46 = pkin(12) * t51 + t90;
t45 = -pkin(10) * t51 - t89;
t40 = t531 * t56;
t37 = t531 * t51;
t28 = t113 * t81 - t156 * t47 + t157 * t48 + t439 * t82;
t27 = t113 * t82 + t156 * t48 + t157 * t47 - t439 * t81;
t26 = 0.2e1 * ((-t233 * t472 + t521 * t586) * t813 + (-t236 * t788 + t562 * t565) * t710 + ((-t177 * t225 + t444) * t813 + (t133 * t177 + t172 * t563 + t562 * t99) * t787) * pkin(5)) * t534 - 0.2e1 * (t562 * t447 + ((t436 - t854) * t813 + (t563 * t700 + (t131 * t177 + t170 * t563 + t562 * t97) * t230) * pkin(4)) * t779) * t533 + 0.2e1 * t845 * t121 + 0.2e1 * t843 * t120;
t22 = -mrSges(6,2) * t33 - mrSges(6,3) * t496;
t21 = mrSges(6,1) * t33 + mrSges(6,3) * t495;
t20 = -t158 * t609 + 0.2e1 * ((-t233 * t475 + t522 * t586) * t813 + (t185 * t565 - t186 * t788) * t710 + ((t128 * t225 + t445) * t813 + (t104 * t177 + t129 * t172 + t185 * t99) * t787) * pkin(5)) * t534 - 0.2e1 * (t185 * t447 + ((-t129 * t225 + t438) * t813 + (t129 * t700 + (t102 * t177 + t129 * t170 + t185 * t97) * t230) * pkin(4)) * t779) * t533 + t462 + t850 * t694 + 0.2e1 * t843 * t94 + 0.2e1 * t845 * t95;
t19 = mrSges(6,1) * t496 - mrSges(6,2) * t495;
t13 = pkin(12) * t18 + t49 * t782;
t12 = -pkin(10) * t18 + t49 * t783;
t8 = -pkin(10) * t11 - t36;
t7 = pkin(12) * t11 + t35;
t4 = mrSges(6,1) * t497 + mrSges(6,2) * t498;
t1 = mrSges(6,1) * t499 + mrSges(6,2) * t500;
t9 = [-0.2e1 * t518 * t154 * Ifges(9,1) + (-0.2e1 * t466 + 0.2e1 * t465 + (mrSges(3,1) * t392 - mrSges(3,2) * t388) * t841 + 0.2e1 * (t388 ^ 2 - t392 ^ 2) * Ifges(3,4) + 0.2e1 * (-Ifges(3,1) + Ifges(3,2)) * t388 * t392 + (-0.2e1 * mrSges(4,1) * t366 - 0.2e1 * mrSges(8,1) * t201 + 0.2e1 * mrSges(4,2) * t365 + mrSges(8,2) * t838 + (m(4) + m(8)) * t860) * t791) * qJD(2) + 0.2e1 * (-t466 + t465) * qJD(3) + (-m(10) * t224 + mrSges(10,1) * t190 + mrSges(10,2) * t520) * t153 * t635 - 0.2e1 * t520 * t125 * Ifges(10,1) + t66 * t24 * t552 * t835 + (-mrSges(8,1) * t126 + mrSges(8,2) * t127) * t860 + (m(5) * t346 + mrSges(5,2) * t88) * t861 + 0.2e1 * t224 * (-mrSges(10,1) * t124 + mrSges(10,2) * t125) + t127 * Ifges(8,1) * t838 + 0.2e1 * t124 * Ifges(10,2) * t190 + 0.2e1 * t153 * Ifges(9,2) * t227 + 0.2e1 * (-t124 * t520 + t125 * t190) * Ifges(10,4) + 0.2e1 * (-t126 * t519 + t127 * t201) * Ifges(8,4) + 0.2e1 * (-t153 * t518 + t154 * t227) * Ifges(9,4) + 0.2e1 * t126 * Ifges(8,2) * t201 + 0.2e1 * t29 * Ifges(11,2) * t83 + (t69 * t840 + t34 * t60 + (qJD(4) * t68 + t21) * t631 + (qJD(4) * t587 + t15 + t768) * t88) * t390 + (t68 * t840 + (-qJD(4) * t69 + t22) * t631 + t587 * t34 + (-t766 - t14 + (-t60 - t772) * qJD(4)) * t88) * t386 + 0.2e1 * t96 * (-mrSges(11,1) * t83 + mrSges(11,2) * t84) + 0.2e1 * t88 * t34 * Ifges(5,1) + 0.2e1 * (t331 * t366 + t365 * t440) * Ifges(4,4) + 0.2e1 * (mrSges(4,2) * t369 + Ifges(4,1) * t365) * t331 + 0.2e1 * (t29 * t84 + t30 * t83) * Ifges(11,4) + 0.2e1 * (-t33 * t88 - t34 * t87) * Ifges(5,4) + 0.2e1 * t84 * t30 * Ifges(11,1) + (mrSges(5,1) * t861 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t33 + t752) * t87 + (-mrSges(9,1) * t153 + mrSges(9,2) * t154) * t841 + 0.2e1 * t346 * (mrSges(5,1) * t33 + mrSges(5,2) * t34) + 0.2e1 * (pkin(15) * (mrSges(7,1) * t245 + mrSges(7,2) * t246) + (-t245 ^ 2 + t246 ^ 2) * Ifges(7,4) + (-Ifges(7,2) + Ifges(7,1)) * t245 * t246) * t160 + 0.2e1 * (m(11) * t96 - mrSges(11,1) * t29 + mrSges(11,2) * t30) * ((-t201 * t383 + t380 * t519) * pkin(4) + t369); (-Ifges(3,5) * t388 - Ifges(3,6) * t392) * qJD(2) + (t127 * t161 - t519 * t85) * Ifges(8,5) + ((t240 * t126 - t241 * t127 + (t201 * t241 - t240 * t519) * t159) * mrSges(8,3) + (-t387 * t331 + t366 * t639 + t391 * t477) * mrSges(4,3)) * pkin(1) + t649 + t647 + t45 * t19 + (t126 * t161 + t201 * t85) * Ifges(8,6) + t846 + (-t33 * t90 - t34 * t89 - t35 * t87 - t36 * t88) * mrSges(5,3) + (t27 * t83 - t28 * t84 + t29 * t80 - t30 * t79) * mrSges(11,3) + t8 * t67 + (t3 * t827 + t7 * t68 + t11 * t831 + t46 * t22 + t39 * t832 + t542 * t51 + (t38 * t828 - t46 * t69 + t51 * t830) * qJD(4)) * t390 + (t11 * t830 - t46 * t21 + t38 * t833 + t2 * t828 - t7 * t69 + t541 * t51 + (t39 * t828 - t46 * t68 + t51 * t540) * qJD(4)) * t386 + (-t245 * t718 + t246 * t86) * Ifges(7,6) + (t245 * t86 + t246 * t718) * Ifges(7,5) + t625 * t829 + t473 * t11 + t474 * t51 + t467 * t20 + t468 * t64; 0.2e1 * t45 * t1 + 0.2e1 * t161 * Ifges(8,3) * t85 + 0.2e1 * t162 * Ifges(7,3) * t86 - 0.2e1 * t8 * t37 + 0.4e1 * (t35 * t90 + t36 * t89) * t836 + (t46 * t552 * t7 + 0.4e1 * t45 * t8) * t835 + 0.2e1 * (m(11) * t79 + mrSges(11,1) * t64) * t28 + 0.2e1 * (m(11) * t80 - mrSges(11,2) * t64) * t27 + (mrSges(10,1) * t313 + mrSges(10,2) * t314) * t254 * t635 + 0.2e1 * t469 * t20 + (t46 * t527 + t732 + t737 - 0.2e1 * t753 + 0.2e1 * t754) * t11 + (t488 + (-t161 * t719 - t240 * t85) * mrSges(8,2) + (-t161 * t720 + t241 * t85) * mrSges(8,1)) * t636 + (0.2e1 * t762 - 0.2e1 * t764 + 0.2e1 * t771 + t759 + t761 + (t731 - t738) * qJD(4) + t7 * t527) * t51; t649 * t815 + t54 * t19 + t12 * t67 + t622 * t829 + (t101 * t190 + t124 * t815) * Ifges(10,6) + t468 * t76 + t474 * t56 + t467 * t26 + t473 * t18 + ((-t49 * t87 - t34) * t115 + (t49 * t88 - t33) * t114) * mrSges(5,3) * pkin(5) + t846 * t188 + (t5 * t828 - t13 * t69 + t18 * t830 - t53 * t21 + t41 * t833 + t541 * t56 + (t42 * t828 - t53 * t68 + t540 * t56) * qJD(4)) * t386 + (t18 * t831 + t53 * t22 + t42 * t832 + t6 * t827 + t13 * t68 + t542 * t56 + (t41 * t828 - t53 * t69 + t56 * t830) * qJD(4)) * t390 + (-t101 * t520 + t125 * t815) * Ifges(10,5) + (-Ifges(9,5) * t518 + Ifges(9,6) * t227 + (-t190 * t313 - t314 * t520) * t820) * t106 + t647; t45 * t4 + t748 - t8 * t40 + (m(6) * t8 + t1) * t54 + (m(6) * t45 - t37) * t12 + pkin(1) * t488 + (mrSges(11,1) * t28 - mrSges(11,2) * t27 + Ifges(11,3) * t20) * t76 + t469 * t26 + (Ifges(9,3) + 0.2e1 * t510) * t106 + (t730 / 0.2e1 + t735 / 0.2e1 + t53 * t577) * t11 + (t732 / 0.2e1 + t737 / 0.2e1 + t754 - t753 + t46 * t577) * t18 + (t532 * t11 + 0.2e1 * ((t49 * t90 + t36) * t115 + (t35 - t863) * t114) * t836) * pkin(5) + ((t313 * t650 + t314 * t566) * mrSges(10,2) + (t313 * t566 - t314 * t650) * mrSges(10,1)) * pkin(2) + (t759 / 0.2e1 + t761 / 0.2e1 + t771 + t762 - t764 + (t731 / 0.2e1 - t738 / 0.2e1) * qJD(4) + t7 * t577) * t56 + (t758 / 0.2e1 + t760 / 0.2e1 + t770 + (t729 / 0.2e1 - t736 / 0.2e1) * qJD(4) + t13 * t577 + t476) * t51 + 0.2e1 * t643 * t835 * (t13 * t46 + t53 * t7); 0.2e1 * t54 * t4 + 0.2e1 * t76 * Ifges(11,3) * t26 - 0.2e1 * t12 * t40 + (t13 * t53 * t552 + 0.4e1 * t12 * t54) * t835 + (0.2e1 * t748 + (-mrSges(10,1) * t314 + mrSges(10,2) * t313) * t106 * t635) * t182 + (t527 * t53 + t532 * t634 + t730 + t735) * t18 + ((0.2e1 * Ifges(9,3) + 0.4e1 * t510) * t106 + ((t101 * t313 + t314 * t707) * mrSges(10,2) + (-t101 * t314 + t313 * t707) * mrSges(10,1)) * t635) * t188 + (0.2e1 * t770 + t760 + t758 + (t729 - t736) * qJD(4) + t13 * t527 + 0.2e1 * t476) * t56; -Ifges(6,6) * t739 + t531 * t24 + ((-mrSges(6,2) * t66 - Ifges(6,6) * t88) * t390 + (-mrSges(6,1) * t66 - Ifges(6,5) * t88) * t386) * qJD(4) + t752; -t530 * t7 + (-t46 * t751 + (mrSges(6,2) * t46 - Ifges(6,6) * t51) * t386) * qJD(4) + t625; -t530 * t13 + (-t53 * t751 + (mrSges(6,2) * t53 - Ifges(6,6) * t56) * t386) * qJD(4) + t622; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
