% Calculate joint inertia matrix for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1TE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_inertiaJ_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1TE_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1TE_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:31:54
% EndTime: 2020-04-12 20:50:35
% DurationCPUTime: 348.94s
% Computational Cost: add. (8706343->865), mult. (13125139->1420), div. (596496->33), fcn. (8315136->28), ass. (0->555)
t455 = pkin(6) ^ 2;
t437 = sin(pkin(20));
t441 = cos(pkin(20));
t443 = sin(qJ(3));
t447 = cos(qJ(3));
t405 = t437 * t447 + t441 * t443;
t628 = pkin(6) * t405;
t586 = 0.2e1 * pkin(1) * t628 + t455;
t682 = pkin(1) ^ 2;
t386 = t586 + t682;
t449 = pkin(13) ^ 2;
t459 = pkin(2) ^ 2;
t380 = t386 - t449 + t459;
t402 = -pkin(1) - t628;
t406 = t437 * t443 - t441 * t447;
t668 = -pkin(2) - pkin(13);
t374 = (pkin(1) - t668) * (pkin(1) + t668) + t586;
t663 = pkin(13) - pkin(2);
t375 = (pkin(1) - t663) * (pkin(1) + t663) + t586;
t594 = t375 * t374;
t462 = sqrt(-t594);
t592 = t406 * t462;
t575 = pkin(6) * t592;
t345 = -t380 * t402 - t575;
t627 = pkin(6) * t406;
t346 = t380 * t627 - t402 * t462;
t384 = 0.1e1 / t386;
t460 = 0.1e1 / pkin(2);
t593 = t384 * t460;
t640 = cos(qJ(2));
t444 = sin(qJ(2));
t648 = -t444 / 0.2e1;
t320 = (-t640 * t346 / 0.2e1 + t345 * t648) * t593;
t689 = (t640 * t345 / 0.2e1 + t346 * t648) * t593;
t271 = rSges(9,1) * t689 + rSges(9,2) * t320;
t623 = m(9) * t271;
t638 = sin(pkin(19));
t641 = cos(pkin(19));
t515 = t444 * t641 - t640 * t638;
t499 = pkin(7) * t515;
t477 = (-0.2e1 * t499 + pkin(1)) * pkin(1);
t681 = pkin(7) ^ 2;
t397 = t477 + t681;
t464 = pkin(3) ^ 2;
t583 = pkin(8) ^ 2 - t464;
t381 = t397 - t583;
t404 = -t499 + pkin(1);
t410 = t444 * t638 + t640 * t641;
t666 = pkin(7) + pkin(8);
t667 = pkin(7) - pkin(8);
t463 = sqrt(-((-pkin(3) + t666) * (pkin(3) + t667) + t477) * ((pkin(3) + t666) * (-pkin(3) + t667) + t477));
t497 = t515 * t463;
t626 = pkin(7) * t410;
t578 = pkin(1) * t626;
t595 = 0.4e1 / t463 * (t666 * t667 - t464 + t477) * t578;
t555 = -t595 / 0.2e1;
t679 = -0.2e1 * pkin(1);
t310 = (t497 + (t404 * t679 - t381 + t555) * t410) * pkin(7);
t579 = -0.2e1 * pkin(7) * t410 ^ 2;
t591 = t410 * t463;
t318 = t404 * t595 / 0.2e1 + (pkin(1) * t579 - t381 * t515 - t591) * pkin(7);
t347 = -pkin(7) * t591 + t404 * t381;
t348 = t381 * t626 + t404 * t463;
t458 = 0.1e1 / pkin(3);
t438 = cos(pkin(23));
t677 = 0.1e1 / t397;
t581 = -t677 / 0.2e1;
t553 = t438 * t581;
t434 = sin(pkin(23));
t580 = t677 / 0.2e1;
t554 = t434 * t580;
t325 = (t347 * t553 + t348 * t554) * t458;
t324 = 0.1e1 / t325 ^ 2;
t552 = t438 * t580;
t326 = (t347 * t554 + t348 * t552) * t458;
t395 = 0.1e1 / t397 ^ 2;
t546 = t395 * t578;
t518 = t348 * t546;
t519 = t347 * t546;
t695 = ((t310 * t554 + t318 * t552 + t434 * t519 + t438 * t518) / t325 - (t310 * t553 + t318 * t554 + t434 * t518 - t438 * t519) * t326 * t324) / (t324 * t326 ^ 2 + 0.1e1) * t458 + 0.1e1;
t473 = t397 + t583;
t482 = pkin(1) * t515 - pkin(7);
t574 = pkin(1) * t591;
t680 = 0.1e1 / pkin(8);
t470 = t680 * (-t473 * t482 - t574);
t468 = t470 * t580;
t472 = pkin(1) * t473;
t471 = t680 * (t410 * t472 - t463 * t482);
t469 = t677 * t471;
t639 = sin(pkin(18));
t642 = cos(pkin(18));
t330 = t642 * t468 - t639 * t469 / 0.2e1;
t331 = t642 * t469 / 0.2e1 + t639 * t468;
t513 = rSges(7,1) * t330 - rSges(7,2) * t331;
t694 = 2 * pkin(4);
t457 = pkin(4) ^ 2;
t456 = pkin(5) ^ 2;
t551 = t443 * t580;
t597 = t348 * t447;
t328 = (t347 * t551 + t580 * t597) * t458;
t550 = t447 * t581;
t329 = (t347 * t550 + t348 * t551) * t458;
t431 = pkin(23) + pkin(22);
t420 = sin(t431);
t421 = cos(t431);
t293 = t328 * t421 - t329 * t420;
t633 = t293 * pkin(5);
t588 = t633 * t694 + t456;
t274 = t457 + t588;
t584 = pkin(9) ^ 2 - pkin(11) ^ 2;
t266 = t274 - t584;
t279 = pkin(4) * t293 + pkin(5);
t665 = (-pkin(9) - pkin(11));
t263 = ((pkin(4) - t665) * (pkin(4) + t665)) + t588;
t664 = (pkin(11) - pkin(9));
t264 = ((pkin(4) - t664) * (pkin(4) + t664)) + t588;
t461 = sqrt(-t264 * t263);
t533 = t328 * t420 + t421 * t329;
t688 = t533 * t461;
t193 = -pkin(4) * t688 + t266 * t279;
t195 = pkin(4) * t266 * t533 + t279 * t461;
t440 = cos(pkin(21));
t272 = 0.1e1 / t274;
t452 = 0.1e1 / pkin(11);
t602 = t272 * t452;
t436 = sin(pkin(21));
t650 = t436 / 0.2e1;
t160 = (t193 * t650 + t195 * t440 / 0.2e1) * t602;
t161 = (-t193 * t440 / 0.2e1 + t195 * t650) * t602;
t408 = t443 * t640 + t444 * t447;
t409 = -t444 * t443 + t447 * t640;
t144 = -t160 * t408 + t161 * t409;
t145 = t160 * t409 + t161 * t408;
t114 = Icges(5,4) * t145 + Icges(5,2) * t144;
t693 = t114 / 0.2e1;
t115 = Icges(5,1) * t145 + Icges(5,4) * t144;
t692 = t115 / 0.2e1;
t691 = t144 / 0.2e1;
t690 = t145 / 0.2e1;
t445 = sin(qJ(1));
t647 = t445 / 0.2e1;
t448 = cos(qJ(1));
t646 = -t448 / 0.2e1;
t644 = t462 / 0.2e1;
t313 = t445 * t320;
t314 = t445 * t689;
t255 = Icges(9,5) * t313 - Icges(9,6) * t314 - Icges(9,3) * t448;
t315 = t448 * t320;
t316 = t448 * t689;
t256 = Icges(9,5) * t315 - Icges(9,6) * t316 + Icges(9,3) * t445;
t257 = Icges(9,4) * t313 - Icges(9,2) * t314 - Icges(9,6) * t448;
t258 = Icges(9,4) * t315 - Icges(9,2) * t316 + Icges(9,6) * t445;
t259 = Icges(9,1) * t313 - Icges(9,4) * t314 - Icges(9,5) * t448;
t260 = Icges(9,1) * t315 - Icges(9,4) * t316 + Icges(9,5) * t445;
t175 = -(-t255 * t448 - t257 * t314 + t259 * t313) * t448 + (-t256 * t448 - t258 * t314 + t260 * t313) * t445;
t344 = 0.1e1 / t345 ^ 2;
t385 = 0.1e1 / t386 ^ 2;
t596 = (t374 + t375) * pkin(1) * t627 / t644;
t556 = -t596 / 0.2e1;
t629 = pkin(6) * t346;
t637 = pkin(1) * t406;
t651 = t384 / 0.2e1;
t208 = 0.2e1 * (((t402 * t556 + (t380 * t405 - t592) * pkin(6)) * t651 + (-t384 * t406 * t455 + t385 * t629) * t637) / t345 - ((-t405 * t462 + (t556 - t380) * t406) * t651 + (t345 * t385 + t384 * t402) * t637) * t344 * t629) * pkin(2) / (t344 * t346 ^ 2 + 0.1e1) * t386 * t460;
t152 = t208 * t175;
t432 = t445 ^ 2;
t433 = t448 ^ 2;
t585 = t432 + t433;
t400 = t409 * t448;
t401 = t408 * t448;
t142 = -t160 * t400 - t161 * t401;
t143 = -t160 * t401 + t161 * t400;
t106 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t445;
t108 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t445;
t113 = Icges(5,5) * t145 + Icges(5,6) * t144;
t687 = t106 * t691 + t108 * t690 + t113 * t647 + t142 * t693 + t143 * t692;
t398 = t409 * t445;
t399 = t408 * t445;
t140 = -t160 * t398 - t161 * t399;
t141 = -t160 * t399 + t161 * t398;
t105 = Icges(5,4) * t141 + Icges(5,2) * t140 - Icges(5,6) * t448;
t107 = Icges(5,1) * t141 + Icges(5,4) * t140 - Icges(5,5) * t448;
t686 = t105 * t691 + t107 * t690 + t113 * t646 + t140 * t693 + t141 * t692;
t494 = -Icges(3,5) * t444 - Icges(3,6) * t640;
t508 = Icges(7,5) * t330 - Icges(7,6) * t331;
t309 = (t497 + (t555 - 0.3e1 * t681 + (0.4e1 * t499 - pkin(1)) * pkin(1) - t583) * t410) * pkin(1);
t317 = -t472 * t515 + t482 * t555 + t579 * t682 - t574;
t327 = 0.1e1 / t330 ^ 2;
t466 = t470 * t578;
t467 = t471 * t578;
t512 = t680 * t642 * t580;
t530 = t680 * t677 * t639;
t184 = ((t317 * t512 + t309 * t530 / 0.2e1 + (t466 * t639 + t467 * t642) * t395) / t330 - (t309 * t512 - t317 * t530 / 0.2e1 + (t466 * t642 - t467 * t639) * t395) * t331 * t327) / (t327 * t331 ^ 2 + 0.1e1);
t678 = t184 ^ 2;
t685 = -(-Icges(7,3) * t448 + t445 * t508) * t678 + Icges(3,3) * t448 - t445 * t494;
t684 = (Icges(7,3) * t445 + t448 * t508) * t678 + Icges(3,3) * t445 + t448 * t494;
t179 = t695 * t445;
t180 = t695 * t448;
t298 = -t444 * t325 - t326 * t640;
t289 = t298 * t445;
t297 = t325 * t640 - t444 * t326;
t290 = t297 * t445;
t291 = t298 * t448;
t292 = t297 * t448;
t683 = t179 * (Icges(8,5) * t291 - Icges(8,6) * t292 + Icges(8,3) * t445) - t180 * (Icges(8,5) * t289 - Icges(8,6) * t290 - Icges(8,3) * t448);
t243 = (t443 * t310 * t581 + t318 * t550 + (-t347 * t443 - t597) * t546) * t458;
t244 = (t310 * t550 + t318 * t551 + (-t347 * t447 + t348 * t443) * t546) * t458;
t205 = t243 * t421 + t244 * t420;
t502 = pkin(5) * (t263 + t264) * t694;
t173 = t205 * t502;
t206 = -t243 * t420 + t244 * t421;
t229 = 0.1e1 / t461;
t654 = -t229 / 0.2e1;
t560 = t533 * t654;
t487 = t173 * t560 - t206 * t461;
t547 = -0.2e1 * pkin(5) * t279 - t266;
t128 = (t205 * t547 + t487) * pkin(4);
t632 = pkin(5) * t533;
t549 = -0.2e1 * t457 * t632;
t561 = t229 * t279 / 0.2e1;
t604 = t205 * t461;
t129 = t173 * t561 + t205 * t549 + (t206 * t266 - t604) * pkin(4);
t156 = 0.1e1 / t161;
t273 = 0.1e1 / t274 ^ 2;
t634 = pkin(5) * t273;
t577 = pkin(4) * t634;
t544 = t440 * t577;
t522 = t205 * t544;
t545 = t436 * t577;
t523 = t205 * t545;
t603 = t272 * t440;
t557 = t603 / 0.2e1;
t558 = -t603 / 0.2e1;
t559 = t272 * t650;
t157 = 0.1e1 / t161 ^ 2;
t605 = t157 * t160;
t606 = 0.1e1 / (t157 * t160 ^ 2 + 0.1e1) * t452;
t662 = ((t128 * t559 + t129 * t557 + t193 * t523 + t195 * t522) * t156 - (t128 * t558 + t129 * t559 - t193 * t522 + t195 * t523) * t605) * t606 + 0.1e1;
t36 = t662 * t445;
t676 = t36 / 0.2e1;
t37 = t662 * t448;
t675 = -t37 / 0.2e1;
t265 = t274 + t584;
t278 = -pkin(4) - t633;
t194 = t265 * t632 - t278 * t461;
t478 = pkin(4) * (-t272 * t456 * t533 + t194 * t634);
t192 = -pkin(5) * t688 - t265 * t278;
t191 = 0.1e1 / t192 ^ 2;
t454 = 0.1e1 / pkin(9);
t531 = pkin(9) / (t191 * t194 ^ 2 + 0.1e1) * t274 * t454;
t488 = pkin(5) * t191 * t194 * t531;
t490 = pkin(4) * (t192 * t273 + t272 * t278);
t516 = 0.1e1 / t192 * t531;
t562 = t278 * t654;
t653 = t272 / 0.2e1;
t69 = 0.2e1 * ((t173 * t562 + (t206 * t265 - t604) * pkin(5)) * t653 + t205 * t478) * t516 - 0.2e1 * ((-t205 * t265 + t487) * t653 + t205 * t490) * t488;
t67 = t445 * t69 + t179;
t670 = t67 / 0.2e1;
t68 = (-t69 - t695) * t448;
t669 = t68 / 0.2e1;
t202 = t533 * t502;
t486 = t202 * t560 - t293 * t461;
t154 = (t533 * t547 + t486) * pkin(4);
t155 = t202 * t561 + t533 * t549 + (t266 * t293 - t688) * pkin(4);
t520 = t533 * t544;
t521 = t533 * t545;
t661 = ((t154 * t559 + t155 * t557 + t193 * t521 + t195 * t520) * t156 - (t154 * t558 + t155 * t559 - t193 * t520 + t195 * t521) * t605) * t606 + 0.1e1;
t660 = -t140 / 0.2e1;
t659 = -t142 / 0.2e1;
t658 = -t144 / 0.2e1;
t655 = t208 / 0.2e1;
t379 = t449 - t455 + t459 + (-pkin(1) - 0.2e1 * t628) * pkin(1);
t652 = -t379 / 0.2e1;
t439 = cos(pkin(22));
t649 = -t439 / 0.2e1;
t645 = -t462 / 0.2e1;
t636 = pkin(1) * t444;
t635 = pkin(2) * t689;
t631 = pkin(5) * t408;
t630 = pkin(5) * t445;
t625 = pkin(10) * t141;
t373 = rSges(4,1) * t408 + rSges(4,2) * t409;
t624 = m(4) * t373;
t442 = sin(qJ(4));
t446 = cos(qJ(4));
t132 = -t141 * t442 - t446 * t448;
t133 = t141 * t446 - t442 * t448;
t74 = Icges(6,5) * t133 + Icges(6,6) * t132 - Icges(6,3) * t140;
t76 = Icges(6,4) * t133 + Icges(6,2) * t132 - Icges(6,6) * t140;
t78 = Icges(6,1) * t133 + Icges(6,4) * t132 - Icges(6,5) * t140;
t29 = -t144 * t74 + (-t442 * t76 + t446 * t78) * t145;
t622 = t29 * t37;
t134 = -t143 * t442 + t445 * t446;
t135 = t143 * t446 + t442 * t445;
t75 = Icges(6,5) * t135 + Icges(6,6) * t134 - Icges(6,3) * t142;
t77 = Icges(6,4) * t135 + Icges(6,2) * t134 - Icges(6,6) * t142;
t79 = Icges(6,1) * t135 + Icges(6,4) * t134 - Icges(6,5) * t142;
t30 = -t144 * t75 + (-t442 * t77 + t446 * t79) * t145;
t621 = t30 * t36;
t396 = t400 * pkin(5);
t82 = -Icges(6,3) * t144 + (Icges(6,5) * t446 - Icges(6,6) * t442) * t145;
t84 = -Icges(6,5) * t144 + (Icges(6,1) * t446 - Icges(6,4) * t442) * t145;
t620 = t145 * t446 * t84 - t144 * t82;
t435 = sin(pkin(22));
t601 = t272 * t454;
t158 = (t435 * t192 / 0.2e1 + t194 * t649) * t601;
t159 = (t192 * t649 - t435 * t194 / 0.2e1) * t601;
t126 = t158 * t298 + t159 * t297;
t127 = -t158 * t297 + t159 * t298;
t100 = rSges(11,1) * t126 + rSges(11,2) * t127;
t619 = m(11) * t100;
t83 = -Icges(6,6) * t144 + (Icges(6,4) * t446 - Icges(6,2) * t442) * t145;
t616 = t442 * t83;
t615 = t445 * pkin(16);
t614 = t448 * rSges(3,3);
t613 = t448 * rSges(7,3);
t514 = -rSges(6,1) * t133 - rSges(6,2) * t132;
t80 = -rSges(6,3) * t140 - t514;
t612 = -pkin(12) * t140 + t625 + t80;
t81 = t135 * rSges(6,1) + t134 * rSges(6,2) - t142 * rSges(6,3);
t611 = -t143 * pkin(10) + pkin(12) * t142 - t81;
t85 = -rSges(6,3) * t144 + (rSges(6,1) * t446 - rSges(6,2) * t442) * t145;
t610 = pkin(10) * t145 - pkin(12) * t144 + t85;
t609 = Icges(3,4) * t444;
t608 = Icges(7,4) * t330;
t607 = Icges(7,4) * t331;
t600 = t292 * t435;
t590 = 0.1e1 / pkin(13) * t460;
t378 = 0.1e1 / t379 ^ 2;
t589 = t208 + (0.1e1 / t379 * t596 / 0.2e1 + t378 * t575 * t679) / (-t378 * t594 + 0.1e1);
t481 = rSges(4,1) * t398 - rSges(4,2) * t399 - rSges(4,3) * t448;
t566 = t400 * rSges(4,1) - t401 * rSges(4,2) + t445 * rSges(4,3);
t349 = t445 * t481 + t448 * t566;
t587 = t448 * t396 + t398 * t630;
t576 = t208 * t635;
t573 = t408 * t630;
t572 = t448 * t631;
t571 = t640 * pkin(1);
t32 = t134 * t83 + t135 * t84 - t142 * t82;
t570 = -t30 / 0.2e1 - t32 / 0.2e1;
t31 = t132 * t83 + t133 * t84 - t140 * t82;
t569 = -t31 / 0.2e1 - t29 / 0.2e1;
t124 = -t158 * t292 + t159 * t291;
t125 = -t158 * t291 - t159 * t292;
t95 = t124 * rSges(11,1) + t125 * rSges(11,2) + t445 * rSges(11,3);
t110 = t143 * rSges(5,1) + t142 * rSges(5,2) + t445 * rSges(5,3);
t249 = (t315 * t652 - t316 * t645) * t590;
t250 = (t315 * t644 - t316 * t652) * t590;
t222 = t249 * rSges(10,1) + t250 * rSges(10,2) + t445 * rSges(10,3);
t568 = t291 * rSges(8,1) - t292 * rSges(8,2) + t445 * rSges(8,3);
t262 = t315 * rSges(9,1) - t316 * rSges(9,2) + t445 * rSges(9,3);
t567 = t445 * rSges(7,3) + t513 * t448;
t565 = Icges(3,4) * t640;
t251 = (t320 * t645 + t652 * t689) * t590;
t252 = (t320 * t652 + t644 * t689) * t590;
t226 = rSges(10,1) * t251 + rSges(10,2) * t252;
t548 = -t226 - t635;
t535 = t445 * t636 - t615;
t529 = t445 * t571;
t528 = t448 * t571;
t363 = Icges(4,4) * t398 - Icges(4,2) * t399 - Icges(4,6) * t448;
t364 = Icges(4,4) * t400 - Icges(4,2) * t401 + Icges(4,6) * t445;
t365 = Icges(4,1) * t398 - Icges(4,4) * t399 - Icges(4,5) * t448;
t366 = Icges(4,1) * t400 - Icges(4,4) * t401 + Icges(4,5) * t445;
t370 = Icges(4,5) * t408 + Icges(4,6) * t409;
t371 = Icges(4,4) * t408 + Icges(4,2) * t409;
t372 = Icges(4,1) * t408 + Icges(4,4) * t409;
t527 = (t364 * t409 + t366 * t408 + t370 * t445 - t371 * t401 + t372 * t400) * t647 + (t363 * t409 + t365 * t408 - t370 * t448 - t371 * t399 + t372 * t398) * t646;
t526 = 0.2e1 * t585;
t525 = t585 * t636;
t430 = t448 * pkin(16);
t524 = -t448 * t636 + t430;
t517 = -t571 - t373;
t267 = t271 ^ 2;
t511 = t267 * t526;
t510 = Icges(7,1) * t330 - t607;
t509 = -Icges(7,2) * t331 + t608;
t176 = -(t255 * t445 - t257 * t316 + t259 * t315) * t448 + (t256 * t445 - t258 * t316 + t260 * t315) * t445;
t261 = rSges(9,1) * t313 - rSges(9,2) * t314 - rSges(9,3) * t448;
t227 = t261 * t445 + t262 * t448;
t505 = t289 * t439 + t290 * t435;
t301 = Icges(7,2) * t330 + t607;
t302 = Icges(7,1) * t331 + t608;
t504 = -t301 * t331 + t302 * t330;
t503 = t396 + t524;
t501 = -t571 - t631;
t500 = -pkin(5) * t398 + t535;
t498 = -rSges(3,1) * t444 - rSges(3,2) * t640;
t496 = -Icges(3,1) * t444 - t565;
t495 = -Icges(3,2) * t640 - t609;
t413 = -Icges(3,2) * t444 + t565;
t414 = Icges(3,1) * t640 - t609;
t491 = -t413 * t640 - t414 * t444;
t489 = (t313 * t445 + t315 * t448) * pkin(2);
t485 = -t525 + t587;
t484 = t501 * t445;
t483 = t501 * t448;
t109 = rSges(5,1) * t141 + rSges(5,2) * t140 - rSges(5,3) * t448;
t480 = rSges(8,1) * t289 - rSges(8,2) * t290 - rSges(8,3) * t448;
t247 = (t313 * t652 - t314 * t645) * t590;
t248 = (t313 * t644 - t314 * t652) * t590;
t221 = rSges(10,1) * t247 + rSges(10,2) * t248 - rSges(10,3) * t448;
t122 = -t158 * t290 + t159 * t289;
t123 = -t158 * t289 - t159 * t290;
t94 = rSges(11,1) * t122 + rSges(11,2) * t123 - rSges(11,3) * t448;
t475 = t445 * rSges(3,3) + t448 * t498;
t417 = rSges(2,1) * t448 - rSges(2,2) * t445;
t416 = rSges(3,1) * t640 - t444 * rSges(3,2);
t415 = -rSges(2,1) * t445 - rSges(2,2) * t448;
t412 = Icges(3,5) * t640 - Icges(3,6) * t444;
t377 = t430 + t475;
t376 = t614 + (-pkin(16) - t498) * t445;
t368 = t517 * t448;
t367 = t517 * t445;
t362 = Icges(4,5) * t400 - Icges(4,6) * t401 + Icges(4,3) * t445;
t361 = Icges(4,5) * t398 - Icges(4,6) * t399 - Icges(4,3) * t448;
t360 = t448 * t475 + (t445 * t498 - t614) * t445;
t354 = t524 + t566;
t353 = -t481 + t535;
t343 = -t525 + t349;
t308 = (-t362 * t448 - t364 * t399 + t366 * t398) * t445 - (-t361 * t448 - t363 * t399 + t365 * t398) * t448;
t307 = t445 * ((t362 * t445 - t364 * t401 + t366 * t400) * t445 - (t361 * t445 - t363 * t401 + t365 * t400) * t448);
t303 = rSges(7,1) * t331 + rSges(7,2) * t330;
t300 = Icges(7,5) * t331 + Icges(7,6) * t330;
t277 = -pkin(15) * t448 + t567;
t276 = t613 + (pkin(15) - t513) * t445;
t275 = t291 * t439 * pkin(4);
t270 = Icges(9,1) * t689 + Icges(9,4) * t320;
t269 = Icges(9,4) * t689 + Icges(9,2) * t320;
t268 = Icges(9,5) * t689 + Icges(9,6) * t320;
t254 = t430 + t262;
t253 = -t261 - t615;
t242 = rSges(8,1) * t297 + rSges(8,2) * t298;
t241 = Icges(8,1) * t297 + Icges(8,4) * t298;
t240 = Icges(8,4) * t297 + Icges(8,2) * t298;
t239 = Icges(8,5) * t297 + Icges(8,6) * t298;
t238 = (t297 * t439 - t298 * t435) * pkin(4);
t237 = Icges(8,1) * t291 - Icges(8,4) * t292 + Icges(8,5) * t445;
t236 = Icges(8,1) * t289 - Icges(8,4) * t290 - Icges(8,5) * t448;
t235 = Icges(8,4) * t291 - Icges(8,2) * t292 + Icges(8,6) * t445;
t234 = Icges(8,4) * t289 - Icges(8,2) * t290 - Icges(8,6) * t448;
t231 = t524 + t568;
t230 = -t480 + t535;
t225 = Icges(10,1) * t251 + Icges(10,4) * t252;
t224 = Icges(10,4) * t251 + Icges(10,2) * t252;
t223 = Icges(10,5) * t251 + Icges(10,6) * t252;
t220 = Icges(10,1) * t249 + Icges(10,4) * t250 + Icges(10,5) * t445;
t219 = Icges(10,1) * t247 + Icges(10,4) * t248 - Icges(10,5) * t448;
t218 = Icges(10,4) * t249 + Icges(10,2) * t250 + Icges(10,6) * t445;
t217 = Icges(10,4) * t247 + Icges(10,2) * t248 - Icges(10,6) * t448;
t216 = Icges(10,5) * t249 + Icges(10,6) * t250 + Icges(10,3) * t445;
t215 = Icges(10,5) * t247 + Icges(10,6) * t248 - Icges(10,3) * t448;
t214 = t548 * t448;
t213 = t548 * t445;
t212 = pkin(2) * t315 + t222 + t430;
t211 = -pkin(2) * t313 - t221 - t615;
t210 = t258 * t320 + t260 * t689;
t209 = t257 * t320 + t259 * t689;
t204 = t589 * t448;
t203 = t589 * t445;
t201 = t268 * t445 - t269 * t316 + t270 * t315;
t200 = -t268 * t448 - t269 * t314 + t270 * t313;
t189 = t221 * t445 + t222 * t448 + t489;
t178 = t218 * t252 + t220 * t251;
t177 = t217 * t252 + t219 * t251;
t174 = t227 * t208;
t172 = -t180 * t242 - t528;
t171 = -t179 * t242 - t529;
t170 = t223 * t445 + t224 * t250 + t225 * t249;
t169 = -t223 * t448 + t224 * t248 + t225 * t247;
t168 = t216 * t445 + t218 * t250 + t220 * t249;
t167 = t215 * t445 + t217 * t250 + t219 * t249;
t166 = -t216 * t448 + t218 * t248 + t220 * t247;
t165 = -t215 * t448 + t217 * t248 + t219 * t247;
t164 = -t204 * t226 - t448 * t576;
t163 = -t203 * t226 - t445 * t576;
t153 = t208 * t176;
t151 = (t448 * t567 + (t445 * t513 - t613) * t445) * t184;
t150 = t179 * t480 + t180 * t568 - t525;
t149 = -t167 * t448 + t168 * t445;
t148 = -t165 * t448 + t166 * t445;
t147 = t203 * t221 + t204 * t222 + t208 * t489;
t119 = -t167 * t204 + t168 * t203;
t118 = -t165 * t204 + t166 * t203;
t116 = rSges(5,1) * t145 + rSges(5,2) * t144;
t104 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t445;
t103 = Icges(5,5) * t141 + Icges(5,6) * t140 - Icges(5,3) * t448;
t102 = t503 + t110;
t101 = -t109 + t500;
t99 = Icges(11,1) * t126 + Icges(11,4) * t127;
t98 = Icges(11,4) * t126 + Icges(11,2) * t127;
t97 = Icges(11,5) * t126 + Icges(11,6) * t127;
t93 = Icges(11,1) * t124 + Icges(11,4) * t125 + Icges(11,5) * t445;
t92 = Icges(11,1) * t122 + Icges(11,4) * t123 - Icges(11,5) * t448;
t91 = Icges(11,4) * t124 + Icges(11,2) * t125 + Icges(11,6) * t445;
t90 = Icges(11,4) * t122 + Icges(11,2) * t123 - Icges(11,6) * t448;
t89 = Icges(11,5) * t124 + Icges(11,6) * t125 + Icges(11,3) * t445;
t88 = Icges(11,5) * t122 + Icges(11,6) * t123 - Icges(11,3) * t448;
t87 = pkin(4) * t600 + t275 + t524 + t95;
t86 = -pkin(4) * t505 + t535 - t94;
t73 = 0.2e1 * ((t202 * t562 + (t293 * t265 - t688) * pkin(5)) * t653 + t533 * t478) * t516 - 0.2e1 * ((-t265 * t533 + t486) * t653 + t533 * t490) * t488;
t66 = t503 - t611;
t65 = -t625 + (pkin(12) + rSges(6,3)) * t140 + t500 + t514;
t62 = t126 * t93 + t127 * t91;
t61 = t126 * t92 + t127 * t90;
t59 = t661 * t448;
t58 = t661 * t445;
t57 = t68 * t100 - t180 * t238 - t528;
t56 = -t67 * t100 - t179 * t238 - t529;
t53 = t142 * t85 - t144 * t81;
t52 = -t140 * t85 + t144 * t80;
t51 = t104 * t445 + t106 * t142 + t108 * t143;
t50 = t103 * t445 + t105 * t142 + t107 * t143;
t49 = -t104 * t448 + t106 * t140 + t108 * t141;
t48 = -t103 * t448 + t105 * t140 + t107 * t141;
t47 = t140 * t81 - t142 * t80;
t46 = t124 * t99 + t125 * t98 + t445 * t97;
t45 = t122 * t99 + t123 * t98 - t448 * t97;
t44 = t124 * t93 + t125 * t91 + t445 * t89;
t43 = t124 * t92 + t125 * t90 + t445 * t88;
t42 = t122 * t93 + t123 * t91 - t448 * t89;
t41 = t122 * t92 + t123 * t90 - t448 * t88;
t40 = -t116 * t59 - t572;
t39 = -t116 * t58 - t573;
t35 = (t445 * t94 + t448 * t95) * t73;
t34 = -t37 * t116 + t483;
t33 = -t36 * t116 + t484;
t28 = t134 * t77 + t135 * t79 - t142 * t75;
t27 = t134 * t76 + t135 * t78 - t142 * t74;
t26 = t132 * t77 + t133 * t79 - t140 * t75;
t25 = t132 * t76 + t133 * t78 - t140 * t74;
t24 = (-t145 * t616 + t620) * t144;
t23 = t180 * t275 + t67 * t94 - t68 * t95 + (t179 * t505 + t180 * t600) * pkin(4) - t525;
t22 = t109 * t58 + t110 * t59 + t587;
t21 = -t59 * t610 - t572;
t20 = -t58 * t610 - t573;
t19 = t109 * t36 + t110 * t37 + t485;
t18 = -t37 * t610 + t483;
t17 = -t36 * t610 + t484;
t16 = (-t43 * t448 + t44 * t445) * t73;
t15 = (-t41 * t448 + t42 * t445) * t73;
t14 = t43 * t68 + t44 * t67;
t13 = t41 * t68 + t42 * t67;
t12 = -t50 * t59 + t51 * t58;
t11 = -t48 * t59 + t49 * t58;
t10 = t36 * t51 - t37 * t50;
t9 = t36 * t49 - t37 * t48;
t8 = -t27 * t59 + t28 * t58;
t7 = -t25 * t59 + t26 * t58;
t6 = t58 * t612 - t59 * t611 + t587;
t5 = -t140 * t27 - t142 * t28 - t144 * t32;
t4 = -t140 * t25 - t142 * t26 - t144 * t31;
t3 = -t27 * t37 + t28 * t36;
t2 = -t25 * t37 + t26 * t36;
t1 = t36 * t612 - t37 * t611 + t485;
t38 = [(t115 - t616) * t145 + (t86 ^ 2 + t87 ^ 2) * m(11) + (t211 ^ 2 + t212 ^ 2) * m(10) + (t253 ^ 2 + t254 ^ 2) * m(9) + (t230 ^ 2 + t231 ^ 2) * m(8) + (t276 ^ 2 + t277 ^ 2) * m(7) + (t65 ^ 2 + t66 ^ 2) * m(6) + (t101 ^ 2 + t102 ^ 2) * m(5) + (t353 ^ 2 + t354 ^ 2) * m(4) + (t376 ^ 2 + t377 ^ 2) * m(3) + m(2) * (t415 ^ 2 + t417 ^ 2) + t330 * t301 + t331 * t302 + t620 + t640 * t414 + t127 * t98 + t298 * t240 + t297 * t241 + t126 * t99 + t144 * t114 + t252 * t224 + t251 * t225 - t444 * t413 + t408 * t372 + t409 * t371 + t689 * t270 + t320 * t269 + Icges(2,3); (-t253 * t448 - t254 * t445) * t623 + (t56 * t87 + t57 * t86) * m(11) + (t17 * t66 + t18 * t65) * m(6) + t527 + (t101 * t34 + t102 * t33) * m(5) + (t353 * t368 + t354 * t367) * m(4) + (t211 * t214 + t212 * t213) * m(10) + (t171 * t231 + t172 * t230) * m(8) + t621 / 0.2e1 - t622 / 0.2e1 + t31 * t675 + t32 * t676 + (-t376 * t448 - t377 * t445) * t416 * m(3) + (t46 + t62) * t670 + (t45 + t61) * t669 + (t235 * t298 + t237 * t297 + t239 * t445 - t240 * t292 + t241 * t291) * t179 / 0.2e1 - (t234 * t298 + t236 * t297 - t239 * t448 - t240 * t290 + t241 * t289) * t180 / 0.2e1 - t686 * t37 + t687 * t36 + (-t444 * (Icges(3,6) * t445 + t448 * t495) + (Icges(3,5) * t445 + t448 * t496) * t640 + t445 * t412 + t448 * t491 + t170 + t178 + t201 + t210) * t647 + (-t444 * (-Icges(3,6) * t448 + t445 * t495) + (-Icges(3,5) * t448 + t445 * t496) * t640 - t448 * t412 + t445 * t491 + t169 + t177 + t200 + t209) * t646 + ((-t276 * t448 - t277 * t445) * t303 * m(7) + ((Icges(7,6) * t445 + t448 * t509) * t330 + (Icges(7,5) * t445 + t448 * t510) * t331 + t300 * t445 + t448 * t504) * t647 + ((-Icges(7,6) * t448 + t445 * t509) * t330 + (-Icges(7,5) * t448 + t445 * t510) * t331 - t300 * t448 + t445 * t504) * t646) * t184; (t23 ^ 2 + t56 ^ 2 + t57 ^ 2) * m(11) + t307 + m(9) * (0.2e1 * t227 ^ 2 + t511) / 0.2e1 + (t1 ^ 2 + t17 ^ 2 + t18 ^ 2) * m(6) + (t19 ^ 2 + t33 ^ 2 + t34 ^ 2) * m(5) + (t343 ^ 2 + t367 ^ 2 + t368 ^ 2) * m(4) + (t189 ^ 2 + t213 ^ 2 + t214 ^ 2) * m(10) + (t150 ^ 2 + t171 ^ 2 + t172 ^ 2) * m(8) + m(7) * t151 ^ 2 + t67 * t14 + t68 * t13 - t180 * ((-t235 * t290 + t237 * t289) * t179 - (-t234 * t290 + t236 * t289) * t180) + t179 * ((-t235 * t292 + t237 * t291) * t179 - (-t234 * t292 + t236 * t291) * t180) + m(3) * t360 ^ 2 + (m(7) * t303 ^ 2 * t678 / 0.2e1 + m(3) * t416 ^ 2 / 0.2e1) * t526 - (t9 + t2) * t37 + (t3 + t10) * t36 + (t179 * t683 + t684 * t432 + t149 + t176) * t445 + (t180 * t683 - t175 - t148 - t308 + t685 * t433 + (t685 * t445 + t684 * t448) * t445) * t448; (t101 * t40 + t102 * t39) * m(5) + (t20 * t66 + t21 * t65) * m(6) + (t163 * t212 + t164 * t211) * m(10) - (t177 / 0.2e1 + t169 / 0.2e1) * t204 + (t170 / 0.2e1 + t178 / 0.2e1) * t203 - (-t569 + t686) * t59 + (-t570 + t687) * t58 + (-t353 * t624 + (-t61 / 0.2e1 - t45 / 0.2e1 - t86 * t619) * t73 + (-t200 / 0.2e1 - t209 / 0.2e1 - t253 * t623) * t208) * t448 + (-t354 * t624 + (t62 / 0.2e1 + t46 / 0.2e1 - t87 * t619) * t73 + (t201 / 0.2e1 + t210 / 0.2e1 - t254 * t623) * t208) * t445 + t527; t307 + t203 * t149 / 0.2e1 - t204 * t148 / 0.2e1 + t16 * t670 + (t19 * t22 + t33 * t39 + t34 * t40) * m(5) + t15 * t669 + m(11) * t23 * t35 + m(4) * t343 * t349 + (t147 * t189 + t163 * t213 + t164 * t214) * m(10) + (t1 * t6 + t17 * t20 + t18 * t21) * m(6) - (t9 / 0.2e1 + t2 / 0.2e1) * t59 + (t3 / 0.2e1 + t10 / 0.2e1) * t58 - (t11 / 0.2e1 + t7 / 0.2e1) * t37 + (t8 / 0.2e1 + t12 / 0.2e1) * t36 + (t176 * t655 - t367 * t624 + t119 / 0.2e1 + t153 / 0.2e1 + (t14 / 0.2e1 - t56 * t619) * t73) * t445 + (-t152 - t308 - t118 / 0.2e1 - t368 * t624 + (-t13 / 0.2e1 - t57 * t619) * t73) * t448 + (t174 * t227 + t511 * t655) * m(9); (t147 ^ 2 + t163 ^ 2 + t164 ^ 2) * m(10) + m(9) * t174 ^ 2 + m(4) * t349 ^ 2 + m(11) * t35 ^ 2 + (t22 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(5) - t204 * t118 + (t20 ^ 2 + t21 ^ 2 + t6 ^ 2) * m(6) + t203 * t119 + t307 - (t7 + t11) * t59 + (t8 + t12) * t58 + (t153 * t208 + t16 * t73) * t445 + (-t15 * t73 - t152 * t208 - t308) * t448 + t585 * (m(11) * t100 ^ 2 * t73 ^ 2 + m(9) * t208 ^ 2 * t267 + m(4) * t373 ^ 2); (t52 * t65 + t53 * t66) * m(6) - t24 + t570 * t142 + t569 * t140; t3 * t659 + t5 * t676 + t4 * t675 + (t1 * t47 + t17 * t53 + t18 * t52) * m(6) + (t621 - t622) * t658 + t2 * t660; (t20 * t53 + t21 * t52 + t47 * t6) * m(6) + t7 * t660 - t59 * t4 / 0.2e1 + t58 * t5 / 0.2e1 + (-t29 * t59 + t30 * t58) * t658 + t8 * t659; -t140 * t4 + (t47 ^ 2 + t52 ^ 2 + t53 ^ 2) * m(6) - t142 * t5 - t144 * (-t29 * t140 - t30 * t142 - t24);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t38(1), t38(2), t38(4), t38(7); t38(2), t38(3), t38(5), t38(8); t38(4), t38(5), t38(6), t38(9); t38(7), t38(8), t38(9), t38(10);];
Mq = res;