% Calculate joint inertia matrix for
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_inertiaJ_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_inertiaJ_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2IC_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2IC_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:46:50
% EndTime: 2020-05-02 23:47:00
% DurationCPUTime: 10.59s
% Computational Cost: add. (24877->665), mult. (22085->986), div. (144->14), fcn. (22205->48), ass. (0->393)
t522 = sin(qJ(1));
t527 = cos(qJ(1));
t695 = t522 * t527;
t519 = sin(qJ(5));
t524 = cos(qJ(5));
t467 = rSges(6,1) * t519 + rSges(6,2) * t524;
t468 = rSges(6,1) * t524 - rSges(6,2) * t519;
t517 = qJ(2) + qJ(3);
t506 = qJ(4) + t517;
t492 = cos(t506);
t649 = t492 * t527;
t490 = sin(t506);
t650 = t490 * t527;
t305 = rSges(6,3) * t650 + t522 * t467 + t468 * t649;
t697 = pkin(9) * t649 + pkin(11) * t650 + t305;
t521 = sin(qJ(2));
t526 = cos(qJ(2));
t455 = rSges(11,1) * t521 + rSges(11,2) * t526;
t456 = rSges(11,1) * t526 - rSges(11,2) * t521;
t511 = -qJ(7) + pkin(19);
t498 = -qJ(10) + t511;
t484 = sin(t498);
t485 = cos(t498);
t317 = t455 * t484 + t456 * t485;
t311 = t317 * t522;
t495 = sin(t511);
t497 = cos(t511);
t682 = pkin(4) * (t495 * t521 + t497 * t526);
t684 = pkin(1) * t526;
t607 = -t682 - t684;
t269 = t522 * t607 + t311;
t312 = t317 * t527;
t270 = t527 * t607 + t312;
t486 = -qJ(2) + t498;
t477 = sin(t486);
t478 = cos(t486);
t403 = -rSges(11,1) * t477 + rSges(11,2) * t478;
t340 = t403 - pkin(4) * sin(qJ(2) - t511);
t685 = pkin(1) * t521;
t322 = t340 - t685;
t554 = -Icges(11,5) * t477 + Icges(11,6) * t478;
t323 = -Icges(11,3) * t527 + t522 * t554;
t324 = Icges(11,3) * t522 + t527 * t554;
t513 = t522 ^ 2;
t633 = Icges(11,4) * t477;
t555 = Icges(11,2) * t478 - t633;
t326 = Icges(11,6) * t522 + t527 * t555;
t632 = Icges(11,4) * t478;
t556 = -Icges(11,1) * t477 + t632;
t328 = Icges(11,5) * t522 + t527 * t556;
t578 = t326 * t478 - t328 * t477;
t325 = -Icges(11,6) * t527 + t522 * t555;
t327 = -Icges(11,5) * t527 + t522 * t556;
t579 = -t325 * t478 + t327 * t477;
t210 = t522 * (t513 * t324 + (t579 * t527 + (-t323 + t578) * t522) * t527);
t514 = t527 ^ 2;
t212 = t514 * t323 + (t578 * t522 + (-t324 + t579) * t527) * t522;
t617 = -t527 * t212 + t210;
t187 = m(11) * (t322 * t403 + (t269 * t522 + t270 * t527) * t317) + t617;
t284 = -t522 * t682 + t311;
t285 = -t527 * t682 + t312;
t188 = m(11) * (t340 * t403 + (t284 * t522 + t285 * t527) * t317) + t617;
t533 = 0.1e1 / pkin(3);
t630 = qJ(7) + pkin(20);
t642 = qJ(6) - qJ(2);
t610 = t630 - t642;
t483 = cos(t610);
t475 = 0.1e1 / t483;
t653 = (-pkin(3) * t483 - pkin(1) * cos(t642)) * t475;
t625 = t533 * t653;
t696 = t188 * t625 + t187;
t501 = sin(t517);
t504 = cos(t517);
t585 = Icges(4,5) * t504 - Icges(4,6) * t501;
t370 = -Icges(4,3) * t527 + t522 * t585;
t371 = Icges(4,3) * t522 + t527 * t585;
t644 = t524 * t527;
t648 = t519 * t522;
t405 = -t492 * t648 - t644;
t645 = t522 * t524;
t647 = t519 * t527;
t406 = t492 * t645 - t647;
t651 = t490 * t522;
t292 = Icges(6,5) * t406 + Icges(6,6) * t405 + Icges(6,3) * t651;
t294 = Icges(6,4) * t406 + Icges(6,2) * t405 + Icges(6,6) * t651;
t296 = Icges(6,1) * t406 + Icges(6,4) * t405 + Icges(6,5) * t651;
t217 = t292 * t651 + t294 * t405 + t296 * t406;
t407 = -t492 * t647 + t645;
t408 = t492 * t644 + t648;
t293 = Icges(6,5) * t408 + Icges(6,6) * t407 + Icges(6,3) * t650;
t295 = Icges(6,4) * t408 + Icges(6,2) * t407 + Icges(6,6) * t650;
t297 = Icges(6,1) * t408 + Icges(6,4) * t407 + Icges(6,5) * t650;
t218 = t293 * t651 + t295 * t405 + t297 * t406;
t203 = -t217 * t527 + t218 * t522;
t584 = Icges(5,5) * t492 - Icges(5,6) * t490;
t347 = -Icges(5,3) * t527 + t522 * t584;
t348 = Icges(5,3) * t522 + t527 * t584;
t669 = Icges(5,4) * t492;
t591 = -Icges(5,2) * t490 + t669;
t352 = Icges(5,6) * t522 + t527 * t591;
t670 = Icges(5,4) * t490;
t598 = Icges(5,1) * t492 - t670;
t356 = Icges(5,5) * t522 + t527 * t598;
t574 = -t352 * t490 + t356 * t492;
t351 = -Icges(5,6) * t527 + t522 * t591;
t355 = -Icges(5,5) * t527 + t522 * t598;
t575 = t351 * t490 - t355 * t492;
t639 = -t514 * t347 - (t574 * t522 + (-t348 + t575) * t527) * t522 - t203;
t219 = t292 * t650 + t294 * t407 + t296 * t408;
t220 = t293 * t650 + t295 * t407 + t297 * t408;
t204 = -t219 * t527 + t220 * t522;
t640 = (t513 * t348 + t204 + (t575 * t527 + (-t347 + t574) * t522) * t527) * t522;
t537 = t522 * (-t370 * t695 + t513 * t371) + t640 + (-t514 * t370 + t371 * t695 + t639) * t527;
t515 = qJ(2) + qJ(8);
t505 = qJ(9) + t515;
t489 = sin(t505);
t491 = cos(t505);
t424 = rSges(10,1) * t489 + rSges(10,2) * t491;
t499 = sin(t515);
t693 = -pkin(2) * t499 + t424;
t631 = qJ(4) + pkin(18);
t608 = pkin(19) + qJ(3) + t631;
t548 = t608 + t642;
t539 = -(2 * qJ(7)) - pkin(20) + t548;
t549 = t608 - t642;
t542 = pkin(20) + t549;
t692 = cos(qJ(10) - t539) + cos(qJ(10) - t542);
t516 = qJ(2) + qJ(7);
t500 = sin(t516);
t503 = cos(t516);
t602 = rSges(8,1) * t500 + rSges(8,2) * t503;
t691 = -rSges(8,3) * t527 - t522 * t602;
t678 = rSges(4,2) * t501;
t679 = rSges(4,1) * t504;
t690 = t522 * (-t678 + t679) - rSges(4,3) * t527;
t687 = t522 / 0.2e1;
t686 = -t527 / 0.2e1;
t520 = sin(qJ(3));
t525 = cos(qJ(3));
t646 = t521 * t525;
t681 = pkin(5) * (t520 * t526 + t646);
t496 = sin(t631);
t680 = pkin(5) * t496;
t676 = rSges(5,3) * t527;
t674 = Icges(3,4) * t521;
t673 = Icges(3,4) * t526;
t672 = Icges(4,4) * t501;
t671 = Icges(4,4) * t504;
t518 = sin(qJ(6));
t668 = Icges(7,4) * t518;
t523 = cos(qJ(6));
t667 = Icges(7,4) * t523;
t666 = Icges(8,4) * t500;
t665 = Icges(8,4) * t503;
t664 = Icges(9,4) * t499;
t502 = cos(t515);
t663 = Icges(9,4) * t502;
t662 = Icges(10,4) * t489;
t661 = Icges(10,4) * t491;
t634 = t513 + t514;
t189 = m(11) * (t317 ^ 2 * t634 + t403 ^ 2) + t617;
t660 = t189 / pkin(8) ^ 2;
t231 = -t292 * t492 + (-t294 * t519 + t296 * t524) * t490;
t659 = t231 * t527;
t232 = -t293 * t492 + (-t295 * t519 + t297 * t524) * t490;
t658 = t232 * t522;
t330 = -Icges(6,6) * t492 + (Icges(6,4) * t524 - Icges(6,2) * t519) * t490;
t657 = t330 * t519;
t332 = 0.1e1 / t692;
t656 = t332 * t533;
t448 = cos(qJ(7) + qJ(10) - t608);
t530 = 0.1e1 / pkin(10);
t655 = (-pkin(10) * t448 - pkin(5) * cos(qJ(3) + t498)) * t530;
t622 = -qJ(8) + pkin(17) + qJ(3);
t482 = cos(qJ(9) - t622);
t654 = (-pkin(12) * t482 + pkin(2) * cos(t622)) / pkin(12);
t652 = t467 * t527;
t580 = Icges(10,5) * t489 + Icges(10,6) * t491;
t345 = -Icges(10,3) * t527 + t522 * t580;
t346 = Icges(10,3) * t522 + t527 * t580;
t587 = Icges(10,2) * t491 + t662;
t350 = Icges(10,6) * t522 + t527 * t587;
t594 = Icges(10,1) * t489 + t661;
t354 = Icges(10,5) * t522 + t527 * t594;
t576 = t350 * t491 + t354 * t489;
t349 = -Icges(10,6) * t527 + t522 * t587;
t353 = -Icges(10,5) * t527 + t522 * t594;
t577 = -t349 * t491 - t353 * t489;
t215 = t514 * t345 + (t576 * t522 + (-t346 + t577) * t527) * t522;
t643 = t527 * t215;
t422 = rSges(10,1) * t491 - t489 * rSges(10,2);
t619 = -pkin(2) * t502 + t422;
t336 = t619 * t522;
t337 = t619 * t527;
t213 = t522 * (t513 * t346 + (t577 * t527 + (-t345 + t576) * t522) * t527);
t616 = t213 - t643;
t197 = m(10) * (t693 * t424 + (t336 * t522 + t337 * t527) * t422) + t616;
t581 = -Icges(9,5) * t499 - Icges(9,6) * t502;
t366 = -Icges(9,3) * t527 + t522 * t581;
t367 = Icges(9,3) * t522 + t527 * t581;
t588 = -Icges(9,2) * t502 - t664;
t373 = Icges(9,6) * t522 + t527 * t588;
t595 = -Icges(9,1) * t499 - t663;
t379 = Icges(9,5) * t522 + t527 * t595;
t572 = -t373 * t502 - t379 * t499;
t372 = -Icges(9,6) * t527 + t522 * t588;
t378 = -Icges(9,5) * t527 + t522 * t595;
t573 = t372 * t502 + t378 * t499;
t226 = t514 * t366 + (t572 * t522 + (-t367 + t573) * t527) * t522;
t512 = 0.1e1 / sin(qJ(9));
t534 = 0.1e1 / pkin(2);
t628 = pkin(6) * t512 * t534;
t605 = t628 * t654;
t441 = rSges(9,1) * t502 - rSges(9,2) * t499;
t443 = -rSges(9,1) * t499 - rSges(9,2) * t502;
t613 = t213 + t522 * (t513 * t367 + (t573 * t527 + (-t366 + t572) * t522) * t527) + m(10) * (t336 ^ 2 + t337 ^ 2 + t693 ^ 2) + m(9) * (t441 ^ 2 * t634 + t443 ^ 2);
t615 = t482 * t628;
t641 = ((-t215 - t226) * t527 + t613) * t615 + t197 * t605;
t553 = rSges(5,1) * t649 - rSges(5,2) * t650 + t522 * rSges(5,3);
t603 = rSges(5,1) * t492 - rSges(5,2) * t490;
t286 = t522 * (t522 * t603 - t676) + t527 * t553;
t335 = -rSges(6,3) * t492 + t468 * t490;
t638 = -pkin(9) * t490 + pkin(11) * t492 - t335;
t507 = t522 * rSges(8,3);
t289 = t522 * t691 + t527 * (-t527 * t602 + t507);
t635 = t522 * rSges(4,3) + t527 * t679;
t290 = t522 * t690 + t527 * (-t527 * t678 + t635);
t493 = pkin(5) * t520 + pkin(1);
t413 = pkin(5) * t646 + t493 * t526;
t421 = rSges(5,1) * t490 + rSges(5,2) * t492;
t637 = -t413 - t421;
t531 = 0.1e1 / pkin(8);
t629 = t531 * t680;
t543 = -qJ(7) + t548;
t544 = -qJ(7) + t549;
t241 = (-(cos(t544) + cos(t543)) * pkin(1) + (-cos(t539) - cos(t542)) * pkin(3)) * pkin(4) + ((cos(qJ(10) - t544) + cos(qJ(10) - t543)) * pkin(1) + t692 * pkin(3)) * pkin(8);
t627 = t241 * t656;
t447 = 0.1e1 / t448;
t626 = t447 * t655;
t623 = -t413 + t638;
t487 = pkin(5) * t504;
t277 = t487 + t286;
t440 = rSges(4,1) * t501 + rSges(4,2) * t504;
t621 = -t440 - t684;
t442 = rSges(8,1) * t503 - rSges(8,2) * t500;
t620 = -t442 - t684;
t618 = -t421 - t681;
t614 = t531 * t627;
t388 = -Icges(11,5) * t478 - Icges(11,6) * t477;
t389 = -Icges(11,2) * t477 - t632;
t390 = -Icges(11,1) * t478 - t633;
t567 = t389 * t478 - t390 * t477;
t612 = (-t326 * t477 - t328 * t478 + t388 * t522 + t527 * t567) * t687 + (-t325 * t477 - t327 * t478 - t388 * t527 + t522 * t567) * t686;
t415 = -Icges(10,5) * t491 + Icges(10,6) * t489;
t417 = Icges(10,2) * t489 - t661;
t419 = -Icges(10,1) * t491 + t662;
t562 = t417 * t491 + t419 * t489;
t611 = (t350 * t489 - t354 * t491 + t415 * t522 + t527 * t562) * t687 + (t349 * t489 - t353 * t491 - t415 * t527 + t522 * t562) * t686;
t304 = -t652 + (rSges(6,3) * t490 + t468 * t492) * t522;
t238 = t522 * t304 + t513 * (pkin(9) * t492 + pkin(11) * t490) + t697 * t527;
t609 = t638 - t681;
t606 = pkin(4) * (pkin(1) * (sin(qJ(10) - t642) + sin(qJ(10) + t642)) + (sin(qJ(10) - t610) + sin(qJ(10) + t610)) * pkin(3)) * t530 * t656;
t464 = rSges(3,1) * t521 + rSges(3,2) * t526;
t466 = rSges(7,1) * t523 - rSges(7,2) * t518;
t329 = -Icges(6,3) * t492 + (Icges(6,5) * t524 - Icges(6,6) * t519) * t490;
t331 = -Icges(6,5) * t492 + (Icges(6,1) * t524 - Icges(6,4) * t519) * t490;
t239 = t329 * t651 + t330 * t405 + t331 * t406;
t195 = -t239 * t492 + (t217 * t522 + t218 * t527) * t490;
t240 = t329 * t650 + t330 * t407 + t331 * t408;
t196 = -t240 * t492 + (t219 * t522 + t220 * t527) * t490;
t601 = t195 * t686 + t196 * t687 + t203 * t651 / 0.2e1 + t204 * t650 / 0.2e1 - t492 * (t658 - t659) / 0.2e1;
t234 = t487 + t238;
t600 = -Icges(3,1) * t521 - t673;
t599 = Icges(4,1) * t504 - t672;
t597 = Icges(7,1) * t523 - t668;
t596 = -Icges(8,1) * t500 - t665;
t593 = -Icges(3,2) * t526 - t674;
t592 = -Icges(4,2) * t501 + t671;
t590 = -Icges(7,2) * t518 + t667;
t589 = -Icges(8,2) * t503 - t666;
t586 = -Icges(3,5) * t521 - Icges(3,6) * t526;
t583 = Icges(7,5) * t523 - Icges(7,6) * t518;
t582 = -Icges(8,5) * t500 - Icges(8,6) * t503;
t374 = -Icges(8,6) * t527 + t522 * t589;
t380 = -Icges(8,5) * t527 + t522 * t596;
t571 = t374 * t503 + t380 * t500;
t375 = Icges(8,6) * t522 + t527 * t589;
t381 = Icges(8,5) * t522 + t527 * t596;
t570 = -t375 * t503 - t381 * t500;
t418 = Icges(5,2) * t492 + t670;
t420 = Icges(5,1) * t490 + t669;
t561 = -t418 * t490 + t420 * t492;
t434 = -Icges(9,2) * t499 + t663;
t437 = Icges(9,1) * t502 - t664;
t560 = -t434 * t502 - t437 * t499;
t435 = -Icges(8,2) * t500 + t665;
t438 = Icges(8,1) * t503 - t666;
t559 = -t435 * t503 - t438 * t500;
t436 = Icges(4,2) * t504 + t672;
t439 = Icges(4,1) * t501 + t671;
t558 = -t436 * t501 + t439 * t504;
t460 = -Icges(3,2) * t521 + t673;
t462 = Icges(3,1) * t526 - t674;
t557 = -t460 * t526 - t462 * t521;
t552 = -pkin(14) + t466;
t551 = t527 * t639 + t640;
t368 = -Icges(8,3) * t527 + t522 * t582;
t369 = Icges(8,3) * t522 + t527 * t582;
t224 = t522 * (t513 * t369 + (t571 * t527 + (-t368 + t570) * t522) * t527);
t227 = t514 * t368 + (t570 * t522 + (-t369 + t571) * t527) * t522;
t550 = t210 + t224 + (-t227 - t212) * t527;
t546 = -t527 * t227 + t224 + t617;
t545 = -pkin(4) * (t495 * t526 - t497 * t521) - t455 * t485 + t456 * t484;
t416 = Icges(5,5) * t490 + Icges(5,6) * t492;
t541 = -t659 / 0.2e1 + t658 / 0.2e1 + (t352 * t492 + t356 * t490 + t416 * t522 + t527 * t561 + t240) * t687 + (t351 * t492 + t355 * t490 - t416 * t527 + t522 * t561 + t239) * t686;
t432 = Icges(8,5) * t503 - Icges(8,6) * t500;
t540 = t612 + (-t375 * t500 + t381 * t503 + t432 * t522 + t527 * t559) * t687 + (-t374 * t500 + t380 * t503 - t432 * t527 + t522 * t559) * t686;
t233 = t234 - t685;
t260 = t623 * t522;
t261 = t623 * t527;
t266 = t277 - t685;
t298 = t638 * t522;
t299 = t638 * t527;
t306 = t637 * t522;
t307 = t637 * t527;
t181 = m(6) * (t233 * t238 + t260 * t298 + t261 * t299) + m(5) * (t266 * t286 + (-t306 * t522 - t307 * t527) * t421) + t551;
t267 = t609 * t522;
t268 = t609 * t527;
t313 = t618 * t522;
t314 = t618 * t527;
t182 = m(6) * (t234 * t238 + t267 * t298 + t268 * t299) + m(5) * (t277 * t286 + (-t313 * t522 - t314 * t527) * t421) + t551;
t308 = rSges(10,3) * t527 + (-pkin(15) - t693) * t522;
t510 = t527 * pkin(15);
t309 = rSges(10,3) * t522 + t527 * t693 + t510;
t338 = rSges(9,3) * t527 + (-pkin(15) - t443) * t522;
t339 = rSges(9,3) * t522 + t443 * t527 + t510;
t431 = Icges(9,5) * t502 - Icges(9,6) * t499;
t538 = m(10) * (t308 * t337 + t309 * t336) + m(9) * (-t338 * t527 - t339 * t522) * t441 + t611 + (-t373 * t499 + t379 * t502 + t431 * t522 + t527 * t560) * t687 + (-t372 * t499 + t378 * t502 - t431 * t527 + t522 * t560) * t686;
t433 = Icges(4,5) * t501 + Icges(4,6) * t504;
t536 = t541 + ((Icges(4,6) * t522 + t527 * t592) * t504 + (Icges(4,5) * t522 + t527 * t599) * t501 + t433 * t522 + t527 * t558) * t687 + ((-Icges(4,6) * t527 + t522 * t592) * t504 + (-Icges(4,5) * t527 + t522 * t599) * t501 - t433 * t527 + t522 * t558) * t686;
t288 = t290 - t685;
t357 = t621 * t522;
t359 = t621 * t527;
t535 = m(6) * (t233 * t234 + t260 * t267 + t261 * t268) + m(5) * (t266 * t277 + t306 * t313 + t307 * t314) + m(4) * (t288 * t290 + (-t357 * t522 - t359 * t527) * t440) + t641 + t537;
t494 = sin(t630);
t488 = -pkin(15) + t685;
t473 = t488 * t522;
t471 = rSges(2,1) * t527 - rSges(2,2) * t522;
t470 = rSges(3,1) * t526 - rSges(3,2) * t521;
t465 = -rSges(2,1) * t522 - rSges(2,2) * t527;
t463 = rSges(7,1) * t518 + rSges(7,2) * t523;
t458 = Icges(3,5) * t526 - Icges(3,6) * t521;
t409 = pkin(5) * t525 * t526 - t493 * t521 + pkin(15);
t394 = Icges(3,3) * t522 + t527 * t586;
t393 = -Icges(3,3) * t527 + t522 * t586;
t392 = Icges(7,3) * t522 + t527 * t583;
t391 = -Icges(7,3) * t527 + t522 * t583;
t387 = t409 * t527;
t384 = rSges(3,3) * t527 + (-pkin(15) + t464) * t522;
t365 = rSges(7,3) * t522 + t527 * t552;
t364 = rSges(3,3) * t522 - t464 * t527 + t510;
t363 = rSges(7,3) * t527 - t522 * t552;
t360 = t620 * t527;
t358 = t620 * t522;
t321 = (-t488 - t678) * t527 + t635;
t320 = t507 + (-t488 - t602) * t527;
t319 = t473 - t690;
t318 = t473 - t691;
t310 = t490 * t524 * t331;
t301 = t387 + t553;
t300 = t676 + (-t409 - t603) * t522;
t287 = t289 - t685;
t259 = rSges(11,3) * t522 + (-t488 - t545) * t527;
t258 = rSges(11,3) * t527 + t522 * t545 + t473;
t251 = t387 + t697;
t250 = t652 + (-t409 + (-pkin(9) - t468) * t492 + (-pkin(11) - rSges(6,3)) * t490) * t522;
t248 = -t305 * t492 - t335 * t650;
t247 = t304 * t492 + t335 * t651;
t243 = (t304 * t527 - t305 * t522) * t490;
t242 = -t492 * t329 - t490 * t657 + t310;
t192 = m(11) * (t258 * t527 + t259 * t522) * t317 + t612;
t185 = -t242 * t492 + m(6) * (t247 * t250 + t248 * t251) + ((t232 / 0.2e1 + t240 / 0.2e1) * t527 + (t231 / 0.2e1 + t239 / 0.2e1) * t522) * t490;
t184 = m(6) * (t250 * t299 + t251 * t298) + m(5) * (-t300 * t527 - t301 * t522) * t421 + t541;
t183 = m(5) * (t421 ^ 2 * t634 + t286 ^ 2) + m(6) * (t238 ^ 2 + t298 ^ 2 + t299 ^ 2) + t551;
t180 = (t189 * t241 * t332 * t531 + t188 * t653) * t533 + t187;
t179 = m(6) * (t238 * t243 + t247 * t299 + t248 * t298) + t601;
t178 = t183 * t626 + t182;
t177 = -t183 * t606 + t181;
t176 = m(6) * (t250 * t268 + t251 * t267) + m(5) * (t300 * t314 + t301 * t313) + t184 * t626 + t536 + t447 * t192 * t629 + t538 * t615 + m(4) * (-t319 * t527 - t321 * t522) * t440 + (m(10) * (t308 * t527 + t309 * t522) * t422 + t611) * t605;
t175 = m(6) * (t234 * t243 + t247 * t268 + t248 * t267) + t179 * t626 + t601;
t174 = m(6) * (t233 * t243 + t247 * t261 + t248 * t260) - t179 * t606 + t601;
t173 = t538 + pkin(1) * t494 / pkin(7) * t475 * ((t523 * (Icges(7,6) * t522 + t527 * t590) + t518 * (Icges(7,5) * t522 + t527 * t597)) * t687 + (t523 * (-Icges(7,6) * t527 + t522 * t590) + t518 * (-Icges(7,5) * t527 + t522 * t597)) * t686 + m(7) * (-t363 * t527 - t365 * t522) * t463 + (t513 / 0.2e1 + t514 / 0.2e1) * (Icges(7,5) * t518 + Icges(7,6) * t523)) + (-t458 * t527 + t522 * t557) * t686 + (t458 * t522 + t527 * t557) * t687 + (m(11) * (t258 * t285 + t259 * t284) + m(8) * (-t318 * t527 - t320 * t522) * t442 + t540) * t625 + t540 - t184 * t606 + t536 + m(3) * (-t364 * t522 - t384 * t527) * t470 + t192 * t614 + m(4) * (t319 * t359 + t321 * t357) + m(8) * (t318 * t360 + t320 * t358) + m(5) * (t300 * t307 + t301 * t306) + m(11) * (t258 * t270 + t259 * t269) + m(6) * (t250 * t261 + t251 * t260) + (-t521 * (-Icges(3,6) * t527 + t522 * t593) + t526 * (-Icges(3,5) * t527 + t522 * t600)) * t686 + (-t521 * (Icges(3,6) * t522 + t527 * t593) + t526 * (Icges(3,5) * t522 + t527 * t600)) * t687;
t1 = [m(2) * (t465 ^ 2 + t471 ^ 2) + m(3) * (t364 ^ 2 + t384 ^ 2) + m(7) * (t363 ^ 2 + t365 ^ 2) + m(9) * (t338 ^ 2 + t339 ^ 2) + m(8) * (t318 ^ 2 + t320 ^ 2) + m(4) * (t319 ^ 2 + t321 ^ 2) + m(10) * (t308 ^ 2 + t309 ^ 2) + m(5) * (t300 ^ 2 + t301 ^ 2) + m(6) * (t250 ^ 2 + t251 ^ 2) + m(11) * (t258 ^ 2 + t259 ^ 2) + (-t329 + t418) * t492 + t310 + (t420 - t657) * t490 + t518 * (Icges(7,1) * t518 + t667) + t523 * (Icges(7,2) * t523 + t668) + t526 * t462 - t521 * t460 + t503 * t438 + t504 * t436 - t499 * t434 - t500 * t435 + t501 * t439 + t502 * t437 - t491 * t419 + t489 * t417 - t477 * t389 - t478 * t390 + Icges(2,3), t173, t176, t185; t173, t546 + t537 + m(4) * (t288 ^ 2 + t357 ^ 2 + t359 ^ 2) + m(8) * (t287 ^ 2 + t358 ^ 2 + t360 ^ 2) + m(11) * (t269 ^ 2 + t270 ^ 2 + t322 ^ 2) + m(5) * (t266 ^ 2 + t306 ^ 2 + t307 ^ 2) + m(6) * (t233 ^ 2 + t260 ^ 2 + t261 ^ 2) + m(3) * (t470 ^ 2 * t634 + t464 ^ 2) + t613 - t527 * (t514 * t393 - t394 * t695) + t522 * (-t393 * t695 + t513 * t394) + pkin(1) ^ 2 * t494 ^ 2 / pkin(7) ^ 2 / t483 ^ 2 * (t522 * (-t391 * t695 + t513 * t392) - t527 * (t514 * t391 - t392 * t695) + m(7) * (t463 ^ 2 * t634 + t466 ^ 2)) - t643 - t527 * t226 + (t546 + 0.2e1 * m(11) * (t269 * t284 + t270 * t285 + t322 * t340) + 0.2e1 * m(8) * (t287 * t289 + (-t358 * t522 - t360 * t527) * t442) + t550 + (m(8) * (t442 ^ 2 * t634 + t289 ^ 2) + m(11) * (t284 ^ 2 + t285 ^ 2 + t340 ^ 2) + t550) * t625) * t625 + (t180 + t696) * t614 + (-t177 - t181) * t606, -t182 * t606 + (t177 * t655 + t180 * t629) * t447 + t535, t174; t176, -t178 * t606 + (t181 * t655 + (t531 * t696 + t627 * t660) * t680) * t447 + t535, m(6) * (t234 ^ 2 + t267 ^ 2 + t268 ^ 2) + m(5) * (t277 ^ 2 + t313 ^ 2 + t314 ^ 2) + m(4) * (t440 ^ 2 * t634 + t290 ^ 2) + t641 * t615 + (t482 * t197 + (m(10) * (t422 ^ 2 * t634 + t424 ^ 2) + t616) * t654) * t534 ^ 2 * t512 ^ 2 * pkin(6) ^ 2 * t654 + pkin(5) ^ 2 * t496 ^ 2 / t448 ^ 2 * t660 + t537 + (t182 + t178) * t626, t175; t185, t174, t175, t492 ^ 2 * t242 + m(6) * (t243 ^ 2 + t247 ^ 2 + t248 ^ 2) + (t527 * t196 + t522 * t195 - t492 * (t231 * t522 + t232 * t527)) * t490;];
Mq = t1;