% Calculate joint inertia matrix for
% picker2Dm1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1DE2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1DE2_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm1DE2_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:25:35
% EndTime: 2020-05-10 20:27:33
% DurationCPUTime: 105.25s
% Computational Cost: add. (2091648->861), mult. (5761406->1455), div. (75736->42), fcn. (1521766->38), ass. (0->661)
t777 = 4 * pkin(1);
t387 = sin(qJ(2));
t340 = pkin(3) * t387;
t621 = 0.8e1 * t340;
t352 = t387 ^ 2;
t351 = t387 * t352;
t434 = pkin(3) ^ 2;
t451 = pkin(3) * t434;
t683 = t351 * t451;
t583 = 0.32e2 * t683;
t429 = pkin(4) ^ 2;
t388 = sin(qJ(1));
t390 = cos(qJ(2));
t674 = t388 * t390;
t316 = pkin(3) * t674;
t515 = pkin(1) * t316;
t440 = pkin(1) ^ 2;
t443 = pkin(7) ^ 2;
t651 = t440 / 0.3e1 + t443;
t250 = -0.4e1 / 0.9e1 * t515 + 0.4e1 / 0.9e1 * t434 - t429 / 0.9e1 + t651;
t367 = -t429 / 0.6e1;
t376 = 0.2e1 / 0.3e1 * t434;
t483 = t443 - t515;
t262 = t367 + t376 + t483;
t375 = 0.4e1 / 0.3e1 * t434;
t347 = t440 + t443;
t369 = -t429 / 0.3e1;
t556 = t369 + t347;
t305 = t375 + t556;
t338 = -t440 / 0.3e1 + t443;
t519 = -0.2e1 * t316;
t391 = cos(qJ(1));
t356 = t391 ^ 2;
t663 = t440 * t356;
t599 = 0.6e1 * t663;
t355 = t391 * t356;
t447 = pkin(1) * t440;
t678 = t355 * t447;
t614 = pkin(7) * t678;
t747 = pkin(7) * t391;
t632 = 0.4e1 * t747;
t668 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t192 = 0.4e1 * t614 + t250 * t599 + t305 * t668 + (t262 * t632 + t338 * t519) * pkin(1);
t410 = 6 * t440;
t497 = -0.4e1 * t515;
t371 = -0.2e1 / 0.3e1 * t429;
t419 = 0.2e1 * t443;
t555 = t371 + t376 + t419;
t452 = t434 ^ 2;
t554 = t371 + t347;
t655 = (t376 + t554) * t347 + t452;
t206 = t305 * t497 + (t410 + t555) * t434 + t655;
t303 = -0.2e1 * t515;
t379 = -t434 / 0.3e1;
t336 = t379 + t443;
t270 = t336 * t303;
t669 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t220 = t305 * t669 + t270;
t358 = 0.10e2 / 0.3e1 * t440;
t221 = (t358 + t555) * t434 + t655;
t346 = -0.3e1 * t434 + t443;
t529 = 0.8e1 * t614;
t289 = t346 * t529;
t552 = t434 + t347;
t314 = -t429 + t552;
t343 = pkin(1) * t391;
t321 = t343 + pkin(7);
t339 = t347 ^ 2;
t438 = t440 ^ 2;
t393 = 15 * t438;
t400 = 0.18e2 * t443;
t401 = -0.2e1 * t429;
t403 = -0.6e1 * t429;
t405 = 0.2e1 * t434;
t442 = t443 ^ 2;
t414 = 0.3e1 * t442;
t431 = t451 ^ 2;
t258 = t303 + t305;
t640 = -t434 + t443;
t774 = 4 * t440;
t326 = t640 * t774;
t625 = 0.4e1 * t343;
t194 = pkin(7) * t258 * t625 + t326 * t356 + t206;
t626 = 0.2e1 * t343;
t325 = pkin(7) * t626;
t349 = -(3 * t440) + t443;
t600 = 0.4e1 * t663;
t269 = t325 + t600 + t349;
t608 = 0.8e1 * t683;
t622 = 0.6e1 * t340;
t472 = t194 * t622 + t269 * t608;
t577 = 0.12e2 * t663;
t681 = t352 * t434;
t581 = 0.12e2 * t681;
t418 = 0.3e1 * t443;
t648 = (15 * t440) + t418;
t776 = 6 * pkin(1);
t159 = t192 * t581 + t289 + t220 * t577 + t431 + (-t429 + t434 + t648) * t452 + (t393 + (t400 + t403 + 0.6e1 * t434) * t440 + t414 + (t401 + t405) * t443) * t434 + t339 * t314 + t472 * t321 + (t206 * t747 - t221 * t316) * t776;
t341 = pkin(3) * t390;
t587 = t440 * t341;
t477 = -t388 * t447 + t587;
t280 = 0.2e1 * t477;
t620 = 0.2e1 * t341;
t323 = pkin(7) * t620;
t412 = 2 * t440;
t345 = t412 + t434;
t348 = -t440 + t443;
t208 = t348 * t341 + t280 * t356 + (t323 * t391 + t345 * t388) * pkin(1);
t411 = 3 * t440;
t644 = t411 + t443;
t330 = t434 + t644;
t288 = t330 * t341;
t331 = 0.3e1 * t434 + t347;
t687 = t331 * t388;
t240 = -pkin(1) * t687 + t288;
t773 = 8 * t440;
t598 = t451 * t773;
t728 = pkin(3) * t438;
t319 = t598 + 0.4e1 * t728;
t565 = t447 * t669;
t243 = t319 * t390 + 0.4e1 * t388 * t565;
t409 = 5 * t438;
t637 = t442 + t452;
t395 = 10 * t440;
t647 = t395 + t419;
t661 = t443 * t440;
t264 = t434 * t647 + t409 + t637 + 0.6e1 * t661;
t404 = 0.5e1 * t452;
t416 = 0.6e1 * t443;
t275 = t404 + (t395 + t416) * t434 + t339;
t569 = t321 * t683;
t507 = -0.8e1 * t569;
t607 = -0.4e1 * t681;
t688 = t321 * t387;
t342 = pkin(1) * t388;
t292 = -t342 + t341;
t505 = t356 * t587;
t633 = 0.2e1 * t747;
t205 = -0.2e1 * t505 + t288 + (t292 * t633 - t687) * pkin(1);
t762 = -0.4e1 * t205;
t172 = t208 * t607 + t243 * t356 + (t688 * t762 + (-t264 + t529) * t390) * pkin(3) + (-0.4e1 * t240 * t747 + (t275 + t507) * t388) * pkin(1);
t282 = t340 + t321;
t759 = 0.2e1 * t387;
t634 = pkin(7) * t759;
t322 = pkin(3) * t634;
t636 = t443 - t429;
t329 = t440 + t636;
t279 = t322 + t329;
t245 = t303 + t279;
t320 = t340 + pkin(7);
t300 = t320 * t391;
t357 = t434 * t774;
t664 = t434 * t443;
t327 = t357 - 0.4e1 * t664;
t420 = -0.2e1 * t443;
t421 = 0.2e1 * pkin(3);
t575 = -0.4e1 * pkin(3) * pkin(7) * t329;
t268 = t322 + t640 + 0.2e1 * t681;
t691 = t268 * t356;
t182 = t327 * t352 + t387 * t575 - t438 - (t443 - (t421 + pkin(4)) * pkin(4)) * (t443 + (t421 - pkin(4)) * pkin(4)) + (t420 + 0.2e1 * t429 - 0.4e1 * t434 - 0.4e1 * t691) * t440 + (-t245 * t300 + t279 * t316) * t777;
t444 = sqrt(t182);
t146 = t159 * t282 + t172 * t444;
t779 = t146 ^ 2;
t778 = 2 * pkin(1);
t422 = 2 * pkin(2);
t428 = t429 ^ 2;
t639 = t438 + t442;
t643 = t419 - t429;
t662 = t443 * t429;
t473 = t643 * t440 + t428 / 0.6e1 + t639 - t662;
t467 = 0.5e1 / 0.6e1 * t452 + t473;
t228 = (t358 + t643) * t434 + t467;
t775 = -0.6e1 * t228;
t263 = 0.4e1 / 0.3e1 * t663 + t325 + t338;
t770 = t418 - t429 - t434;
t312 = t770 * t395;
t402 = -0.5e1 * t429;
t645 = t402 - 0.5e1 * t434;
t350 = t352 ^ 2;
t684 = t350 * t452;
t686 = t339 * (-t434 + t329);
t772 = 0.7e1 * t431 + ((35 * t440) + 0.15e2 * t443 + t645) * t452 + ((21 * t438) + t312 + 0.9e1 * t442 + (t403 - 0.6e1 * t434) * t443) * t434 + t686 - 0.24e2 * t263 * t684;
t368 = -t429 / 0.4e1;
t771 = t368 + t434 / 0.2e1;
t372 = -0.3e1 / 0.2e1 * t429;
t652 = t428 / 0.2e1 - t452 / 0.2e1;
t495 = -0.3e1 * t662 + t414 + t652;
t658 = t347 * ((t372 + t419) * t440 - 0.3e1 / 0.2e1 * t662 + t639 + t652) + t431;
t189 = t515 * t775 + (t393 + (t400 - 0.9e1 * t429) * t440 + t495) * t434 + (t372 + t648) * t452 + t658;
t370 = -t429 / 0.2e1;
t310 = t370 + t552;
t520 = -0.4e1 * t316;
t487 = t310 * t520;
t484 = pkin(1) * t487;
t209 = t484 + (t410 + t643) * t434 + t467;
t222 = t310 * t669 + t270;
t628 = pkin(7) * t343;
t573 = 0.6e1 * t628;
t176 = t209 * t573 + t222 * t577 + t189 + t289;
t415 = 0.8e1 * t443;
t246 = t447 * t520 + t357 + (4 * t438) + (t401 + t415) * t440;
t257 = -t440 + t483 + t771;
t193 = t529 + t246 * t356 + t310 * t349 + (t257 * t632 + t519 * t668) * pkin(1);
t769 = t176 * t621 + t193 * t583;
t768 = 0.2e1 * pkin(7);
t767 = -2 * pkin(1);
t277 = t405 + t279;
t595 = t268 * t343;
t210 = t277 * t320 + 0.2e1 * t595;
t212 = t277 * t391 + (0.4e1 * t356 - 0.2e1) * t320 * pkin(1);
t653 = -t316 + t300;
t251 = pkin(1) + t653;
t171 = t210 * t388 + t212 * t341 + t251 * t444;
t385 = cos(pkin(8));
t721 = sin(pkin(8));
t266 = -t385 * t391 - t388 * t721;
t726 = pkin(5) * t266;
t256 = t726 * t767;
t426 = pkin(5) ^ 2;
t247 = t256 + 0.2e1 * t426;
t254 = -pkin(1) + t726;
t265 = t385 * t388 - t391 * t721;
t657 = t256 + t426;
t233 = -(t778 + pkin(5)) * pkin(5) + t657;
t234 = pkin(5) * (t778 - pkin(5)) + t657;
t446 = sqrt(-t233 * t234);
t187 = -pkin(5) * t265 * t247 - t254 * t446;
t703 = t171 * t187;
t313 = t405 + t411 + t636;
t241 = t313 + t322 + t497;
t301 = t316 - pkin(1);
t672 = t390 * t391;
t588 = pkin(3) * t672;
t689 = t320 * t388;
t253 = t588 + t689;
t694 = t253 * t444;
t170 = -t241 * t300 + t694 + (t301 * t634 + t313 * t674) * pkin(3) + (-0.2e1 * t691 + (0.2e1 * t352 - 0.4e1) * t434 - t329) * pkin(1);
t693 = t265 * t446;
t185 = pkin(5) * t693 - t247 * t254;
t705 = t170 * t185;
t475 = t705 / 0.4e1 + t703 / 0.4e1;
t593 = pkin(1) * t300;
t218 = t303 + t322 + t552 + 0.2e1 * t593;
t216 = 0.1e1 / t218;
t435 = 0.1e1 / pkin(3);
t242 = t440 + t657;
t238 = 0.1e1 / t242;
t697 = t238 / pkin(5);
t570 = t435 * t697;
t500 = t216 * t570;
t141 = t475 * t500;
t702 = t185 * t171;
t704 = t170 * t187;
t474 = t704 / 0.4e1 - t702 / 0.4e1;
t142 = t474 * t500;
t389 = sin(pkin(9));
t392 = cos(pkin(9));
t133 = t141 * t392 - t142 * t389;
t765 = 0.1e1 / t133;
t764 = 0.1e1 / t185;
t763 = 0.1e1 / t133 ^ 2;
t761 = 0.8e1 * t282;
t758 = -0.2e1 * t391;
t424 = (pkin(6) ^ 2);
t757 = 2 * t424;
t756 = -0.6e1 * t444;
t755 = pkin(1) * pkin(7);
t725 = pkin(6) * t133;
t128 = t725 * t422;
t659 = t128 + t424;
t121 = -(t422 + pkin(6)) * pkin(6) + t659;
t122 = pkin(6) * (t422 - pkin(6)) + t659;
t709 = t121 * t122;
t445 = sqrt(-t709);
t96 = 0.1e1 / t445;
t754 = -t96 / 0.2e1;
t297 = t443 + t434 / 0.4e1 + t440 / 0.4e1 - t429 / 0.8e1;
t649 = 0.4e1 / 0.7e1 * t443 - t429 / 0.7e1;
t203 = -0.32e2 / 0.21e2 * t297 * t515 + 0.5e1 / 0.42e2 * t452 + (0.16e2 / 0.21e2 * t440 + t649) * t434 + t438 / 0.7e1 + t649 * t440 + t442 - 0.3e1 / 0.7e1 * t662 + t428 / 0.42e2;
t299 = t651 + t771;
t381 = 0.4e1 / 0.3e1 * t440;
t204 = -0.8e1 / 0.3e1 * t299 * t515 + 0.5e1 / 0.18e2 * t452 + (t381 + t369) * t434 + t442 - t438 / 0.3e1 + t428 / 0.18e2 + (t375 + 0.2e1 / 0.3e1 * t440 + t371) * t443;
t259 = -t452 / 0.6e1 + t473;
t383 = t440 / 0.2e1;
t650 = t383 + t443;
t261 = -0.2e1 / 0.3e1 * t515 + t368 + t650;
t417 = 0.4e1 * t443;
t333 = (t417 + t429) * t440;
t380 = -0.2e1 / 0.3e1 * t434;
t337 = t380 + t443;
t523 = 0.16e2 * t614;
t354 = t356 ^ 2;
t679 = t354 * t438;
t605 = 0.8e1 * t679;
t175 = t337 * t605 + t261 * t523 + 0.14e2 * t203 * t663 + t348 * t452 + (t333 - 0.10e2 / 0.3e1 * t438 + 0.2e1 * t442 - t662) * t434 + t259 * t668 + (0.6e1 * t204 * t747 + t338 * t487) * pkin(1);
t638 = t442 - t438;
t195 = t336 * t484 - t431 + (-t358 - t636) * t452 + (t333 + t452 / 0.6e1 - t428 / 0.6e1 + t638) * t434 + t259 * t443;
t408 = 7 * t438;
t199 = (t372 + t418 + (7 * t440)) * t452 + (t408 + (t402 + 0.10e2 * t443) * t440 + t495) * t434 + t658;
t456 = pkin(7) * t443;
t328 = -0.12e2 * pkin(7) * t447 + t456 * t777;
t335 = -(8 * t438) + 0.12e2 * t661;
t213 = t328 * t391 + t335 * t356 + t523 + t605 + t639 - 0.6e1 * t661;
t224 = t303 * t669 + t310 * t346;
t290 = 0.16e2 * (t637 - 0.6e1 * t664) * t438;
t344 = -0.30e2 * t429 + 0.60e2 * t443;
t642 = t428 - t452;
t478 = 0.6e1 * t442 + t642 - 0.6e1 * t662;
t140 = t290 * t354 + 0.32e2 * t224 * t614 + 0.24e2 * t195 * t663 + (t401 + t417 + (28 * t440)) * t431 + t314 * t686 + (0.24e2 * t175 * t352 + t344 * t438 + t442 * t403 + t478 * t410 + t642 * t419 + (28 * t447 ^ 2) + 0.4e1 * t456 ^ 2) * t434 + t769 * t321 + 0.8e1 * (t189 * t747 - t199 * t316) * pkin(1) + (0.16e2 * t213 * t350 + t344 * t440 + (70 * t438) + t452 + t478) * t452;
t295 = 0.7e1 / 0.6e1 * t434 + t367 + t650;
t377 = t434 / 0.3e1;
t557 = t367 + t377 + t443;
t298 = t381 + t557;
t219 = -t295 * t342 + t298 * t341;
t307 = t440 + t557;
t332 = t405 + t348;
t230 = t307 * t341 - t332 * t342 / 0.2e1;
t558 = t429 / 0.3e1 + t377 + t419;
t494 = -0.8e1 / 0.3e1 * t679 + t434 * t348 - 0.5e1 / 0.3e1 * t438 + t558 * t440 + t443 * (t369 + t336);
t746 = pkin(7) * t438;
t631 = -0.4e1 * t746;
t359 = -0.20e2 / 0.3e1 * t440;
t559 = 0.2e1 / 0.3e1 * t429 + t376 + t417;
t560 = 0.4e1 / 0.3e1 * t429 + t375 + t420;
t736 = t452 / 0.2e1 - (t359 + t559) * t434 / 0.2e1 + 0.3e1 / 0.2e1 * t438 - t560 * t440 / 0.2e1 - t442 / 0.2e1;
t177 = t388 * t355 * t631 + t219 * t600 + t494 * t341 + (t230 * t632 + t388 * t736) * pkin(1);
t646 = t401 - 0.2e1 * t434;
t553 = t416 + t646;
t215 = t452 + (t371 + t380 + t647) * t434 + t409 + t553 * t440 + t443 * (t371 + t337);
t211 = t215 * t341;
t304 = 0.8e1 / 0.3e1 * t434 + t556;
t306 = t369 + t376 + t644;
t223 = -t304 * t342 + t306 * t341;
t311 = 0.5e1 / 0.6e1 * t434 + t383 + t367;
t236 = t311 * t620 + t342 * t669;
t586 = t447 * t341;
t504 = t355 * t586;
t602 = -0.4e1 * t663;
t229 = t404 + (t395 + t553) * t434 + (t380 + t554) * t347;
t698 = t229 * t388;
t178 = -0.8e1 * pkin(7) * t504 + t236 * t602 + t211 + (t223 * t632 - t698) * pkin(1);
t196 = -pkin(1) * t698 + t211;
t227 = -0.3e1 * t452 + (t359 + t560) * t434 + t559 * t440 + t638;
t231 = -0.5e1 / 0.3e1 * t452 + (-t440 + t558) * t434 + t443 * (t379 + t556);
t627 = -0.2e1 * t342;
t197 = t227 * t341 + t231 * t627;
t281 = 0.4e1 * t477;
t291 = 0.2e1 * t342 + t341;
t309 = t370 + t330;
t200 = t349 * t341 + t281 * t356 + (t291 * t747 + t309 * t388) * t778;
t296 = t443 + 0.5e1 / 0.2e1 * t434 + 0.3e1 / 0.2e1 * t440 + t370;
t685 = t346 * t388;
t232 = t296 * t341 + pkin(1) * t685 / 0.2e1;
t486 = 0.24e2 * t336 * t679 - t431 - ((21 * t440) + t770) * t452 - (t443 * t646 + t312 + t414 + (35 * t438)) * t434 - (t408 + (t415 + t645) * t440 + t443 * (-t434 + t636)) * t347;
t582 = -0.12e2 * t681;
t147 = t232 * t523 + t200 * t507 + t177 * t582 - 0.6e1 * t197 * t663 + (-0.6e1 * t178 * t688 + t390 * t486) * pkin(3) + (-0.6e1 * t196 * t747 + t772 * t388) * pkin(1);
t118 = t140 * t282 + t147 * t444;
t144 = 0.1e1 / t146;
t737 = t216 / 0.2e1;
t537 = t435 * t737;
t152 = qJ(1) + atan2(t171 * t537, t170 * t537);
t430 = 0.1e1 / pkin(4);
t666 = t430 * t435;
t530 = t666 / 0.2e1;
t102 = atan2(t444 * t530, t118 * t144 * t530) + t152;
t98 = sin(t102);
t753 = pkin(4) * t98;
t99 = cos(t102);
t752 = pkin(4) * t99;
t450 = pkin(2) ^ 2;
t125 = t450 + t659;
t425 = 0.1e1 / pkin(6);
t550 = 0.1e1 / t125 * t425 / 0.2e1;
t126 = t128 + t757;
t127 = -pkin(2) - t725;
t135 = t141 * t389 + t142 * t392;
t724 = pkin(6) * t135;
t64 = -t126 * t127 - t445 * t724;
t65 = t126 * t724 - t127 * t445;
t45 = atan2(t65 * t550, t64 * t550) + t152;
t43 = sin(t45);
t751 = pkin(6) * t43;
t44 = cos(t45);
t750 = pkin(6) * t44;
t749 = pkin(7) * t356;
t748 = pkin(7) * t388;
t745 = t144 / 0.2e1;
t744 = t144 / 0.4e1;
t673 = t388 * t391;
t690 = t320 * t356;
t190 = (t268 * t673 + t341 * t690) * t773 + (t245 * t689 + t279 * t588) * t777;
t353 = t388 ^ 2;
t181 = 0.1e1 / t444;
t739 = t181 / 0.2e1;
t540 = t251 * t739;
t156 = -t694 + t190 * t540 + t268 * t353 * t767 + t210 * t391 + (-t277 - 0.8e1 * t593) * t316;
t743 = t156 / 0.4e1;
t742 = -t170 / 0.4e1;
t741 = t170 / 0.4e1;
t740 = t171 / 0.4e1;
t738 = t185 / 0.4e1;
t735 = -t444 / 0.4e1;
t734 = t444 / 0.4e1;
t733 = t445 / 0.2e1;
t732 = 0.2e1 * t265 ^ 2;
t731 = pkin(1) * t338;
t730 = pkin(1) * t356;
t729 = pkin(2) / t125 ^ 2;
t727 = pkin(3) * t444;
t723 = t199 * pkin(1);
t547 = t144 * t740;
t469 = t118 * t547 + t170 * t734;
t665 = t430 / pkin(3) ^ 2;
t571 = t216 * t665;
t104 = t469 * t571;
t677 = t387 * t388;
t271 = t672 + t677;
t100 = t271 * t104;
t548 = t144 * t741;
t468 = t118 * t548 + t171 * t735;
t103 = t468 * t571;
t676 = t387 * t391;
t272 = t674 - t676;
t73 = -t103 * t272 - t100;
t71 = 0.1e1 / t73 ^ 2;
t101 = t272 * t104;
t72 = -t103 * t271 + t101;
t722 = t71 * t72;
t150 = sin(t152);
t169 = 0.1e1 / t170 ^ 2;
t506 = pkin(3) / (t169 * t171 ^ 2 + 0.1e1) * t218 * t435;
t476 = -0.2e1 * t169 * t171 * t506;
t480 = 0.2e1 / t170 * t506;
t217 = 0.1e1 / t218 ^ 2;
t517 = pkin(1) * t588;
t590 = pkin(3) * t677;
t516 = pkin(1) * t590;
t654 = 0.2e1 * t516 + t323;
t237 = 0.2e1 * t517 + t654;
t700 = t217 * t237;
t536 = -t700 / 0.2e1;
t671 = t390 * t434;
t566 = t387 * t671;
t274 = t323 + 0.4e1 * t566;
t611 = t434 * t748;
t179 = t274 * t602 + (t327 * t759 + t575) * t390 + (-t654 * t300 + 0.2e1 * t390 ^ 2 * t611 + (-t245 * t672 - t279 * t677) * pkin(3)) * t777;
t539 = t253 * t739;
t589 = pkin(3) * t676;
t675 = t387 * t434;
t155 = (t316 - t589) * t444 + t179 * t539 - 0.2e1 * t274 * t730 - (t323 + 0.4e1 * t516) * t300 - 0.2e1 * t352 * t611 - t313 * t590 + (t675 * t777 + (-t241 * t391 + t301 * t768) * pkin(3)) * t390;
t706 = t155 * t216;
t154 = t179 * t540 + t274 * t388 * t626 + ((t388 * t444 - t212) * t387 + (t391 * t444 + (t320 * t768 + t277) * t388 + (-pkin(1) + 0.2e1 * t730 + t747) * t620) * t390) * pkin(3);
t707 = t154 * t216;
t93 = (t707 / 0.2e1 + t171 * t536) * t480 + (t706 / 0.2e1 + t170 * t536) * t476;
t720 = t150 * t93;
t158 = t653 * t444 + t190 * t539 + (t241 * t320 + 0.4e1 * t595) * t388 + (t633 * t675 + (t313 * t391 + t690 * t777) * pkin(3)) * t390;
t244 = t253 * t778;
t699 = t217 * t244;
t534 = t699 / 0.2e1;
t94 = 0.1e1 + (t156 * t737 + t171 * t534) * t480 + (t158 * t737 + t170 * t534) * t476;
t719 = t150 * t94;
t151 = cos(t152);
t718 = t151 * t93;
t717 = t151 * t94;
t716 = -0.6e1 * t755;
t715 = -0.4e1 * t755;
t488 = pkin(7) * t505;
t294 = -0.4e1 * t488;
t324 = pkin(7) * t627;
t489 = t586 * t749;
t498 = t355 * t565;
t528 = pkin(7) * t600;
t541 = t190 * t739;
t670 = t391 * t440;
t564 = t388 * t670;
t578 = -0.24e2 * t670;
t606 = 0.8e1 * t681;
t609 = -0.8e1 * t683;
t660 = t447 * t356;
t612 = pkin(7) * t660;
t623 = -0.4e1 * t340;
t521 = -0.24e2 * t612;
t656 = -0.24e2 * t336 * t504 + t521 * t685;
t119 = (t353 * t351 * t598 + (t343 * t345 - 0.2e1 * t678) * t607 + 0.4e1 * t498 + t331 * t528 + t275 * t343) * t444 + t172 * t541 + ((-0.8e1 / 0.3e1 * t504 + t294 - 0.2e1 * t338 * t517) * t581 - 0.24e2 * t305 * t488 - 0.6e1 * t221 * t517 + t656) * t282 + ((t280 * t391 * t606 + t243 * t758 - 0.24e2 * t489) * t444 + (0.12e2 * (-t250 * t670 - t612) * t581 + t220 * t578) * t282 + (0.4e1 * t205 * t387 * t727 - t472 * t282 - t159 + ((t341 * t606 + 0.4e1 * t240) * t444 + (-0.48e2 * t262 * t681 - 0.6e1 * t206) * t282) * pkin(7)) * pkin(1)) * t388 + ((t609 * t343 + ((0.4e1 * t316 * t391 - 0.2e1 * t749) * t440 + (-0.2e1 * t292 * t748 - t331 * t391) * pkin(1)) * t623) * t444 + ((t324 - 0.8e1 * t564) * t608 + (-0.8e1 * t488 - 0.2e1 * t326 * t673 + (-t258 * t748 - t305 * t588) * t777) * t622) * t282) * t321;
t465 = t217 * t469;
t544 = t170 * t181 / 0.8e1;
t145 = 0.1e1 / t779;
t545 = -t145 * t171 / 0.4e1;
t527 = 0.32e2 / 0.3e1 * t438;
t490 = t355 * t527;
t491 = 0.64e2 / 0.3e1 * t297 * t447;
t502 = 0.32e2 * t569;
t509 = t669 * t746;
t522 = -0.48e2 * t612;
t579 = -0.32e2 * t355 * t438;
t580 = 0.24e2 * t681;
t584 = -0.96e2 * t336 * t355;
t592 = pkin(1) * t668;
t594 = t310 * t731;
t601 = -0.2e1 * t663;
t603 = 0.8e1 * t670;
t604 = -0.4e1 * t678;
t613 = pkin(7) * t663;
t624 = -0.6e1 * t340;
t81 = ((t309 * t626 + t528 + t604) * t507 + (0.12e2 * t353 * t356 * t746 + t295 * t604 + t354 * t631) * t582 + 0.12e2 * t231 * t678 + (t736 * t582 + t772) * t343 + (t332 * t601 * t582 + t229 * t599 + t346 * t605) * pkin(7)) * t444 + t147 * t541 + t294 * t282 * t502 + (((pkin(7) * t304 * t602 - t229 * t343 - 0.4e1 * t498) * t756 + t656 * t761) * t688 + ((-pkin(7) * t354 * t527 - 0.16e2 * t299 * t613 - t355 * t491 - 0.4e1 * t391 * t594) * t580 - 0.64e2 * t354 * t509 + t447 * t310 * t584 - 0.48e2 * t228 * t613 - 0.8e1 * t391 * t723 + ((t592 * t758 + t604) * t583 + (-0.24e2 * t310 * t613 + t343 * t775) * t621) * t321) * t282 * t390) * pkin(3) + ((0.16e2 * (t335 * t758 - t328 + t522 + t579) * t684 + (t246 * t758 + t257 * t715 + t521) * t502 + (-0.28e2 * t203 * t670 + t204 * t716 + t261 * t522 + t337 * t579) * t580 + t321 * (t209 * t716 + t222 * t578) * t621 - 0.4e1 * t290 * t355 - 0.96e2 * t224 * t612 - 0.48e2 * t195 * t670 - 0.8e1 * t189 * t755) * t282 + ((-0.8e1 * t219 * t670 + t341 * t490) * t582 + t390 * t584 * t728 + t232 * t522 + 0.12e2 * t197 * t670 + (0.2e1 * (-t281 * t391 - t291 * t755) * t609 + (t223 * t715 + t236 * t603 + 0.24e2 * t489) * t624) * t321) * t444 + (-t769 * t282 - t140 + (t178 * t622 + t200 * t608 + (-0.24e2 * t324 + 0.64e2 * t564) * t684 + (0.48e2 * t230 * t681 + 0.6e1 * t196) * pkin(7)) * t444) * pkin(1)) * t388;
t714 = t103 + (t244 * t465 + (t158 * t734 + t190 * t544 + t81 * t547 + (t119 * t545 + t144 * t743) * t118) * t216) * t665;
t567 = t336 * t660;
t260 = 0.24e2 * t567 * t590;
t525 = pkin(7) * t589;
t481 = t440 * t388 * t525;
t278 = 0.4e1 * t481;
t293 = t525 * t767;
t542 = t179 * t739;
t680 = t352 * t451;
t568 = t321 * t680;
t120 = (t293 * t607 + (-0.8e1 * t208 * t671 - t319 * t356 + ((-t348 + t601) * t607 + t264 + (t330 * t625 - 0.8e1 * t678) * pkin(7)) * pkin(3)) * t387 + ((t293 + (-t330 + 0.2e1 * t663) * t340) * t623 + (pkin(3) * t762 - 0.24e2 * t342 * t680) * t390) * t321) * t444 + t172 * t542 + (0.24e2 * t269 * t390 * t568 + (t278 + (0.8e1 / 0.3e1 * t660 + 0.2e1 * t731) * t590) * t581 + 0.24e2 * t192 * t566 + t260 + 0.24e2 * t305 * t481 + 0.6e1 * t221 * t516 + 0.6e1 * ((pkin(7) * t603 + t305 * t777) * t388 * t681 + t194 * t341) * t321) * t282 + t159 * t341;
t574 = -0.4e1 * t628;
t610 = pkin(7) * t670;
t682 = t351 * t452;
t91 = t147 * t542 + (0.32e2 * t278 * t282 - 0.8e1 * t293 * t444) * t569 + (0.24e2 * (-0.4e1 * t263 * t342 * t682 - t177 * t675 - t200 * t568) * t444 + (0.48e2 * t175 * t675 + 0.96e2 * t193 * t568 + 0.64e2 * t213 * t682) * t282) * t390 + ((t140 + (t176 * t761 + t178 * t756) * t321) * t390 + (t321 * t260 * t761 + ((t298 * t602 + t307 * t574 - t494) * t582 - 0.16e2 * t296 * t614 + t227 * t599 + t215 * t573 + ((-t349 + t602) * t609 + (t306 * t574 + 0.8e1 * t311 * t663 - t215 + t529) * t624) * t321 - t486) * t444 + ((pkin(7) * t490 + 0.16e2 * t299 * t610 + t356 * t491 + 0.4e1 * t594) * t580 + 0.64e2 * t355 * t509 + 0.96e2 * t310 * t567 + 0.48e2 * t228 * t610 + 0.8e1 * t723 + ((0.2e1 * t592 + 0.4e1 * t660) * t583 + (t228 * t776 + 0.24e2 * t310 * t610) * t621) * t321) * t282 * t388) * t387) * pkin(3);
t713 = t103 - (-t237 * t465 + (t155 * t734 + t179 * t544 + t91 * t547 + (t120 * t545 + t154 * t744) * t118) * t216) * t665;
t712 = pkin(3) * t719 + t342;
t711 = pkin(2) * t719 + t342;
t710 = 0.1e1 / t118 ^ 2 * t779;
t708 = t763 * t135;
t526 = t265 * t778;
t701 = 0.1e1 / t446 * (t233 + t234) * pkin(5) * t526;
t696 = t238 / pkin(1);
t239 = 0.1e1 / t242 ^ 2;
t695 = t239 * t265;
t692 = t266 * t446;
t667 = t425 / pkin(2);
t538 = t701 / 0.2e1;
t160 = (-t692 + (t254 * t778 - t247 + t538) * t265) * pkin(5);
t161 = -t254 * t701 / 0.2e1 + t426 * pkin(1) * t732 + (t247 * t266 - t693) * pkin(5);
t499 = t697 * t699;
t596 = pkin(1) * t695;
t82 = (t475 * t499 + ((t705 / 0.2e1 + t703 / 0.2e1) * t596 + (t158 * t738 + t160 * t741 + t161 * t740 + t187 * t743) * t697) * t216) * t435;
t83 = (-t474 * t499 + ((-t704 / 0.2e1 + t702 / 0.2e1) * t596 + (-t158 * t187 / 0.4e1 + t161 * t742 + t160 * t740 + t156 * t738) * t697) * t216) * t435;
t60 = t389 * t83 + t392 * t82;
t619 = t60 * t729;
t535 = -t700 / 0.4e1;
t470 = t706 / 0.4e1 + t170 * t535;
t471 = t707 / 0.4e1 + t171 * t535;
t109 = (t185 * t470 + t187 * t471) * t570;
t110 = (t185 * t471 - t187 * t470) * t570;
t77 = t109 * t392 + t110 * t389;
t618 = t77 * t729;
t617 = pkin(2) * t718;
t616 = pkin(3) * t718;
t63 = 0.1e1 / t64 ^ 2;
t615 = pkin(6) * t125 / (t63 * t65 ^ 2 + 0.1e1);
t585 = pkin(5) * t695;
t576 = 0.1e1 / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t450 * t763 * t709 / 0.4e1) * t667;
t563 = t127 * t754;
t562 = -t765 * t96 / 0.4e1;
t561 = t135 * t754;
t551 = -t118 * t145 / 0.2e1;
t549 = t763 * t733;
t546 = t145 * t742;
t543 = -t171 * t181 / 0.8e1;
t533 = t697 / 0.2e1;
t532 = -t696 / 0.2e1;
t531 = t696 / 0.2e1;
t524 = pkin(2) * t135 * t757;
t518 = -0.2e1 * pkin(2) * t127 + t126;
t485 = -0.2e1 * t63 * t65 * t615;
t496 = 0.2e1 / t64 * t615;
t501 = pkin(6) * t550;
t479 = pkin(6) * (-t121 - t122) * t422;
t59 = t77 * t479;
t78 = t109 * t389 - t110 * t392;
t16 = ((t59 * t563 + t77 * t524 + (t126 * t78 + t445 * t77) * pkin(6)) * t550 - t65 * t618) * t496 + (-t64 * t618 + (-t445 * t78 + t518 * t77 + t561 * t59) * t501) * t485 + t93;
t106 = 0.1e1 / (t182 * t710 + 0.1e1);
t466 = -0.2e1 * pkin(4) * t106 * t666 * t710 * t727;
t482 = t106 / t118 * t146 * t739;
t42 = t179 * t482 + (t120 * t551 + t745 * t91) * t466 + t93;
t493 = -pkin(2) * t717 - t343;
t492 = -pkin(3) * t717 - t343;
t46 = t60 * t479;
t61 = t389 * t82 - t392 * t83;
t8 = ((t46 * t563 + t60 * t524 + (t126 * t61 + t445 * t60) * pkin(6)) * t550 - t65 * t619) * t496 + (-t64 * t619 + (-t445 * t61 + t46 * t561 + t518 * t60) * t501) * t485 + t94;
t41 = t190 * t482 + (t119 * t551 + t745 * t81) * t466 + t94;
t464 = t217 * t468;
t287 = -rSges(2,1) * t391 + rSges(2,2) * t388;
t286 = rSges(8,1) * t390 - rSges(8,2) * t387;
t285 = rSges(2,1) * t388 + rSges(2,2) * t391;
t284 = rSges(8,1) * t387 + rSges(8,2) * t390;
t255 = -pkin(1) * t266 + pkin(5);
t248 = t256 + t412;
t188 = -pkin(1) * t265 * t248 + t255 * t446;
t186 = pkin(1) * t693 + t248 * t255;
t184 = 0.1e1 / t186 ^ 2;
t183 = 0.1e1 / t185 ^ 2;
t167 = qJ(1) + atan2(t187 * t533, t185 * t533);
t166 = pkin(8) + atan2(t188 * t531, t186 * t532);
t165 = cos(t167);
t164 = sin(t167);
t163 = cos(t166);
t162 = sin(t166);
t149 = rSges(6,1) * t163 - rSges(6,2) * t162;
t148 = rSges(6,1) * t162 + rSges(6,2) * t163;
t139 = -rSges(3,1) * t151 + rSges(3,2) * t150;
t138 = -rSges(3,1) * t150 - rSges(3,2) * t151;
t117 = (-((t255 * t538 + t440 * pkin(5) * t732 + (t248 * t266 - t693) * pkin(1)) * t531 + t188 * t585) / t186 - (-t186 * t585 + (-t692 + (-0.2e1 * t255 * pkin(5) - t248 + t538) * t265) * pkin(1) * t532) * t188 * t184) / (t184 * t188 ^ 2 + 0.1e1) * t242 * t778;
t113 = 0.1e1 + (t161 * t764 * t697 + (-t160 * t183 * t697 + (-t185 * t183 + t764) * t239 * t526) * t187) * t242 / (t183 * t187 ^ 2 + 0.1e1) * pkin(5);
t108 = -t343 + t113 * (rSges(9,1) * t165 - rSges(9,2) * t164);
t107 = t342 - t113 * (rSges(9,1) * t164 + rSges(9,2) * t165);
t97 = 0.1e1 / (t135 ^ 2 * t763 + 0.1e1);
t90 = atan2(t135, t133) + t152;
t89 = cos(t90);
t88 = sin(t90);
t85 = pkin(2) * t720;
t84 = pkin(3) * t720;
t80 = t139 * t94 - t343;
t79 = -t138 * t94 + t342;
t70 = 0.1e1 / t73;
t69 = -rSges(5,1) * t99 + rSges(5,2) * t98;
t68 = -rSges(5,1) * t98 - rSges(5,2) * t99;
t67 = rSges(7,1) * t89 - rSges(7,2) * t88;
t66 = rSges(7,1) * t88 + rSges(7,2) * t89;
t58 = 0.1e1 / (t71 * t72 ^ 2 + 0.1e1);
t55 = atan2(t72, t73) + t102;
t54 = cos(t55);
t53 = sin(t55);
t50 = (-t237 * t464 + (t91 * t548 + t154 * t735 + t179 * t543 + (t120 * t546 + t155 * t744) * t118) * t216) * t665;
t48 = (t244 * t464 + (t81 * t548 + t156 * t735 + t190 * t543 + (t119 * t546 + t158 * t744) * t118) * t216) * t665;
t40 = (-t708 * t77 + t765 * t78) * t97 + t93;
t38 = atan2(t667 * t733, -t133) + t45;
t37 = cos(t38);
t36 = sin(t38);
t35 = (-t60 * t708 + t61 * t765) * t97 + t94;
t34 = rSges(11,1) * t54 - rSges(11,2) * t53;
t33 = rSges(11,1) * t53 + rSges(11,2) * t54;
t32 = t42 * t69 - t616;
t31 = -t42 * t68 + t84;
t30 = t41 * t69 + t492;
t29 = -t41 * t68 + t712;
t28 = t35 * t67 - t343;
t27 = -t35 * t66 + t342;
t26 = rSges(4,1) * t44 - rSges(4,2) * t43;
t25 = rSges(4,1) * t43 + rSges(4,2) * t44;
t24 = -rSges(10,1) * t37 + rSges(10,2) * t36;
t23 = -rSges(10,1) * t36 - rSges(10,2) * t37;
t18 = ((-t271 * t50 - t272 * t713 - t100) * t70 - ((-t104 - t50) * t272 + t713 * t271) * t722) * t58 + t42;
t17 = ((t714 * t272 + (t104 - t48) * t271) * t70 - (-t271 * t714 - t272 * t48 + t101) * t722) * t58 + t41;
t15 = (t549 * t77 + t562 * t59) * t576 + t16;
t14 = t18 * t34 - t42 * t752 - t616;
t13 = -t18 * t33 + t42 * t753 + t84;
t12 = t17 * t34 - t41 * t752 + t492;
t11 = -t17 * t33 + t41 * t753 + t712;
t10 = t16 * t26 - t617;
t9 = -t16 * t25 + t85;
t7 = t26 * t8 + t493;
t6 = -t25 * t8 + t711;
t5 = (t46 * t562 + t549 * t60) * t576 + t8;
t4 = t15 * t24 + t16 * t750 - t617;
t3 = -t15 * t23 - t16 * t751 + t85;
t2 = t24 * t5 + t8 * t750 + t493;
t1 = -t23 * t5 - t8 * t751 + t711;
t19 = [t94 ^ 2 * Icges(3,3) + t35 ^ 2 * Icges(7,3) + t5 ^ 2 * Icges(10,3) + t113 ^ 2 * Icges(9,3) + t41 ^ 2 * Icges(5,3) + t8 ^ 2 * Icges(4,3) + t17 ^ 2 * Icges(11,3) + (t1 ^ 2 + t2 ^ 2) * m(10) + (t27 ^ 2 + t28 ^ 2) * m(7) + (t11 ^ 2 + t12 ^ 2) * m(11) + (t29 ^ 2 + t30 ^ 2) * m(5) + m(9) * (t107 ^ 2 + t108 ^ 2) + (t6 ^ 2 + t7 ^ 2) * m(4) + (t79 ^ 2 + t80 ^ 2) * m(3) + m(2) * (t285 ^ 2 + t287 ^ 2) + Icges(2,3) + (Icges(6,3) + m(6) * (t148 ^ 2 + t149 ^ 2)) * t117 ^ 2; t5 * Icges(10,3) * t15 + (t10 * t7 + t6 * t9) * m(4) + t41 * Icges(5,3) * t42 + t8 * Icges(4,3) * t16 + (t1 * t3 + t2 * t4) * m(10) + t17 * Icges(11,3) * t18 + (t11 * t13 + t12 * t14) * m(11) + (t29 * t31 + t30 * t32) * m(5) + (t94 * Icges(3,3) + (-t138 * t79 + t139 * t80) * m(3)) * t93 + (t35 * Icges(7,3) + (-t27 * t66 + t28 * t67) * m(7)) * t40; t18 ^ 2 * Icges(11,3) + t42 ^ 2 * Icges(5,3) + (t3 ^ 2 + t4 ^ 2) * m(10) + t16 ^ 2 * Icges(4,3) + (t13 ^ 2 + t14 ^ 2) * m(11) + (t31 ^ 2 + t32 ^ 2) * m(5) + (t10 ^ 2 + t9 ^ 2) * m(4) + t15 ^ 2 * Icges(10,3) + m(8) * (t284 ^ 2 + t286 ^ 2) + Icges(8,3) + (Icges(3,3) + (t138 ^ 2 + t139 ^ 2) * m(3)) * t93 ^ 2 + ((t66 ^ 2 + t67 ^ 2) * m(7) + Icges(7,3)) * t40 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t19(1), t19(2); t19(2), t19(3);];
Mq = res;
