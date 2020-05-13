% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m2OL_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2OL_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 02:46:34
% EndTime: 2020-05-03 02:47:10
% DurationCPUTime: 15.91s
% Computational Cost: add. (3829->579), mult. (4845->845), div. (0->0), fcn. (1289->98), ass. (0->442)
t637 = pkin(3) * m(7);
t246 = m(6) * rSges(6,1) + t637;
t656 = m(7) * rSges(7,3);
t659 = m(6) * rSges(6,2);
t254 = t656 + t659;
t748 = t254 * cos(qJ(5)) + t246 * sin(qJ(5));
t361 = qJ(4) + qJ(5);
t747 = t254 * cos(t361) + t246 * sin(t361);
t746 = t748 * pkin(5) * qJD(5);
t330 = qJ(3) + t361;
t745 = t254 * cos(t330) + t246 * sin(t330);
t397 = 2 * qJ(2);
t363 = qJ(3) + t397;
t329 = qJ(4) + t363;
t289 = qJ(5) + t329;
t720 = t254 * cos(t289) + t246 * sin(t289);
t394 = 0.2e1 * qJ(5);
t579 = 0.2e1 * qJ(3) + t397;
t516 = 0.2e1 * qJ(4) + t579;
t230 = qJ(6) + t394 + t516;
t370 = rSges(7,2) * rSges(7,3);
t624 = rSges(7,1) * pkin(3);
t682 = -t624 + t370;
t741 = t682 * m(7) - Icges(7,6);
t371 = rSges(7,1) * rSges(7,3);
t625 = pkin(3) * rSges(7,2);
t681 = -t625 - t371;
t742 = t681 * m(7) + Icges(7,5);
t744 = t742 * cos(t230) + t741 * sin(t230);
t285 = qJ(5) + t516;
t238 = qJ(6) + t285;
t180 = sin(t238);
t470 = -qJ(6) + t516;
t239 = qJ(5) + t470;
t181 = sin(t239);
t743 = pkin(5) * (t181 + t180);
t683 = -t624 - t370;
t740 = t683 * m(7) + Icges(7,6);
t357 = qJ(6) + qJ(2);
t331 = qJ(3) + t357;
t283 = qJ(4) + t331;
t240 = qJ(5) + t283;
t159 = 0.2e1 * t240;
t139 = sin(t159);
t358 = -qJ(6) + qJ(2);
t332 = qJ(3) + t358;
t284 = qJ(4) + t332;
t241 = qJ(5) + t284;
t160 = 0.2e1 * t241;
t140 = sin(t160);
t399 = rSges(7,2) ^ 2;
t400 = rSges(7,1) ^ 2;
t729 = Icges(7,1) - Icges(7,2);
t164 = m(7) * (-t399 + t400) - t729;
t290 = qJ(2) + t330;
t393 = 0.2e1 * qJ(6);
t236 = t393 + t290;
t178 = sin(t236);
t237 = -0.2e1 * qJ(6) + t290;
t179 = sin(t237);
t223 = cos(t290);
t390 = qJD(1) ^ 2;
t548 = qJD(1) * qJD(6);
t457 = t164 * t548 / 0.2e1;
t578 = pkin(1) * t390;
t198 = 0.2e1 * t290;
t610 = (rSges(6,1) * t659 + rSges(7,3) * t637 - Icges(6,4)) * cos(t198);
t231 = t394 + t470;
t615 = t740 * sin(t231);
t650 = t390 / 0.2e1;
t739 = (t179 + t178) * t457 + t254 * t223 * t578 + (-t615 - t744) * t650 + (t610 + (t140 + t139) * t164 / 0.8e1) * t390;
t347 = qJD(3) + qJD(4);
t299 = qJD(2) + t347;
t738 = t299 ^ 2;
t350 = qJD(2) + qJD(3);
t737 = (t350 ^ 2);
t710 = -2 * Icges(7,4);
t736 = t370 * m(7) / 0.2e1;
t182 = sin(t240);
t344 = qJD(5) + qJD(6);
t294 = qJD(4) + t344;
t247 = qJD(3) + t294;
t202 = qJD(2) + t247;
t734 = t182 * t202;
t183 = sin(t241);
t345 = qJD(5) - qJD(6);
t295 = qJD(4) + t345;
t248 = qJD(3) + t295;
t203 = qJD(2) + t248;
t733 = t203 * t183;
t384 = rSges(5,1) * m(5);
t667 = pkin(5) * m(7);
t196 = m(6) * pkin(5) + t384 + t667;
t362 = qJ(3) + qJ(4);
t333 = qJ(2) + t362;
t267 = sin(t333);
t732 = t267 * pkin(1) * t196;
t388 = m(5) + m(6);
t638 = pkin(2) * m(7);
t143 = t388 * pkin(2) + m(4) * rSges(4,1) + t638;
t364 = qJ(2) + qJ(3);
t314 = sin(t364);
t731 = t314 * pkin(1) * t143;
t229 = 0.2e1 * t333;
t619 = (pkin(5) ^ 2 * (m(6) + m(7)) + (rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - Icges(5,1) + Icges(5,2)) * sin(t229);
t277 = 0.2e1 * t364;
t730 = (pkin(2) ^ 2 * (m(7) + t388) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2)) * sin(t277);
t245 = -t371 + t625;
t728 = m(7) * t245 + Icges(7,5);
t574 = pkin(2) * t390;
t328 = qJ(4) + t579;
t288 = qJ(5) + t328;
t686 = t254 * cos(t288);
t100 = t574 * t686;
t191 = cos(t237);
t228 = t390 + (2 * t737);
t359 = qJ(5) + qJ(6);
t326 = qJ(4) + t359;
t260 = sin(t326);
t360 = qJ(5) - qJ(6);
t327 = qJ(4) + t360;
t261 = sin(t327);
t623 = rSges(7,3) * qJD(6);
t342 = rSges(7,2) * t623;
t352 = sin(t393);
t372 = sin(qJ(6));
t391 = pkin(3) ^ 2;
t398 = rSges(7,3) ^ 2;
t508 = t400 / 0.2e1 + t399 / 0.2e1;
t643 = (rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6);
t409 = (t398 - t391 - t508) * m(7) - t643 + Icges(6,1) + Icges(7,1) / 0.2e1 - Icges(6,2) + Icges(7,2) / 0.2e1 - Icges(7,3);
t234 = -qJ(6) + t288;
t188 = cos(t234);
t536 = rSges(7,2) * t638;
t486 = t390 * t536;
t444 = t188 * t486;
t232 = qJ(6) + t288;
t186 = cos(t232);
t445 = t186 * t486;
t645 = m(7) * t390;
t525 = t645 / 0.4e1;
t482 = rSges(7,1) * t525;
t269 = cos(t327);
t484 = t269 * t536;
t657 = m(7) * rSges(7,2);
t278 = rSges(7,1) * t657 - Icges(7,4);
t507 = t278 * t548;
t346 = qJD(4) + qJD(5);
t296 = qJD(3) + t346;
t251 = qJD(2) + t296;
t571 = pkin(3) * t251;
t521 = rSges(7,1) * t571;
t553 = qJD(6) * t251;
t569 = Icges(7,6) * qJD(6);
t577 = pkin(2) * t228;
t153 = sin(t198);
t611 = t153 * t390;
t343 = rSges(7,1) * t623;
t378 = cos(qJ(6));
t570 = Icges(7,5) * qJD(6);
t671 = -0.2e1 * rSges(7,2);
t626 = t378 * (m(7) * (t571 * t671 + t343) - t570);
t669 = m(7) / 0.4e1;
t233 = qJ(6) + t289;
t175 = sin(t233);
t235 = -qJ(6) + t289;
t177 = sin(t235);
t189 = cos(t235);
t286 = qJ(3) + t326;
t211 = sin(t286);
t287 = qJ(3) + t327;
t212 = sin(t287);
t220 = cos(t287);
t293 = 0.2e1 * qJD(2) ^ 2 + t390;
t527 = pkin(4) * t650;
t648 = pkin(4) * t293;
t528 = t648 / 0.2e1;
t547 = pkin(4) * t657;
t219 = cos(t286);
t697 = rSges(7,2) * t219;
t187 = cos(t233);
t699 = rSges(7,2) * t187;
t673 = t745 * t528 + t720 * t527 - (t390 * t189 + t220 * t293) * t547 / 0.4e1 + (t697 + (t211 + t212) * rSges(7,1)) * t648 * t669 + (t525 * t699 + (t175 + t177) * t482) * pkin(4);
t268 = cos(t326);
t696 = rSges(7,2) * t268;
t174 = sin(t232);
t176 = sin(t234);
t714 = pkin(2) * (t176 + t174);
t690 = t246 * sin(t288);
t99 = t574 * t690;
t726 = t99 / 0.2e1 + t100 / 0.2e1 - t228 * t484 / 0.4e1 - t191 * t507 + t673 - t164 * t352 * t553 - t409 * t611 / 0.2e1 - t444 / 0.4e1 + t445 / 0.4e1 + t482 * t714 + qJD(6) * (-t372 * ((t342 + 0.2e1 * t521) * m(7) - t569) + t626) + t739 + (t747 / 0.2e1 + (t696 + (t261 + t260) * rSges(7,1)) * t669) * t577;
t519 = -t638 / 0.2e1;
t725 = (rSges(7,1) * t260 + t696) * t294 * t519;
t308 = sin(t358);
t316 = cos(t358);
t349 = qJD(2) - qJD(6);
t668 = pkin(4) * m(7);
t489 = t349 * t668 / 0.2e1;
t658 = m(7) * rSges(7,1);
t534 = -t658 / 0.2e1;
t724 = pkin(4) * t349 * t308 * t534 + rSges(7,2) * t316 * t489;
t723 = -(t211 * rSges(7,1) + t697) * t247 * t668 / 0.2e1;
t380 = cos(qJ(4));
t660 = m(5) * rSges(5,2);
t537 = pkin(2) * t660;
t374 = sin(qJ(4));
t604 = t196 * t374;
t722 = (pkin(2) * t604 + t380 * t537) * qJD(4);
t721 = t747 * pkin(2) * t346;
t271 = cos(t329);
t606 = t196 * sin(t329);
t719 = t271 * t660 + t606;
t313 = sin(t363);
t321 = cos(t363);
t662 = m(4) * rSges(4,2);
t718 = t143 * t313 + t321 * t662;
t157 = -t372 * rSges(7,1) - t378 * rSges(7,2);
t717 = t296 * t745;
t320 = cos(t362);
t605 = t196 * sin(t362);
t716 = t347 * (-t320 * t660 - t605);
t307 = sin(t357);
t315 = cos(t357);
t348 = qJD(2) + qJD(6);
t715 = t348 * (rSges(7,1) * t307 + rSges(7,2) * t315);
t375 = sin(qJ(3));
t381 = cos(qJ(3));
t713 = qJD(3) * (-t143 * t375 - t381 * t662);
t141 = cos(t159);
t142 = cos(t160);
t190 = cos(t236);
t215 = sin(t290);
t551 = pkin(3) * qJD(1);
t301 = pkin(1) * t551;
t336 = t399 + t400;
t354 = cos(t393);
t351 = pkin(1) * qJD(1);
t514 = rSges(6,1) * t351;
t568 = Icges(7,3) * qJD(6);
t633 = m(7) * qJD(1);
t712 = ((qJD(6) * t336 + t301) * t633 + qJD(1) * (m(6) * t514 + t568)) * t215 + t190 * t507 + (-0.2e1 * t354 * t553 + (t141 / 0.4e1 - t142 / 0.4e1) * t390) * t278;
t709 = 2 * Icges(7,4);
t367 = qJD(5) / 0.2e1;
t224 = t367 + t299;
t707 = -0.2e1 * t224;
t479 = qJD(1) / 0.2e1;
t208 = sin(t283);
t702 = rSges(7,1) * t208;
t265 = sin(t331);
t701 = rSges(7,1) * t265;
t216 = cos(t283);
t698 = rSges(7,2) * t216;
t273 = cos(t331);
t695 = rSges(7,2) * t273;
t317 = cos(t359);
t693 = rSges(7,2) * t317;
t692 = t196 * sin(t328);
t691 = t246 * sin(t285);
t687 = t254 * cos(t285);
t613 = t728 * cos(t231);
t368 = qJD(4) / 0.2e1;
t369 = qJD(3) / 0.2e1;
t641 = qJD(2) / 0.2e1;
t679 = (t641 + t369 + t368 + t367) * pkin(3);
t253 = pkin(2) * t381 + pkin(4);
t376 = sin(qJ(2));
t382 = cos(qJ(2));
t572 = pkin(3) * t223;
t575 = pkin(2) * t375;
t678 = (pkin(5) * (t382 * (-t374 * t375 + t380 * t381) - t376 * (t374 * t381 + t375 * t380)) + t572 + t253 * t382 - t376 * t575 + pkin(1)) * qJD(1);
t535 = t660 / 0.2e1;
t651 = -t390 / 0.2e1;
t675 = t527 * t606 + t528 * t605 - t651 * t619 + (t271 * t390 + t293 * t320) * pkin(4) * t535;
t136 = qJD(6) * t223 + qJD(1);
t383 = cos(qJ(1));
t558 = qJD(6) * t215;
t377 = sin(qJ(1));
t517 = t299 * t377;
t72 = qJD(5) * t377 + t517;
t26 = -t383 * t558 + t72;
t73 = t251 * t383;
t27 = -t377 * t558 - t73;
t602 = t215 * t383;
t603 = t215 * t377;
t583 = t378 * t383;
t585 = t377 * t372;
t74 = -t223 * t585 + t583;
t584 = t377 * t378;
t587 = t372 * t383;
t75 = t223 * t584 + t587;
t76 = -t223 * t587 - t584;
t77 = t223 * t583 - t585;
t408 = (Icges(7,5) * t223 + (t372 * t710 + t729 * t378) * t215) * t136 + (-Icges(7,5) * t602 + t76 * t709 + t729 * t77) * t26 + (-Icges(7,5) * t603 + t74 * t709 + t729 * t75) * t27;
t646 = pkin(5) * t390;
t102 = t646 * t691;
t104 = t646 * t687;
t366 = -qJD(6) / 0.2e1;
t163 = t366 + t224;
t209 = sin(t284);
t275 = cos(t333);
t318 = cos(t360);
t546 = pkin(5) * t657;
t505 = t318 * t546;
t463 = t345 * t505;
t193 = cos(t239);
t506 = t390 * t546;
t464 = t193 * t506;
t192 = cos(t238);
t465 = t192 * t506;
t483 = rSges(7,1) * t645 / 0.2e1;
t520 = m(7) * t548;
t672 = t390 * t732 + t712 + t104 + t275 * t578 * t660 + t163 * t463 + t102 - t464 / 0.2e1 + t465 / 0.2e1 + t746 * t707 + t483 * t743 + (t698 + (t209 + t208) * rSges(7,1)) * pkin(5) * t520;
t670 = -0.2e1 * qJD(6);
t664 = pkin(5) * ((t374 * t382 + t376 * t380) * t381 + t375 * (-t374 * t376 + t380 * t382));
t663 = m(3) * rSges(3,1);
t661 = m(4) * rSges(4,3);
t386 = m(5) * rSges(5,3);
t385 = m(6) * rSges(6,3);
t655 = rSges(3,2) * m(3);
t654 = -t202 / 0.2e1;
t653 = -t203 / 0.2e1;
t652 = t293 / 0.2e1;
t644 = rSges(6,1) * rSges(6,3);
t636 = pkin(4) * qJD(1);
t632 = rSges(7,3) * t251;
t629 = t409 * t153;
t356 = m(4) + t388;
t620 = ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + (m(7) + t356) * pkin(4) ^ 2) * sin(t397);
t341 = qJD(2) + t369;
t252 = t368 + t341;
t170 = t367 + t252;
t365 = qJD(6) / 0.2e1;
t144 = t365 + t170;
t616 = t144 * t247;
t291 = t368 + t350;
t205 = t367 + t291;
t151 = t365 + t205;
t612 = t151 * t294;
t162 = t365 + t224;
t609 = t162 * t344;
t608 = t164 * t251;
t599 = t248 * t212;
t598 = t248 * t220;
t595 = t261 * t295;
t594 = t275 * t299;
t592 = (rSges(5,2) * t384 - Icges(5,4)) * cos(t229);
t591 = (rSges(4,1) * t662 - Icges(4,4)) * cos(t277);
t590 = (rSges(3,1) * t655 - Icges(3,4)) * cos(t397);
t322 = cos(t364);
t589 = t322 * t350;
t310 = sin(t360);
t588 = t345 * t310;
t586 = t377 * t157;
t582 = t383 * t157;
t550 = pkin(3) * qJD(6);
t581 = t550 / 0.2e1 + t351;
t580 = t385 + t386;
t573 = pkin(3) * t215;
t194 = cos(t240);
t565 = qJD(1) * t194;
t195 = cos(t241);
t564 = qJD(1) * t195;
t563 = qJD(1) * t251;
t562 = qJD(1) * t278;
t559 = qJD(2) * t382;
t557 = qJD(6) * t683;
t556 = qJD(6) * t682;
t555 = qJD(6) * t681;
t554 = qJD(6) * t245;
t552 = pkin(2) * qJD(1);
t544 = pkin(4) * t633;
t543 = pkin(5) * t633;
t538 = 0.16e2 * t633;
t530 = 0.2e1 * t351;
t524 = -rSges(6,2) * rSges(6,3) / 0.2e1;
t522 = m(7) * t552;
t515 = rSges(7,1) * t550;
t324 = rSges(7,1) * t351;
t323 = rSges(7,2) * t351;
t513 = rSges(7,2) * t550;
t512 = t611 / 0.4e1;
t509 = t251 * t278 / 0.2e1;
t298 = qJD(3) + t349;
t297 = qJD(3) + t348;
t495 = t144 * t544;
t145 = t366 + t170;
t494 = t145 * t544;
t493 = t162 * t543;
t492 = t163 * t543;
t487 = -t344 * t667 / 0.2e1;
t475 = rSges(7,2) * t633 / 0.2e1;
t474 = t151 * t522;
t152 = t366 + t205;
t473 = t152 * t522;
t462 = pkin(5) * t588 * t658;
t461 = pkin(4) * t599 * t658;
t460 = t547 * t598;
t448 = pkin(5) * t475;
t249 = qJD(4) + t297;
t250 = qJD(4) + t298;
t447 = t249 * t250 * t667 / 0.2e1;
t443 = rSges(7,1) * t595 * t638;
t442 = t295 * t484;
t158 = rSges(7,1) * t378 - rSges(7,2) * t372;
t431 = pkin(2) * t475;
t430 = t297 * t298 * t638 / 0.2e1;
t429 = -rSges(7,3) * t215 + t158 * t223;
t425 = pkin(5) * t209 * t250 * t534;
t309 = sin(t359);
t420 = rSges(7,1) * t309 * t487;
t419 = t487 * t693;
t266 = sin(t332);
t417 = rSges(7,1) * t266 * t298 * t519;
t413 = -pkin(2) * (t375 * t382 + t376 * t381) * qJD(3) - (t253 * t376 + t382 * t575) * qJD(2);
t412 = t136 * (-Icges(7,5) * t372 - Icges(7,6) * t378) * t215 + t26 * (Icges(7,5) * t76 - Icges(7,6) * t77) + t27 * (Icges(7,5) * t74 - Icges(7,6) * t75);
t410 = t215 * t412;
t407 = (Icges(7,6) * t602 + t77 * t710 + t729 * t76) * t26 + (Icges(7,6) * t603 + t75 * t710 + t729 * t74) * t27 + (-Icges(7,6) * t223 + (-t729 * t372 + t378 * t710) * t215) * t136;
t217 = cos(t284);
t406 = ((-t217 * t548 - t317 * t609) * rSges(7,2) + (-t163 * t588 - t309 * t609) * rSges(7,1)) * pkin(5);
t274 = cos(t332);
t404 = t406 + ((-t268 * t612 - t274 * t548) * rSges(7,2) + (-t152 * t595 - t260 * t612) * rSges(7,1)) * pkin(2);
t270 = cos(t328);
t137 = t390 * t270 * t537;
t255 = t323 / 0.2e1;
t256 = -t323 / 0.2e1;
t302 = Icges(7,6) * t548;
t57 = t613 * t650;
t92 = t574 * t692;
t403 = t152 * t442 + t99 + t100 + t92 + ((t256 - t554) * m(7) - t570) * t564 + ((t255 - t555) * m(7) - t570) * t565 + ((t324 - 0.2e1 * t556) * t538 + 0.32e2 * t302) * t182 / 0.32e2 + ((t324 - 0.2e1 * t557) * t538 - 0.32e2 * t302) * t183 / 0.32e2 + (-t191 * t562 - t352 * t608 + t626) * qJD(6) + ((t266 + t265) * rSges(7,1) + t695) * pkin(2) * t520 - t57 + t739 + t322 * t578 * t662 - 0.2e1 * t721 * t205 - 0.2e1 * t722 * t291 + (t591 + t592 + t731) * t390 + ((t521 + t342 / 0.2e1) * m(7) - t569 / 0.2e1) * t372 * t670 - t651 * t730 + t672 - t444 / 0.2e1 + t445 / 0.2e1 + t483 * t714 + t137;
t292 = -0.32e2 * Icges(7,5) * t548;
t258 = -t324 / 0.2e1;
t257 = t324 / 0.2e1;
t197 = pkin(4) * t356 + t663;
t44 = t157 * t215;
t34 = t581 - t679;
t33 = t581 + t679;
t25 = t223 * rSges(7,3) + t158 * t215;
t20 = -t158 * t377 + t223 * t582;
t19 = t158 * t383 + t223 * t586;
t6 = t429 * t383 + t586;
t5 = t429 * t377 - t582;
t1 = [(-(t530 * t656 + (rSges(6,2) * t530 + qJD(5) * t644) * m(6) - Icges(6,5) * qJD(5) + t299 * (m(6) * t644 - Icges(6,5))) * t223 - 0.2e1 * ((t508 * qJD(6) + t301) * m(7) + (qJD(5) * t524 + t514) * m(6) + Icges(6,6) * t367 + t568 / 0.2e1 + t299 * (m(6) * t524 + Icges(6,6) / 0.2e1)) * t215) * t251 + t295 * t269 * t431 - 0.2e1 * ((t270 * t660 + t692) * t291 + (t686 + t690) * t205) * t552 - 0.2e1 * (t720 * t170 + t719 * t252 + t718 * t341) * t636 + ((t217 + t216) * t447 - t186 * t474 + t188 * t473 + t189 * t494 - t192 * t493 + t193 * t492 + (t274 + t273) * t430) * rSges(7,2) - (t530 * t662 + t350 * (t580 * pkin(2) + rSges(4,1) * t661 - Icges(4,5))) * t589 - (((t580 + t661) * pkin(4) + rSges(3,3) * t663 - Icges(3,5)) * qJD(2) + t530 * t655) * t559 + (t203 * t142 / 0.2e1 + t141 * t654) * t562 - 0.2e1 * t737 * (-rSges(4,2) * t661 / 0.2e1 + Icges(4,6) / 0.2e1) * t314 - 0.2e1 * t738 * (-rSges(5,2) * t386 / 0.2e1 + Icges(5,6) / 0.2e1) * t267 + (t178 * t608 / 0.4e1 + t190 * t509) * (0.2e1 * qJD(6) + t251) + (-t179 * t608 / 0.4e1 + t191 * t509) * (t670 + t251) + t724 * t348 + t297 * t417 + (pkin(4) * t713 + t419 + t420 - t721 - t722 + t723 + t725 + (t687 + t691) * pkin(5) * t707 - t462 / 0.2e1 - t461 / 0.2e1 - t443 / 0.2e1 + (t613 + t615) * (t366 + t251) + t744 * (t365 + t251) - (t202 * t139 + t203 * t140) * t164 / 0.4e1 + (-0.2e1 * t590 - t620) * qJD(2) + (-0.2e1 * t731 - 0.2e1 * t591 - t730) * t350 + (-0.2e1 * t732 - 0.2e1 * t592 - t619) * t299 - t746) * qJD(1) + (t716 - t717) * t636 - 0.2e1 * qJD(2) * ((-rSges(3,3) * t655 / 0.2e1 + Icges(3,6) / 0.2e1) * qJD(2) + (t197 + t668) * t351) * t376 + t202 * ((t203 * t371 + t34 * t671) * m(7) - Icges(7,5) * t203) * t194 / 0.2e1 + ((t202 * t371 + t33 * t671) * m(7) - Icges(7,5) * t202) * t195 * t653 + t352 * t457 + t157 * m(7) * qJD(6) * (pkin(1) * qJD(6) + t551) + t354 * t507 - (Icges(7,6) * t653 + t203 * t736) * t734 - (Icges(7,6) * t654 + t202 * t736) * t733 + t345 * t318 * t448 - (t530 * t660 + t299 * (pkin(5) * t385 + rSges(5,3) * t384 - Icges(5,5))) * t594 + (-0.2e1 * t610 + t629) * t563 + t489 * t715 + t430 * t701 + t447 * t702 + t460 * t479 + t249 * t425 + (-t174 * t474 - t175 * t495 - t176 * t473 - t177 * t494 - t180 * t493 - t181 * t492 + (-t33 * t733 - t34 * t734) * m(7)) * rSges(7,1) - t495 * t699; t404 * m(7) + t403 + (-t629 / 0.2e1 + t619 / 0.2e1 + t620 / 0.2e1 + t590 + (t197 * t376 + t382 * t655) * pkin(1)) * t390 + ((t718 + t719 + t720) * t390 + (t376 * t578 + (-t219 * t616 + t145 * t598 + (-t189 / 0.2e1 + t187 / 0.2e1) * t390 + (t315 - t316) * t548) * rSges(7,2) + (-t211 * t616 - t145 * t599 + (t177 / 0.2e1 + t175 / 0.2e1) * t390 + (t307 + t308) * t548) * rSges(7,1)) * m(7)) * pkin(4) + 0.2e1 * (-t170 * t717 + t252 * t716 + t341 * t713) * pkin(4); (-0.2e1 * Icges(6,1) - Icges(7,1) + 0.2e1 * Icges(6,2) - Icges(7,2) + 0.2e1 * Icges(7,3) + 0.2e1 * t643) * t512 + t673 + t403 + ((t313 * t650 + t375 * t652) * t143 + (t321 * t650 + t381 * t652) * t662) * pkin(4) + ((0.2e1 * t391 - 0.2e1 * t398 + t336) * t512 + t404) * m(7) + t675; m(7) * t406 + t92 / 0.2e1 + (t604 / 0.2e1 + t380 * t535) * t577 + (-(t258 + t556) * t633 + t302) * t182 + ((t257 - t557) * t633 - t302) * t183 + (t592 - t613 / 0.2e1) * t390 + t672 + t137 / 0.2e1 + t675 + t726 + ((t323 - 0.2e1 * t555) * t538 + t292) * t194 / 0.32e2 + (-0.16e2 * (t323 + 0.2e1 * t554) * t633 + t292) * t195 / 0.32e2; t712 + t482 * t743 + ((t255 + t343 + t513) * m(7) - t570) * t565 + ((t256 + t343 - t513) * m(7) - t570) * t564 + t104 / 0.2e1 - t57 + (-(t258 + t342 - t515) * t633 + t302) * t182 + ((t257 + t342 + t515) * t633 - t302) * t183 + t102 / 0.2e1 - t464 / 0.4e1 + t465 / 0.4e1 + t726 + (-t505 / 0.4e1 + (t748 / 0.2e1 + (t693 + (t309 + t310) * rSges(7,1)) * t669) * pkin(5)) * (t390 + 0.2e1 * t738); t298 * t274 * t431 + t299 * t462 / 0.2e1 - t299 * t463 / 0.2e1 - m(7) * ((-t73 * t573 - t136 * t5 + t27 * t25 - t377 * t678 + (-t299 * t664 + t413) * t383) * (-t136 * t19 + t27 * t44) + (t136 * t6 - t26 * t25 + t413 * t377 + t383 * t678 - t517 * t664 - t72 * t573) * (t136 * t20 - t26 * t44) + (pkin(2) * t589 + pkin(4) * t559 + pkin(5) * t594 + t26 * t5 - t27 * t6 + (t377 * t72 + t383 * t73) * t572) * (t19 * t26 - t20 * t27)) - t136 * (t412 * t223 + (-t408 * t372 + t407 * t378) * t215) / 0.2e1 - t26 * (-t383 * t410 + t407 * t77 + t408 * t76) / 0.2e1 - t27 * (-t377 * t410 + t407 * t75 + t408 * t74) / 0.2e1 + t299 * t419 + t299 * t420 - (m(7) * t336 + Icges(7,3)) * t215 * t563 + t461 * t641 + t250 * t217 * t448 + (((-rSges(7,2) * t632 - t324) * m(7) + Icges(7,6) * t251) * t372 + ((rSges(7,1) * t632 - t323) * m(7) - Icges(7,5) * t251) * t378) * qJD(6) + (-t460 / 0.2e1 + t723) * qJD(2) + (-t442 / 0.2e1 + t443 / 0.2e1 + t725) * t350 + ((t741 * t182 + t742 * t194) * t202 + (t740 * t183 + t728 * t195) * t203) * t479 - (pkin(4) * t715 + (t698 + t702) * pkin(5) * t249 + (t695 + t701) * pkin(2) * t297) * t633 / 0.2e1 + (t417 + t425 + t724) * qJD(1);];
tauc = t1(:);
