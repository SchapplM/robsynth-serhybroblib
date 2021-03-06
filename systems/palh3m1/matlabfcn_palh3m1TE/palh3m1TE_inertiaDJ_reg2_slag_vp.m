% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh3m1TE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_inertiaDJ_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 00:38:25
% EndTime: 2020-04-18 01:56:14
% DurationCPUTime: 1090.01s
% Computational Cost: add. (33087348->599), mult. (51619444->1335), div. (1944862->28), fcn. (32258958->22), ass. (0->523)
t272 = sin(qJ(2));
t276 = cos(qJ(2));
t593 = sin(pkin(16));
t594 = cos(pkin(16));
t249 = t272 * t594 + t276 * t593;
t644 = -0.4e1 * t249;
t284 = pkin(5) ^ 2;
t289 = pkin(1) ^ 2;
t493 = 2 * pkin(1);
t343 = -t272 * t593 + t276 * t594;
t583 = pkin(5) * t343;
t504 = t493 * t583 + t289;
t239 = t284 + t504;
t236 = 0.1e1 / t239;
t237 = 0.1e1 / t239 ^ 2;
t238 = t236 * t237;
t243 = t249 * qJD(2);
t515 = t243 * t249;
t466 = t289 * t515;
t420 = t284 * t466;
t395 = t238 * t420;
t643 = 0.4e1 * t395;
t270 = sin(qJ(4));
t262 = t270 ^ 2;
t274 = cos(qJ(4));
t263 = t274 ^ 2;
t503 = t262 - t263;
t429 = qJD(4) * t503;
t381 = 0.2e1 * t429;
t623 = t238 * t284 * t289 * t644;
t288 = 0.1e1 / pkin(2);
t500 = pkin(2) ^ 2 - pkin(6) ^ 2;
t234 = t239 + t500;
t240 = pkin(1) + t583;
t618 = pkin(5) + pkin(6);
t619 = pkin(5) - pkin(6);
t231 = (pkin(2) + t618) * (-pkin(2) + t619) + t504;
t232 = (-pkin(2) + t618) * (pkin(2) + t619) + t504;
t519 = t232 * t231;
t291 = sqrt(-t519);
t221 = pkin(5) * t249 * t234 + t240 * t291;
t275 = cos(qJ(3));
t522 = t221 * t275;
t511 = t249 * t291;
t220 = -pkin(5) * t511 + t234 * t240;
t271 = sin(qJ(3));
t527 = t220 * t271;
t357 = t522 / 0.2e1 + t527 / 0.2e1;
t488 = pkin(1) * pkin(5) * t237;
t442 = t243 * t488;
t370 = (t231 + t232) * pkin(5) * t493;
t223 = t243 * t370;
t244 = t343 * qJD(2);
t592 = pkin(1) * t284;
t426 = t515 * t592;
t226 = 0.1e1 / t291;
t608 = t226 / 0.2e1;
t455 = t240 * t608;
t514 = t243 * t291;
t189 = t223 * t455 - 0.2e1 * t426 + (t244 * t234 - t514) * pkin(5);
t540 = t189 * t271;
t607 = -t249 / 0.2e1;
t637 = t223 * t226;
t415 = t607 * t637;
t210 = pkin(5) * t415;
t580 = t240 * pkin(1);
t447 = t234 + 0.2e1 * t580;
t513 = t244 * t291;
t188 = t210 + (-t243 * t447 - t513) * pkin(5);
t541 = t188 * t275;
t523 = t221 * t271;
t526 = t220 * t275;
t630 = t523 - t526;
t142 = ((-t540 / 0.2e1 + t541 / 0.2e1 - t357 * qJD(3)) * t236 - t630 * t442) * t288;
t384 = t522 + t527;
t350 = t243 * t384;
t358 = -t523 / 0.2e1 + t526 / 0.2e1;
t539 = t189 * t275;
t542 = t188 * t271;
t143 = (t350 * t488 + (t539 / 0.2e1 + t542 / 0.2e1 + t358 * qJD(3)) * t236) * t288;
t261 = pkin(18) + pkin(19);
t259 = sin(t261);
t260 = cos(t261);
t137 = t142 * t260 - t143 * t259;
t285 = pkin(4) ^ 2;
t516 = t236 * t288;
t206 = t358 * t516;
t207 = t357 * t516;
t174 = -t206 * t260 + t207 * t259;
t585 = pkin(4) * t174;
t642 = 2 * pkin(3);
t506 = t585 * t642 + t285;
t617 = (-pkin(8) - pkin(10));
t162 = ((pkin(3) - t617) * (pkin(3) + t617)) + t506;
t616 = (pkin(10) - pkin(8));
t163 = ((pkin(3) - t616) * (pkin(3) + t616)) + t506;
t550 = t163 * t162;
t290 = sqrt(-t550);
t555 = t137 * t290;
t257 = -pkin(1) * t275 + pkin(4);
t286 = pkin(3) ^ 2;
t169 = t286 + t506;
t166 = 0.1e1 / t169;
t279 = 0.1e1 / pkin(10);
t167 = 0.1e1 / t169 ^ 2;
t590 = pkin(3) * t167;
t487 = pkin(4) * t590;
t501 = pkin(8) ^ 2 - pkin(10) ^ 2;
t165 = t169 - t501;
t171 = pkin(3) * t174 + pkin(4);
t390 = -t206 * t259 - t207 * t260;
t589 = pkin(3) * t390;
t129 = t165 * t589 + t171 * t290;
t264 = sin(pkin(17));
t559 = t129 * t264;
t636 = t390 * t290;
t127 = -pkin(3) * t636 + t165 * t171;
t265 = cos(pkin(17));
t560 = t127 * t265;
t333 = (t559 - t560) * t487;
t605 = -t265 / 0.2e1;
t606 = t264 / 0.2e1;
t145 = 0.1e1 / t290;
t609 = -t390 / 0.2e1;
t458 = t145 * t609;
t391 = -t142 * t259 - t143 * t260;
t635 = t391 * t290;
t492 = 0.2e1 * pkin(4);
t369 = pkin(3) * (t162 + t163) * t492;
t99 = t137 * t369;
t353 = t458 * t99 - t635;
t445 = -0.2e1 * pkin(4) * t171 - t165;
t64 = (t137 * t445 + t353) * pkin(3);
t586 = pkin(4) * t390;
t448 = -0.2e1 * t286 * t586;
t611 = t145 / 0.2e1;
t459 = t171 * t611;
t66 = t99 * t459 + t137 * t448 + (t165 * t391 - t555) * pkin(3);
t43 = ((t64 * t605 + t66 * t606) * t166 + t137 * t333) * t279;
t558 = t129 * t265;
t561 = t127 * t264;
t332 = (t558 + t561) * t487;
t604 = t265 / 0.2e1;
t44 = ((t66 * t604 + t64 * t606) * t166 + t137 * t332) * t279;
t548 = t166 * t279;
t95 = (-t560 / 0.2e1 + t559 / 0.2e1) * t548;
t567 = t271 * t95;
t96 = (t558 / 0.2e1 + t561 / 0.2e1) * t548;
t29 = pkin(1) * (qJD(3) * (t275 * t96 + t567) + t271 * t44) + t257 * t43;
t92 = 0.1e1 / t95;
t224 = t249 * t370;
t621 = -0.2e1 * t249 ^ 2;
t194 = t224 * t455 + t592 * t621 + (t234 * t343 - t511) * pkin(5);
t535 = t194 * t275;
t456 = t224 * t608;
t413 = t456 + t234;
t512 = t343 * t291;
t192 = (-t512 + (-t413 - 0.2e1 * t580) * t249) * pkin(5);
t538 = t192 * t271;
t359 = t535 / 0.2e1 + t538 / 0.2e1;
t441 = t249 * t488;
t160 = (t236 * t359 + t384 * t441) * t288;
t536 = t194 * t271;
t537 = t192 * t275;
t360 = -t536 / 0.2e1 + t537 / 0.2e1;
t161 = (-t236 * t360 + t441 * t630) * t288;
t141 = -t160 * t259 - t161 * t260;
t106 = t141 * t369;
t140 = -t160 * t260 + t161 * t259;
t348 = t106 * t458 - t140 * t290;
t70 = (t141 * t445 + t348) * pkin(3);
t554 = t141 * t290;
t72 = t106 * t459 + t141 * t448 + (t140 * t165 - t554) * pkin(3);
t322 = (t70 * t605 + t72 * t606) * t166 + t141 * t333;
t566 = t279 * t96;
t313 = t322 * t566;
t317 = t279 * ((t72 * t604 + t70 * t606) * t166 + t141 * t332);
t93 = 0.1e1 / t95 ^ 2;
t633 = t93 * t313 - t92 * t317;
t94 = t96 ^ 2;
t80 = t93 * t94 + 0.1e1;
t78 = 0.1e1 / t80;
t21 = -t633 * t78 + 0.1e1;
t77 = pkin(1) * t271 * t96 + t257 * t95;
t19 = -pkin(9) * t21 - t77;
t489 = 0.2e1 * t19;
t486 = 0.2e1 * t270;
t502 = t262 + t263;
t417 = 0.2e1 * t502;
t565 = t43 * t92 * t93;
t574 = t93 * t96;
t641 = 0.2e1 * (t44 * t574 - t94 * t565) / t80 ^ 2;
t164 = t169 + t501;
t170 = -pkin(3) - t585;
t126 = -pkin(4) * t636 - t164 * t170;
t123 = 0.1e1 / t126 ^ 2;
t128 = t164 * t586 - t170 * t290;
t125 = t128 ^ 2;
t103 = t123 * t125 + 0.1e1;
t562 = t123 * t128;
t122 = 0.1e1 / t126;
t446 = t170 * t642 - t164;
t63 = (t137 * t446 + t353) * pkin(4);
t573 = t122 * t123 * t63;
t450 = -0.2e1 * t285 * t589;
t460 = -t145 * t170 / 0.2e1;
t638 = t164 * t391;
t65 = t99 * t460 + t137 * t450 + (-t555 + t638) * pkin(4);
t640 = 0.2e1 / t103 ^ 2 * (-t125 * t573 + t65 * t562);
t575 = t78 * t93;
t639 = t43 * t575;
t138 = t390 * t369;
t347 = t138 * t458 - t174 * t290;
t87 = (t390 * t445 + t347) * pkin(3);
t89 = t138 * t459 + t390 * t448 + (t165 * t174 - t636) * pkin(3);
t321 = (t87 * t605 + t89 * t606) * t166 + t390 * t333;
t312 = t321 * t566;
t316 = t279 * ((t89 * t604 + t87 * t606) * t166 + t390 * t332);
t634 = t93 * t312 - t92 * t316;
t168 = t166 * t167;
t556 = t137 * t286;
t469 = t285 * t556;
t396 = t168 * t279 * t469;
t355 = 0.4e1 * t264 * t396;
t440 = t279 * t487;
t410 = t264 * t440;
t632 = t390 * t355 + t391 * t410;
t218 = t244 * t370 - 0.8e1 * t420;
t467 = 0.1e1 / t519 * t224 * t637;
t344 = t467 / 0.4e1 + t218 * t608;
t383 = -0.2e1 * t243 * t343 + t244 * t644;
t153 = t210 + t344 * t240 + t383 * t592 + (-t243 * t413 - t513) * pkin(5);
t318 = (-t343 * t223 / 0.2e1 - t244 * t224 / 0.2e1 + t218 * t607) * t226 + t514 - t249 * t467 / 0.4e1;
t155 = 0.4e1 * t426 + (-t244 * t447 + t318) * pkin(5);
t598 = t271 / 0.2e1;
t107 = (t630 * t643 + (-t155 * t275 / 0.2e1 + t153 * t598 + t359 * qJD(3)) * t236 + (t630 * t244 - (-t536 + t537) * t243 + (qJD(3) * t384 + t540 - t541) * t249) * t488) * t288;
t108 = (-t350 * t623 + (t153 * t275 / 0.2e1 + t155 * t598 + t360 * qJD(3)) * t236 + (t384 * t244 - (-t535 - t538) * t243 + (-qJD(3) * t630 + t539 + t542) * t249) * t488) * t288;
t84 = -t107 * t260 - t108 * t259;
t631 = t141 * t355 + t84 * t410;
t629 = qJD(2) + qJD(3);
t628 = 0.2e1 * t565 * t78;
t101 = 0.1e1 / t103;
t416 = pkin(3) * pkin(4) * pkin(8) * t101 * t137;
t581 = pkin(8) * t169;
t480 = t101 * t581;
t428 = t123 * t480;
t479 = t122 * t581;
t627 = -0.2e1 * t122 * t416 - t63 * t428 - t479 * t640;
t478 = t128 * t581;
t427 = t123 * t478;
t624 = 0.2e1 * t101 * t478 * t573 + 0.2e1 * t416 * t562 + t427 * t640 - t65 * t428;
t268 = cos(pkin(19));
t524 = t221 * t268;
t266 = sin(pkin(19));
t529 = t220 * t266;
t386 = t524 + t529;
t351 = t243 * t386;
t600 = t268 / 0.2e1;
t603 = t266 / 0.2e1;
t148 = ((t188 * t603 + t189 * t600) * t236 + t351 * t488) * t288;
t525 = t221 * t266;
t528 = t220 * t268;
t387 = -t525 + t528;
t352 = t243 * t387;
t601 = -t268 / 0.2e1;
t149 = ((t188 * t601 + t189 * t603) * t236 - t352 * t488) * t288;
t156 = ((t192 * t603 + t194 * t600) * t236 + t386 * t441) * t288;
t157 = ((t192 * t601 + t194 * t603) * t236 - t387 * t441) * t288;
t203 = (-t528 / 0.2e1 + t525 / 0.2e1) * t516;
t196 = 0.1e1 / t203 ^ 2;
t204 = (t524 / 0.2e1 + t529 / 0.2e1) * t516;
t198 = t204 ^ 2;
t181 = t196 * t198 + 0.1e1;
t179 = 0.1e1 / t181;
t195 = 0.1e1 / t203;
t335 = t189 * t249 + t194 * t243 + t221 * t244;
t336 = t188 * t249 + t192 * t243 + t220 * t244;
t534 = t196 * t204;
t468 = t179 * t534;
t543 = t179 * t196;
t544 = t179 * t195;
t553 = t149 * t195 * t196;
t564 = 0.2e1 * (t148 * t534 - t198 * t553) / t181 ^ 2;
t51 = ((-t351 * t623 + (t153 * t600 + t155 * t603) * t236 + (t266 * t336 + t268 * t335) * t488) * t544 - (t352 * t623 + (t153 * t603 + t155 * t601) * t236 + (t266 * t335 - t268 * t336) * t488) * t468) * t288 + (-t149 * t543 - t195 * t564) * t156 + ((0.2e1 * t179 * t553 + t196 * t564) * t204 - t148 * t543) * t157;
t36 = -t634 * t78 + 0.1e1;
t622 = 0.2e1 * t36;
t283 = 0.1e1 / pkin(6);
t517 = t236 * t283;
t233 = t239 - t500;
t241 = -pkin(1) * t343 - pkin(5);
t222 = pkin(1) * t249 * t233 - t241 * t291;
t273 = sin(pkin(15));
t521 = t222 * t273;
t219 = -pkin(1) * t511 - t233 * t241;
t277 = cos(pkin(15));
t530 = t219 * t277;
t205 = (t530 / 0.2e1 + t521 / 0.2e1) * t517;
t199 = 0.1e1 / t205;
t200 = 0.1e1 / t205 ^ 2;
t620 = -t99 / 0.2e1;
t613 = -t106 / 0.2e1;
t324 = (t137 * t613 + t141 * t620) * t145 - t84 * t290;
t572 = t145 / t550 * t99;
t465 = t572 / 0.4e1;
t431 = -0.8e1 * t469;
t54 = t141 * t431 + t369 * t84;
t346 = t106 * t465 + t54 * t611;
t394 = -t141 * t391 - t390 * t84;
t451 = t286 * t492;
t457 = t548 / 0.2e1;
t83 = t107 * t259 - t108 * t260;
t307 = (t346 * t171 + (-t137 * t140 + t394) * t451 + (t83 * t165 + t324) * pkin(3)) * t457;
t419 = -t390 * t572 / 0.4e1;
t320 = (t140 * t620 + t391 * t613 + t54 * t609) * t145 - t83 * t290 + t106 * t419;
t449 = 0.4e1 * pkin(4) * t556;
t309 = (t141 * t449 + (t445 * t84 + t320) * pkin(3)) * t548;
t364 = t265 * t396;
t341 = 0.4e1 * t129 * t364;
t342 = -0.4e1 * t127 * t364;
t409 = t265 * t440;
t373 = t141 * t409;
t374 = t141 * t410;
t375 = t137 * t409;
t376 = t137 * t410;
t379 = t84 * t409;
t435 = t279 * t44 * t575;
t483 = t78 * t574;
t576 = t78 * t92;
t5 = (t127 * t631 + t129 * t379 + t141 * t341 + t265 * t307 + t309 * t606 + t66 * t373 + t64 * t374 + t72 * t375 + t70 * t376) * t576 - t317 * t639 - (-t127 * t379 + t129 * t631 + t141 * t342 + t264 * t307 + t309 * t605 - t64 * t373 + t66 * t374 - t70 * t375 + t72 * t376) * t483 - t322 * t435 + t313 * t628 + t633 * t641;
t4 = -pkin(9) * t5 - t29;
t615 = t21 * t4;
t614 = t21 * t5;
t612 = -t138 / 0.2e1;
t610 = t166 / 0.2e1;
t267 = sin(pkin(18));
t602 = t267 / 0.2e1;
t269 = cos(pkin(18));
t599 = -t269 / 0.2e1;
t597 = -t273 / 0.2e1;
t596 = t273 / 0.2e1;
t595 = t277 / 0.2e1;
t591 = pkin(3) * t166;
t588 = pkin(3) * t267;
t587 = pkin(3) * t269;
t281 = 0.1e1 / pkin(8);
t584 = pkin(4) * t281;
t582 = pkin(5) * t289;
t579 = t270 * t5;
t578 = t274 * t5;
t497 = qJD(3) * t271;
t499 = qJD(2) * t272;
t507 = t275 * t276;
t229 = -t271 * t499 - t272 * t497 + t507 * t629;
t382 = t271 * t276 + t272 * t275;
t230 = t629 * t382;
t246 = t271 * t272 - t507;
t26 = -t229 * t95 + t230 * t96 + t246 * t44 - t382 * t43;
t74 = t246 * t96 - t382 * t95;
t577 = t74 * t26;
t570 = t262 * t36;
t569 = t263 * t36;
t209 = pkin(1) * t415;
t491 = 0.2e1 * t241 * pkin(5);
t444 = -t233 + t491;
t187 = t209 + (t243 * t444 - t513) * pkin(1);
t425 = pkin(5) * t466;
t454 = -t226 * t241 / 0.2e1;
t190 = t223 * t454 - 0.2e1 * t425 + (t244 * t233 - t514) * pkin(1);
t389 = t521 + t530;
t150 = ((t187 * t595 + t190 * t596) * t236 + t389 * t442) * t283;
t520 = t222 * t277;
t531 = t219 * t273;
t208 = (t520 / 0.2e1 - t531 / 0.2e1) * t517;
t202 = t208 ^ 2;
t184 = t200 * t202 + 0.1e1;
t201 = t199 * t200;
t388 = -t520 + t531;
t151 = ((t187 * t597 + t190 * t595) * t236 - t388 * t442) * t283;
t532 = t208 * t151;
t563 = 0.2e1 * (-t150 * t201 * t202 + t200 * t532) / t184 ^ 2;
t557 = t137 * t285;
t552 = t150 * t208;
t549 = t166 * t170;
t547 = t166 * t281;
t533 = t200 * t208;
t518 = t236 * t273;
t510 = t267 * t128;
t509 = t269 * t126;
t508 = t270 * t274;
t498 = qJD(2) * t276;
t496 = qJD(3) * t275;
t495 = qJD(4) * t270;
t494 = qJD(4) * t274;
t27 = t229 * t96 + t230 * t95 + t246 * t43 + t382 * t44;
t75 = t246 * t95 + t382 * t96;
t490 = 0.2e1 * t75 * t27;
t485 = 0.2e1 * t274;
t484 = 0.2e1 * t285;
t482 = -0.2e1 * pkin(13) * qJD(2);
t481 = 0.4e1 * t168 * t286;
t258 = pkin(1) * t499;
t477 = pkin(1) * t497;
t476 = pkin(1) * t496;
t475 = t36 * t508;
t474 = t74 * t508;
t472 = t75 * t495;
t471 = t75 * t494;
t462 = t270 * t494;
t461 = t272 * t498;
t453 = t236 * t595;
t255 = -pkin(1) * t276 - pkin(13);
t452 = pkin(3) * t484;
t228 = -pkin(4) * t230 + t258;
t443 = t503 * t74;
t439 = t126 * t487;
t438 = t128 * t487;
t32 = -pkin(4) * t95 - pkin(9) * t36;
t437 = qJD(4) * t32 * t622;
t436 = -0.4e1 * t474;
t434 = 0.2e1 * t472;
t433 = 0.2e1 * t471;
t432 = t128 * t481;
t20 = t21 ^ 2;
t424 = t20 * t462;
t35 = t36 ^ 2;
t423 = t35 * t462;
t422 = t36 * t462;
t73 = t74 ^ 2;
t421 = t73 * t462;
t76 = -pkin(1) * t567 + t96 * t257;
t414 = t456 + t233;
t412 = 0.2e1 * t480;
t323 = (t137 * t612 + t390 * t620) * t145 - t635;
t85 = t369 * t391 + t390 * t431;
t345 = t138 * t465 + t85 * t611;
t392 = -0.2e1 * t390 * t391;
t306 = (t345 * t171 + (-t137 * t174 + t392) * t451 + (-t137 * t165 + t323) * pkin(3)) * t457;
t319 = (t174 * t620 + t391 * t612 + t85 * t609) * t145 + t555 + t138 * t419;
t308 = (t390 * t449 + (t391 * t445 + t319) * pkin(3)) * t548;
t371 = t390 * t409;
t372 = t390 * t410;
t377 = t391 * t409;
t8 = (t127 * t632 + t129 * t377 + t265 * t306 + t308 * t606 + t390 * t341 + t66 * t371 + t64 * t372 + t89 * t375 + t87 * t376) * t576 - t316 * t639 - (-t127 * t377 + t129 * t632 + t264 * t306 + t308 * t605 + t390 * t342 - t64 * t371 + t66 * t372 - t87 * t375 + t89 * t376) * t483 - t321 * t435 + t312 * t628 + t634 * t641;
t6 = -pkin(4) * t43 - pkin(9) * t8;
t411 = t32 * t8 + t36 * t6;
t408 = t27 * t36 + t75 * t8;
t12 = -pkin(9) * t27 - pkin(11) * t26 + t228;
t235 = -pkin(4) * t246 + t255;
t53 = -pkin(9) * t75 - pkin(11) * t74 + t235;
t404 = t12 * t75 + t27 * t53;
t18 = pkin(11) * t21 + t76;
t403 = t18 * t75 + t19 * t74;
t402 = t26 * t75 + t27 * t74;
t33 = pkin(4) * t96 + pkin(11) * t36;
t401 = t32 * t74 + t33 * t75;
t400 = t101 * t281 * t479;
t399 = t101 * t427;
t397 = t18 * t417;
t178 = t203 * t276 - t204 * t272;
t177 = t203 * t272 + t204 * t276;
t113 = t156 * t544 - t157 * t468 + 0.1e1;
t366 = t21 * t494 + t579;
t365 = t21 * t495 - t578;
t363 = t273 * t395;
t361 = qJD(4) * (t19 * t36 + t21 * t32);
t354 = t277 * t643;
t349 = 0.2e1 * t402;
t191 = (-t512 + (-t414 + t491) * t249) * pkin(1);
t337 = t187 * t249 + t191 * t243 + t219 * t244;
t193 = t224 * t454 + t582 * t621 + (t233 * t343 - t511) * pkin(1);
t334 = t190 * t249 + t193 * t243 + t222 * t244;
t331 = -t26 * t508 + t429 * t74;
t330 = qJD(4) * t436 - t503 * t26;
t328 = t19 * t8 + t21 * t6 + t32 * t5 + t36 * t4;
t28 = t44 * t257 + t96 * t477 + (-t271 * t43 - t95 * t496) * pkin(1);
t3 = pkin(11) * t5 + t28;
t327 = t18 * t27 + t19 * t26 + t3 * t75 + t4 * t74;
t7 = pkin(4) * t44 + pkin(11) * t8;
t326 = t26 * t32 + t27 * t33 + t6 * t74 + t7 * t75;
t325 = t281 * (t126 * t481 + 0.2e1 * t591) * t557;
t250 = 0.2e1 * t255 * t258;
t182 = 0.1e1 / t184;
t159 = ((t191 * t595 + t193 * t596) * t236 + t389 * t441) * t283;
t158 = ((t191 * t597 + t193 * t595) * t236 - t388 * t441) * t283;
t154 = 0.4e1 * t425 + (t244 * t444 + t318) * pkin(1);
t152 = t209 - t344 * t241 + t383 * t582 + (-t243 * t414 - t513) * pkin(1);
t147 = (-t177 * t267 - t178 * t269) * pkin(3) + t255;
t133 = -qJD(2) * t177 - t148 * t272 + t149 * t276;
t132 = qJD(2) * t178 + t148 * t276 + t149 * t272;
t114 = (t158 * t199 - t159 * t533) * t182;
t111 = pkin(1) * t204 + t113 * t588;
t110 = pkin(1) * t203 + t113 * t587;
t104 = t258 + (-t132 * t267 - t133 * t269) * pkin(3);
t98 = (-t509 / 0.2e1 - t510 / 0.2e1) * t547;
t97 = (t126 * t602 + t128 * t599) * t547;
t88 = t138 * t460 + t390 * t450 + (t164 * t174 - t636) * pkin(4);
t86 = (t390 * t446 + t347) * pkin(4);
t82 = (t390 * t438 + t88 * t610) * t281;
t81 = (t390 * t439 + t86 * t610) * t281;
t71 = t106 * t460 + t141 * t450 + (t140 * t164 - t554) * pkin(4);
t69 = (t141 * t446 + t348) * pkin(4);
t68 = -t177 * t97 + t178 * t98;
t67 = t177 * t98 + t178 * t97;
t58 = (t141 * t438 + t71 * t610) * t281;
t57 = (t141 * t439 + t69 * t610) * t281;
t56 = t110 * t98 - t111 * t97;
t55 = t110 * t97 + t111 * t98;
t52 = (-t150 * t182 * t200 - t199 * t563) * t158 + (t533 * t563 + (-t151 * t200 + 0.2e1 * t201 * t552) * t182) * t159 + ((t152 * t453 + t222 * t354 - t154 * t518 / 0.2e1 - 0.4e1 * t219 * t363) * t199 - (t154 * t453 + t219 * t354 + t152 * t518 / 0.2e1 + 0.4e1 * t222 * t363) * t533 + ((t199 * t334 - t337 * t533) * t277 + (-t199 * t337 - t334 * t533) * t273) * t488) * t182 * t283;
t50 = pkin(1) * t149 + t51 * t587;
t49 = pkin(1) * t148 + t51 * t588;
t48 = (t122 * t82 - t81 * t562) * t412;
t46 = ((t63 * t599 - t267 * t65 / 0.2e1) * t166 + (-t509 - t510) * t137 * t487) * t281;
t45 = ((t65 * t599 + t63 * t602) * t166 + (t267 * t439 - t269 * t438) * t137) * t281;
t40 = (t122 * t58 - t57 * t562) * t412 + t113;
t24 = -t132 * t97 + t133 * t98 - t177 * t45 + t178 * t46;
t23 = t132 * t98 + t133 * t97 + t177 * t46 + t178 * t45;
t16 = t110 * t46 - t111 * t45 - t49 * t97 + t50 * t98;
t15 = t110 * t45 + t111 * t46 + t49 * t98 + t50 * t97;
t10 = 0.2e1 * ((-t170 * t345 + t392 * t452) * t610 + (-t174 * t591 + t390 * t432) * t557 + ((-t137 * t164 + t323) * t610 + (t128 * t391 + t137 * t88 + t390 * t65) * t590) * pkin(4)) * t400 - 0.2e1 * (t390 * t325 + ((t319 - t638) * t610 + (t391 * t549 + (t126 * t391 + t137 * t86 + t390 * t63) * t167) * pkin(3)) * t584) * t399 + 0.2e1 * t627 * t82 + 0.2e1 * t624 * t81;
t9 = 0.2e1 * ((-t170 * t346 + t394 * t452) * t610 + (-t140 * t591 + t141 * t432) * t557 + ((t83 * t164 + t324) * t610 + (t84 * t128 + t137 * t71 + t141 * t65) * t590) * pkin(4)) * t400 - 0.2e1 * (t141 * t325 + ((-t84 * t164 + t320) * t610 + (t84 * t549 + (t126 * t84 + t137 * t69 + t141 * t63) * t167) * pkin(3)) * t584) * t399 + t51 + 0.2e1 * t57 * t624 + 0.2e1 * t58 * t627;
t2 = t331 * t36 - t8 * t474;
t1 = t21 * t331 - t5 * t474;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t461, 0.2e1 * (-t272 ^ 2 + t276 ^ 2) * qJD(2), 0, -0.2e1 * t461, 0, 0, t272 * t482, t276 * t482, 0, 0, 0.2e1 * t382 * t229, -0.2e1 * t229 * t246 - 0.2e1 * t230 * t382, 0, 0.2e1 * t246 * t230, 0, 0, -0.2e1 * t230 * t255 - 0.2e1 * t246 * t258, -0.2e1 * t229 * t255 - 0.2e1 * t258 * t382, 0, t250, 0.2e1 * t577, t349, 0, t490, 0, 0, -0.2e1 * t228 * t75 - 0.2e1 * t235 * t27, 0.2e1 * t228 * t74 + 0.2e1 * t235 * t26, 0, 0.2e1 * t235 * t228, 0.2e1 * t263 * t577 - 0.2e1 * t421, t26 * t436 + t381 * t73, -t402 * t485 + t434 * t74, 0.2e1 * t262 * t577 + 0.2e1 * t421, t270 * t349 + t433 * t74, t490, -t404 * t485 + t434 * t53, t404 * t486 + t433 * t53, -0.2e1 * (t12 * t74 + t26 * t53) * t502, t53 * t12 * t417, 0.2e1 * t532, 0.2e1 * t151 * t205 + 0.2e1 * t552, 0, 0.2e1 * t205 * t150, 0, 0, -0.2e1 * pkin(7) * t150, 0.2e1 * pkin(7) * t151, 0, 0, 0.2e1 * t177 * t132, 0.2e1 * t132 * t178 + 0.2e1 * t133 * t177, 0, 0.2e1 * t178 * t133, 0, 0, -0.2e1 * t133 * t255 - 0.2e1 * t178 * t258, 0.2e1 * t132 * t255 + 0.2e1 * t177 * t258, 0, t250, 0.2e1 * t67 * t23, 0.2e1 * t23 * t68 + 0.2e1 * t24 * t67, 0, 0.2e1 * t68 * t24, 0, 0, -0.2e1 * t104 * t68 - 0.2e1 * t147 * t24, 0.2e1 * t104 * t67 + 0.2e1 * t147 * t23, 0, 0.2e1 * t147 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t498, 0, -t499, 0, 0, 0, 0, 0, 0, 0, -t229, 0, t230, 0, 0, 0, (-t229 * t275 - t230 * t271 + (-t246 * t275 + t271 * t382) * qJD(3)) * pkin(1), 0, 0, 0, t21 * t26 + t5 * t74, 0, t21 * t27 + t5 * t75, 0, 0, 0, -t26 * t77 + t27 * t76 + t28 * t75 - t29 * t74, 0, -t1, t21 * t330 - t443 * t5, -t75 * t579 + (-t27 * t270 - t471) * t21, t1, -t75 * t578 + (-t27 * t274 + t472) * t21, 0, t270 * t327 + t403 * t494, t274 * t327 - t403 * t495, 0, 0, 0, 0, t114 * t151 + t208 * t52, 0, t114 * t150 + t205 * t52, 0, 0, 0, 0, 0, 0, 0, t113 * t132 + t177 * t51, 0, t113 * t133 + t178 * t51, 0, 0, 0, (-t132 * t203 + t133 * t204 + t148 * t178 - t149 * t177) * pkin(1), 0, 0, 0, t23 * t40 + t67 * t9, 0, t24 * t40 + t68 * t9, 0, 0, 0, t15 * t68 - t16 * t67 - t23 * t56 + t24 * t55, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t477, 0.2e1 * t476, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t614, 0.2e1 * t21 * t29 + 0.2e1 * t5 * t77, -0.2e1 * t21 * t28 - 0.2e1 * t5 * t76, 0, 0.2e1 * t28 * t76 + 0.2e1 * t29 * t77, 0.2e1 * t262 * t614 + 0.2e1 * t424, -t20 * t381 + 0.4e1 * t508 * t614, 0, 0.2e1 * t263 * t614 - 0.2e1 * t424, 0, 0, -0.2e1 * t274 * t615 + t365 * t489, t366 * t489 + t486 * t615, t21 * t3 * t417 + t397 * t5, t3 * t397 + t4 * t489, 0, 0, 0, 0, 0, 0.2e1 * t114 * t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113 * t51, (t113 * t149 + t203 * t51) * t493, (-t113 * t148 - t204 * t51) * t493, 0, 0.2e1 * (t148 * t204 + t149 * t203) * t289, 0, 0, 0, 0, 0, 0.2e1 * t40 * t9, 0.2e1 * t16 * t40 + 0.2e1 * t56 * t9, -0.2e1 * t15 * t40 - 0.2e1 * t55 * t9, 0, 0.2e1 * t15 * t55 + 0.2e1 * t16 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, 0, t230, 0, 0, 0, 0, 0, 0, 0, t26 * t36 + t74 * t8, 0, t408, 0, 0, 0, (-t26 * t95 + t27 * t96 - t43 * t74 + t44 * t75) * pkin(4), 0, -t2, t330 * t36 - t443 * t8, -t270 * t408 - t36 * t471, t2, -t274 * t408 + t36 * t472, 0, t270 * t326 + t401 * t494, t274 * t326 - t401 * t495, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t67 + t23 * t48, 0, t10 * t68 + t24 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t477, t476, 0, 0, 0, 0, 0, 0, 0, t21 * t8 + t36 * t5, t29 * t36 + t77 * t8 + (t21 * t43 + t5 * t95) * pkin(4), -t28 * t36 - t76 * t8 + (-t21 * t44 - t5 * t96) * pkin(4), 0, (t28 * t96 + t29 * t95 + t43 * t77 + t44 * t76) * pkin(4), t5 * t570 + (t262 * t8 + 0.2e1 * t422) * t21, 0.2e1 * t5 * t475 + (-t36 * t381 + 0.2e1 * t508 * t8) * t21, 0, t5 * t569 + (t263 * t8 - 0.2e1 * t422) * t21, 0, 0, t270 * t361 - t274 * t328, t270 * t328 + t274 * t361, t502 * (t18 * t8 + t21 * t7 + t3 * t36 + t33 * t5), t19 * t6 + t32 * t4 + (t18 * t7 + t3 * t33) * t502, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t40 + t48 * t9, t10 * t56 + t16 * t48, -t10 * t55 - t15 * t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t622, (t36 * t43 + t8 * t95) * t492, (-t36 * t44 - t8 * t96) * t492, 0, (t43 * t95 + t44 * t96) * t484, 0.2e1 * t570 * t8 + 0.2e1 * t423, -t35 * t381 + 0.4e1 * t475 * t8, 0, 0.2e1 * t569 * t8 - 0.2e1 * t423, 0, 0, t270 * t437 - t411 * t485, t274 * t437 + t411 * t486, (t33 * t8 + t36 * t7) * t417, t33 * t417 * t7 + 0.2e1 * t32 * t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t48 * t10, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t274 - t495 * t74, 0, -t26 * t270 - t494 * t74, -t27, t12 * t274 - t495 * t53, -t12 * t270 - t494 * t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t366, 0, -t365, 0, -t18 * t494 - t270 * t3, t18 * t495 - t274 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270 * t8 + t36 * t494, 0, t274 * t8 - t36 * t495, 0, -t270 * t7 - t33 * t494, -t274 * t7 + t33 * t495, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
