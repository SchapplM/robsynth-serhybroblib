% Calculate inertial parameters regressor of gravitation load for
% picker2Dm1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = picker2Dm1TE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1TE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t295 = cos(pkin(9));
t291 = sin(pkin(8));
t292 = cos(pkin(8));
t736 = sin(qJ(1));
t738 = cos(qJ(1));
t195 = -t291 * t736 - t292 * t738;
t729 = pkin(5) * t195;
t764 = -2 * pkin(1);
t185 = t729 * t764;
t325 = pkin(5) ^ 2;
t687 = t185 + t325;
t806 = 2 * pkin(1);
t171 = -(t806 + pkin(5)) * pkin(5) + t687;
t172 = pkin(5) * (t806 - pkin(5)) + t687;
t336 = sqrt(-t171 * t172);
t194 = -t291 * t738 + t292 * t736;
t730 = pkin(5) * t194;
t597 = 0.1e1 / t336 * (t171 + t172) * t806 * t730 / 0.2e1;
t609 = -pkin(1) + t729;
t759 = 0.2e1 * t325;
t686 = t185 + t759;
t703 = t195 * t336;
t113 = -pkin(5) * t703 + (t609 * t806 + t597 - t686) * t730;
t590 = pkin(5) * t686;
t734 = pkin(1) * t194;
t115 = -t609 * t597 + t195 * t590 + (-pkin(5) * t336 + t734 * t759) * t194;
t294 = cos(qJ(2));
t354 = t738 ^ 2;
t348 = (pkin(4) ^ 2);
t350 = (pkin(7) ^ 2);
t343 = (pkin(3) ^ 2);
t754 = 2 * t343;
t666 = t350 + t754;
t610 = -t348 + t666;
t339 = pkin(1) ^ 2;
t757 = 3 * t339;
t555 = t757 + t610;
t604 = t294 * t736;
t570 = pkin(1) * t604;
t735 = sin(qJ(2));
t648 = t735 * pkin(7);
t586 = 0.2e1 * t648;
t441 = (-0.4e1 * t570 + t586) * pkin(3) + t555;
t635 = pkin(3) * t735;
t557 = t635 + pkin(7);
t437 = t557 * t441;
t517 = pkin(3) * t586 + t350;
t357 = t735 ^ 2;
t691 = t343 * t357;
t651 = 0.2e1 * t691;
t463 = -t343 + t517 + t651;
t455 = t463 * t764;
t565 = pkin(3) * t604;
t766 = 0.1e1 / pkin(3);
t714 = t766 / 0.2e1;
t803 = 0.1e1 / t714;
t712 = pkin(7) * t803;
t470 = (t565 - pkin(1)) * t712;
t513 = pkin(3) * t555;
t472 = t736 * t513;
t668 = t348 - t350;
t612 = t339 - t668;
t753 = 4 * t343;
t341 = t339 ^ 2;
t438 = pkin(1) * (0.2e1 * (-t570 + t648) * pkin(3) + t612);
t434 = t557 * t438;
t495 = t339 + t517;
t477 = -t348 + t495;
t460 = pkin(1) * pkin(3) * t477;
t454 = t736 * t460;
t715 = -t766 / 0.4e1;
t509 = pkin(7) * t612 / t715;
t670 = t339 - t350;
t516 = t670 * t753;
t693 = t339 * t354;
t656 = -0.4e1 * t693;
t338 = sqrt(t463 * t656 - 0.4e1 * t738 * t434 + 0.4e1 * t294 * t454 + t357 * t516 + t735 * t509 - t341 - (2 * t610 * t339) - (t350 - (t803 + pkin(4)) * pkin(4)) * (t350 + (t803 - pkin(4)) * pkin(4)));
t483 = t736 * t557;
t605 = t294 * t738;
t566 = pkin(3) * t605;
t452 = t566 + t483;
t792 = t452 * t338;
t422 = t294 * t472 + t792 + t354 * t455 - t738 * t437 + t735 * t470 + (t651 - t753 - t612) * pkin(1);
t484 = t557 * t738;
t451 = -t565 + t484;
t433 = t451 * t806 + t343 + t495;
t431 = 0.1e1 / t433;
t414 = t431 * t422;
t410 = t766 * t414;
t177 = t339 + t687;
t761 = 0.1e1 / t177;
t765 = 0.1e1 / pkin(5);
t731 = t765 * t761;
t404 = t410 * t731;
t453 = t463 * t738;
t466 = t754 + t477;
t429 = t453 * t806 + t466 * t557;
t777 = (0.4e1 * t354 - 0.2e1) * pkin(1);
t430 = pkin(3) * (t466 * t738 + t557 * t777);
t447 = pkin(1) + t451;
t424 = t294 * t430 + t447 * t338 + t429 * t736;
t420 = t431 * t424;
t418 = t420 / 0.4e1;
t413 = t766 * t418 * t731;
t704 = t194 * t336;
t456 = pkin(5) * t704 - t609 * t686;
t457 = -t194 * t590 - t336 * t609;
t175 = 0.1e1 / t177 ^ 2;
t705 = t175 * t194;
t750 = pkin(1) / 0.2e1;
t130 = 0.1e1 / t338;
t756 = 8 * t339;
t251 = pkin(3) * t294;
t492 = t557 * t251;
t778 = t354 * t492 + t736 * t453;
t428 = 0.4e1 * t736 * t434 + 0.4e1 * t460 * t605 + t756 * t778;
t427 = t130 * t428;
t426 = t427 / 0.2e1;
t650 = pkin(7) * t738;
t589 = 0.2e1 * t650;
t542 = t294 * t589;
t599 = t735 * t343;
t763 = 4 * pkin(1);
t371 = t431 * (t451 * t338 + t452 * t426 + t736 * t437 + t513 * t605 + t542 * t599 + t763 * t778);
t370 = t766 * t371;
t432 = 0.1e1 / t433 ^ 2;
t415 = t432 * t422;
t791 = t452 * t806;
t389 = t415 * t791;
t386 = t766 * t389;
t632 = t731 / 0.4e1;
t776 = (t370 + t386) * t632;
t419 = t766 * t420;
t779 = t456 * t410 + t457 * t419;
t253 = pkin(1) * t738;
t665 = t736 ^ 2;
t367 = t431 * (-t792 + t447 * t426 + t665 * t455 + t429 * t738 + (-0.8e1 * t483 * t253 - t466 * t736) * t251);
t366 = t766 * t367;
t421 = t432 * t424;
t398 = t421 * t791;
t394 = t766 * t398;
t633 = -t731 / 0.4e1;
t785 = t366 * t632 - t394 * t633;
t360 = t113 * t404 / 0.4e1 + t115 * t413 + t779 * t750 * t705 + t785 * t457 + t776 * t456;
t396 = t457 * t410;
t399 = t456 * t419;
t58 = -t115 * t404 / 0.4e1 + t113 * t413 + (-pkin(1) * t396 / 0.2e1 + t399 * t750) * t705 - t776 * t457 + t785 * t456;
t737 = sin(pkin(9));
t41 = t295 * t360 + t737 * t58;
t321 = 2 * pkin(2);
t323 = pkin(6) ^ 2;
t100 = t396 * t633 + t399 * t632;
t362 = t779 * t632;
t92 = t737 * t100 + t295 * t362;
t746 = pkin(6) * t92;
t91 = t746 * t321;
t711 = t323 + t91;
t84 = -(t321 + pkin(6)) * pkin(6) + t711;
t85 = pkin(6) * (t321 - pkin(6)) + t711;
t511 = pkin(6) * (-t84 - t85) * t321;
t23 = t41 * t511;
t337 = sqrt(-t84 * t85);
t42 = t295 * t58 - t737 * t360;
t762 = -2 * pkin(2);
t89 = t91 + 0.2e1 * t323;
t90 = -pkin(2) - t746;
t608 = t90 * t762 + t89;
t67 = 0.1e1 / t337;
t93 = t100 * t295 - t737 * t362;
t623 = t67 * t93 / 0.2e1;
t10 = (t23 * t623 + t337 * t42 + t41 * t608) * pkin(6);
t407 = -t410 / 0.2e1;
t639 = -t738 / 0.2e1;
t107 = t736 * t407 + t419 * t639;
t636 = t736 / 0.2e1;
t108 = t738 * t407 + t419 * t636;
t593 = t323 * t93 * t762;
t624 = -t67 * t90 / 0.2e1;
t11 = t23 * t624 + t41 * t593 + (t337 * t41 - t42 * t89) * pkin(6);
t745 = pkin(6) * t93;
t43 = t337 * t745 - t89 * t90;
t324 = 0.1e1 / pkin(6);
t88 = (pkin(2) ^ 2) + t711;
t717 = t324 / t88;
t44 = -t337 * t90 - t89 * t745;
t749 = -t44 / 0.2e1;
t638 = t738 / 0.2e1;
t784 = t366 / 0.2e1 + t394 / 0.2e1;
t77 = t370 * t639 - t386 * t638 + t736 * t784 - t107;
t637 = -t736 / 0.2e1;
t78 = t370 * t636 - t386 * t637 + t738 * t784 - t108;
t718 = t108 * t44;
t720 = t107 * t43;
t747 = pkin(2) / t88 ^ 2;
t788 = (t718 + t720) * t747;
t3 = t41 * t788 + (t78 * t43 / 0.2e1 - t107 * t10 / 0.2e1 + t77 * t749 - t108 * t11 / 0.2e1) * t717;
t789 = (t718 / 0.2e1 + t720 / 0.2e1) * t717;
t817 = t3 * t92 - t41 * t789;
t545 = t735 * t738;
t200 = -t545 + t604;
t476 = t294 * (0.4e1 * t599 + t712);
t465 = t354 * t476;
t473 = pkin(3) * t484;
t544 = t735 * t736;
t519 = pkin(1) * t544;
t504 = 0.4e1 * t519;
t553 = t294 * t599;
t649 = pkin(7) * t736;
t576 = t343 * t649;
t744 = pkin(7) * t294;
t758 = -4 * t339;
t760 = t294 ^ 2;
t805 = 8 * pkin(1);
t129 = t465 * t758 - 0.4e1 * t438 * t566 - 0.4e1 * (t519 + t744) * t473 * t806 + t760 * t576 * t805 - 0.4e1 * t735 * t454 + (0.2e1 * t735 * t516 + t509) * t294;
t706 = t129 * t130;
t598 = t706 / 0.2e1;
t662 = 0.2e1 * t744;
t732 = pkin(3) * t338;
t377 = t431 * (t200 * t732 + t452 * t598 + t465 * t764 - t441 * t566 - (t504 + t662) * t473 + t553 * t763 - 0.2e1 * t357 * t576 + t294 * t470 - t735 * t472);
t375 = t766 * t377;
t199 = t544 + t605;
t440 = t199 * t806 + t662;
t387 = t440 * t415;
t775 = t375 * t632 + t387 * t633;
t582 = 0.2e1 * t253;
t380 = t431 * (t199 * t732 + t447 * t598 + (0.2e1 * pkin(7) * t492 + t466 * t251 + t476 * t582) * t736 - t735 * t430 + t343 * t760 * (t589 + t777));
t379 = t766 * t380;
t397 = t440 * t421;
t780 = t379 * t632 + t397 * t633;
t361 = t456 * t775 + t457 * t780;
t72 = t456 * t780 - t457 * t775;
t53 = t295 * t361 + t737 * t72;
t36 = t53 * t511;
t54 = t295 * t72 - t737 * t361;
t17 = (t337 * t54 + t36 * t623 + t53 * t608) * pkin(6);
t18 = t36 * t624 + t53 * t593 + (t337 * t53 - t54 * t89) * pkin(6);
t743 = t107 / 0.2e1;
t781 = -t375 / 0.2e1 + t387 / 0.2e1;
t81 = t379 * t639 + t397 * t638 + t736 * t781;
t748 = t81 / 0.2e1;
t82 = t379 * t636 + t397 * t637 + t738 * t781;
t5 = t53 * t788 - (t82 * t44 / 0.2e1 + t108 * t18 / 0.2e1 + t43 * t748 + t17 * t743) * t717;
t816 = t5 * t92 - t53 * t789;
t327 = 0.1e1 / pkin(4);
t694 = t327 / t343;
t378 = t380 * t694;
t391 = pkin(3) * t397 * t694;
t406 = t414 * t694;
t402 = t130 * t406 / 0.8e1;
t417 = t420 * t694;
t536 = pkin(1) * t565;
t682 = t339 / 0.3e1 + t350;
t182 = -0.4e1 / 0.9e1 * t536 + 0.4e1 / 0.9e1 * t343 - t348 / 0.9e1 + t682;
t273 = -t348 / 0.6e1;
t282 = 0.2e1 / 0.3e1 * t343;
t493 = t350 - t536;
t191 = t282 + t273 + t493;
t281 = 0.4e1 / 0.3e1 * t343;
t257 = t339 + t350;
t275 = -t348 / 0.3e1;
t616 = t275 + t257;
t223 = t281 + t616;
t249 = -t339 / 0.3e1 + t350;
t538 = -0.2e1 * t565;
t588 = 0.4e1 * t650;
t654 = 0.6e1 * t693;
t263 = t738 * t354;
t340 = pkin(1) * t339;
t698 = t263 * t340;
t658 = pkin(7) * t698;
t696 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t135 = 0.4e1 * t658 + t182 * t654 + t223 * t696 + (t191 * t588 + t249 * t538) * pkin(1);
t312 = 6 * t339;
t277 = -0.2e1 / 0.3e1 * t348;
t320 = 2 * t350;
t615 = t277 + t282 + t320;
t345 = t343 ^ 2;
t614 = t277 + t257;
t684 = (t282 + t614) * t257 + t345;
t149 = -0.4e1 * t223 * t536 + (t312 + t615) * t343 + t684;
t285 = -t343 / 0.3e1;
t247 = t285 + t350;
t515 = -0.2e1 * t536;
t198 = t247 * t515;
t697 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t158 = t223 * t697 + t198;
t264 = 0.10e2 / 0.3e1 * t339;
t159 = (t264 + t615) * t343 + t684;
t238 = pkin(7) * t582;
t259 = -3 * t339 + t350;
t655 = 0.4e1 * t693;
t197 = t238 + t655 + t259;
t256 = -3 * t343 + t350;
t591 = 0.8e1 * t658;
t211 = t256 * t591;
t231 = t343 + t612;
t236 = t253 + pkin(7);
t250 = t257 ^ 2;
t296 = 15 * t341;
t303 = 18 * t350;
t304 = -2 * t348;
t306 = -6 * t348;
t335 = t350 ^ 2;
t315 = 3 * t335;
t344 = pkin(3) * t343;
t328 = t344 ^ 2;
t640 = 0.12e2 * t691;
t642 = 0.12e2 * t693;
t660 = 6 * pkin(1);
t319 = 3 * t350;
t679 = 15 * t339 + t319;
t261 = t735 * t357;
t699 = t261 * t344;
t187 = t515 + t223;
t239 = 4 * (-t343 + t350) * t339;
t584 = pkin(7) * t253;
t137 = 0.4e1 * t187 * t584 + t239 * t354 + t149;
t801 = 0.6e1 * t137;
t111 = t135 * t640 + t211 + t158 * t642 + t328 + ((t343 - t348 + t679) * t345) + ((t296 + (t303 + t306 + 6 * t343) * t339 + t315 + (t304 + t754) * t350) * t343) + (t250 * t231) + (0.8e1 * t197 * t699 + t635 * t801) * t236 + (t149 * t650 - t159 * t565) * t660;
t600 = t736 * t340;
t644 = t339 * t251;
t478 = -t600 + t644;
t206 = 0.2e1 * t478;
t313 = 2 * t339;
t255 = t313 + t343;
t151 = -t670 * t251 + t206 * t354 + (pkin(3) * t542 + t736 * t255) * pkin(1);
t667 = t350 + t757;
t241 = t343 + t667;
t210 = t241 * t251;
t242 = 3 * t343 + t257;
t601 = t736 * t242;
t176 = -pkin(1) * t601 + t210;
t657 = t344 * t756;
t804 = 0.4e1 / t766;
t235 = t341 * t804 + t657;
t552 = t736 * t697;
t178 = t235 * t294 + 0.4e1 * t340 * t552;
t311 = 5 * t341;
t671 = t335 + t345;
t298 = 10 * t339;
t678 = t298 + t320;
t690 = t350 * t339;
t193 = t343 * t678 + t311 + t671 + 6 * t690;
t307 = 5 * t345;
t317 = 6 * t350;
t202 = t307 + (t298 + t317) * t343 + t250;
t554 = t736 * t699;
t607 = t236 * t735;
t652 = -0.4e1 * t691;
t252 = pkin(1) * t736;
t214 = -t252 + t251;
t573 = t354 * t644;
t148 = -0.2e1 * t573 + t210 + (t214 * t589 - t601) * pkin(1);
t800 = -0.4e1 * t148;
t124 = t151 * t652 + t178 * t354 + (t607 * t800 + (-t193 + t591) * t294) * pkin(3) + (-0.4e1 * t176 * t650 + t736 * t202 - 0.8e1 * t236 * t554) * pkin(1);
t208 = t253 + t557;
t103 = t111 * t208 + t124 * t338;
t102 = 0.1e1 / t103 ^ 2;
t289 = t339 / 0.2e1;
t681 = t289 + t350;
t217 = 0.7e1 / 0.6e1 * t343 + t273 + t681;
t287 = 0.4e1 / 0.3e1 * t339;
t283 = t343 / 0.3e1;
t617 = t273 + t283 + t350;
t220 = t287 + t617;
t157 = -t217 * t252 + t220 * t251;
t618 = t348 / 0.3e1 + t283 + t320;
t163 = -(t343 * t670) - 0.5e1 / 0.3e1 * t341 + t618 * t339 + t350 * (t275 + t247);
t265 = -0.20e2 / 0.3e1 * t339;
t318 = 4 * t350;
t619 = 0.2e1 / 0.3e1 * t348 + t282 + t318;
t620 = 0.4e1 / 0.3e1 * t348 + t281 - (2 * t350);
t164 = -t345 + (t265 + t619) * t343 - (3 * t341) + t620 * t339 + t335;
t225 = t339 + t617;
t243 = -t339 + t666;
t168 = t225 * t251 - t243 * t252 / 0.2e1;
t606 = t263 * t736;
t550 = t341 * t606;
t262 = t354 ^ 2;
t692 = t341 * t262;
t127 = -0.4e1 * pkin(7) * t550 + t157 * t655 + (-0.8e1 / 0.3e1 * t692 + t163) * t251 + (t164 * t637 + t168 * t588) * pkin(1);
t286 = -0.2e1 / 0.3e1 * t343;
t248 = t286 + t350;
t677 = t304 - 2 * t343;
t613 = t317 + t677;
t156 = t345 + (t277 + t286 + t678) * t343 + t311 + (t613 * t339) + t350 * (t277 + t248);
t153 = t156 * t251;
t167 = t307 + ((t298 + t613) * t343) + (t286 + t614) * t257;
t602 = t736 * t167;
t139 = -pkin(1) * t602 + t153;
t673 = t335 - t341;
t165 = -(3 * t345) + (t265 + t620) * t343 + t619 * t339 + t673;
t169 = -0.5e1 / 0.3e1 * t345 + (-t339 + t618) * t343 + t350 * (t285 + t616);
t579 = -0.2e1 * t252;
t140 = t165 * t251 + t169 * t579;
t786 = t319 - t343 - t348;
t230 = t786 * t298;
t310 = 7 * t341;
t316 = 8 * t350;
t611 = -t343 - t668;
t305 = -5 * t348;
t676 = t305 - 5 * t343;
t141 = t328 + ((21 * t339 + t786) * t345) + ((t350 * t677 + t230 + t315 + 35 * t341) * t343) + ((t310 + (t316 + t676) * t339 + t350 * t611) * t257);
t207 = 0.4e1 * t478;
t213 = 0.2e1 * t252 + t251;
t276 = -t348 / 0.2e1;
t227 = t276 + t241;
t143 = t259 * t251 + t207 * t354 + (t213 * t650 + t736 * t227) * t806;
t701 = t250 * (t339 + t611);
t755 = -6 * t343;
t150 = 0.7e1 * t328 + ((35 * t339 + 15 * t350 + t676) * t345) + ((21 * t341 + t230 + 9 * t335 + (t306 + t755) * t350) * t343) + t701;
t218 = t350 + 0.5e1 / 0.2e1 * t343 + 0.3e1 / 0.2e1 * t339 + t276;
t170 = t218 * t251 + t256 * t252 / 0.2e1;
t192 = 0.4e1 / 0.3e1 * t693 + t238 + t249;
t260 = t357 ^ 2;
t700 = t260 * t345;
t797 = -0.24e2 * t700;
t525 = t736 * t797;
t622 = t236 * t699;
t577 = -0.8e1 * t622;
t585 = 0.16e2 * t658;
t641 = -0.12e2 * t691;
t222 = 0.8e1 / 0.3e1 * t343 + t616;
t224 = t275 + t282 + t667;
t161 = -t222 * t252 + t224 * t251;
t229 = 0.5e1 / 0.6e1 * t343 + t289 + t273;
t174 = t229 * t294 * t803 + pkin(1) * t552;
t646 = t263 * t251;
t574 = t340 * t646;
t128 = -0.8e1 * pkin(7) * t574 + t174 * t656 + t153 + (t161 * t588 - t602) * pkin(1);
t802 = -0.6e1 * t128;
t104 = t170 * t585 + t143 * t577 + t127 * t641 - 0.6e1 * t140 * t693 + (t607 * t802 + (0.24e2 * t247 * t692 - t141) * t294) * pkin(3) + (-0.6e1 * t139 * t650 + t736 * t150 + t192 * t525) * pkin(1);
t219 = t350 + t343 / 0.4e1 + t339 / 0.4e1 - t348 / 0.8e1;
t326 = t348 ^ 2;
t680 = 0.4e1 / 0.7e1 * t350 - t348 / 0.7e1;
t689 = t350 * t348;
t146 = -0.32e2 / 0.21e2 * t219 * t536 + 0.5e1 / 0.42e2 * t345 + (0.16e2 / 0.21e2 * t339 + t680) * t343 + t341 / 0.7e1 + t680 * t339 + t335 - 0.3e1 / 0.7e1 * t689 + t326 / 0.42e2;
t274 = -t348 / 0.4e1;
t787 = t274 + t343 / 0.2e1;
t221 = t682 + t787;
t147 = -0.8e1 / 0.3e1 * t221 * t536 + 0.5e1 / 0.18e2 * t345 + (t287 + t275) * t343 + t335 - t341 / 0.3e1 + t326 / 0.18e2 + (t281 + 0.2e1 / 0.3e1 * t339 + t277) * t350;
t672 = t335 + t341;
t675 = t320 - t348;
t485 = (t675 * t339) + t326 / 0.6e1 + t672 - t689;
t188 = -t345 / 0.6e1 + t485;
t190 = -0.2e1 / 0.3e1 * t536 + t274 + t681;
t244 = (t318 + t348) * t339;
t228 = t276 + t343 + t257;
t798 = -0.4e1 * t228;
t508 = t565 * t798;
t653 = 0.8e1 * t692;
t125 = t248 * t653 + t190 * t585 + 0.14e2 * t146 * t693 - (t670 * t345) + (t244 + (2 * t335) - 0.10e2 / 0.3e1 * t341 - t689) * t343 + t188 * t696 + (0.6e1 * t147 * t650 + t249 * t508) * pkin(1);
t469 = 0.5e1 / 0.6e1 * t345 + t485;
t166 = (t264 + t675) * t343 + t469;
t278 = -0.3e1 / 0.2e1 * t348;
t683 = t326 / 0.2e1 - t345 / 0.2e1;
t556 = -(3 * t689) + t315 + t683;
t688 = t257 * ((t278 + t320) * t339 - 0.3e1 / 0.2e1 * t689 + t672 + t683) + t328;
t133 = -0.6e1 * t166 * t536 + (t296 + ((t303 - 9 * t348) * t339) + t556) * t343 + (t278 + t679) * t345 + t688;
t494 = pkin(1) * t508;
t152 = t494 + ((t312 + t675) * t343) + t469;
t160 = t228 * t697 + t198;
t126 = 0.6e1 * t152 * t584 + t160 * t642 + t133 + t211;
t528 = t340 * t565;
t179 = -0.4e1 * t528 + (4 * t341) + ((t753 + t304 + t316) * t339);
t186 = -t339 + t493 + t787;
t136 = t591 + t179 * t354 + t259 * t228 + (t186 * t588 + t538 * t696) * pkin(1);
t138 = t247 * t494 - t328 + (-t264 + t668) * t345 + (t244 + t345 / 0.6e1 - t326 / 0.6e1 + t673) * t343 + t188 * t350;
t142 = (t278 + t319 + (7 * t339)) * t345 + (t310 + ((t305 + 10 * t350) * t339) + t556) * t343 + t688;
t351 = pkin(7) * t350;
t240 = -0.12e2 * pkin(7) * t340 + t351 * t763;
t246 = -8 * t341 + 12 * t690;
t154 = t240 * t738 + t246 * t354 + t585 + t653 + t672 - (6 * t690);
t162 = t256 * t228 + t515 * t697;
t212 = 16 * (t350 * t755 + t671) * t341;
t254 = -30 * t348 + 60 * t350;
t674 = t326 - t345;
t503 = 6 * t335 + t674 - 6 * t689;
t578 = 0.8e1 * t635;
t643 = 0.32e2 * t699;
t98 = t212 * t262 + 0.32e2 * t162 * t658 + 0.24e2 * t138 * t693 + (t304 + t318 + 28 * t339) * t328 + (t231 * t701) + (0.24e2 * t125 * t357 + (t254 * t341) + (t335 * t306) + (t503 * t312) + (t674 * t320) + (28 * t340 ^ 2) + 0.4e1 * t351 ^ 2) * t343 + (t126 * t578 + t136 * t643) * t236 + (t133 * t650 - t142 * t565) * t805 + (0.16e2 * t154 * t260 + (t254 * t339) + (70 * t341) + t345 + t503) * t345;
t79 = t104 * t338 + t98 * t208;
t796 = -t102 * t79 / 0.4e1;
t408 = t417 * t796;
t101 = 0.1e1 / t103;
t416 = t101 * t418;
t412 = t416 * t694;
t490 = t340 * t354 * t544;
t733 = pkin(3) * t247;
t467 = t490 * t733;
t189 = 0.24e2 * t467;
t546 = t738 * t736;
t512 = t339 * t546;
t471 = t735 * t512;
t468 = pkin(7) * t471;
t462 = pkin(3) * t468;
t204 = 0.4e1 * t462;
t520 = pkin(1) * t545;
t751 = pkin(3) * pkin(7);
t489 = t520 * t751;
t215 = -0.2e1 * t489;
t461 = 0.24e2 * t468;
t491 = t341 * t263 * t544;
t497 = t249 * t519;
t502 = pkin(3) * t519;
t505 = -0.4e1 * t520;
t549 = t735 * t698;
t510 = t549 * t751;
t702 = t236 * t294;
t560 = t344 * t357 * t702;
t530 = -0.24e2 * t560;
t603 = t354 * t735;
t548 = t339 * t603;
t533 = -0.4e1 * t548;
t547 = t735 * t692;
t567 = t343 * t607;
t569 = pkin(3) * t607;
t647 = pkin(3) * t702;
t774 = 0.64e2 / 0.3e1 * t219;
t795 = t261 * t345;
t799 = 0.6e1 * t166;
t59 = (-0.96e2 * t192 * t570 * t795 + t215 * t577 + t143 * t530 - 0.24e2 * t127 * t553 - 0.6e1 * (0.8e1 * t229 * t548 - t735 * t156 + (t224 * t505 + 0.8e1 * t549) * pkin(7)) * t567 + t647 * t802 - 0.24e2 * t547 * t733 - 0.16e2 * t218 * t510 + 0.6e1 * t156 * t489 + t141 * t635 + ((-t735 * t259 + t533) * t577 + (0.8e1 / 0.3e1 * t547 + t220 * t533 + t225 * pkin(7) * t505 - t163 * t735) * t641 + 0.6e1 * t165 * t548) * pkin(3)) * t338 + t104 * t598 + 0.8e1 * (0.8e1 * t154 * t294 * t795 + 0.4e1 * (t204 + (0.2e1 * t519 * t696 + 0.4e1 * t490) * pkin(3)) * t622 + 0.12e2 * t136 * t560 + 0.3e1 * (t490 * t774 + 0.4e1 * t228 * t497 + (0.16e2 * t221 * t471 + 0.32e2 / 0.3e1 * t491) * pkin(7)) * pkin(3) * t691 + 0.6e1 * t125 * t553 + (t189 + (t228 * t461 + t519 * t799) * pkin(3)) * t569 + t126 * t647 + 0.8e1 * t491 * t697 * t751 + 0.12e2 * t228 * t467 + t462 * t799 + t142 * t502) * t208 + t98 * t251;
t409 = t414 / 0.4e1;
t405 = t101 * t409;
t740 = -t338 / 0.4e1;
t68 = (t79 * t405 + t420 * t740) * t694;
t722 = t101 * t79;
t594 = t694 / 0.4e1;
t595 = -t694 / 0.4e1;
t782 = pkin(3) * t387 * t595 + t377 * t594;
t488 = 0.6e1 * t502;
t540 = -0.4e1 * t569;
t83 = (t530 * t252 + (t215 + (t670 * t735 - 0.2e1 * t548) * pkin(3)) * t652 - 0.8e1 * t151 * t553 + (t215 + (-t735 * t241 + 0.2e1 * t548) * pkin(3)) * t540 + t647 * t800 - 0.8e1 * t510 - t235 * t603 + 0.4e1 * t241 * t489 + t193 * t635) * t338 + t124 * t598 + (0.24e2 * t197 * t560 + (t204 + (0.8e1 / 0.3e1 * t490 + 0.2e1 * t497) * pkin(3)) * t640 + 0.24e2 * t135 * t553 + 0.6e1 * (t223 * t504 + 0.8e1 * t468) * t567 + t647 * t801 + t189 + pkin(3) * t223 * t461 + t159 * t488) * t208 + t111 * t251;
t723 = t129 * t402 + t83 * t408 + t59 * t412 + (t378 / 0.4e1 - t391 / 0.4e1) * t722 + t782 * t338 - t68;
t401 = t406 * t796;
t403 = t405 * t694;
t411 = -t130 * t417 / 0.8e1;
t69 = (t338 * t409 + t416 * t79) * t694;
t739 = t338 / 0.4e1;
t724 = t129 * t411 + t378 * t740 + t391 * t739 + t83 * t401 + t59 * t403 + t722 * t782 + t69;
t14 = -t724 * t199 + t723 * t200;
t15 = t723 * t199 + t724 * t200;
t716 = -t766 / 0.2e1;
t571 = t327 * t79 * t716;
t531 = t102 * t571;
t501 = t108 * t531;
t564 = t101 * t327 * t714;
t522 = t108 * t564;
t532 = t79 * t564;
t563 = t706 * t715;
t628 = t338 * t716;
t444 = t82 * t532 + t59 * t522 + t83 * t501 + (t107 * t563 + t628 * t81) * t327;
t45 = -t199 * t68 + t200 * t69;
t46 = t199 * t69 + t200 * t68;
t500 = t107 * t327 * t628 + t79 * t522;
t631 = t107 * t714;
t572 = t79 * t631;
t627 = t338 * t714;
t790 = (t101 * t572 + t108 * t627) * t327;
t811 = (-t102 * t572 * t83 + t101 * (t81 * t79 * t714 + t59 * t631) - t108 * t563 - t628 * t82) * t327;
t815 = t14 * t790 + t15 * t500 + t444 * t46 + t45 * t811;
t814 = -t14 * t500 + t15 * t790 - t444 * t45 + t46 * t811;
t364 = t367 * t694;
t392 = t398 * t694;
t543 = pkin(7) * t573;
t216 = -0.4e1 * t543;
t237 = pkin(7) * t579;
t498 = -0.8e1 * t512;
t551 = t354 * t600;
t524 = -0.48e2 * t551;
t507 = pkin(7) * t524;
t523 = -0.32e2 * t550;
t526 = -0.2e1 * t546;
t527 = t341 * t565;
t529 = pkin(1) * t554;
t534 = t247 * t574;
t535 = pkin(7) * t551;
t537 = pkin(1) * t566;
t539 = -0.2e1 * t566;
t559 = t697 * t698;
t575 = t354 * t649;
t580 = -0.4e1 * t252;
t581 = -0.6e1 * t252;
t583 = pkin(7) * t252;
t587 = -0.4e1 * t649;
t592 = pkin(7) * t655;
t659 = pkin(7) * t251;
t506 = -0.24e2 * t535;
t685 = t256 * t506 - 0.24e2 * t534;
t56 = (pkin(1) * (-0.8e1 / 0.3e1 * t512 + t237) * t525 + 0.8e1 * t143 * t529 + (t157 * t498 + (-0.4e1 * t340 * t217 + 0.32e2 / 0.3e1 * t527) * t263) * t641 + t128 * t488 - 0.96e2 * t247 * t263 * t527 + t170 * t507 + 0.12e2 * t169 * t698 + 0.12e2 * t140 * t512 + 0.6e1 * t139 * t583 + (t192 * t797 - t164 * t641 / 0.2e1 + t150) * t253 + ((t168 * t580 - 0.2e1 * t243 * t693 + (0.12e2 * t354 * t665 - 0.4e1 * t262) * t341) * t641 + t256 * t653 + t167 * t654) * pkin(7)) * t338 + t104 * t426 + (0.16e2 * (-t240 * t736 + t246 * t526 + t507 + t523) * t700 - 0.32e2 * t136 * t529 + 0.24e2 * (-0.28e2 * t146 * t512 + t248 * t523 + (t147 * t581 + t190 * t524) * pkin(7) + (-t698 * t774 + t249 * t253 * t798 + (-0.16e2 * t221 * t693 - 0.32e2 / 0.3e1 * t692) * pkin(7)) * t251) * t691 - 0.8e1 * t126 * t502 - 0.4e1 * t212 * t606 - 0.64e2 * t659 * t692 * t697 - 0.96e2 * t162 * t535 - 0.96e2 * t228 * t534 - 0.48e2 * t138 * t512 - 0.48e2 * t166 * t543 - 0.8e1 * t133 * t583 - 0.8e1 * t142 * t537) * t208 - t98 * t252 + ((-0.8e1 * (t207 * t526 + t592 - 0.4e1 * t698 + (-t213 * t649 + t738 * t227) * t806) * t699 - 0.6e1 * (0.8e1 * t174 * t512 - 0.4e1 * t559 - t167 * t253 + (t161 * t580 + (t222 * t758 + 0.24e2 * t528) * t354) * pkin(7)) * t635) * t338 + ((t179 * t526 + t216 + (-0.24e2 * t575 - 0.4e1 * t646) * t340 + (t186 * t587 + t539 * t696) * pkin(1)) * t643 + (0.24e2 * (-t228 * t354 * t659 - t160 * t546) * t339 + (-t152 * t649 - t166 * t566) * t660 + t685) * t578) * t208) * t236;
t783 = t371 * t594 - t389 * t595;
t80 = (t665 * t261 * t657 + (t206 * t526 - 0.2e1 * t698 + (pkin(7) * t538 + t738 * t255) * pkin(1)) * t652 + 0.4e1 * t148 * t502 + ((t294 * t546 * t804 - 0.2e1 * pkin(7) * t354) * t339 + (-0.2e1 * t214 * t649 - t738 * t242) * pkin(1)) * t540 + t506 * t251 + 0.4e1 * t559 + t178 * t526 + t242 * t592 + 0.4e1 * t176 * t583 + (t577 + t202) * t253) * t338 + t124 * t426 + (-0.8e1 * t197 * t529 + 0.8e1 * (t498 + t237) * t622 + (-0.12e2 * t182 * t512 + t216 + (-0.12e2 * t575 - 0.8e1 / 0.3e1 * t646) * t340 + (t191 * t587 + t249 * t539) * pkin(1)) * t640 - 0.6e1 * t137 * t502 + 0.6e1 * (-0.8e1 * t543 + t239 * t526 + (-t187 * t649 - t223 * t566) * t763) * t569 - 0.24e2 * t158 * t512 - 0.24e2 * t223 * t543 + t149 * pkin(7) * t581 - 0.6e1 * t159 * t537 + t685) * t208 - t111 * t252;
t725 = t428 * t402 + t80 * t408 + t56 * t412 + (t364 / 0.4e1 + t392 / 0.4e1) * t722 + t783 * t338 + t68;
t726 = t364 * t740 - t392 * t739 + t80 * t401 + t56 * t403 + t428 * t411 + t722 * t783 - t69;
t12 = -t726 * t199 + t725 * t200;
t13 = t725 * t199 + t726 * t200;
t425 = t327 * t766 * t427;
t561 = t327 * t627;
t25 = t78 * t101 * t571 + t77 * t561 + t108 * t425 / 0.4e1 + (t80 * t531 + t56 * t564) * t107;
t423 = t78 * t561 - t107 * t425 / 0.4e1 + t77 * t532 + t56 * t522 + t80 * t501;
t808 = t12 * t790 + t13 * t500 + t25 * t45 + t423 * t46;
t807 = -t12 * t500 + t13 * t790 + t25 * t46 - t423 * t45;
t719 = t107 * t44;
t481 = (t108 * t43 - t719) * t747;
t621 = -t717 / 0.2e1;
t558 = t108 * t621;
t562 = t43 * t621;
t436 = t77 * t562 + t10 * t558 + (t11 * t743 + t78 * t749) * t717 + t41 * t481;
t499 = t717 * t719 / 0.2e1 + t43 * t558;
t768 = t41 * t499 + t436 * t92;
t435 = t82 * t562 + t17 * t558 + (t18 * t743 + t44 * t748) * t717 + t53 * t481;
t767 = t435 * t92 + t499 * t53;
t742 = -t337 / 0.2e1;
t741 = t337 / 0.2e1;
t728 = t23 * t67;
t727 = t36 * t67;
t710 = t77 * pkin(3) + t252;
t709 = t77 * pkin(2) + t252;
t708 = -t78 * pkin(3) - t253;
t707 = -t78 * pkin(2) - t253;
t695 = t324 / pkin(2);
t664 = -t761 / 0.2e1;
t663 = t761 / 0.2e1;
t645 = pkin(5) * t705;
t634 = t765 * t664;
t630 = -t728 / 0.4e1;
t629 = t728 / 0.4e1;
t626 = -t727 / 0.4e1;
t625 = t727 / 0.4e1;
t596 = t291 * t664;
t568 = t738 * t731;
t541 = -g(1) * t82 - g(2) * t81;
t518 = t736 * t634;
t487 = t107 * t42 + t108 * t41 + t77 * t92 - t78 * t93;
t486 = t107 * t54 + t108 * t53 + t81 * t93 + t82 * t92;
t482 = -g(1) * t736 + g(2) * t738;
t459 = t107 * t41 - t108 * t42 - t77 * t93 - t78 * t92;
t458 = t107 * t53 - t108 * t54 + t81 * t92 - t82 * t93;
t448 = t456 * t734;
t446 = t457 * t738;
t445 = t457 * t736;
t442 = t456 * t634;
t334 = 1 / pkin(1);
t205 = t482 * pkin(1);
t184 = -pkin(1) * t195 + pkin(5);
t180 = t185 + t313;
t132 = -t180 * t734 + t184 * t336;
t131 = pkin(1) * t704 + t180 * t184;
t116 = t184 * t597 + t194 ^ 2 * pkin(5) * t313 + (t180 * t195 - t704) * pkin(1);
t114 = (-t703 + (-0.2e1 * t184 * pkin(5) - t180 + t597) * t194) * pkin(1);
t97 = (t114 * t292 * t664 + t116 * t596) * t334 + (-t131 * t292 - t132 * t291) * t645;
t96 = (t116 * t292 * t663 + t114 * t596) * t334 + (-t131 * t291 + t132 * t292) * t645;
t95 = t765 * t445 * t663 - t115 * t568 / 0.2e1 + t738 * t442 + t113 * t518 + (-t446 * t734 - t736 * t448) * t175;
t94 = t736 * t442 + t113 * t568 / 0.2e1 + t446 * t634 + t115 * t518 + (-t445 * t734 + t738 * t448) * t175;
t1 = [0, 0, 0, 0, 0, 0, t482, -g(1) * t738 - g(2) * t736, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t77 + g(2) * t78, -g(1) * t78 - g(2) * t77, 0, t205, 0, 0, 0, 0, 0, 0, -g(1) * t436 - g(2) * t3, g(1) * t3 - g(2) * t436, 0, -g(1) * t709 - g(2) * t707, 0, 0, 0, 0, 0, 0, -g(1) * t423 - g(2) * t25, g(1) * t25 - g(2) * t423, 0, -g(1) * t710 - g(2) * t708, 0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t96, g(1) * t96 - g(2) * t97, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t487 + g(2) * t459, -g(1) * t459 + g(2) * t487, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t94 + g(2) * t95, -g(1) * t95 - g(2) * t94, 0, t205, 0, 0, 0, 0, 0, 0, -g(2) * t817 - g(1) * t768 + (-g(2) * (t436 * t742 + t499 * t630) - g(1) * (-t3 * t742 + t630 * t789)) * t695, -g(2) * t768 + g(1) * t817 + (-g(2) * (t3 * t741 - t629 * t789) - g(1) * (t436 * t741 + t499 * t629)) * t695, 0, -g(2) * (pkin(6) * t3 + t707) - g(1) * (pkin(6) * t436 + t709), 0, 0, 0, 0, 0, 0, -g(1) * t808 - g(2) * t807, g(1) * t807 - g(2) * t808, 0, -g(2) * (pkin(4) * t25 + t708) - g(1) * (pkin(4) * t423 + t710); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t541, g(1) * t81 - g(2) * t82, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t435 - g(2) * t5, g(1) * t5 - g(2) * t435, 0, t541 * pkin(2), 0, 0, 0, 0, 0, 0, -g(1) * t444 - g(2) * t811, g(1) * t811 - g(2) * t444, 0, t541 * pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t486 + g(2) * t458, -g(1) * t458 + g(2) * t486, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t294 - g(2) * t735, g(1) * t735 - g(2) * t294, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t816 - g(1) * t767 + (-g(2) * (t435 * t742 + t499 * t626) - g(1) * (-t5 * t742 + t626 * t789)) * t695, -g(2) * t767 + g(1) * t816 + (-g(2) * (t5 * t741 - t625 * t789) - g(1) * (t435 * t741 + t499 * t625)) * t695, 0, -g(2) * (pkin(2) * t81 + pkin(6) * t5) - g(1) * (pkin(2) * t82 + pkin(6) * t435), 0, 0, 0, 0, 0, 0, -g(1) * t815 - g(2) * t814, g(1) * t814 - g(2) * t815, 0, -g(2) * (pkin(3) * t81 + pkin(4) * t811) - g(1) * (pkin(3) * t82 + pkin(4) * t444);];
taug_reg = t1;