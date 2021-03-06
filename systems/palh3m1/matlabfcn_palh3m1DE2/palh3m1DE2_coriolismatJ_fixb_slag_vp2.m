% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m1DE2
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
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m1DE2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1DE2_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 21:25:41
% EndTime: 2020-04-20 01:17:33
% DurationCPUTime: 6557.61s
% Computational Cost: add. (171181980->1041), mult. (259117204->1885), div. (10993188->28), fcn. (164483540->34), ass. (0->767)
t1001 = (-pkin(8) - pkin(10));
t1045 = 2 * pkin(3);
t517 = pkin(4) ^ 2;
t516 = pkin(5) ^ 2;
t521 = pkin(1) ^ 2;
t504 = sin(qJ(2));
t505 = sin(pkin(16));
t946 = cos(qJ(2));
t947 = cos(pkin(16));
t482 = t504 * t505 - t946 * t947;
t930 = pkin(5) * t482;
t764 = -0.2e1 * pkin(1) * t930 + t521;
t472 = t516 + t764;
t466 = 0.1e1 / t472;
t520 = 0.1e1 / pkin(2);
t793 = t466 * t520;
t759 = pkin(2) ^ 2 - pkin(6) ^ 2;
t464 = t472 + t759;
t473 = pkin(1) - t930;
t483 = t504 * t947 + t946 * t505;
t1002 = pkin(5) + pkin(6);
t1003 = pkin(5) - pkin(6);
t461 = (pkin(2) + t1002) * (-pkin(2) + t1003) + t764;
t462 = (-pkin(2) + t1002) * (pkin(2) + t1003) + t764;
t796 = t462 * t461;
t523 = sqrt(-t796);
t448 = pkin(5) * t483 * t464 + t473 * t523;
t503 = sin(qJ(3));
t801 = t448 * t503;
t787 = t483 * t523;
t447 = -pkin(5) * t787 + t464 * t473;
t508 = cos(qJ(3));
t804 = t447 * t508;
t438 = (-t804 / 0.2e1 + t801 / 0.2e1) * t793;
t800 = t448 * t508;
t805 = t447 * t503;
t439 = (t800 / 0.2e1 + t805 / 0.2e1) * t793;
t493 = pkin(18) + pkin(19);
t489 = sin(t493);
t490 = cos(t493);
t599 = t438 * t490 + t439 * t489;
t934 = pkin(4) * t599;
t766 = t934 * t1045 + t517;
t390 = ((pkin(3) - t1001) * (pkin(3) + t1001)) + t766;
t1000 = (pkin(10) - pkin(8));
t391 = ((pkin(3) - t1000) * (pkin(3) + t1000)) + t766;
t830 = t391 * t390;
t522 = sqrt(-t830);
t1040 = t522 * t599;
t752 = 0.2e1 * pkin(4);
t591 = pkin(3) * (t390 + t391) * t752;
t598 = t438 * t489 - t439 * t490;
t366 = t598 * t591;
t376 = 0.1e1 / t522;
t959 = -t598 / 0.2e1;
t686 = t376 * t959;
t565 = t366 * t686 - t1040;
t518 = pkin(3) ^ 2;
t397 = t518 + t766;
t760 = pkin(8) ^ 2 - pkin(10) ^ 2;
t393 = t397 - t760;
t399 = pkin(3) * t599 + pkin(4);
t661 = -0.2e1 * pkin(4) * t399 - t393;
t322 = (t598 * t661 + t565) * pkin(3);
t1039 = t598 * t522;
t1017 = -0.2e1 * t598;
t931 = pkin(4) * t518;
t663 = t931 * t1017;
t965 = t376 / 0.2e1;
t687 = t399 * t965;
t324 = t366 * t687 + t598 * t663 + (t393 * t599 - t1039) * pkin(3);
t394 = 0.1e1 / t397;
t511 = 0.1e1 / pkin(10);
t1049 = t393 * t598;
t359 = pkin(3) * t1049 + t399 * t522;
t496 = sin(pkin(17));
t837 = t359 * t496;
t357 = -pkin(3) * t1039 + t393 * t399;
t497 = cos(pkin(17));
t838 = t357 * t497;
t603 = t837 - t838;
t395 = 0.1e1 / t397 ^ 2;
t940 = pkin(3) * t395;
t749 = pkin(4) * t940;
t557 = t603 * t749;
t957 = -t497 / 0.2e1;
t958 = t496 / 0.2e1;
t290 = ((t322 * t957 + t324 * t958) * t394 + t598 * t557) * t511;
t836 = t359 * t497;
t839 = t357 * t496;
t602 = t836 + t839;
t556 = t602 * t749;
t956 = t497 / 0.2e1;
t291 = ((t322 * t958 + t324 * t956) * t394 + t598 * t556) * t511;
t827 = t394 * t511;
t331 = (-t838 / 0.2e1 + t837 / 0.2e1) * t827;
t328 = 0.1e1 / t331 ^ 2;
t332 = (t836 / 0.2e1 + t839 / 0.2e1) * t827;
t330 = t332 ^ 2;
t314 = t328 * t330 + 0.1e1;
t312 = 0.1e1 / t314;
t849 = t328 * t332;
t720 = t312 * t849;
t327 = 0.1e1 / t331;
t852 = t312 * t327;
t240 = -t290 * t720 + t291 * t852;
t239 = t240 + 0.1e1;
t507 = cos(qJ(4));
t502 = sin(qJ(4));
t918 = Ifges(6,4) * t502;
t614 = Ifges(6,2) * t507 + t918;
t200 = t614 * t239;
t1058 = t200 / 0.4e1;
t868 = t239 * t507;
t238 = Ifges(6,4) * t868;
t869 = t239 * t502;
t202 = Ifges(6,1) * t869 + t238;
t1057 = -t202 / 0.4e1;
t911 = Ifges(6,6) * t507;
t915 = Ifges(6,5) * t502;
t1053 = -t915 / 0.2e1 - t911 / 0.2e1;
t329 = t327 * t328;
t848 = t329 * t330;
t247 = -0.2e1 * t290 * t848 + 0.2e1 * t291 * t849;
t595 = -t801 + t804;
t1012 = pkin(1) * pkin(5);
t467 = 0.1e1 / t472 ^ 2;
t750 = t467 * t1012;
t654 = t483 * t750;
t479 = t483 ^ 2;
t942 = pkin(1) * t516;
t732 = t479 * t942;
t1037 = 0.2e1 * t1012 * (t461 + t462);
t452 = t483 * t1037;
t454 = 0.1e1 / t523;
t797 = t454 * t452;
t426 = t473 * t797 / 0.2e1 - 0.2e1 * t732 + (-t482 * t464 - t787) * pkin(5);
t813 = t426 * t503;
t679 = -t797 / 0.2e1;
t789 = t482 * t523;
t928 = t473 * pkin(1);
t424 = (t789 + (-t464 + t679 - 0.2e1 * t928) * t483) * pkin(5);
t816 = t424 * t508;
t384 = ((-t813 / 0.2e1 + t816 / 0.2e1) * t466 + t595 * t654) * t520;
t594 = t800 + t805;
t812 = t426 * t508;
t817 = t424 * t503;
t385 = ((t812 / 0.2e1 + t817 / 0.2e1) * t466 + t594 * t654) * t520;
t373 = t384 * t490 - t385 * t489;
t600 = -t384 * t489 - t385 * t490;
t1038 = t600 * t522;
t342 = t373 * t591;
t566 = t342 * t686 - t1038;
t300 = (t661 * t373 + t566) * pkin(3);
t833 = t373 * t522;
t302 = t342 * t687 + t373 * t663 + (t393 * t600 - t833) * pkin(3);
t273 = ((t300 * t957 + t302 * t958) * t394 + t373 * t557) * t511;
t966 = -t366 / 0.2e1;
t967 = -t342 / 0.2e1;
t548 = (t373 * t966 + t598 * t967) * t376 - t1038;
t1043 = t342 * t366;
t835 = t373 * t598;
t717 = t518 * t835;
t636 = t517 * t717;
t325 = t591 * t600 - 0.8e1 * t636;
t377 = t376 / t830;
t963 = t377 / 0.4e1;
t563 = t963 * t1043 + t325 * t965;
t601 = t600 * t1017;
t272 = t563 * t399 + (-t373 * t599 + t601) * t518 * t752 + (-t373 * t393 + t548) * pkin(3);
t274 = ((t300 * t958 + t302 * t956) * t394 + t373 * t556) * t511;
t964 = -t377 / 0.4e1;
t684 = t598 * t964;
t540 = (t325 * t959 + t599 * t967 + t600 * t966) * t376 + t833 + t684 * t1043;
t276 = 0.4e1 * pkin(4) * t717 + (t600 * t661 + t540) * pkin(3);
t558 = t302 * t598 + t324 * t373 + t359 * t600;
t559 = t300 * t598 + t322 * t373 + t357 * t600;
t396 = t394 * t395;
t1044 = 0.4e1 * t396;
t567 = t636 * t1044;
t552 = (-((t272 * t958 + t276 * t957) * t394 + t603 * t567 + (t496 * t558 - t497 * t559) * t749) * t511 * t332 - t273 * t291 - t274 * t290) * t312;
t313 = 0.1e1 / t314 ^ 2;
t851 = t313 * t327;
t723 = t274 * t851;
t847 = t329 * t332;
t722 = t290 * t847;
t748 = 0.2e1 * t312;
t769 = ((t272 * t956 + t276 * t958) * t394 + t602 * t567 + (t496 * t559 + t497 * t558) * t749) * t511 * t852 + t273 * t722 * t748;
t850 = t313 * t332;
t71 = -t247 * t723 + (t247 * t273 * t850 + t552) * t328 + t769;
t1056 = t1053 * t71;
t229 = -t273 * t720 + t274 * t852;
t228 = t229 + 0.1e1;
t989 = t228 / 0.2e1;
t953 = -t502 / 0.4e1;
t1055 = t507 / 0.4e1;
t168 = t614 * t228;
t872 = t228 * t507;
t227 = Ifges(6,4) * t872;
t873 = t228 * t502;
t170 = Ifges(6,1) * t873 + t227;
t1047 = t170 * t1055 + t168 * t953;
t1054 = Ifges(5,5) * t989 + t1047;
t234 = -0.2e1 * t273 * t848 + 0.2e1 * t274 * t849;
t721 = t291 * t851;
t64 = -t234 * t721 + (t234 * t290 * t850 + t552) * t328 + t769;
t21 = t614 * t64;
t917 = Ifges(6,4) * t507;
t616 = Ifges(6,1) * t502 + t917;
t22 = t616 * t64;
t27 = t614 * t71;
t28 = t616 * t71;
t315 = atan2(t332, t331);
t310 = sin(t315);
t311 = cos(t315);
t582 = pkin(4) * (mrSges(5,1) * t310 + mrSges(5,2) * t311);
t564 = t229 * t582;
t944 = pkin(1) * t503;
t307 = t310 * t944;
t943 = pkin(1) * t508;
t486 = pkin(4) - t943;
t285 = t311 * t486 + t307;
t177 = t285 * t240 - t311 * t943 + t307;
t306 = t310 * t486;
t853 = t311 * t503;
t178 = (t239 * t853 + t310 * t508) * pkin(1) - t240 * t306;
t622 = t178 * mrSges(5,1) - t177 * mrSges(5,2);
t935 = pkin(4) * t311;
t51 = pkin(11) * t64 + t229 * t935;
t494 = t502 ^ 2;
t495 = t507 ^ 2;
t761 = t494 + t495;
t655 = t761 * t51;
t40 = pkin(11) * t71 + t177;
t656 = t761 * t40;
t925 = t64 - t71;
t1052 = (t656 - t655) * mrSges(6,3) - t925 * Ifges(5,3) - (t22 / 0.2e1 - t28 / 0.2e1) * t502 - (-t27 / 0.2e1 + t21 / 0.2e1) * t507 + t622 + t564;
t500 = cos(pkin(19));
t802 = t448 * t500;
t498 = sin(pkin(19));
t807 = t447 * t498;
t596 = t802 + t807;
t814 = t426 * t500;
t819 = t424 * t498;
t380 = ((t814 / 0.2e1 + t819 / 0.2e1) * t466 + t596 * t654) * t520;
t803 = t448 * t498;
t806 = t447 * t500;
t597 = -t803 + t806;
t815 = t426 * t498;
t818 = t424 * t500;
t381 = ((-t818 / 0.2e1 + t815 / 0.2e1) * t466 - t597 * t654) * t520;
t435 = (-t806 / 0.2e1 + t803 / 0.2e1) * t793;
t428 = 0.1e1 / t435 ^ 2;
t436 = (t802 / 0.2e1 + t807 / 0.2e1) * t793;
t430 = t436 ^ 2;
t415 = t428 * t430 + 0.1e1;
t411 = 0.1e1 / t415;
t811 = t428 * t436;
t715 = t411 * t811;
t427 = 0.1e1 / t435;
t820 = t411 * t427;
t346 = t380 * t820 - t381 * t715;
t345 = t346 + 0.1e1;
t1051 = t381 * t427 * t428;
t392 = t397 + t760;
t1050 = t392 * t598;
t612 = t911 + t915;
t746 = 0.2e1 * t483;
t792 = t479 * t521;
t714 = t516 * t792;
t634 = t466 * t467 * t714;
t590 = 0.4e1 * t634;
t990 = -t228 / 0.4e1;
t398 = -pkin(3) - t934;
t356 = -pkin(4) * t1039 - t392 * t398;
t358 = pkin(4) * t1050 - t398 * t522;
t513 = 0.1e1 / pkin(8);
t960 = t394 / 0.2e1;
t680 = t513 * t960;
t334 = atan2(t358 * t680, t356 * t680);
t333 = cos(t334);
t499 = sin(pkin(18));
t501 = cos(pkin(18));
t945 = sin(t334);
t309 = -t501 * t333 - t499 * t945;
t418 = atan2(t436, t435);
t409 = sin(t418);
t410 = cos(t418);
t378 = -t409 * t504 + t410 * t946;
t379 = t409 * t946 + t410 * t504;
t573 = -t499 * t333 + t501 * t945;
t280 = t309 * t378 + t379 * t573;
t975 = t280 / 0.2e1;
t281 = t309 * t379 - t378 * t573;
t974 = t281 / 0.2e1;
t480 = t503 * t504 - t508 * t946;
t481 = t503 * t946 + t508 * t504;
t283 = -t310 * t481 - t311 * t480;
t972 = t283 / 0.2e1;
t1035 = t310 * t480 - t311 * t481;
t970 = -t1035 / 0.2e1;
t968 = t1035 / 0.2e1;
t584 = pkin(1) * t853 - t306;
t218 = pkin(11) * t228 - t584;
t648 = t218 * t761;
t936 = pkin(4) * t310;
t231 = pkin(11) * t239 + t936;
t647 = t231 * t761;
t1042 = t392 * t599;
t1041 = t392 * t600;
t159 = t584 * t229;
t726 = t71 / 0.4e1 - t64 / 0.4e1;
t1036 = t726 * t283;
t860 = t1035 * t507;
t261 = mrSges(6,1) * t283 - mrSges(6,3) * t860;
t773 = t507 * t261;
t861 = t1035 * t502;
t260 = -mrSges(6,2) * t283 - mrSges(6,3) * t861;
t780 = t502 * t260;
t1034 = -t780 / 0.2e1 - t773 / 0.2e1;
t889 = t507 * mrSges(6,1);
t891 = t502 * mrSges(6,2);
t1033 = -t891 / 0.2e1 + t889 / 0.2e1;
t901 = t281 * Ifges(9,5);
t902 = t280 * Ifges(9,6);
t253 = t901 + t902;
t1032 = t902 / 0.2e1 + t253 / 0.2e1;
t912 = Ifges(6,6) * t502;
t914 = Ifges(6,5) * t507;
t1031 = t914 / 0.2e1 - t912 / 0.2e1;
t888 = t507 * mrSges(6,2);
t892 = t502 * mrSges(6,1);
t619 = t888 + t892;
t167 = t619 * t228;
t198 = t619 * t239;
t217 = -pkin(9) * t228 - t285;
t232 = -pkin(9) * t239 - t935;
t1028 = -t232 * t167 / 0.2e1 - t217 * t198 / 0.2e1;
t201 = -Ifges(6,2) * t869 + t238;
t617 = Ifges(6,1) * t507 - t918;
t203 = t617 * t239;
t1026 = (t203 / 0.2e1 - t200 / 0.2e1) * t502 + (t202 / 0.2e1 + t201 / 0.2e1) * t507;
t169 = -Ifges(6,2) * t873 + t227;
t171 = t617 * t228;
t993 = t170 / 0.2e1;
t995 = -t168 / 0.2e1;
t1025 = (t171 / 0.2e1 + t995) * t502 + (t993 + t169 / 0.2e1) * t507;
t352 = 0.1e1 / t356;
t515 = 0.1e1 / pkin(6);
t794 = t466 * t515;
t463 = t472 - t759;
t474 = pkin(1) * t482 - pkin(5);
t449 = pkin(1) * t483 * t463 - t474 * t523;
t506 = sin(pkin(15));
t799 = t449 * t506;
t446 = -pkin(1) * t787 - t463 * t474;
t509 = cos(pkin(15));
t808 = t446 * t509;
t437 = (t808 / 0.2e1 + t799 / 0.2e1) * t794;
t431 = 0.1e1 / t437;
t353 = 0.1e1 / t356 ^ 2;
t432 = 0.1e1 / t437 ^ 2;
t1022 = 0.4e1 * t217;
t1021 = 0.4e1 * t232;
t1020 = -0.4e1 * t373;
t1019 = 0.2e1 * t373;
t1018 = -0.2e1 * t397;
t1016 = 0.2e1 * t598;
t1015 = 0.2e1 * t496;
t1014 = m(5) / 0.4e1;
t1013 = m(6) / 0.4e1;
t1011 = pkin(3) * pkin(4);
t1007 = t22 / 0.4e1;
t1006 = t51 / 0.2e1;
t1005 = -t64 / 0.2e1;
t574 = (0.1e1 / t796 * t452 ^ 2 / 0.4e1 - t482 * t1037 / 0.2e1 - 0.4e1 * t714) * t454;
t641 = -t464 - t797;
t739 = 0.6e1 * t482 * t483;
t389 = t739 * t942 + t574 * t473 + (t641 * t483 + t789) * pkin(5);
t951 = t503 / 0.2e1;
t553 = (t523 - t574) * t483;
t387 = 0.4e1 * t732 + (t553 + (-t641 + 0.2e1 * t928) * t482) * pkin(5);
t962 = -t387 / 0.2e1;
t363 = (-t595 * t590 + (t389 * t951 + t508 * t962) * t466 + ((t813 - t816) * t746 + t595 * t482) * t750) * t520;
t961 = t389 / 0.2e1;
t364 = (t594 * t590 + (t387 * t951 + t508 * t961) * t466 + ((t812 + t817) * t746 - t594 * t482) * t750) * t520;
t339 = -t363 * t490 - t364 * t489;
t369 = t373 ^ 2;
t772 = t517 * t518;
t737 = -0.8e1 * t772;
t316 = t339 * t591 + t369 * t737;
t338 = t363 * t489 - t364 * t490;
t341 = t342 ^ 2;
t576 = -t342 * t373 * t376 - t339 * t522;
t744 = t600 * t1020;
t268 = (t316 * t965 + t341 * t963) * t399 + (t339 * t1017 + t744) * t931 + (t338 * t393 + t576) * pkin(3);
t544 = (t316 * t959 - t342 * t600) * t376 - t338 * t522 + t341 * t684;
t751 = 0.4e1 * t931;
t269 = t369 * t751 + (t661 * t339 + t544) * pkin(3);
t716 = t396 * t772;
t644 = 0.4e1 * t716;
t627 = t496 * t644;
t587 = t369 * t627;
t635 = t497 * t716;
t588 = 0.4e1 * t359 * t635;
t589 = -0.4e1 * t357 * t635;
t653 = t312 * t749;
t828 = t394 * t497;
t681 = t828 / 0.2e1;
t682 = -t828 / 0.2e1;
t683 = t394 * t958;
t719 = t313 * t849;
t743 = t373 * t1015;
t834 = t373 * t497;
t74 = -t234 * t723 + (t234 * t719 + (t273 * t847 - t274 * t328) * t748) * t273 + (((t268 * t681 + t269 * t683 + t357 * t587 + t369 * t588) * t327 - (t268 * t683 + t269 * t682 + t359 * t587 + t369 * t589) * t849) * t312 + ((t300 * t743 + 0.2e1 * t302 * t834) * t327 - (-0.2e1 * t300 * t834 + t302 * t743) * t849 + (t602 * t327 - t603 * t849) * t339) * t653) * t511;
t1004 = t74 / 0.2e1;
t151 = t228 * t283;
t999 = t151 / 0.2e1;
t998 = -t151 / 0.2e1;
t152 = t228 * t1035;
t997 = -t152 / 0.2e1;
t996 = t152 / 0.2e1;
t994 = t169 / 0.4e1;
t175 = t239 * t283;
t992 = t175 / 0.2e1;
t991 = -t175 / 0.2e1;
t988 = -t231 / 0.2e1;
t987 = t232 / 0.2e1;
t986 = t239 / 0.2e1;
t985 = t239 / 0.4e1;
t982 = Ifges(9,4) * t974 + Ifges(9,2) * t975;
t981 = Ifges(9,1) * t974 + Ifges(9,4) * t975;
t979 = t260 / 0.2e1;
t978 = -t261 / 0.2e1;
t662 = t398 * t1045 - t392;
t321 = (t598 * t662 + t565) * pkin(4);
t651 = t598 * t749;
t317 = (t321 * t960 + t356 * t651) * t513;
t937 = pkin(3) * t517;
t664 = t937 * t1017;
t688 = -t376 * t398 / 0.2e1;
t323 = t366 * t688 + t598 * t664 + (-t1039 + t1042) * pkin(4);
t318 = (t323 * t960 + t358 * t651) * t513;
t842 = t353 * t358;
t577 = t317 * t842 - t318 * t352;
t747 = 0.2e1 * t397;
t355 = t358 ^ 2;
t337 = t353 * t355 + 0.1e1;
t335 = 0.1e1 / t337;
t929 = pkin(8) * t335;
t270 = t577 * t747 * t929;
t976 = -t270 / 0.4e1;
t973 = -t283 / 0.2e1;
t971 = t283 / 0.4e1;
t969 = -t1035 / 0.4e1;
t955 = t498 / 0.2e1;
t954 = -t502 / 0.2e1;
t952 = t502 / 0.2e1;
t950 = -t507 / 0.4e1;
t949 = t507 / 0.2e1;
t948 = t509 / 0.2e1;
t941 = pkin(3) * t394;
t939 = pkin(3) * t499;
t938 = pkin(3) * t501;
t933 = pkin(4) * t481;
t932 = pkin(4) * t513;
t927 = t474 * pkin(5);
t926 = -Ifges(4,1) + Ifges(4,2);
t924 = t612 * t74;
t923 = pkin(1) * qJD(2);
t922 = pkin(4) * qJD(3);
t921 = mrSges(6,3) * t218;
t920 = mrSges(6,3) * t502;
t919 = mrSges(6,3) * t507;
t176 = t239 * t1035;
t913 = Ifges(5,6) * t176;
t235 = t309 * t270;
t236 = t573 * t270;
t191 = t235 * t379 - t236 * t378;
t910 = Ifges(9,6) * t191;
t909 = Ifges(6,3) * t176;
t908 = Ifges(9,3) * t270;
t299 = (t662 * t373 + t566) * pkin(4);
t652 = t373 * t749;
t294 = (t299 * t960 + t356 * t652) * t513;
t301 = t342 * t688 + t373 * t664 + (-t833 + t1041) * pkin(4);
t295 = (t301 * t960 + t358 * t652) * t513;
t841 = t353 * t397;
t734 = pkin(8) * t841;
t637 = t358 * t734;
t625 = t335 * t637;
t592 = -0.2e1 * t625;
t735 = pkin(8) * t352 * t397;
t640 = t335 * t735;
t251 = t294 * t592 + 0.2e1 * t295 * t640;
t225 = t309 * t251;
t226 = t573 * t251;
t319 = t345 * t378;
t320 = t345 * t379;
t146 = t225 * t378 + t226 * t379 + t309 * t319 + t320 * t573;
t907 = t146 * Ifges(9,5);
t147 = -t225 * t379 + t226 * t378 - t309 * t320 + t319 * t573;
t906 = t147 * Ifges(9,6);
t905 = t151 * Ifges(5,5);
t904 = t152 * Ifges(5,6);
t615 = -Ifges(6,2) * t502 + t917;
t898 = t283 * Ifges(6,6);
t245 = t1035 * t615 + t898;
t246 = t283 * Ifges(6,5) + t1035 * t617;
t257 = t619 * t1035;
t897 = t1035 * Ifges(5,5);
t899 = t283 * Ifges(5,6);
t263 = t897 - t899;
t41 = -pkin(9) * t71 - t178;
t733 = t229 * t936;
t50 = -pkin(9) * t64 + t733;
t874 = t228 * t176;
t695 = -t874 / 0.4e1;
t555 = t695 - t1036;
t570 = t217 * t619;
t199 = t612 * t239;
t95 = t619 * t151;
t579 = t152 * t199 / 0.4e1 - t95 * t987;
t96 = -mrSges(6,2) * t152 + t151 * t920;
t580 = t151 * t1057 + t231 * t96 / 0.2e1;
t97 = mrSges(6,1) * t152 + t151 * t919;
t581 = t151 * t1058 + t97 * t988;
t98 = -t904 - t905;
t629 = -t905 / 0.4e1 + t98 / 0.4e1;
t67 = Ifges(6,6) * t152 - t615 * t151;
t354 = t352 * t353;
t840 = t354 * t355;
t287 = -0.2e1 * t321 * t840 + 0.2e1 * t323 * t842;
t718 = t354 * t358 * t397;
t645 = 0.4e1 * t718;
t336 = 0.1e1 / t337 ^ 2;
t646 = t336 * t747;
t698 = 0.4e1 * t358 * t1011;
t699 = -0.4e1 * t352 * t1011;
t745 = -0.2e1 * t841;
t742 = t518 * t1044;
t554 = t513 * t517 * (t356 * t742 + 0.2e1 * t941);
t626 = t513 * t640;
t593 = 0.2e1 * t626;
t829 = t394 * t398;
t768 = (t554 * t835 + ((t540 - t1041) * t960 + (t600 * t829 + (t299 * t598 + t321 * t373 + t356 * t600) * t395) * pkin(3)) * t932) * t592 + ((-t398 * t563 + 0.2e1 * t601 * t937) * t960 + (t358 * t598 * t742 - t599 * t941) * t517 * t373 + ((-t373 * t392 + t548) * t960 + (t301 * t598 + t323 * t373 + t358 * t600) * t940) * pkin(4)) * t593;
t123 = ((t294 * t842 - t295 * t352) * t287 * t646 + ((t321 * t745 + t598 * t699) * t295 + (t321 * t645 + (t323 * t1018 + t598 * t698) * t353) * t294) * t335) * pkin(8) + t768;
t275 = -0.2e1 * t299 * t840 + 0.2e1 * t301 * t842;
t124 = (t577 * t275 * t646 + ((t299 * t745 + t373 * t699) * t318 + (t299 * t645 + (t301 * t1018 + t373 * t698) * t353) * t317) * t335) * pkin(8) + t768;
t671 = t124 / 0.4e1 - t123 / 0.4e1;
t68 = Ifges(6,5) * t152 - t617 * t151;
t689 = t310 * t997;
t728 = t1006 - t40 / 0.2e1;
t77 = t176 * Ifges(6,6) - t615 * t175;
t78 = t176 * Ifges(6,5) - t617 * t175;
t859 = t584 * t176;
t343 = pkin(1) * t409 + t345 * t939;
t344 = pkin(1) * t410 + t345 * t938;
t278 = t309 * t343 - t344 * t573;
t862 = t278 * t191;
t190 = -t235 * t378 - t236 * t379;
t277 = t309 * t344 + t343 * t573;
t863 = t277 * t190;
t250 = t345 + t251;
t864 = t250 * t191;
t865 = t250 * t190;
t119 = mrSges(6,1) * t176 + t175 * t919;
t875 = t218 * t119;
t118 = -mrSges(6,2) * t176 + t175 * t920;
t876 = t218 * t118;
t189 = t235 * t343 - t236 * t344;
t879 = t189 * t281;
t188 = -t235 * t344 - t236 * t343;
t880 = t188 * t280;
t881 = t177 * t283;
t896 = t285 * mrSges(5,3);
t90 = t906 + t907;
t2 = t90 * t976 - t726 * t263 + (t50 / 0.2e1 - t41 / 0.2e1) * t257 + t671 * t253 + t629 * t239 + (-t864 / 0.2e1 + t147 * t976 + t671 * t280) * Ifges(9,6) + (-t152 * t239 / 0.4e1 + t874 / 0.2e1 + t1036) * Ifges(5,6) + (t879 / 0.2e1 + t863 / 0.2e1 - t862 / 0.2e1 - t880 / 0.2e1) * mrSges(9,3) + (-t859 / 0.2e1 + t881 / 0.2e1 + (t689 + (t229 * t973 + t999) * t311) * pkin(4)) * mrSges(5,3) + (-t876 / 0.2e1 + t77 * t990 + t67 * t985 + t728 * t260 - t726 * t245 + t555 * Ifges(6,6) + t580) * t507 + (t78 * t990 + t68 * t985 + t875 / 0.2e1 - t728 * t261 - t726 * t246 + t555 * Ifges(6,5) + t581) * t502 + ((t1007 - t28 / 0.4e1) * t507 + (t27 / 0.4e1 - t21 / 0.4e1) * t502 - t726 * Ifges(5,5) + (t733 / 0.2e1 + t178 / 0.2e1) * mrSges(5,3)) * t1035 + (-t865 / 0.2e1 + t146 * t976 + t671 * t281) * Ifges(9,5) - (-t570 / 0.2e1 + t896 / 0.2e1 - t1054) * t175 + t579;
t903 = t2 * qJD(1);
t484 = -t946 * pkin(1) - pkin(13);
t465 = -t480 * pkin(4) + t484;
t256 = t283 * pkin(9) - pkin(11) * t1035 + t465;
t423 = (t789 + (-t463 + t679 + 0.2e1 * t927) * t483) * pkin(1);
t731 = pkin(5) * t792;
t425 = t474 * t679 - 0.2e1 * t731 + (-t463 * t482 - t787) * pkin(1);
t383 = ((t423 * t948 + t425 * t506 / 0.2e1) * t466 + (t799 + t808) * t654) * t515;
t798 = t449 * t509;
t809 = t446 * t506;
t440 = (t798 / 0.2e1 - t809 / 0.2e1) * t794;
t434 = t440 ^ 2;
t419 = t432 * t434 + 0.1e1;
t416 = 0.1e1 / t419;
t810 = t432 * t440;
t382 = ((t425 * t948 - t423 * t506 / 0.2e1) * t466 + (t798 - t809) * t654) * t515;
t832 = t382 * t431;
t348 = (-t383 * t810 + t832) * t416;
t367 = (-t378 * t501 - t379 * t499) * pkin(3) + t484;
t420 = atan2(t440, t437);
t413 = sin(t420);
t414 = cos(t420);
t460 = -mrSges(4,1) * t481 + mrSges(4,2) * t480;
t478 = Ifges(4,4) * t480;
t492 = t504 * pkin(1);
t613 = -t912 + t914;
t569 = t283 * t613;
t775 = t507 * t246;
t782 = t502 * t245;
t539 = t465 * mrSges(5,2) + t569 / 0.2e1 + Ifges(5,1) * t968 - t283 * Ifges(5,4) + t775 / 0.2e1 - t782 / 0.2e1;
t549 = t465 * mrSges(5,1) + Ifges(5,4) * t970 + t613 * t968 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t283 + (Ifges(5,2) + Ifges(6,3)) * t972;
t632 = 0.4e1 * t761;
t583 = t632 * t1013;
t606 = t773 + t780;
t659 = m(5) * t465 + mrSges(5,1) * t283 + mrSges(5,2) * t1035;
t660 = t492 - t933;
t94 = pkin(9) * t152 + pkin(11) * t151 + t660;
t5 = -Ifges(4,4) * t481 ^ 2 + t367 * (-mrSges(9,1) * t147 + mrSges(9,2) * t146) - t320 * Ifges(8,2) * t378 + t319 * t379 * Ifges(8,1) + (Ifges(9,1) * t146 + Ifges(9,4) * t147) * t974 + t146 * t981 + (Ifges(9,4) * t146 + Ifges(9,2) * t147) * t975 + t147 * t982 + t606 * t94 + (mrSges(8,1) * t320 + mrSges(8,2) * t319 + t460) * t484 + t659 * t660 + (m(9) * t367 - mrSges(9,1) * t280 + mrSges(9,2) * t281) * (t492 + (-t319 * t499 + t320 * t501) * pkin(3)) + t549 * t152 - t539 * t151 + (Ifges(5,1) * t998 + Ifges(5,4) * t997 + t67 * t954 + t68 * t949) * t1035 + (t502 * t96 + t507 * t97 + t583 * t94) * t256 + (t926 * t481 + t478) * t480 + (t319 * t378 - t320 * t379) * Ifges(8,4) + (-pkin(13) * mrSges(3,1) - Ifges(3,4) * t504 + (-mrSges(4,1) * t480 - mrSges(8,1) * t378 - mrSges(4,2) * t481 + mrSges(8,2) * t379 + (m(4) + m(8)) * t484) * pkin(1)) * t504 + (-pkin(13) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t504 + Ifges(3,4) * t946) * t946 + ((pkin(7) * mrSges(7,2) + Ifges(7,4) * t414) * t414 + (pkin(7) * mrSges(7,1) - Ifges(7,4) * t413 + (Ifges(7,1) - Ifges(7,2)) * t414) * t413) * t348;
t893 = t5 * qJD(1);
t890 = t502 * t68;
t887 = t507 * t67;
t117 = pkin(9) * t176 + pkin(11) * t175 - t933;
t9 = t484 * t460 + t480 * t478 + t367 * (-mrSges(9,1) * t191 + mrSges(9,2) * t190) + (Ifges(9,1) * t190 + Ifges(9,4) * t191) * t974 + t190 * t981 + (Ifges(9,4) * t190 + Ifges(9,2) * t191) * t975 + t191 * t982 + t549 * t176 - t539 * t175 + (t78 * t949 + t77 * t954 + Ifges(5,1) * t991 - Ifges(5,4) * t176 / 0.2e1) * t1035 + (t118 * t502 + t119 * t507) * t256 + (t256 * t583 + t606) * t117 + (-t659 * pkin(4) - Ifges(4,4) * t481 + t926 * t480) * t481;
t886 = t9 * qJD(1);
t401 = t598 ^ 2;
t360 = t401 * t737 + t591 * t599;
t365 = t366 ^ 2;
t575 = -t366 * t376 * t598 - t1040;
t741 = -0.4e1 * t599 * t598;
t297 = (t360 * t965 + t365 * t963) * t399 + (t1017 * t599 + t741) * t931 + (t575 - t1049) * pkin(3);
t543 = (t360 * t959 - t366 * t599) * t376 + t1039 + t365 * t684;
t298 = t401 * t751 + (t599 * t661 + t543) * pkin(3);
t586 = t401 * t627;
t740 = t598 * t1015;
t823 = t598 * t497;
t824 = t599 * t497;
t825 = t599 * t496;
t101 = -t247 * t721 + (t247 * t719 + (-t291 * t328 + t722) * t748) * t290 + (((t297 * t681 + t298 * t683 + t357 * t586 + t401 * t588) * t327 - (t297 * t683 + t298 * t682 + t359 * t586 + t401 * t589) * t849) * t312 + ((t322 * t740 + 0.2e1 * t324 * t823 + t357 * t825 + t359 * t824) * t327 - (-0.2e1 * t322 * t823 + t324 * t740 - t357 * t824 + t359 * t825) * t849) * t653) * t511;
t885 = t612 * t101;
t884 = qJD(2) * t74;
t258 = t614 * t1035;
t259 = t616 * t1035;
t620 = t889 - t891;
t571 = t217 * t620;
t665 = t495 / 0.2e1 + t494 / 0.2e1;
t709 = Ifges(6,3) * t996;
t779 = t507 * t168;
t784 = t502 * t170;
t10 = t709 + (Ifges(6,5) * t998 + t94 * mrSges(6,1) / 0.2e1 + t218 * t261 / 0.2e1) * t507 + (Ifges(6,6) * t999 - t94 * mrSges(6,2) / 0.2e1 + t218 * t979) * t502 + (-t569 / 0.4e1 + t782 / 0.4e1 - t258 * t950 - t775 / 0.4e1 - t259 * t953) * t228 + (-t571 / 0.2e1 + t779 / 0.4e1 + t502 * t994 + t784 / 0.4e1 + t171 * t950 + t665 * t921) * t1035;
t883 = t10 * qJD(1);
t237 = Ifges(6,5) * t868;
t777 = t507 * t200;
t783 = t502 * t202;
t538 = (t620 * t987 - t777 / 0.4e1 + t201 * t953 - t783 / 0.4e1 + t203 * t1055 - t665 * t231 * mrSges(6,3)) * t1035 + t237 * t971;
t578 = -t898 / 0.4e1 - t245 / 0.4e1 - t259 / 0.4e1;
t666 = -t258 / 0.4e1 + t246 / 0.4e1;
t13 = -t909 / 0.2e1 + (t260 * t988 + Ifges(6,6) * t991 + t117 * mrSges(6,2) / 0.2e1 + t578 * t239) * t502 + (t231 * t978 + Ifges(6,5) * t992 - t117 * mrSges(6,1) / 0.2e1 + t666 * t239) * t507 + t538;
t882 = t13 * qJD(1);
t878 = t217 * t167;
t870 = t232 * t198;
t867 = t240 * t310;
t866 = t240 * t311;
t858 = t310 * t1035;
t855 = t311 * t283;
t846 = t346 * t409;
t845 = t346 * t410;
t844 = 0.2e1 * (-t430 * t1051 + t380 * t811) / t415 ^ 2;
t831 = t383 * t431 * t432;
t843 = 0.2e1 * (t382 * t810 - t434 * t831) / t419 ^ 2;
t795 = t466 * t506;
t791 = t482 * t506;
t790 = t482 * t509;
t788 = t483 * t506;
t781 = t502 * t246;
t776 = t507 * t245;
t774 = t507 * t260;
t572 = t781 / 0.2e1 + t776 / 0.2e1;
t84 = (-t502 * t261 + t774) * t256 + (-t258 * t954 - t259 * t949 + t612 * t973 - t572) * t1035;
t771 = t84 * qJD(1);
t770 = t123 - t124;
t158 = t285 * t229;
t763 = Ifges(4,5) * t480 + Ifges(4,6) * t481;
t762 = mrSges(4,1) * t944 + mrSges(4,2) * t943;
t292 = -t380 * t427 * t844 + (-(-t597 * t590 + (t389 * t955 + t500 * t962) * t466 + ((t815 - t818) * t746 + t597 * t482) * t750) * t715 + (t596 * t590 + (t387 * t955 + t500 * t961) * t466 + ((t814 + t819) * t746 - t596 * t482) * t750) * t820) * t520 + (t811 * t844 + 0.2e1 * (t436 * t1051 - t380 * t428) * t411) * t381;
t696 = t929 * t1011;
t607 = t696 * t842;
t623 = t336 * t637;
t624 = t718 * t929;
t628 = t358 * t644;
t639 = t335 * t734;
t630 = -0.2e1 * t639;
t631 = t352 * t696;
t638 = t336 * t735;
t685 = t398 * t964;
t103 = ((t316 * t688 + t339 * t664 + t341 * t685 + t744 * t937) * t960 + t369 * t628 + ((t338 * t392 + t576) * t960 + (t301 * t1019 + t339 * t358) * t940) * pkin(4)) * t593 + (t369 * t554 + ((-t339 * t392 + t544) * t960 + (t299 * t395 * t1019 + (t356 * t395 + t829) * t339) * pkin(3)) * t932) * t592 + t292 + (t631 * t1020 - 0.2e1 * t275 * t638 + t299 * t630) * t295 + (0.2e1 * t275 * t623 + 0.4e1 * t299 * t624 + t301 * t630 + 0.4e1 * t373 * t607) * t294;
t758 = qJD(2) * t103;
t757 = qJD(2) * t228;
t756 = qJD(3) * t101;
t755 = qJD(3) * t239;
t754 = qJD(3) * t502;
t753 = qJD(3) * t507;
t738 = t509 * t746;
t727 = t1005 + t71 / 0.2e1;
t702 = t40 * t954;
t700 = -t888 / 0.2e1;
t694 = -t873 / 0.4e1;
t693 = t873 / 0.4e1;
t692 = t872 / 0.4e1;
t691 = -t861 / 0.4e1;
t690 = t860 / 0.4e1;
t678 = t466 * t948;
t36 = pkin(11) * t74 + t158;
t658 = t36 * t761;
t86 = pkin(4) * t866 + pkin(11) * t101;
t657 = t86 * t761;
t650 = mrSges(6,1) * t702 + t40 * t700 - t1056;
t649 = -t1053 * t64 + (-t892 / 0.2e1 + t700) * t51;
t642 = -t463 - t797;
t618 = t189 * mrSges(9,1) - t188 * mrSges(9,2);
t197 = t620 * t239;
t20 = t620 * t64;
t611 = -t50 * t197 - t232 * t20;
t610 = t506 * t634;
t609 = t218 * t632;
t608 = t231 * t632;
t605 = -t285 * t310 - t311 * t584;
t585 = t509 * t590;
t568 = t613 * t228;
t562 = mrSges(9,1) * t277 - mrSges(9,2) * t278 + Ifges(9,3) * t250;
t560 = t263 / 0.2e1 + t897 / 0.2e1 - t899 / 0.2e1;
t166 = t620 * t228;
t26 = t620 * t71;
t551 = -m(9) * (t188 * t278 + t189 * t277) + t217 * t26 + t41 * t166 - t762;
t550 = t171 * t869 / 0.4e1 + t868 * t994 + t200 * t694 + t203 * t693 - t1028 + (t201 + t202) * t692 + t1047 * t239;
t547 = t779 / 0.2e1 + t784 / 0.2e1 + t285 * mrSges(5,1) + t584 * mrSges(5,2);
t546 = t777 / 0.2e1 + t783 / 0.2e1 + (mrSges(5,1) * t311 - mrSges(5,2) * t310) * pkin(4);
t30 = t614 * t74;
t31 = t616 * t74;
t545 = t159 * mrSges(5,1) - t158 * mrSges(5,2) + Ifges(5,3) * t74 + t30 * t949 + t31 * t952;
t542 = mrSges(6,3) * t648 + t547;
t541 = mrSges(6,3) * t647 + t546;
t388 = t521 * pkin(5) * t739 - t574 * t474 + (t642 * t483 + t789) * pkin(1);
t386 = 0.4e1 * t731 + (t553 + (-t642 - 0.2e1 * t927) * t482) * pkin(1);
t293 = -t832 * t843 + (t810 * t843 + 0.2e1 * (-t382 * t432 + t440 * t831) * t416) * t383 + ((t388 * t678 + t449 * t585 - t386 * t795 / 0.2e1 - 0.4e1 * t446 * t610) * t431 - (t386 * t678 + t446 * t585 + t388 * t795 / 0.2e1 + 0.4e1 * t449 * t610) * t810 + ((-0.2e1 * t423 * t788 + t425 * t738 + t446 * t791 - t449 * t790) * t431 - (t423 * t738 + 0.2e1 * t425 * t788 - t446 * t790 - t449 * t791) * t810) * t750) * t416 * t515;
t289 = pkin(1) * t845 + t292 * t939;
t288 = -pkin(1) * t846 + t292 * t938;
t139 = 0.2e1 * ((t360 * t688 + t365 * t685 + t599 * t664 + t741 * t937) * t960 + t401 * t628 + ((t575 - t1050) * t960 + (t323 * t1016 + t358 * t599) * t940) * pkin(4)) * t626 - 0.2e1 * (t401 * t554 + ((t543 - t1042) * t960 + (t599 * t829 + (t321 * t1016 + t356 * t599) * t395) * pkin(3)) * t932) * t625 + 0.2e1 * (t631 * t1017 - t287 * t638 - t321 * t639) * t318 + 0.2e1 * (t607 * t1016 + t287 * t623 + 0.2e1 * t321 * t624 - t323 * t639) * t317;
t133 = -t225 * t343 + t226 * t344 + t288 * t309 + t289 * t573;
t132 = t225 * t344 + t226 * t343 - t288 * t573 + t289 * t309;
t102 = -0.4e1 * t177 * t584 + 0.4e1 * t178 * t285;
t85 = pkin(4) * t867 - pkin(9) * t101;
t47 = t616 * t101;
t46 = t614 * t101;
t45 = t620 * t101;
t37 = -pkin(9) * t74 - t159;
t32 = t1026 * t239 + t870;
t29 = t620 * t74;
t16 = t1025 * t228 + t878;
t15 = t50 * t1021 + t51 * t608;
t14 = t41 * t1022 + t40 * t609;
t12 = t909 / 0.2e1 + t1034 * t231 + (t578 * t502 + t666 * t507) * t239 + t538 - t1031 * t175 + t1033 * t117;
t11 = t571 * t968 + t568 * t971 + t245 * t694 - t259 * t693 + t779 * t969 + t171 * t690 + t709 + t1033 * t94 - t1031 * t151 + (-t258 + t246) * t692 + (t169 + t170) * t691 + t761 * t921 * t970 + t1034 * t218;
t8 = (mrSges(6,2) * t1006 + Ifges(6,6) * t1005) * t507 + (mrSges(6,1) * t1006 + Ifges(6,5) * t1005) * t502 + t550 + t650;
t7 = t1056 + (t888 / 0.2e1 + t892 / 0.2e1) * t40 + t550 + t649;
t6 = ((-t170 / 0.4e1 - t169 / 0.4e1) * t507 + (-t171 / 0.4e1 + t168 / 0.4e1) * t502) * t239 + ((t1057 - t201 / 0.4e1) * t507 + (-t203 / 0.4e1 + t1058) * t502) * t228 + t649 + t650 + t1028;
t4 = t15 * t1013 - t124 * t908 + t541 * t64 + (mrSges(6,3) * t655 + Ifges(5,3) * t64 + t21 * t949 + t22 * t952 - t564) * t239 + t611;
t3 = t102 * t1014 + t14 * t1013 + t618 * t250 + t562 * t123 + t542 * t71 + (mrSges(6,3) * t656 + t71 * Ifges(5,3) + t27 * t949 + t28 * t952 + t622) * t228 - t551;
t1 = (t50 + t41) * t257 / 0.2e1 + (t880 + t862) * mrSges(9,3) / 0.2e1 - (t879 + t863) * mrSges(9,3) / 0.2e1 + t27 * t691 + t77 * t692 + t78 * t693 + Ifges(5,6) * t695 + t28 * t690 + t579 - (t881 - t859) * mrSges(5,3) / 0.2e1 + (t1007 * t1035 + t51 * t979 + t580) * t507 + (t64 + t71) * (t612 * t971 + t781 / 0.4e1 + t897 / 0.4e1 - t899 / 0.4e1 + t776 / 0.4e1 + t263 / 0.4e1) + (t124 + t123) * (t901 / 0.4e1 + t902 / 0.4e1 + t253 / 0.4e1) + (t21 * t969 + t51 * t978 + t581) * t502 - (t90 / 0.4e1 + t906 / 0.4e1 + t907 / 0.4e1) * t270 + t250 * (Ifges(9,5) * t190 + t910) / 0.4e1 + (t887 / 0.4e1 + t890 / 0.4e1 - t904 / 0.4e1 + t629) * t239 + t612 * t874 / 0.4e1 + Ifges(9,6) * t864 / 0.4e1 + Ifges(9,5) * t865 / 0.4e1 + ((t311 * t999 + t689 + (t858 / 0.2e1 - t855 / 0.2e1) * t229) * pkin(4) + t178 * t970) * mrSges(5,3) + t570 * t991 + t896 * t992 + t40 * t774 / 0.2e1 + t876 * t949 + t875 * t954 - t1054 * t175 + t763 + t261 * t702 + t913 * t990;
t17 = [qJD(2) * t5 + qJD(3) * t9 + qJD(4) * t84, t893 + t1 * qJD(3) + t11 * qJD(4) + t560 * t884 + (t901 / 0.2e1 + t1032) * t758 + ((t480 * t508 - t481 * t503) * mrSges(4,3) + (-t319 * t410 - t320 * t409 + (t378 * t410 + t379 * t409) * t346) * mrSges(8,3)) * t923 + (t612 * t996 + t887 / 0.2e1 + t890 / 0.2e1 - t905 / 0.2e1 - t904 / 0.2e1 + t98 / 0.2e1) * t757 + (Ifges(3,5) * t946 - Ifges(3,6) * t504 - t217 * t95 + t37 * t257 + t924 * t972 + t763 + (t378 * t292 - t320 * t345) * Ifges(8,6) + t348 ^ 2 * (t414 * Ifges(7,5) - t413 * Ifges(7,6)) + (t90 / 0.2e1 + t907 / 0.2e1 + t906 / 0.2e1) * t250 + (t132 * t280 - t133 * t281 - t146 * t277 + t147 * t278) * mrSges(9,3) + (-t1035 * t159 + t151 * t285 + t152 * t584 - t158 * t283) * mrSges(5,3) + (t245 * t1004 - t151 * t993 + t218 * t96 + t36 * t260 + t31 * t968) * t507 + (t246 * t1004 - t151 * t995 - t218 * t97 - t36 * t261 + t30 * t970) * t502 + t293 * (Ifges(7,5) * t413 + Ifges(7,6) * t414) + (t292 * t379 + t319 * t345) * Ifges(8,5)) * qJD(2), t886 + t1 * qJD(2) + t12 * qJD(4) + (-t176 * t310 + (-t855 + t858) * t240) * mrSges(5,3) * t922 + (t231 * t118 + t86 * t260 + t47 * t968 + t77 * t986) * t753 + (-t231 * t119 - t86 * t261 + t46 * t970 + t78 * t986) * t754 + (t560 + t572) * t756 + (-t270 * t910 - t239 * t913 + t176 * t199 / 0.2e1 + t85 * t257 + t885 * t972 + t763 + t1032 * t139 + (t139 * t974 - t270 * t190) * Ifges(9,5) - (-mrSges(5,3) * t935 + t239 * Ifges(5,5) + t200 * t954 + t202 * t949 + t232 * t619) * t175) * qJD(3), t771 + t11 * qJD(2) + t12 * qJD(3) + (-t1035 * t612 - t256 * t619) * qJD(4); -qJD(3) * t2 - qJD(4) * t10 - t893, t3 * qJD(3) + t16 * qJD(4) + t562 * t758 + ((-t292 * t409 - t345 * t845) * mrSges(8,2) + (t292 * t410 - t345 * t846) * mrSges(8,1)) * t923 + t542 * t884 + (mrSges(6,3) * t658 + t545) * t757 + ((t37 * t1022 + t36 * t609) * t1013 - t217 * t29 + t348 * Ifges(7,3) * t293 + t345 * Ifges(8,3) * t292 + 0.4e1 * (-t158 * t584 + t159 * t285) * t1014 - t37 * t166 + (m(9) * t278 - mrSges(9,2) * t250) * t132 + (m(9) * t277 + t250 * mrSges(9,1)) * t133) * qJD(2), -t903 + t3 * qJD(2) + t8 * qJD(4) + (t727 * t200 + t46 * t989) * t753 + (t727 * t202 + t47 * t989) * t754 + (t228 * Ifges(5,3) + t547) * t756 + (m(5) * (t177 * t310 + t178 * t311 + t240 * t605) + (-t228 * t866 + t925 * t310) * mrSges(5,2) + (-t228 * t867 - t925 * t311) * mrSges(5,1)) * t922 + t1052 * t755 + (-t85 * t166 - t41 * t197 - t217 * t45 - t232 * t26 - t611 + t762 + t562 * t139 - (Ifges(9,3) * t770 + t618) * t270 + (t101 * t648 + t228 * t657 - t925 * t647) * mrSges(6,3) + (t85 * t217 + t232 * t41 - t15 / 0.4e1 + t40 * t647 + t86 * t648) * m(6)) * qJD(3), -t883 + t16 * qJD(2) + t8 * qJD(3) + (-t218 * t620 + t568) * qJD(4); qJD(2) * t2 + qJD(4) * t13 - t886, t903 + t4 * qJD(3) + t7 * qJD(4) + t546 * t884 - t1052 * t757 + (t925 * t547 - t50 * t166 - t37 * t197 - t217 * t20 - t232 * t29 + t551 + (-t103 * t270 - t770 * t250) * Ifges(9,3) + (t132 * t270 + t188 * t250 + t770 * t278) * mrSges(9,2) + (-t133 * t270 - t189 * t250 - t770 * t277) * mrSges(9,1) + t545 * t239 + (-t102 / 0.4e1 + (t158 * t310 + t159 * t311 + t229 * t605) * pkin(4)) * m(5) + (t239 * t658 + t74 * t647 + t925 * t648) * mrSges(6,3) + (t217 * t50 + t232 * t37 - t14 / 0.4e1 + t36 * t647 + t51 * t648) * m(6)) * qJD(2), t4 * qJD(2) + (-t139 * t908 - t85 * t197 - t232 * t45 + (t85 * t1021 + t608 * t86) * t1013) * qJD(3) + t32 * qJD(4) + (mrSges(6,3) * t657 - t240 * t582 + t46 * t949 + t47 * t952) * t755 + (t239 * Ifges(5,3) + t541) * t756, t882 + t7 * qJD(2) + t32 * qJD(3) + (-Ifges(6,6) * t869 - t231 * t620 + t237) * qJD(4); qJD(2) * t10 - qJD(3) * t13 - t771, t883 + (-t619 * t36 - t878 + t924) * qJD(2) + t6 * qJD(3) - t1025 * t757, -t882 + t6 * qJD(2) + (-t619 * t86 - t870 + t885) * qJD(3) - t1026 * t755, 0;];
Cq = t17;
