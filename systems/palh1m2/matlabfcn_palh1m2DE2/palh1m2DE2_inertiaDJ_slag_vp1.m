% Calculate time derivative of joint inertia matrix for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE2_inertiaDJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:31
% EndTime: 2020-05-02 21:00:54
% DurationCPUTime: 52.60s
% Computational Cost: add. (149957->1110), mult. (252405->1617), div. (6038->8), fcn. (363791->74), ass. (0->706)
t572 = -qJD(3) - qJD(2);
t598 = cos(qJ(1));
t1062 = t572 * t598;
t590 = sin(qJ(3));
t539 = pkin(5) * t590 + pkin(1);
t591 = sin(qJ(2));
t596 = cos(qJ(3));
t597 = cos(qJ(2));
t827 = t596 * t597;
t1047 = pkin(5) * t827 - t539 * t591;
t830 = t591 * t596;
t447 = t590 * t597 + t830;
t628 = t447 * qJD(3);
t801 = qJD(2) * t597;
t802 = qJD(2) * t591;
t1051 = t539 * t801 - (-t596 * t802 - t628) * pkin(5);
t1061 = t1047 * t1051;
t571 = pkin(22) + pkin(21);
t738 = pkin(18) - t571;
t528 = -pkin(20) + t738;
t653 = qJ(4) + t528;
t496 = qJ(1) + t653;
t467 = sin(t496);
t654 = -qJ(4) + t528;
t499 = -qJ(1) + t654;
t470 = sin(t499);
t581 = qJ(1) + qJ(4);
t554 = sin(t581);
t1055 = t554 / 0.2e1 - t470 / 0.4e1 + t467 / 0.4e1;
t1060 = t1055 * rSges(6,2);
t559 = cos(t581);
t820 = cos(t499) + cos(t496);
t1059 = t559 / 0.2e1 + t820 / 0.4e1;
t592 = sin(qJ(1));
t586 = sin(pkin(19));
t588 = cos(pkin(19));
t435 = t586 * t596 + t588 * t590;
t438 = t586 * t590 - t588 * t596;
t348 = qJ(2) + atan2(t438, t435);
t342 = sin(t348);
t343 = cos(t348);
t686 = Icges(9,5) * t342 + Icges(9,6) * t343;
t995 = Icges(9,3) * t598 + t592 * t686;
t593 = sin(pkin(18));
t599 = cos(pkin(18));
t449 = -t590 * t599 + t593 * t596;
t450 = t590 * t593 + t596 * t599;
t337 = t449 * t591 + t450 * t597;
t546 = sin(t571);
t547 = cos(t571);
t662 = t449 * t597 - t450 * t591;
t226 = t337 * t546 - t547 * t662;
t227 = t337 * t547 + t546 * t662;
t585 = sin(pkin(22));
t917 = cos(pkin(22));
t439 = t593 * t585 + t599 * t917;
t440 = t585 * t599 - t593 * t917;
t155 = -qJ(2) - atan2(t439 * t597 - t440 * t591, t439 * t591 + t440 * t597) + pkin(21) - atan2(t226, t227);
t153 = sin(t155);
t154 = cos(t155);
t657 = Icges(11,5) * t153 - Icges(11,6) * t154;
t998 = Icges(11,3) * t598 + t592 * t657;
t1057 = t998 + t995;
t108 = Icges(11,3) * t592 - t598 * t657;
t229 = Icges(9,3) * t592 - t598 * t686;
t580 = qJ(2) + qJ(3);
t553 = sin(t580);
t558 = cos(t580);
t690 = Icges(4,5) * t558 - Icges(4,6) * t553;
t366 = Icges(4,3) * t592 + t598 * t690;
t1056 = t108 + t229 + t366;
t113 = -Icges(11,5) * t154 - Icges(11,6) * t153;
t255 = Icges(9,5) * t343 - Icges(9,6) * t342;
t419 = Icges(4,5) * t553 + Icges(4,6) * t558;
t1035 = t113 + t255 + t419;
t814 = Icges(11,4) * t154;
t114 = -Icges(11,2) * t153 - t814;
t815 = Icges(11,4) * t153;
t115 = -Icges(11,1) * t154 - t815;
t1004 = -t114 * t154 + t115 * t153;
t902 = Icges(9,4) * t343;
t256 = -Icges(9,2) * t342 + t902;
t903 = Icges(9,4) * t342;
t257 = Icges(9,1) * t343 - t903;
t1005 = -t256 * t343 - t257 * t342;
t911 = Icges(4,4) * t553;
t420 = Icges(4,2) * t558 + t911;
t910 = Icges(4,4) * t558;
t421 = Icges(4,1) * t553 + t910;
t1032 = -t420 * t553 + t421 * t558 - t1004 + t1005;
t587 = cos(pkin(20));
t916 = sin(pkin(20));
t437 = t587 * t593 - t599 * t916;
t441 = t599 * t587 + t593 * t916;
t327 = t437 * t596 - t441 * t590;
t331 = t437 * t590 + t441 * t596;
t1049 = -t327 * t591 - t331 * t597;
t1050 = -t327 * t597 + t331 * t591;
t184 = t1049 * t547 + t1050 * t546;
t942 = rSges(4,2) * t553;
t943 = rSges(4,1) * t558;
t724 = -t942 + t943;
t378 = t724 * t572;
t1054 = m(4) * t378;
t574 = qJD(1) - qJD(4);
t497 = qJ(1) + t654;
t472 = cos(t497);
t498 = -qJ(1) + t653;
t473 = cos(t498);
t582 = qJ(1) - qJ(4);
t560 = cos(t582);
t725 = t560 / 0.2e1 - t473 / 0.4e1 - t472 / 0.4e1;
t468 = sin(t497);
t469 = sin(t498);
t555 = sin(t582);
t726 = t555 / 0.2e1 + t469 / 0.4e1 - t468 / 0.4e1;
t611 = rSges(6,1) * t725 + rSges(6,2) * t726;
t1053 = t574 * t611;
t612 = rSges(6,1) * t726 - rSges(6,2) * t725;
t1052 = t574 * t612;
t1039 = -t1047 + rSges(5,1) * cos(t528) - rSges(5,2) * sin(t528) - pkin(15);
t294 = rSges(5,3) * t598 + t1039 * t592;
t1046 = -t592 * rSges(5,3) + t1039 * t598;
t183 = -t546 * t1049 + t1050 * t547;
t476 = rSges(11,1) * t597 - rSges(11,2) * t591;
t480 = rSges(11,1) * t591 + rSges(11,2) * t597;
t520 = t738 - t580;
t501 = sin(t520);
t502 = cos(t520);
t577 = pkin(18) - pkin(22);
t552 = -qJ(2) + t577;
t796 = pkin(21) - atan2(cos(t552), -sin(t552));
t317 = -atan2(-t501, t502) + t796;
t315 = cos(t317);
t853 = t315 * t598;
t314 = sin(t317);
t855 = t314 * t598;
t208 = t476 * t853 + t480 * t855;
t839 = t501 ^ 2 / t502 ^ 2;
t363 = t572 / (0.1e1 + t839);
t204 = t363 * t839 + qJD(2) + t363;
t1044 = (-t314 * t476 + t315 * t480) * t204;
t135 = atan2(t183, t184) + t580;
t133 = sin(t135);
t595 = cos(qJ(4));
t174 = 0.1e1 / t184 ^ 2;
t318 = t331 * qJD(3);
t319 = t327 * qJD(3);
t800 = qJD(3) * t597;
t831 = t591 * t319;
t832 = t591 * t318;
t62 = -t572 + (((t331 * t800 + t831) * t547 + t546 * (t319 * t597 - t832) - t184 * qJD(2)) / t184 - ((-t327 * t800 + t832) * t547 - (-t318 * t597 - t831) * t546 + t183 * qJD(2)) * t183 * t174) / (t174 * t183 ^ 2 + 0.1e1);
t589 = sin(qJ(4));
t906 = Icges(6,4) * t595;
t697 = -Icges(6,2) * t589 + t906;
t907 = Icges(6,4) * t589;
t134 = cos(t135);
t930 = t134 * t62;
t22 = t697 * t930 + (Icges(6,6) * t62 + (-Icges(6,2) * t595 - t907) * qJD(4)) * t133;
t76 = -Icges(6,6) * t134 + t133 * t697;
t706 = Icges(6,1) * t595 - t907;
t77 = -Icges(6,5) * t134 + t133 * t706;
t1043 = -t589 * t22 + (-t589 * t77 - t595 * t76) * qJD(4);
t566 = t592 * rSges(4,3);
t1042 = qJD(1) * (t598 * t724 + t566);
t694 = Icges(9,2) * t343 + t903;
t231 = Icges(9,6) * t592 - t598 * t694;
t641 = t572 * t256;
t121 = qJD(1) * t231 + t592 * t641;
t703 = Icges(9,1) * t342 + t902;
t233 = Icges(9,5) * t592 - t598 * t703;
t642 = t572 * t257;
t123 = qJD(1) * t233 + t592 * t642;
t991 = Icges(9,5) * t598 + t592 * t703;
t993 = Icges(9,6) * t598 + t592 * t694;
t674 = t342 * t993 - t343 * t991;
t1041 = t121 * t343 + t123 * t342 - t572 * t674;
t120 = qJD(1) * t993 + t598 * t641;
t122 = qJD(1) * t991 + t598 * t642;
t672 = t231 * t342 - t233 * t343;
t1040 = -t120 * t343 - t122 * t342 - t572 * t672;
t659 = Icges(11,1) * t153 - t814;
t996 = Icges(11,5) * t598 + t592 * t659;
t658 = -Icges(11,2) * t154 + t815;
t997 = Icges(11,6) * t598 + t592 * t658;
t684 = t153 * t996 - t154 * t997;
t1012 = t598 * t684;
t675 = -t342 * t991 - t343 * t993;
t1015 = t598 * t675;
t698 = -Icges(4,2) * t553 + t910;
t367 = -Icges(4,6) * t598 + t592 * t698;
t707 = Icges(4,1) * t558 - t911;
t369 = -Icges(4,5) * t598 + t592 * t707;
t667 = t367 * t553 - t369 * t558;
t1016 = t598 * t667;
t365 = -Icges(4,3) * t598 + t592 * t690;
t1037 = -t1012 + t1015 + t1016 + (-t365 + t1057) * t592;
t368 = Icges(4,6) * t592 + t598 * t698;
t370 = Icges(4,5) * t592 + t598 * t707;
t666 = t368 * t553 - t370 * t558;
t673 = t231 * t343 + t233 * t342;
t110 = Icges(11,6) * t592 - t598 * t658;
t112 = Icges(11,5) * t592 - t598 * t659;
t683 = t110 * t154 - t112 * t153;
t1036 = (-t666 - t673 + t683) * t598 + t1056 * t592;
t1034 = t1055 * rSges(6,1);
t131 = t694 * t572;
t132 = t703 * t572;
t375 = t698 * t572;
t376 = t707 * t572;
t225 = 0.1e1 / t227 ^ 2;
t408 = t450 * qJD(3);
t409 = t449 * qJD(3);
t241 = qJD(2) * t337 + t408 * t597 + t409 * t591;
t242 = qJD(2) * t662 - t408 * t591 + t409 * t597;
t68 = (-(t241 * t547 + t242 * t546) / t227 + (-t241 * t546 + t242 * t547) * t226 * t225) / (t225 * t226 ^ 2 + 0.1e1);
t51 = t658 * t68;
t52 = t659 * t68;
t631 = t114 * t68;
t632 = t115 * t68;
t638 = t420 * t572;
t639 = t421 * t572;
t1033 = -t131 * t343 - t132 * t342 - (t256 * t342 - t257 * t343) * t572 + (-t376 + t638) * t558 + (t375 + t639) * t553 + (t51 - t632) * t154 + (-t52 - t631) * t153 + t1035 * qJD(1);
t573 = qJD(1) + qJD(4);
t1031 = t1059 * t573;
t1030 = t1032 * qJD(1) - t657 * t68 + (-t686 + t690) * t572;
t489 = -rSges(7,1) * t593 + rSges(7,2) * t599;
t490 = rSges(7,1) * t599 + rSges(7,2) * t593;
t594 = sin(pkin(17));
t600 = cos(pkin(17));
t356 = -t489 * t600 - t490 * t594;
t358 = t489 * t594 - t490 * t600;
t266 = t356 * t597 + t358 * t591;
t1029 = -t356 * t591 + t358 * t597;
t273 = atan2(t438, -t435) + t348;
t268 = sin(t273);
t1028 = t268 / 0.2e1;
t269 = cos(t273);
t1027 = -t269 / 0.2e1;
t689 = Icges(6,5) * t595 - Icges(6,6) * t589;
t21 = t689 * t930 + (Icges(6,3) * t62 + (-Icges(6,5) * t589 - Icges(6,6) * t595) * qJD(4)) * t133;
t926 = t589 * t76;
t1026 = -t62 * t926 - t21;
t913 = Icges(3,4) * t591;
t699 = Icges(3,2) * t597 + t913;
t383 = Icges(3,6) * t592 - t598 * t699;
t912 = Icges(3,4) * t597;
t708 = Icges(3,1) * t591 + t912;
t385 = Icges(3,5) * t592 - t598 * t708;
t664 = t383 * t597 + t385 * t591;
t1023 = t592 * t664;
t1022 = t592 * t666;
t1021 = t592 * t673;
t448 = t593 * t600 - t594 * t599;
t452 = t593 * t594 + t599 * t600;
t338 = t448 * t591 + t452 * t597;
t733 = t448 * t597 - t452 * t591;
t695 = Icges(7,4) * t733 - Icges(7,2) * t338;
t221 = Icges(7,6) * t592 + t598 * t695;
t904 = Icges(7,4) * t338;
t704 = Icges(7,1) * t733 - t904;
t223 = Icges(7,5) * t592 + t598 * t704;
t677 = t221 * t338 - t223 * t733;
t1020 = t592 * t677;
t901 = Icges(10,4) * t268;
t692 = Icges(10,2) * t269 + t901;
t192 = Icges(10,6) * t592 + t598 * t692;
t900 = Icges(10,4) * t269;
t701 = Icges(10,1) * t268 + t900;
t194 = Icges(10,5) * t592 + t598 * t701;
t678 = t192 * t269 + t194 * t268;
t1019 = t592 * t678;
t1018 = t592 * t683;
t990 = Icges(3,5) * t598 + t592 * t708;
t992 = Icges(3,6) * t598 + t592 * t699;
t665 = -t591 * t990 - t597 * t992;
t1017 = t598 * t665;
t222 = -Icges(7,6) * t598 + t592 * t695;
t224 = -Icges(7,5) * t598 + t592 * t704;
t676 = t222 * t338 - t224 * t733;
t1014 = t598 * t676;
t191 = -Icges(10,6) * t598 + t592 * t692;
t193 = -Icges(10,5) * t598 + t592 * t701;
t679 = t191 * t269 + t193 * t268;
t1013 = t598 * t679;
t316 = -qJ(2) + t317;
t312 = sin(t316);
t313 = cos(t316);
t142 = (rSges(11,1) * t313 + rSges(11,2) * t312) * (-qJD(2) + t204);
t1008 = Icges(5,1) - Icges(5,2);
t216 = -rSges(11,1) * t312 + rSges(11,2) * t313;
t483 = rSges(3,1) * t591 + rSges(3,2) * t597;
t918 = -pkin(15) + t483;
t373 = rSges(3,3) * t598 + t918 * t592;
t779 = t592 * t930;
t803 = qJD(1) * t598;
t626 = t133 * t803 + t779;
t778 = t598 * t930;
t804 = qJD(1) * t592;
t625 = -t133 * t804 + t778;
t1006 = t1046 * t592 - t294 * t598;
t1003 = -t110 * t153 - t112 * t154;
t1002 = -t153 * t997 - t154 * t996;
t1001 = qJD(1) * t365;
t152 = -t572 - qJD(3);
t719 = rSges(10,1) * t269 - rSges(10,2) * t268;
t869 = t572 * t343;
t84 = pkin(2) * t869 + t152 * t719;
t685 = Icges(10,5) * t268 + Icges(10,6) * t269;
t189 = -Icges(10,3) * t598 + t592 * t685;
t687 = Icges(7,5) * t733 - Icges(7,6) * t338;
t220 = -Icges(7,3) * t598 + t592 * t687;
t691 = Icges(3,5) * t591 + Icges(3,6) * t597;
t994 = Icges(3,3) * t598 + t592 * t691;
t371 = -rSges(4,3) * t598 + t592 * t724;
t731 = qJD(1) * t134 - qJD(4);
t876 = t133 * t598;
t989 = t592 * t731 + t62 * t876;
t515 = -qJ(1) + t528;
t492 = sin(t515);
t514 = qJ(1) + t528;
t493 = cos(t514);
t568 = t598 * pkin(15);
t601 = pkin(11) + rSges(6,3);
t494 = cos(t515);
t964 = t494 / 0.2e1;
t491 = sin(t514);
t965 = t491 / 0.2e1;
t985 = (t493 / 0.2e1 + t964) * pkin(9) + (t965 + t492 / 0.2e1) * t601 - t568;
t578 = t592 ^ 2;
t579 = t598 ^ 2;
t650 = t1031 * rSges(6,1) - t573 * t1060;
t166 = t650 + t1053;
t983 = 0.2e1 * t166;
t649 = rSges(6,2) * t1031 + t1034 * t573;
t167 = t649 + t1052;
t982 = 0.2e1 * t167;
t972 = -rSges(6,1) / 0.2e1;
t648 = t559 * t972 - t820 * rSges(6,1) / 0.4e1 + t1060;
t196 = -t611 + t648;
t981 = 0.2e1 * t196;
t647 = t1059 * rSges(6,2) + t1034;
t197 = t612 + t647;
t980 = 0.2e1 * t197;
t979 = -0.2e1 * t592;
t978 = -0.2e1 * t598;
t977 = -pkin(5) / 0.2e1;
t976 = m(4) / 0.2e1;
t975 = m(5) / 0.2e1;
t974 = m(6) / 0.2e1;
t973 = m(8) * pkin(1);
t969 = t338 / 0.2e1;
t968 = t733 / 0.2e1;
t959 = t592 / 0.2e1;
t958 = -t598 / 0.2e1;
t583 = qJ(1) + qJ(2);
t957 = pkin(1) * sin(t583);
t584 = qJ(1) - qJ(2);
t557 = sin(t584);
t956 = pkin(1) * t557;
t562 = cos(t584);
t955 = pkin(1) * t562;
t954 = pkin(1) * (qJD(1) - qJD(2));
t953 = pkin(1) * t591;
t952 = pkin(1) * t597;
t951 = pkin(1) * t598;
t950 = pkin(2) * t342;
t564 = qJ(1) - t580;
t536 = sin(t564);
t949 = pkin(5) * t536;
t538 = cos(t564);
t948 = pkin(5) * t538;
t947 = pkin(5) * (qJD(1) + t572);
t484 = rSges(3,1) * t597 - rSges(3,2) * t591;
t946 = m(3) * t484;
t402 = pkin(5) * t830 + t539 * t597;
t945 = m(5) * t402;
t944 = m(7) * t1029;
t941 = rSges(4,2) * t558;
t939 = rSges(3,3) * t592;
t935 = rSges(7,3) * t598;
t934 = rSges(8,3) * t592;
t933 = rSges(9,3) * t592;
t932 = rSges(10,3) * t592;
t931 = t133 * t62;
t828 = t595 * t598;
t834 = t589 * t592;
t126 = -t134 * t834 - t828;
t829 = t592 * t595;
t833 = t589 * t598;
t127 = t134 * t829 - t833;
t878 = t133 * t592;
t69 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t878;
t929 = t134 * t69;
t128 = -t134 * t833 + t829;
t129 = t134 * t828 + t834;
t70 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t876;
t928 = t134 * t70;
t924 = t592 * rSges(6,1);
t567 = t592 * rSges(6,2);
t565 = t592 * rSges(7,3);
t921 = t595 * t77;
t920 = t598 * rSges(6,1);
t919 = pkin(14) - t266;
t914 = rSges(11,3) * t592;
t324 = t733 * qJD(2);
t905 = Icges(7,4) * t324;
t877 = t133 * t595;
t868 = t572 * t592;
t209 = qJD(1) * t294 - t1051 * t598;
t866 = t209 * t592;
t210 = qJD(1) * t1046 + t1051 * t592;
t865 = t210 * t598;
t856 = t314 * t592;
t854 = t315 * t592;
t846 = t378 * t592;
t845 = t992 * t591;
t844 = t383 * t591;
t434 = rSges(4,1) * t553 + t941;
t841 = t434 * t592;
t837 = t553 * t572;
t836 = t558 * t572;
t264 = t1029 * qJD(2);
t826 = rSges(7,3) * t803 + t264 * t598;
t254 = t266 * t598 + t565;
t207 = t476 * t854 + t480 * t856;
t818 = t598 * t943 + t566;
t372 = -t598 * t942 + t818;
t281 = t592 * t371 + t598 * t372;
t563 = qJ(1) + t580;
t535 = sin(t563);
t774 = (qJD(1) - t572) * t977;
t422 = t535 * t774;
t772 = t947 / 0.2e1;
t354 = t536 * t772 + t422;
t537 = cos(t563);
t423 = t537 * t774;
t355 = t538 * t772 + t423;
t533 = -pkin(15) + t953;
t784 = pkin(1) * t801;
t824 = t533 * t803 + t592 * t784;
t823 = t469 - t468;
t822 = t470 - t467;
t821 = t473 + t472;
t819 = rSges(4,3) * t803 + t804 * t942;
t504 = t535 * t977;
t400 = t949 / 0.2e1 + t504;
t506 = pkin(5) * t537 / 0.2e1;
t401 = -t948 / 0.2e1 + t506;
t783 = qJD(2) * t951;
t785 = qJD(1) * t952;
t817 = t591 * t783 + t592 * t785;
t816 = t578 + t579;
t810 = qJD(1) * t108;
t190 = Icges(10,3) * t592 + t598 * t685;
t809 = qJD(1) * t190;
t219 = Icges(7,3) * t592 + t598 * t687;
t808 = qJD(1) * t219;
t807 = qJD(1) * t229;
t806 = qJD(1) * t366;
t381 = Icges(3,3) * t592 - t598 * t691;
t805 = qJD(1) * t381;
t640 = t572 * t255;
t118 = qJD(1) * t995 + t598 * t640;
t119 = t592 * t640 + t807;
t143 = -t592 * t675 + t598 * t995;
t144 = -t229 * t598 - t1021;
t26 = (t598 * t119 + (t144 + t1015) * qJD(1)) * t598 + (t143 * qJD(1) + (t807 + t1040) * t592 + (-t118 + (t995 - t673) * qJD(1) + t1041) * t598) * t592;
t630 = t68 * t113;
t38 = qJD(1) * t998 + t598 * t630;
t39 = t592 * t630 + t810;
t40 = qJD(1) * t997 + t598 * t631;
t42 = qJD(1) * t996 + t598 * t632;
t54 = t592 * t684 + t598 * t998;
t55 = -t108 * t598 + t1018;
t41 = qJD(1) * t110 + t592 * t631;
t750 = -t68 * t996 - t41;
t43 = qJD(1) * t112 + t592 * t632;
t752 = -t68 * t997 + t43;
t4 = (t598 * t39 + (t55 - t1012) * qJD(1)) * t598 + (t54 * qJD(1) + (t1003 * t68 - t153 * t42 + t154 * t40 + t810) * t592 + (-t38 + t750 * t154 + t752 * t153 + (t998 + t683) * qJD(1)) * t598) * t592;
t211 = -t365 * t598 - t592 * t667;
t212 = -t366 * t598 - t1022;
t637 = t572 * t419;
t282 = t598 * t637 - t1001;
t283 = t592 * t637 + t806;
t284 = -qJD(1) * t367 + t598 * t638;
t286 = -qJD(1) * t369 + t598 * t639;
t287 = qJD(1) * t370 + t592 * t639;
t734 = t367 * t572 + t287;
t285 = qJD(1) * t368 + t592 * t638;
t736 = -t369 * t572 + t285;
t83 = (t598 * t283 + (t212 + t1016) * qJD(1)) * t598 + (t211 * qJD(1) + (-t284 * t553 + t286 * t558 + t368 * t836 + t370 * t837 + t806) * t592 + (-t282 - t734 * t558 + t736 * t553 + (-t365 - t666) * qJD(1)) * t598) * t592;
t795 = -t4 - t83 - t26;
t794 = -0.2e1 * t572;
t793 = pkin(2) * m(10) / 0.2e1;
t23 = t706 * t930 + (Icges(6,5) * t62 + (-Icges(6,1) * t589 - t906) * qJD(4)) * t133;
t75 = -Icges(6,3) * t134 + t133 * t689;
t792 = t23 * t877 + t75 * t931 + t921 * t930;
t791 = 0.2e1 * qJD(1);
t788 = t62 * t929;
t787 = t62 * t928;
t782 = t69 * t878;
t780 = t70 * t876;
t541 = -t924 / 0.2e1;
t569 = t598 * rSges(6,2);
t461 = t569 / 0.2e1 + t541;
t462 = -t569 / 0.2e1 + t541;
t767 = t920 / 0.2e1;
t463 = t767 + t567 / 0.2e1;
t464 = t767 - t567 / 0.2e1;
t485 = t569 + t924;
t486 = t569 - t924;
t487 = t567 + t920;
t488 = t567 - t920;
t147 = t485 * t560 + t486 * t559 + t487 * t554 + t488 * t555 - t823 * t464 - t822 * t463 + t821 * t462 + t820 * t461 + ((-t493 + t494) * t598 + (-t491 - t492) * t592) * rSges(6,3);
t777 = t147 * t930;
t522 = -pkin(1) * cos(t583) / 0.2e1;
t776 = -t954 / 0.2e1;
t775 = t954 / 0.2e1;
t773 = -t947 / 0.2e1;
t72 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t876;
t74 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t876;
t711 = -t589 * t72 + t595 * t74;
t32 = t133 * t711 - t928;
t34 = t128 * t76 + t129 * t77 + t75 * t876;
t771 = t32 / 0.2e1 + t34 / 0.2e1;
t71 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t878;
t73 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t878;
t712 = -t589 * t71 + t595 * t73;
t31 = t133 * t712 - t929;
t33 = t126 * t76 + t127 * t77 + t75 * t878;
t770 = t33 / 0.2e1 + t31 / 0.2e1;
t764 = t434 * t804;
t763 = t558 * t804;
t755 = t804 / 0.2e1;
t753 = -t434 - t952;
t751 = -t110 * t68 - t42;
t749 = t112 * t68 - t40;
t748 = t919 * t592;
t745 = 0.2e1 * m(4);
t743 = 0.2e1 * m(6);
t742 = 2 * m(9);
t741 = 0.2e1 * m(10);
t739 = 2 * m(11);
t737 = t370 * t572 - t284;
t735 = t368 * t572 + t286;
t732 = -qJD(4) * t134 + qJD(1);
t729 = 0.4e1 * t816;
t728 = rSges(6,1) * t755;
t362 = t753 * t598;
t722 = rSges(8,1) * cos(t577) - rSges(8,2) * sin(t577);
t721 = rSges(9,1) * t343 - rSges(9,2) * t342;
t720 = -rSges(9,1) * t342 - rSges(9,2) * t343;
t718 = -rSges(10,1) * t268 - rSges(10,2) * t269;
t27 = t126 * t71 + t127 * t73 + t782;
t28 = t126 * t72 + t127 * t74 + t70 * t878;
t715 = t27 * t592 + t28 * t598;
t29 = t128 * t71 + t129 * t73 + t69 * t876;
t30 = t128 * t72 + t129 * t74 + t780;
t714 = t29 * t592 + t30 * t598;
t713 = t31 * t592 + t32 * t598;
t709 = Icges(3,1) * t597 - t913;
t629 = t338 * qJD(2);
t705 = -Icges(7,1) * t629 - t905;
t702 = Icges(10,1) * t269 - t901;
t700 = -Icges(3,2) * t591 + t912;
t696 = -Icges(7,4) * t629 - Icges(7,2) * t324;
t693 = -Icges(10,2) * t268 + t900;
t688 = -Icges(7,5) * t629 - Icges(7,6) * t324;
t451 = -t590 * t591 + t827;
t661 = t595 * t732;
t660 = t732 * t589;
t656 = -pkin(15) - t720;
t655 = 0.2e1 - 0.2e1 * t816;
t646 = t434 * t572;
t246 = t721 * t598;
t453 = t484 * qJD(2);
t645 = t152 * t702;
t644 = t152 * t693;
t643 = t152 * (Icges(10,5) * t269 - Icges(10,6) * t268);
t171 = -t372 * t804 + t592 * (t592 * t646 + t1042) + t598 * (t1062 * t941 + (t1062 * t553 - t763) * rSges(4,1) + t819) + t371 * t803;
t635 = qJD(2) * t709;
t634 = qJD(2) * t700;
t633 = qJD(2) * (-Icges(3,5) * t597 + Icges(3,6) * t591);
t627 = t451 * qJD(3);
t186 = -t718 - t950;
t624 = pkin(2) * t343 - t719;
t623 = pkin(4) * sin(qJ(2) - t796) - t216;
t622 = -pkin(15) - t186;
t173 = t624 * t598;
t124 = (-t1062 * t343 - t342 * t804) * rSges(9,2) + (-t1062 * t342 + t343 * t804) * rSges(9,1);
t125 = (t342 * t803 - t343 * t868) * rSges(9,2) + (-t342 * t868 - t343 * t803) * rSges(9,1);
t245 = t721 * t592;
t621 = (t572 * t720 * t721 - t246 * t124 - t245 * t125) * t742 + ((-t143 - t211 - t54) * t804 + t1037 * t803) * t598 + (((t118 + t38 + t282) * t592 + (t1021 - t1018 + t1022 - t1037) * qJD(1)) * t592 + (t212 + t55 + t144) * t804 + t1036 * t803 + ((t751 * t153 - t749 * t154 + t737 * t553 + t735 * t558 + t1040 - t119 - t283 - t39) * t592 + (t1002 * t68 + t153 * t43 - t154 * t41 + t285 * t553 - t287 * t558 - t367 * t836 - t369 * t837 - t1001 + t1041) * t598 + ((-t675 + t684 - t667 + t1056) * t592 + t1057 * t598 + t1036) * qJD(1)) * t598) * t592;
t620 = t623 * t592;
t619 = rSges(8,3) * t598 + t592 * t722;
t442 = t480 * qJD(2);
t446 = t476 * qJD(2);
t99 = qJD(1) * t208 + t1044 * t592 - t442 * t854 + t446 * t856;
t217 = rSges(9,3) * t598 + t592 * t656;
t618 = t152 * t718 - t572 * t950;
t137 = -t142 - t784;
t168 = rSges(10,3) * t598 + t592 * t622;
t614 = -t592 * pkin(15) + (-t493 / 0.2e1 + t964) * t601 + (t965 - t492 / 0.2e1) * pkin(9);
t100 = -t442 * t853 + t446 * t855 + (-t314 * t480 - t315 * t476) * t804 + t1044 * t598;
t116 = qJD(1) * t217 + t246 * t572;
t117 = -t721 * t868 + (t598 * t656 - t933) * qJD(1);
t218 = t598 * t720 + t568 + t933;
t606 = m(9) * (-t116 * t245 - t117 * t246 + t124 * t217 + t125 * t218) - t672 * t803 / 0.2e1 + (t1032 * t598 + t1035 * t592 + t368 * t558 + t370 * t553 + t1003) * t803 / 0.2e1 + (-t1030 * t592 + t1033 * t598 - t120 * t342 + t122 * t343 + t153 * t749 + t154 * t751 + t553 * t735 - t558 * t737 + t572 * t673) * t959 + (t1030 * t598 + t1033 * t592 - t121 * t342 + t123 * t343 + t153 * t750 - t154 * t752 + t553 * t734 + t558 * t736 + t572 * t675) * t958 + (t1032 * t592 - t1035 * t598 + t367 * t558 + t369 * t553 - t1002 + t674) * t755;
t602 = pkin(5) ^ 2;
t575 = qJD(1) + qJD(2);
t545 = rSges(6,2) * t803;
t524 = t803 * t972;
t521 = -t957 / 0.2e1;
t510 = pkin(1) * t592 * t802;
t495 = t533 * t592;
t465 = t533 * t804;
t460 = t575 * t522;
t459 = t575 * t957 / 0.2e1;
t454 = t483 * qJD(2);
t364 = -t483 * t598 + t568 + t939;
t361 = t753 * t592;
t347 = qJD(2) * t451 + t627;
t346 = -qJD(2) * t447 - t628;
t345 = -t956 / 0.2e1 + t521 + t401;
t344 = t522 - t955 / 0.2e1 + t400;
t341 = -t371 + t495;
t340 = (-t533 - t942) * t598 + t818;
t334 = t495 + t619;
t333 = t934 + (-t533 - t722) * t598;
t322 = -t539 * t802 + (t596 * t801 + t627) * pkin(5);
t307 = t592 * t633 + t805;
t306 = qJD(1) * t994 + t598 * t633;
t303 = t453 * t592 + (t598 * t918 - t939) * qJD(1);
t302 = 0.4e1 * t602 * t451 * t346;
t301 = qJD(1) * t373 - t598 * t453;
t300 = qJD(1) * t619 - t597 * t783 + t465;
t299 = (t598 * t722 - t934) * qJD(1) + t824;
t297 = t557 * t775 + t355 + t459;
t296 = t562 * t776 + t354 + t460;
t274 = t281 - t953;
t271 = t378 * t598 + t764 + t817;
t270 = qJD(1) * t362 + t510 + t846;
t265 = t266 * qJD(2);
t253 = t266 * t592 - t935;
t248 = -rSges(4,1) * t763 + t465 + (t646 - t784) * t598 + t819;
t247 = -t572 * t841 - t1042 + t824;
t244 = -pkin(14) * t598 + t254;
t243 = t748 + t935;
t240 = t381 * t592 - t598 * t664;
t239 = -t592 * t994 - t1017;
t238 = -t381 * t598 - t1023;
t237 = -t592 * t665 + t598 * t994;
t215 = t216 - t953;
t206 = -t597 * t951 + t208;
t205 = -t592 * t952 + t207;
t199 = t914 + (-t533 - t623) * t598;
t198 = rSges(11,3) * t598 + t495 + t620;
t188 = -t264 * t592 + (t598 * t919 - t565) * qJD(1);
t187 = qJD(1) * t748 + t826;
t178 = t592 * t688 + t808;
t175 = -qJD(1) * t220 + t688 * t598;
t172 = t624 * t592;
t170 = t171 - t784;
t169 = t186 * t598 + t568 + t932;
t159 = -t949 / 0.2e1 + t504 + t522 + t955 / 0.2e1 - t612 + t614 + t647;
t158 = t948 / 0.2e1 + t506 + t956 / 0.2e1 + t521 + t611 + t648 - t985;
t151 = qJD(1) * t614 + t536 * t773 + t562 * t775 - t1052 + t422 + t460 + t649;
t150 = qJD(1) * t985 + t538 * t773 + t557 * t776 - t1053 + t423 + t459 + t650;
t141 = -t220 * t598 - t592 * t676;
t140 = -t219 * t598 - t1020;
t139 = t220 * t592 - t1014;
t138 = t219 * t592 - t598 * t677;
t106 = t142 * t592 + (t598 * t623 - t914) * qJD(1) + t824;
t105 = t465 + qJD(1) * t620 + (rSges(11,3) * qJD(1) + t137) * t598;
t104 = t190 * t592 + t598 * t678;
t103 = t189 * t592 + t1013;
t102 = -t190 * t598 + t1019;
t101 = -t189 * t598 + t592 * t679;
t98 = (t555 + t554) * t545 + t822 * (t728 - t545 / 0.2e1) + t823 * (t728 + t545 / 0.2e1) + t820 * (-rSges(6,2) * t804 / 0.2e1 + t524) + t821 * (rSges(6,2) * t755 + t524) + ((-t559 - t560) * t567 + ((-t559 + t560) * t598 + (-t554 + t555) * t592) * rSges(6,1)) * qJD(1) + (t462 * t823 + t464 * t821 - t485 * t555 + t488 * t560) * t574 + (t461 * t822 + t463 * t820 - t486 * t554 + t559 * t487) * t573;
t97 = t100 + t817;
t96 = -t598 * t785 + t510 + t99;
t86 = t592 * t643 + t809;
t85 = -qJD(1) * t189 + t598 * t643;
t81 = -qJD(1) * t173 + t592 * t618;
t80 = t598 * t618 + t624 * t804;
t79 = -t84 * t592 + (t598 * t622 - t932) * qJD(1);
t78 = qJD(1) * t168 + t598 * t84;
t49 = t731 * t828 + (-t62 * t877 + t660) * t592;
t48 = t592 * t661 + (-t598 * t731 + t62 * t878) * t589;
t47 = -t595 * t989 + t598 * t660;
t46 = t589 * t989 + t598 * t661;
t35 = -t134 * t75 + (t921 - t926) * t133;
t24 = t35 * t931;
t19 = Icges(6,1) * t49 + Icges(6,4) * t48 + Icges(6,5) * t626;
t18 = Icges(6,1) * t47 + Icges(6,4) * t46 + Icges(6,5) * t625;
t17 = Icges(6,4) * t49 + Icges(6,2) * t48 + Icges(6,6) * t626;
t16 = Icges(6,4) * t47 + Icges(6,2) * t46 + Icges(6,6) * t625;
t15 = Icges(6,5) * t49 + Icges(6,6) * t48 + Icges(6,3) * t626;
t14 = Icges(6,5) * t47 + Icges(6,6) * t46 + Icges(6,3) * t625;
t7 = t1026 * t134 + t1043 * t133 + t792;
t6 = t75 * t779 + t126 * t22 + t127 * t23 + t48 * t76 + t49 * t77 + (t21 * t592 + t75 * t803) * t133;
t5 = t75 * t778 + t128 * t22 + t129 * t23 + t46 * t76 + t47 * t77 + (t21 * t598 - t75 * t804) * t133;
t2 = (t62 * t711 - t14) * t134 + (-t16 * t589 + t18 * t595 + t62 * t70 + (-t589 * t74 - t595 * t72) * qJD(4)) * t133;
t1 = (t62 * t712 - t15) * t134 + (-t17 * t589 + t19 * t595 + t62 * t69 + (-t589 * t73 - t595 * t71) * qJD(4)) * t133;
t3 = [0.2e1 * m(5) * (-t1046 * t209 + t210 * t294) + (t324 * t733 - t338 * t629) * Icges(7,4) - (0.2e1 * Icges(7,2) * t733 + t904) * t629 + t733 * t905 + (t247 * t341 + t248 * t340) * t745 + (t150 * t159 + t151 * t158) * t743 + (t168 * t79 + t169 * t78) * t741 + (t116 * t218 + t117 * t217) * t742 + ((-t701 - t693) * t269 + (-t702 + t692) * t268) * t152 + (0.2e1 * Icges(5,4) * t930 + t1008 * t931 + t1026) * t134 + t1004 * t68 + (-t709 + t699) * t802 + 0.2e1 * m(8) * (t299 * t334 + t300 * t333) + 0.2e1 * m(7) * (t187 * t244 + t188 * t243) + 0.2e1 * m(3) * (t301 * t364 + t303 * t373) + 0.2e1 * t338 * t324 * Icges(7,1) + (t105 * t199 + t106 * t198) * t739 - t553 * t376 - t558 * t375 + t343 * t132 - t342 * t131 - t153 * t51 - t154 * t52 + (-t708 - t700) * t801 + (-0.2e1 * Icges(5,4) * t931 + t1008 * t930 + t1043) * t133 + t792 - t421 * t836 - t1005 * t572 + t420 * t837; (t247 * t362 + t248 * t361 + t270 * t340 + t271 * t341) * m(4) + ((-t299 * t598 - t300 * t592) * t597 + (t333 * t592 + t334 * t598) * t802) * t973 + (t105 * t205 + t106 * t206 + t198 * t97 + t199 * t96) * m(11) + (t150 * t344 + t151 * t345 + t158 * t296 + t159 * t297) * m(6) + t606 + ((t221 * t968 + t223 * t969 + t192 * t1028 + t194 * t1027 - t844 / 0.2e1 + t244 * t944 - t364 * t946 + t1046 * t945 + (t385 / 0.2e1 - t333 * t973) * t597) * t598 + (t845 / 0.2e1 + t222 * t968 + t224 * t969 + t191 * t1028 + t193 * t1027 - t243 * t944 + t373 * t946 + t294 * t945 + (-t990 / 0.2e1 + t334 * t973) * t597) * t592) * qJD(1) + (-qJD(2) * t664 + t152 * t678 + (-qJD(1) * t222 + t696 * t598) * t733 + (-qJD(1) * t224 + t705 * t598) * t338 - t221 * t629 + t223 * t324 + t268 * (-qJD(1) * t191 + t598 * t644) - t269 * (-qJD(1) * t193 + t598 * t645) - (qJD(1) * t992 - t598 * t634) * t591 + (qJD(1) * t990 - t598 * t635) * t597) * t959 + (-qJD(2) * t665 + t152 * t679 + (qJD(1) * t221 + t592 * t696) * t733 + (qJD(1) * t223 + t592 * t705) * t338 - t222 * t629 + t224 * t324 + t268 * (qJD(1) * t192 + t592 * t644) - t269 * (qJD(1) * t194 + t592 * t645) - (qJD(1) * t383 - t592 * t634) * t591 + (qJD(1) * t385 - t592 * t635) * t597) * t958 + (Icges(7,5) * t324 - Icges(7,6) * t629 - qJD(2) * t691 + t152 * t685) * (t578 / 0.2e1 + t579 / 0.2e1) + m(10) * (t168 * t80 + t169 * t81 - t172 * t78 - t173 * t79) + m(7) * (-(-t187 * t592 - t188 * t598) * t1029 + (-t243 * t598 - t244 * t592) * t265) + m(5) * ((-t865 - t866) * t402 + t1006 * t322) + m(3) * ((-t301 * t592 - t303 * t598) * t484 - (-t364 * t592 - t373 * t598) * t454); (t137 * t215 + t205 * t96 + t206 * t97) * t739 + (t170 * t274 + t270 * t361 + t271 * t362) * t745 + (-t172 * t81 - t173 * t80 + t186 * t84) * t741 + (t296 * t345 + t297 * t344 - t1061) * t743 + t592 * ((t592 * t85 + (t103 - t1019) * qJD(1)) * t592 + (t104 * qJD(1) + (-t86 + (t190 + t679) * qJD(1)) * t592) * t598) - t598 * ((t598 * t86 + (t102 - t1013) * qJD(1)) * t598 + (t101 * qJD(1) + t809 * t592 + (t678 * qJD(1) - t85) * t598) * t592) + t621 - t598 * t4 - t598 * t83 - t598 * t26 + m(8) * pkin(1) ^ 2 * t591 * t655 * t801 + t592 * ((t592 * t306 + (t239 + t1023) * qJD(1)) * t592 + (t240 * qJD(1) + (-t801 * t990 + t802 * t992) * t598 + (-t307 + (-t385 * t597 + t844) * qJD(2) + (t381 - t665) * qJD(1)) * t592) * t598) + m(3) * (-t454 * t484 * t729 + 0.4e1 * t453 * t483) / 0.2e1 + (t322 * t402 * t729 - 0.4e1 * t1061) * t975 + m(7) * (0.4e1 * (t253 * t592 + t254 * t598) * ((qJD(1) * t253 + t826) * t598 + (-qJD(1) * t254 + (rSges(7,3) * qJD(1) + t264) * t592) * t592) - t1029 * t265 * t729) / 0.2e1 - t598 * ((t598 * t307 + (t238 + t1017) * qJD(1)) * t598 + (t237 * qJD(1) + (t383 * t802 - t385 * t801 + t805) * t592 + (-t306 + (-t597 * t990 + t845) * qJD(2) - t664 * qJD(1)) * t598) * t592) + t592 * ((t592 * t175 + (t139 + t1020) * qJD(1)) * t592 + (t138 * qJD(1) + (-t178 + (t219 - t676) * qJD(1)) * t592) * t598) - t598 * ((t598 * t178 + (t140 + t1014) * qJD(1)) * t598 + (t141 * qJD(1) + t808 * t592 + (-t677 * qJD(1) - t175) * t598) * t592) + ((-t101 - t141 - t237) * t598 + (t102 + t140 + t238) * t592) * t804 + ((-t103 - t139 - t239) * t598 + (t104 + t138 + t240) * t592) * t803; (t150 * t400 + t151 * t401 + t158 * t354 + t159 * t355) * m(6) + (0.2e1 * t1006 * t347 + (-0.2e1 * t866 - 0.2e1 * t865 + (t1046 * t598 + t294 * t592) * t791) * t447) * pkin(5) * t975 + (t247 * t978 + t248 * t979 + (-t340 * t598 + t341 * t592) * t791) * t434 * t976 + (t100 * t198 + t105 * t207 + t106 * t208 + t199 * t99) * m(11) + ((t168 * t598 + t169 * t592) * t342 * t794 + (t78 * t979 + t79 * t978 + (t168 * t592 - t169 * t598) * t791) * t343) * t793 - (-t340 * t592 - t341 * t598) * t1054 + t606; (t296 * t401 + t297 * t400 + t344 * t355 + t345 * t354) * m(6) + (t170 * t281 + t171 * t274 - t270 * t841 + t361 * t846 + t362 * t764) * m(4) + (t100 * t206 + t137 * t216 - t142 * t215 + t205 * t99 + t207 * t96 + t208 * t97) * m(11) + ((-0.2e1 * t84 + (-t172 * t592 - t173 * t598) * t794) * t342 + (0.2e1 * t186 * t572 + t81 * t979 + t80 * t978 + (t172 * t598 - t173 * t592) * t791) * t343) * t793 + t621 + (0.2e1 * (t322 * t447 + t347 * t402) * t975 * t816 + (t346 * t1047 - t451 * t1051) * (m(5) + m(6))) * pkin(5) + (t362 * t1054 + (-qJD(1) * t361 - t271) * t434 * m(4) + t795) * t598; (t100 * t208 - t142 * t216 + t207 * t99) * t739 + t621 + (0.4e1 * t354 * t401 + 0.4e1 * t355 * t400 + t302) * t974 - m(10) * pkin(2) ^ 2 * t342 * t655 * t869 + t795 * t598 + (-t378 * t434 * t729 + 0.4e1 * t171 * t281) * t976 + (t347 * t447 * t602 * t729 + t302) * t975; (t150 * t197 + t151 * t196 + t158 * t167 + t159 * t166) * m(6) + t24 + (-t7 + (t592 * t770 + t598 * t771) * t62) * t134 + ((t2 / 0.2e1 + t5 / 0.2e1) * t598 + (t6 / 0.2e1 + t1 / 0.2e1) * t592 + (-t592 * t771 + t598 * t770) * qJD(1)) * t133; (-t1047 * t777 + t344 * t983 + t345 * t982 + t296 * t981 + t297 * t980 + (-t1047 * t98 + t1051 * t147) * t133) * t974; (t400 * t983 + t401 * t982 + t354 * t981 + t355 * t980 + (-t451 * t777 + (-t147 * t346 - t451 * t98) * t133) * pkin(5)) * t974; (t133 * t713 - t134 * t35) * t931 - t134 * (t24 + (t62 * t713 - t7) * t134 + (t1 * t592 + t2 * t598 + (t31 * t598 - t32 * t592) * qJD(1)) * t133) + ((t62 * t715 - t6) * t134 + (t33 * t62 + (t126 * t16 + t127 * t18 + t48 * t72 + t49 * t74 + (t27 + t780) * qJD(1)) * t598 + (-t28 * qJD(1) + t126 * t17 + t127 * t19 + t48 * t71 + t49 * t73 + (t133 * t15 + t788) * t592 + (t787 + (qJD(1) * t69 + t14) * t133) * t598) * t592) * t133) * t878 + (0.4e1 * t166 * t197 + 0.4e1 * t167 * t196 + (t133 * t98 + t777) * t147 * t133) * t974 + ((t62 * t714 - t5) * t134 + (t34 * t62 + (t128 * t17 + t129 * t19 + t46 * t71 + t47 * t73 + (-t30 - t782) * qJD(1)) * t592 + (t128 * t16 + t129 * t18 + t46 * t72 + t47 * t74 + t29 * qJD(1) + (t133 * t14 + t787) * t598 + (t788 + (-qJD(1) * t70 + t15) * t133) * t592) * t598) * t133) * t876 + t625 * (t133 * t714 - t134 * t34) + t626 * (t133 * t715 - t134 * t33);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;