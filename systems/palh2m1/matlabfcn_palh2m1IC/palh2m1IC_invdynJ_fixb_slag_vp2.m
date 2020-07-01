% Calculate vector of inverse dynamics joint torques with ic for
% palh2m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1IC_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:18
% EndTime: 2020-05-03 01:03:47
% DurationCPUTime: 29.72s
% Computational Cost: add. (5752->1023), mult. (10437->1418), div. (0->0), fcn. (1734->10), ass. (0->555)
t394 = cos(qJ(5));
t902 = mrSges(6,2) * t394;
t407 = m(4) + m(5);
t218 = mrSges(3,1) + (m(6) + t407) * pkin(2);
t389 = sin(qJ(5));
t419 = qJD(1) ^ 2;
t711 = t389 * t419;
t317 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t390 = sin(qJ(4));
t768 = mrSges(6,1) * t390;
t901 = pkin(3) * t768 + t317;
t415 = qJD(5) ^ 2;
t417 = qJD(3) ^ 2;
t418 = qJD(2) ^ 2;
t692 = (t417 + t418);
t572 = (-qJD(4) ^ 2 - t692);
t888 = t415 + t572;
t900 = 2 * qJD(3);
t814 = m(5) + m(6);
t894 = t814 * pkin(3);
t899 = mrSges(4,1) + t894;
t338 = qJD(3) + qJD(4);
t298 = qJD(2) + t338;
t479 = qJD(5) * t298;
t516 = t479 * t902;
t412 = pkin(4) * m(6);
t766 = mrSges(6,1) * t394;
t548 = mrSges(5,1) + t766;
t476 = t548 + t412;
t761 = mrSges(6,2) * t389;
t858 = t476 - t761;
t129 = t858 * g(2);
t398 = cos(qJ(1));
t130 = t858 * g(1);
t336 = qJDD(2) + qJDD(3);
t276 = qJDD(4) + t336;
t240 = Ifges(6,6) * t276;
t242 = Ifges(6,5) * t276;
t274 = t298 ^ 2;
t349 = qJDD(5) * Ifges(6,3);
t381 = (Ifges(6,4) - Ifges(5,5));
t386 = (mrSges(5,2) - mrSges(6,3));
t393 = sin(qJ(1));
t669 = qJD(1) * qJD(5);
t600 = 2 * t669;
t787 = (pkin(6) * t276);
t468 = pkin(1) * t600 + t787;
t670 = qJD(1) * qJD(4);
t592 = pkin(1) * t670;
t358 = mrSges(6,1) * qJDD(1);
t610 = pkin(1) * t358;
t683 = qJD(1) * t298;
t618 = m(6) * t683;
t738 = pkin(1) * qJDD(1);
t648 = -0.2e1 * t738;
t755 = pkin(1) * qJD(1);
t820 = 2 * qJD(4);
t698 = t381 * t820 - 2 * t386 * t755;
t384 = Ifges(6,1) - Ifges(6,2);
t712 = t384 * t389;
t753 = t274 * Ifges(6,4);
t359 = t394 ^ 2;
t830 = 0.2e1 * t359;
t821 = 2 * qJD(2);
t442 = qJD(3) * t820 + t338 * t821 - t572;
t357 = mrSges(6,1) * qJDD(5);
t694 = -mrSges(6,2) * t415 + t357;
t85 = mrSges(6,2) * t442 + t694;
t356 = mrSges(6,2) * qJDD(5);
t237 = mrSges(6,1) * t415 + t356;
t86 = mrSges(6,1) * t442 - t237;
t882 = -2 * pkin(6);
t744 = m(6) * qJDD(1);
t330 = pkin(4) * t744;
t690 = mrSges(5,1) * qJDD(1);
t890 = t330 + t690;
t427 = -((t381 * t900 + t698) * qJD(2)) - (t698 * qJD(3)) + (mrSges(6,1) * t468 - mrSges(6,2) * t648 - pkin(4) * t86 - t242) * t389 + (mrSges(6,2) * t468 - pkin(4) * t85 + t274 * t712 - t240 - 0.2e1 * t610) * t394 + Ifges(5,6) * t276 - t130 * t393 + (2 * t386 * t592) + t753 * t830 - t349 + (t572 * t381) + ((t618 * t882) - 0.2e1 * t890) * pkin(1);
t898 = t129 * t398 + t427;
t353 = Ifges(6,6) * qJDD(1);
t757 = Ifges(6,4) * t389;
t895 = mrSges(6,1) * pkin(4);
t226 = t757 + t895;
t464 = t226 * t298;
t564 = t384 * t683;
t510 = t359 * t564;
t421 = pkin(4) ^ 2;
t523 = t421 * t618;
t420 = pkin(6) ^ 2;
t524 = t420 * t618;
t774 = Ifges(6,3) - Ifges(6,2);
t544 = (Ifges(5,2) - Ifges(5,1) + t774);
t555 = t544 * t670;
t593 = mrSges(6,3) * t683;
t355 = mrSges(6,3) * qJDD(1);
t613 = pkin(4) * t355;
t790 = pkin(4) * t298;
t652 = mrSges(6,2) * t790;
t685 = Ifges(6,5) * qJDD(1);
t739 = Ifges(6,6) * qJD(5);
t740 = Ifges(6,5) * qJD(5);
t567 = mrSges(6,2) * t669;
t216 = -t358 + t567;
t788 = pkin(6) * t216;
t568 = mrSges(6,1) * t669;
t687 = mrSges(6,2) * qJDD(1);
t215 = t568 + t687;
t789 = pkin(6) * t215;
t448 = (0.2e1 * t330 + (4 * t593)) * pkin(6) + 0.2e1 * t510 + (2 * t524) - 0.2e1 * t523 - (2 * t555) + 0.2e1 * t613 + (-0.2e1 * t788 - 0.2e1 * t685 + (0.2e1 * t739 - 0.4e1 * t464) * qJD(1)) * t394 + (-0.2e1 * t789 + 0.2e1 * t353 + (0.4e1 * t652 + 0.2e1 * t740) * qJD(1)) * t389;
t671 = qJD(1) * qJD(3);
t422 = pkin(3) ^ 2;
t799 = m(6) * t422;
t517 = t671 * t799;
t243 = Ifges(4,1) - Ifges(4,2) + t544;
t777 = t422 * m(5);
t450 = -t243 + t777;
t556 = t450 * t671;
t339 = qJD(2) + qJD(3);
t387 = qJD(4) / 0.2e1;
t271 = t387 + t339;
t551 = t687 / 0.2e1;
t560 = mrSges(6,1) * qJD(5) / 0.2e1;
t187 = qJD(1) * t560 + t551;
t188 = t358 / 0.2e1 - t567 / 0.2e1;
t545 = -t670 / 0.2e1;
t481 = -t187 * t389 + t188 * t394 + t386 * t545 + t690 / 0.2e1 + t330 / 0.2e1;
t819 = m(6) * pkin(6);
t36 = (t271 * t819 - (t339 * t386)) * qJD(1) + t481;
t793 = pkin(3) * t390;
t651 = t36 * t793;
t879 = 2 * t683;
t541 = mrSges(6,2) * t879;
t195 = pkin(4) * t541;
t595 = mrSges(6,1) * t683;
t533 = pkin(4) * t595;
t197 = -0.2e1 * t533;
t395 = cos(qJ(4));
t360 = t395 ^ 2;
t869 = Ifges(6,6) * t394;
t461 = (Ifges(6,5) * t389 + t869) * qJD(5);
t484 = t298 * t544;
t263 = t384 * t359;
t806 = pkin(6) * mrSges(6,3);
t496 = t263 + (2 * t806);
t547 = -t353 + t789;
t756 = Ifges(6,5) * t394;
t403 = pkin(4) * mrSges(6,3);
t868 = Ifges(5,4) + t403;
t691 = t420 - t421;
t870 = t691 * m(6);
t754 = ((t197 - t788) * t394 + (t195 - t547) * t389 + t330 * pkin(6) + (-t756 + t868) * qJDD(1) + (t461 - t484 + (-0.2e1 * t394 * t757 + t496 + t870) * t298) * qJD(1)) * t360;
t897 = 0.2e1 * t556 + 0.2e1 * t517 + 0.4e1 * t651 + t448 - 0.4e1 * t754;
t181 = Ifges(6,4) * qJDD(1) + t384 * t669;
t296 = qJDD(1) * t384;
t565 = Ifges(6,4) * t669;
t529 = 4 * t565;
t183 = -t296 + t529;
t318 = mrSges(5,1) + t412;
t261 = t318 * t419;
t597 = mrSges(6,3) + t819;
t280 = -mrSges(5,2) + t597;
t321 = pkin(4) * t687;
t324 = pkin(4) * t358;
t372 = mrSges(6,1) * t419;
t328 = pkin(1) * t372;
t333 = (Ifges(6,1) + t774);
t538 = pkin(4) * t600;
t798 = pkin(1) * t419;
t21 = t129 * t393 + (mrSges(6,1) * t538 + t183 * t389 + t321 + t328) * t394 - (-t324 + (t538 + t798) * mrSges(6,2)) * t389 - t181 * t830 + pkin(1) * t261 + qJDD(1) * t381 + (t333 * t669) + (t280 * g(3));
t392 = sin(qJ(2));
t893 = t21 * t392;
t391 = sin(qJ(3));
t204 = t280 * g(2);
t236 = t280 * t419;
t316 = mrSges(6,2) * pkin(6) - Ifges(6,6);
t400 = g(3) * mrSges(6,1);
t526 = -pkin(6) * t358 + t685;
t805 = g(3) * mrSges(6,2);
t856 = pkin(1) * t236 + (pkin(6) * t687 + t317 * t600 - t353 - t400) * t394 + t204 * t393 - (t316 * t600 + t526 - t805) * t389 + Ifges(5,6) * qJDD(1) - t318 * g(3);
t30 = t856 * t392;
t760 = mrSges(6,2) * t418;
t765 = mrSges(6,1) * t418;
t889 = -t386 + t819;
t76 = (qJDD(2) * t889) + t318 * t418 - t389 * t760 + t394 * t765;
t773 = pkin(2) * t76 + t30;
t891 = t773 * t391;
t176 = t226 * t394;
t649 = pkin(4) * t761;
t845 = -t496 + 0.2e1 * t176 - 0.2e1 * t649;
t886 = -Ifges(6,5) * t888 - qJDD(5) * Ifges(6,6) - pkin(6) * t86 - t276 * t895;
t763 = mrSges(4,2) * t419;
t764 = mrSges(4,2) * t393;
t885 = qJDD(1) * Ifges(4,6) - pkin(1) * t763 - g(2) * t764 - g(3) * t899 - t21 * t390;
t220 = t236 * pkin(3);
t881 = 0.2e1 * t220;
t880 = 0.2e1 * t394;
t315 = mrSges(5,3) * pkin(3) - Ifges(4,5);
t878 = t316 * t389;
t876 = (t551 + t568) * t880 + (t358 - 0.2e1 * t567) * t389;
t686 = Ifges(5,4) * qJDD(1);
t808 = pkin(3) * t36;
t874 = (0.8e1 * (t685 + t788) * t394 + 0.8e1 * t547 * t389 - (0.8e1 * t330 + (16 * t593)) * pkin(6) - 0.8e1 * t613 - 0.8e1 * t686) * t390 - 0.8e1 * t808;
t872 = 4 * Ifges(6,4);
t657 = 8 * t298;
t871 = -t412 / 0.2e1;
t786 = g(1) * t398;
t473 = t317 * t683;
t707 = t391 * t398;
t706 = t391 * t419;
t725 = t181 * t389;
t488 = pkin(4) * t216 - t725;
t866 = t488 * t394;
t811 = mrSges(6,2) * pkin(4);
t463 = (-t712 + t811) * t395;
t231 = t597 * pkin(4) + Ifges(5,4);
t863 = t317 * t394 - t878;
t101 = t231 + t863;
t802 = 0.8e1 * t360 - 0.4e1;
t865 = t802 * t101;
t702 = t395 * t419;
t716 = t316 * t390;
t326 = t419 * t811;
t169 = t384 * t711 - t326;
t727 = t169 * t360;
t812 = mrSges(6,2) * pkin(3);
t864 = 0.2e1 * t727 + (-0.2e1 * t716 - t812) * t702;
t862 = (-t420 - t421) * t276;
t414 = m(5) * pkin(3);
t319 = mrSges(4,1) + t414;
t413 = pkin(3) * m(6);
t252 = t413 + t319;
t205 = t280 * g(1);
t861 = -t205 * t398 - t856;
t860 = 0.2e1 * t463 + t812;
t705 = t392 * t398;
t629 = g(1) * t705;
t737 = pkin(2) * qJDD(2);
t467 = -t629 - t737;
t859 = pkin(3) * t706 + t467;
t170 = -t814 * t422 + t243;
t492 = t761 - t766;
t779 = t393 * g(2);
t855 = -t779 - t798;
t854 = qJD(1) * t339;
t673 = qJD(2) + qJD(3) / 0.2e1;
t586 = t389 * t702;
t758 = Ifges(6,4) * t359;
t608 = (t360 - 0.1e1 / 0.2e1) * t419 * t758;
t715 = t317 * t390;
t369 = Ifges(6,4) * t419;
t615 = mrSges(6,1) * t711;
t542 = pkin(4) * t615;
t179 = t369 + t542;
t726 = t179 * t360;
t813 = mrSges(6,1) * pkin(3);
t852 = (-0.2e1 * t715 - t813) * t586 + 0.4e1 * t608 - 0.2e1 * t726;
t340 = (m(3) + t407);
t850 = pkin(1) * (m(6) + t340) + mrSges(2,1);
t641 = 0.16e2 * t394;
t847 = t226 * t641 - 0.8e1 * t263 - 0.16e2 * t649;
t186 = t899 * g(1);
t668 = qJD(2) * qJD(3);
t486 = 2 * t668 + t692;
t574 = t899 * g(2);
t846 = (t204 * t390 + t574) * t398 + Ifges(4,6) * t336 - t186 * t393 - t692 * t315 + (-(mrSges(6,1) * t486 - t237) * t389 - (mrSges(6,2) * t486 + t694) * t394) * pkin(3);
t178 = 0.2e1 * mrSges(5,1) * t755 - qJD(5) * (-Ifges(6,3) + t384);
t844 = t178 * qJD(4) - t276 * t381 + t572 * Ifges(5,6) - (0.2e1 * Ifges(5,6) * t338 - t178) * qJD(2) - (Ifges(5,6) * t820 - t178) * qJD(3) - t205 * t393;
t329 = t421 * t744;
t532 = pkin(4) * t593;
t543 = pkin(4) * t618;
t566 = Ifges(5,4) * t670;
t578 = t420 * t744;
t843 = -0.2e1 * t329 - 0.8e1 * t532 - 0.8e1 * t566 + 0.2e1 * t578 - (-0.2e1 * t296 + (8 * t565)) * t359 - (-0.4e1 * t355 + 0.8e1 * t543) * pkin(6);
t841 = -t317 * t641 - 0.16e2 * Ifges(5,4) - 0.16e2 * t403 + 0.16e2 * t878;
t839 = 2 * m(6);
t838 = 0.2e1 * pkin(2);
t837 = 0.4e1 * pkin(2);
t836 = 0.2e1 * pkin(3);
t835 = 0.2e1 * pkin(4);
t834 = -0.2e1 * mrSges(6,1);
t833 = 0.2e1 * mrSges(4,2);
t832 = 0.2e1 * mrSges(6,2);
t158 = -t492 + t318;
t728 = t158 * t390;
t627 = pkin(3) * t728;
t342 = pkin(6) * t412;
t598 = t342 + t868;
t98 = t598 + t863;
t95 = -Ifges(4,4) + t98;
t831 = -0.4e1 * (t95 + t627) * t419;
t829 = -0.2e1 * t360;
t828 = 0.2e1 * t360;
t826 = -0.4e1 * t390;
t825 = 0.2e1 * t391;
t824 = 0.8e1 * t391;
t823 = 0.4e1 * t395;
t822 = -0.2e1 * t419;
t815 = -mrSges(5,1) / 0.2e1;
t696 = t418 * mrSges(6,3) + qJDD(2) * t761;
t809 = pkin(2) * ((-mrSges(5,2) + t819) * t418 - t476 * qJDD(2) + t696);
t457 = t492 - t412;
t639 = pkin(6) * t744;
t575 = -qJDD(1) * t386 / 0.2e1 + mrSges(5,1) * t545 + t639 / 0.2e1;
t43 = (-mrSges(5,1) * t339 + t271 * t457) * qJD(1) + t575;
t807 = pkin(3) * t43;
t804 = t828 - 0.1e1;
t803 = 0.4e1 * t360 - 0.2e1;
t447 = (-0.4e1 * t330 - (8 * t593)) * pkin(6) - 0.4e1 * t510 - (4 * t524) + 0.4e1 * t523 + (4 * t555) - 0.4e1 * t613 + (0.4e1 * t788 + 0.4e1 * t685 + (0.8e1 * t464 - 0.4e1 * t739) * qJD(1)) * t394 + (0.4e1 * t789 - 0.4e1 * t353 + (-0.8e1 * t652 - 0.4e1 * t740) * qJD(1)) * t389;
t801 = ((4 * t544 * t854) + t447 - 0.4e1 * t686) * t390 - 0.4e1 * t808;
t472 = t316 * t683;
t487 = 0.4e1 * t488;
t792 = pkin(4) * t215;
t441 = -(4 * t565) + (0.8e1 * t473 - t487) * t394 + (-0.8e1 * t472 - 0.4e1 * t792) * t389 - t843;
t494 = Ifges(5,4) * t854;
t667 = qJDD(1) * t544;
t800 = (t441 + 0.8e1 * t494 + 0.2e1 * t667) * t390 - 0.4e1 * t807;
t797 = pkin(2) * t391;
t796 = pkin(2) * t418;
t795 = pkin(3) * t158;
t785 = g(2) * t398;
t314 = -t763 / 0.2e1;
t709 = t390 * t419;
t580 = t158 * t709;
t102 = -t580 / 0.2e1 + t314;
t784 = t102 * pkin(2);
t313 = qJDD(2) + qJDD(3) / 0.2e1;
t376 = qJDD(4) / 0.2e1;
t249 = t376 + t313;
t259 = t387 + t673;
t718 = t259 * t338;
t783 = (pkin(4) * t718 - pkin(6) * t249) * m(6);
t213 = t252 * t419;
t590 = t280 * t709;
t782 = (t213 + t590) * pkin(2);
t325 = t339 ^ 2;
t244 = t325 + t419;
t780 = t244 * pkin(3);
t485 = t280 * pkin(3);
t778 = t420 * m(6);
t775 = Ifges(5,4) - Ifges(4,4);
t380 = Ifges(4,3) + Ifges(5,3);
t701 = t419 * t360;
t581 = t101 * t701;
t772 = t101 * t822 + 0.4e1 * t581;
t133 = t419 * t795;
t605 = t98 * t709;
t60 = t133 + 0.2e1 * t605;
t99 = t130 * t707;
t771 = -0.2e1 * t392 * t60 + t99;
t769 = mrSges(6,1) * t389;
t767 = mrSges(6,1) * t393;
t762 = mrSges(5,2) * t249;
t430 = t176 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t359 + Ifges(6,3) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - t649 - t806;
t717 = t280 * t390;
t625 = pkin(3) * t717;
t693 = t421 / 0.2e1 - t420 / 0.2e1;
t506 = (-t422 / 0.2e1 + t693) * m(6) - t777 / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + t430 - t625;
t752 = t506 * t706;
t65 = m(6) * t693 + t430;
t751 = t65 * t390;
t750 = (Ifges(6,4) * t276 + t384 * t479) * t359;
t749 = t98 * t360;
t748 = t130 * t705 + t893;
t531 = pkin(3) * t580;
t602 = t98 * t701;
t747 = -0.2e1 * t531 + 0.4e1 * t602;
t603 = t65 * t709;
t746 = -0.4e1 * t603 + t881;
t66 = -t870 + t544 + t845;
t604 = t66 * t709;
t745 = -0.2e1 * t604 + t220;
t743 = mrSges(5,1) * qJD(2);
t742 = mrSges(4,2) * qJD(1);
t741 = Ifges(6,4) * qJD(5);
t736 = t101 * t390;
t735 = (t526 + 0.2e1 * t533 + t805) * t392;
t118 = -qJDD(1) * t316 + t195 - t400;
t734 = t118 * t392;
t731 = (-mrSges(6,1) * t479 - mrSges(6,2) * t249) * t389;
t729 = t130 * t390;
t150 = t183 * t359;
t722 = t205 * t390;
t721 = t215 * t389;
t719 = (qJD(1) + t298) * (qJD(1) - t298);
t714 = t359 * t395;
t396 = cos(qJ(3));
t361 = t396 ^ 2;
t713 = t361 * t392;
t710 = t390 * t391;
t708 = t391 * t392;
t704 = t392 * t419;
t703 = t394 * t419;
t279 = t422 * t336;
t700 = -0.2e1 * t532 - 0.2e1 * t566;
t281 = -qJDD(2) / 0.2e1 - qJDD(3) / 0.2e1;
t689 = mrSges(6,1) * qJDD(2);
t688 = mrSges(4,2) * qJDD(1);
t684 = Ifges(6,3) * qJDD(1);
t681 = qJD(1) * t390;
t223 = mrSges(6,2) * t755 - Ifges(6,5) * qJD(4);
t680 = qJD(2) * (-Ifges(6,5) * qJD(3) + t223);
t638 = mrSges(6,1) * t755;
t679 = qJD(2) * (Ifges(6,6) * t338 + t638);
t678 = qJD(3) * t223;
t677 = qJD(3) * (Ifges(6,6) * qJD(4) + t638);
t675 = t271 * qJD(4);
t511 = t671 * t815 + t575;
t40 = (t259 * t457 - t743) * qJD(1) + t511;
t674 = t40 * t837;
t672 = qJD(1) * qJD(2);
t666 = qJDD(1) * t775;
t665 = qJDD(1) * t391;
t664 = qJDD(1) * t422;
t423 = pkin(2) ^ 2;
t663 = qJDD(2) * t423;
t662 = -0.4e1 * t797;
t661 = 0.2e1 * t795;
t660 = t281 * t836;
t655 = 0.2e1 * t390;
t650 = t43 * t793;
t647 = 0.2e1 * t738;
t645 = qJDD(4) + (2 * qJDD(2));
t103 = t186 + t722;
t432 = -g(3) * mrSges(4,2) + pkin(1) * t213 + pkin(3) * t876 + qJDD(1) * t315 + t856 * t390 + t574 * t393;
t606 = t65 * t701;
t518 = t391 * t606;
t644 = t103 * t705 + t432 * t392 - 0.4e1 * t518;
t637 = pkin(1) * t742;
t635 = pkin(2) * t236;
t634 = pkin(2) * t710;
t633 = pkin(2) * t706;
t632 = pkin(2) * t704;
t631 = pkin(3) * t708;
t628 = t391 * t782;
t626 = pkin(4) * t721;
t623 = -0.2e1 * t675;
t59 = 0.2e1 * t133 + 0.4e1 * t605;
t621 = -t59 * t391 + t748;
t620 = t95 * t822 + t747;
t583 = t391 * t705;
t125 = t205 * t583;
t262 = t376 + t336;
t475 = t819 / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1;
t619 = t125 + 0.2e1 * t603 - 0.4e1 * pkin(3) * (t516 / 0.2e1 + (t298 * t560 + mrSges(6,2) * t262 / 0.2e1) * t389 + (-mrSges(5,2) / 0.4e1 + mrSges(6,3) / 0.4e1 + t819 / 0.4e1) * t419 + t476 * (-qJDD(4) / 0.4e1 + t281) - t475 * t675);
t616 = pkin(3) * t358;
t614 = mrSges(6,1) * t704;
t609 = pkin(3) * t854;
t607 = mrSges(5,1) * t718;
t601 = -2 * t672;
t591 = pkin(2) * t672;
t589 = t280 * t704;
t588 = t360 * t708;
t587 = t392 * t701;
t585 = t390 * t706;
t584 = t390 * t702;
t582 = t391 * t704;
t579 = t316 * t706;
t577 = m(5) * t664;
t576 = t279 - t862;
t571 = t422 + t691;
t569 = t804 * t316;
t563 = t391 * t681;
t562 = t333 * t683;
t561 = t392 * t684;
t557 = t170 * t672;
t332 = Ifges(3,4) + t775;
t554 = t332 * t672;
t553 = t775 * t672;
t552 = t775 * t671;
t540 = -4 * Ifges(6,4) * t683;
t539 = 0.2e1 * t591;
t537 = -0.2e1 * t582;
t536 = 0.8e1 * t582;
t535 = Ifges(5,4) * t339;
t534 = qJD(1) * t657;
t177 = pkin(2) * t589;
t530 = pkin(3) * t590;
t528 = 0.2e1 * t564;
t527 = -Ifges(4,4) + t598;
t525 = -t415 + t623;
t522 = Ifges(6,4) * t582;
t521 = Ifges(6,4) * t390 * t704;
t520 = t392 * t609;
t519 = Ifges(6,4) * t585;
t514 = t392 * t591;
t513 = t317 * t389 * t701;
t509 = t614 * t835;
t507 = t413 / 0.2e1 + t319 / 0.2e1;
t504 = -m(6) * (-pkin(4) * qJDD(2) + pkin(6) * t418) - t696;
t502 = t582 * t826;
t501 = -0.2e1 * t513;
t500 = t779 + t786;
t499 = -g(1) * t393 + t785;
t498 = mrSges(4,2) * g(1) + t729;
t495 = -t324 - 0.2e1 * t473;
t493 = t584 * t758;
t491 = t769 + t902;
t489 = -0.4e1 * t736 - t795;
t483 = t391 * t513;
t478 = (t544 - t778) * qJDD(1);
t477 = t298 * t392 * t563;
t58 = t419 * (m(6) * t571 + t450 - t845);
t44 = -t58 - 0.2e1 * t530;
t471 = mrSges(5,3) + t491;
t470 = t316 * t879 + t321 - t328;
t466 = (-0.8e1 * t355 + 0.16e2 * t543) * pkin(6) + (-0.4e1 * t296 + (16 * t565)) * t359 + 0.4e1 * t329;
t462 = t280 * t718;
t460 = t585 * t835 + t796;
t458 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + t471;
t454 = -pkin(2) * qJDD(1) + t609 * t825;
t453 = -t577 - t578;
t282 = -2 * t565;
t449 = (-0.2e1 * t355 + 0.4e1 * t543) * pkin(6) + t150 + 0.4e1 * t532 + t282 + 0.4e1 * t566 + t329 + (0.2e1 * t473 - t488) * t880 + (-0.4e1 * t472 - 0.2e1 * t792) * t389;
t446 = (-t317 * t534 + t487) * t394 + (t316 * t534 + 0.4e1 * t792) * t389 + t843;
t444 = t407 * t423 - Ifges(3,1) + Ifges(3,2) + t243;
t440 = (-t751 + t485 / 0.2e1) * t536 + t861;
t437 = t444 - t777;
t436 = t498 * t398 + t536 * t749 - t885;
t435 = (-0.2e1 * t667 + (-0.8e1 * t535 + (4 * t741)) * qJD(1) + t446) * t360 + 0.4e1 * t552 - 0.4e1 * t650 + t449;
t434 = -qJDD(5) * Ifges(6,5) + pkin(6) * t85 + t276 * t811 + t389 * (-t276 * t384 + (t479 * t872));
t431 = -0.8e1 * t552 - 0.8e1 * t553 + 0.8e1 * t650;
t429 = t177 + t21;
t428 = (Ifges(6,6) * t888 + t592 * t834 + t434 - 0.2e1 * t677 - 0.2e1 * t679) * t394 + (t592 * t832 + 0.2e1 * t678 + 0.2e1 * t680 - t886) * t389 + mrSges(6,3) * t647 + mrSges(5,2) * t648 - 0.2e1 * t750 - t844 + 0.2e1 * (-t543 + t639) * pkin(1);
t121 = t158 * t633;
t122 = t158 * t632;
t397 = cos(qJ(2));
t362 = t397 ^ 2;
t425 = -(0.2e1 * (t133 * t395 + (t804 * t66 + (t101 * t823 + t485) * t390) * t419) * t713 + ((t122 + t861) * t395 + (t130 * t398 + t429) * t390 + ((t66 * t826 + 0.2e1 * t485) * t395 - 0.2e1 * t627 + t865) * t582) * t396 - 0.2e1 * t66 * t587 + (t391 * t429 + t489 * t704 + t99) * t395 + t707 * t722 + ((-t122 + t856) * t391 - pkin(3) * t589) * t390 + t66 * t704) * t397 - (((-t280 * t395 + t728) * pkin(2) + ((t661 + 0.8e1 * t736) * t395 + 0.2e1 * t625 + t803 * t66) * t391) * t419 * t396 + (t121 + t745) * t395 + (t280 * t633 - t133) * t390 + t772) * t362 - t231 * t419 - ((t748 - t809) * t395 + (t205 * t705 + t773) * t390 + (t395 * t489 + t66 * t829 - t625 + t66) * t706) * t396 - ((t809 - t893) * t391 + (-t244 * t492 + t318 * t325 + t261) * pkin(3) + t336 * t485) * t390 - (t125 + t604 + t891 - t220 - t325 * t485 - (-mrSges(5,1) + t457) * pkin(3) * t336) * t395 + t583 * t729 + ((mrSges(6,1) * t525 - t356 + t372) * t394 + (-t357 + (-t419 - t525) * mrSges(6,2)) * t389 + t261 + (-pkin(4) * t675 + pkin(6) * t262) * t839 + (0.2e1 * qJDD(3) + t645) * mrSges(6,3) + mrSges(5,1) * t623 - 0.2e1 * mrSges(5,2) * t262) * t793 + m(5) * t279 + t380 * qJDD(3) - 0.2e1 * t602 + 0.2e1 * t581 + t862 * m(6) + (t531 - t772 - ((t66 * t823 + t661) * t390 - t865) * t419 * t362 + t620 + (t362 * t881 - t745 + t746) * t395) * t361 + (-t276 + qJDD(4)) * Ifges(5,3);
t424 = pkin(1) ^ 2;
t385 = mrSges(4,3) + mrSges(5,3);
t275 = qJDD(3) + t645;
t207 = -0.4e1 * t519;
t196 = t595 * t882;
t194 = pkin(6) * t541;
t172 = t218 * g(1);
t171 = t218 * g(2);
t146 = t415 / 0.2e1 + t718;
t144 = pkin(4) * t274 - t787;
t94 = -g(2) * t767 + t470;
t93 = mrSges(6,2) * t855 + t495;
t92 = -pkin(2) * t614 + t118;
t90 = t197 + (-g(3) - t632) * mrSges(6,2) - t526;
t79 = t94 * t392;
t78 = t93 * t392;
t53 = (0.4e1 * t133 + 0.8e1 * t605) * t391;
t48 = -0.4e1 * t65 * t587;
t46 = 0.8e1 * t518;
t35 = (t259 * t819 - t386 * t673) * qJD(1) + t481;
t27 = qJDD(1) * t507 + t40 * t390 - t673 * t742;
t25 = t35 * t390 + t688 / 0.2e1 + (t252 * qJD(2) + qJD(3) * t507) * qJD(1);
t16 = (-0.8e1 * t749 + (0.8e1 * t751 - 0.4e1 * t485) * t395 + 0.4e1 * t627 + 0.4e1 * t95) * t419 * t361;
t13 = t885 * t392;
t12 = (t282 + t296 / 0.2e1) * t359 + (t196 + t488) * t394 + (t194 + t792) * t389 + t355 * pkin(6) - t329 / 0.2e1 + (t778 / 0.2e1 - t544 / 0.2e1) * qJDD(1) + (t741 - 0.2e1 * t535 + 0.2e1 * (-Ifges(6,6) * t389 - t342 + t756) * t298) * qJD(1) + t700;
t5 = -t204 * t398 + t428;
t4 = (mrSges(4,2) * g(2) + t129 * t390) * t398 + t427 * t390 - g(1) * t764 + t486 * Ifges(4,6) + (-0.2e1 * t688 - (t900 + t821) * qJD(1) * t252) * pkin(1) + (pkin(3) * t471 - Ifges(4,5)) * t336;
t1 = [0.2e1 * (t355 - t543) * pkin(6) + (t252 + t717) * t391 * t539 + (t801 * t395 + t435 + ((-t420 - t422) * m(6) - t450) * qJDD(1) + 0.4e1 * t553) * t361 + ((-0.2e1 * t280 * t591 + t5 * t392) * t395 + t4 * t392 + (mrSges(4,2) + t728) * t539 + (t800 * t395 - 0.2e1 * t557 + 0.2e1 * t666 + t897) * t391) * t396 + 0.4e1 * ((-t743 / 0.2e1 + (-t766 / 0.2e1 + t761 / 0.2e1 + t871) * t298) * qJD(1) + t511) * t793 + (0.8e1 * (t754 + (t12 * t390 + t807) * t395 - t651 - t510 / 0.2e1 + (-t188 * pkin(6) + t685 / 0.2e1 + (-t739 / 0.2e1 + t464) * qJD(1)) * t394 + (t187 * pkin(6) - t353 / 0.2e1 + (-t652 - t740 / 0.2e1) * qJD(1)) * t389 - t524 / 0.2e1 + (-t330 / 0.2e1 - t593) * pkin(6) + t523 / 0.2e1 - t613 / 0.2e1 + t557 / 0.2e1 - t517 / 0.2e1 - t556 / 0.2e1 - t544 * t545 - t666 / 0.2e1) * t713 + (0.8e1 * t12 * t588 - t898 * t395 + t428 * t390 + (qJD(3) * t315 - t637) * t821 - 0.2e1 * qJD(3) * t637 - t846) * t396 - t171 * t398 + t172 * t393 - Ifges(3,5) * t418 + (mrSges(3,2) * pkin(1) * t601) - (qJDD(2) * Ifges(3,6)) + (t5 * t395 + t4) * t391 + ((t694 + t760) * t394 + (-t237 + t765) * t389 + t385 * t418) * pkin(2) + (((t874 * t391 + t674 + ((t544 * t657) - 0.8e1 * t461 + (-0.8e1 * t870 + t847) * t298) * t563) * t395 + (-t431 + t441) * t391 - 0.4e1 * t25 * pkin(2) + 0.2e1 * t170 * t665) * t396 + (t35 * t662 + t800) * t395 + 0.2e1 * qJDD(1) * t332 + t27 * t662 + ((-t422 + t423) * m(6) + t437) * t601 + t897) * t392 + (t396 * t899 + t218) * t647) * t397 + (-Ifges(6,6) * t879 - pkin(1) * t237 + t194) * t389 + (Ifges(6,5) * t879 + pkin(1) * t694 + t196 - 0.2e1 * t725) * t394 + (t898 * t708 + (t448 + 0.2e1 * t686) * t390 + (t158 * t797 * t821 - t339 * t544 * t655) * qJD(1) + (-t216 * t394 + (t683 * t889) - t721 + t890) * t836) * t395 + ((t754 * t824 + (((-t466 + 0.4e1 * t578 - 0.4e1 * t667) * t390 + 0.8e1 * t807) * t391 + t35 * t837 + 0.8e1 * (t626 + t866) * t710 + (t298 * t841 + (8 * t741)) * t563) * t395 + (-0.8e1 * t651 - 0.4e1 * t666 + 0.4e1 * (qJD(2) * t170 + (-t450 - t799) * qJD(3)) * qJD(1) + t447) * t391 + t27 * t837) * t396 + ((-t683 * t841 + t466 + 0.4e1 * t478 - (8 * t565) - 0.8e1 * t626 - 0.8e1 * t866) * t360 + ((-t298 * t847 + t657 * t870 + 0.8e1 * t461 - (8 * t484)) * t681 - t874) * t395 + t664 * t839 + 0.2e1 * t577 + t529 - 0.2e1 * qJDD(1) * t243 + t431 + t446) * t361 + (t391 * t674 + t801) * t395 + t435 + t444 * qJDD(1) + (pkin(2) - pkin(3)) * (pkin(2) + pkin(3)) * t744 + 0.4e1 * t554 + t25 * t662 + t453) * t362 + ((t846 + (t389 * t886 - t434 * t394 - t869 * t888 + t844) * t390) * t391 + Ifges(3,6) * t418 + t499 * mrSges(3,2) + (-Ifges(3,5) + (t385 + t491) * pkin(2)) * qJDD(2) + ((t750 + (t677 + t679) * t394 + (-t678 - t680) * t389) * t390 - t315 * t668) * t825 + ((((m(6) * t790 - qJD(4) * t492) * t655 + mrSges(4,2) * t821 + qJD(3) * t833) * t391 - 0.2e1 * t218 * qJD(2)) * qJD(1) + ((-t280 * t655 - 0.2e1 * mrSges(4,1) - 0.2e1 * t413 - 0.2e1 * t414) * t391 - (2 * mrSges(3,2))) * qJDD(1)) * pkin(1)) * t392 + (g(1) * t850 + t458 * g(2)) * t393 + (t458 * g(1) - g(2) * t850) * t398 + (t449 + t478 + 0.4e1 * t494) * t360 + ((t340 * t424) + Ifges(3,1) + Ifges(5,1) + Ifges(4,2) + Ifges(6,2) + Ifges(2,3)) * qJDD(1) - t150 - 0.2e1 * t554 + t700 + (t422 + t424) * t744 - 0.2e1 * t552 + (2 * t565) - t453; (t16 + (t46 + (t53 - 0.2e1 * t635) * t395 - 0.4e1 * t752 - 0.4e1 * t784) * t396 + 0.2e1 * (t121 + (-0.2e1 * t751 + t485) * t419) * t395 + 0.2e1 * t628 + (Ifges(3,4) + t95) * t822 + t747) * t362 + (((t249 * t476 + t419 * t475 + t462 - t516 + t731) * t838 + t621) * t395 + 0.2e1 * t752 + (((-t372 / 0.2e1 - t146 * mrSges(6,1) - t356 / 0.2e1) * t394 + (-t357 / 0.2e1 + (t419 / 0.2e1 + t146) * mrSges(6,2)) * t389 + (t815 + t871) * t419 - t783 + t249 * mrSges(6,3) - t607 - t762) * t390 + t314 + (-t417 / 0.2e1 - t668) * mrSges(4,2) + t899 * t313) * t838 + t644) * t396 + ((t30 - pkin(2) * ((t372 + (t415 + 0.2e1 * t718) * mrSges(6,1) + t356) * t394 + (t357 + (-0.2e1 * t146 - t419) * mrSges(6,2)) * t389 + t261 + 0.2e1 * t783 - t275 * mrSges(6,3) + 0.2e1 * t607 + 0.2e1 * t762)) * t391 + t619) * t395 + (t218 * t798 - (mrSges(3,2) * g(3)) - Ifges(3,5) * qJDD(1) + t171 * t393 + (qJDD(1) * t385 + t876) * pkin(2)) * t392 + (0.4e1 * (t395 * t60 + (t65 * t828 - t506) * t419) * t713 + ((0.2e1 * t122 + t440) * t395 + (t391 * t831 + 0.2e1 * t782) * t392 + t436) * t396 + t48 + ((t21 + 0.2e1 * t177) * t391 + t771) * t395 + t103 * t707 + (0.4e1 * t392 * t784 + t432) * t391 + (-0.2e1 * t625 + (t423 - t571) * m(6) + t437 + t845) * t704 - qJDD(1) * Ifges(3,6) + (-t855 + t786) * mrSges(3,2) + t218 * g(3)) * t397 + (t576 + t663) * m(6) + (Ifges(3,4) + t527) * t419 + (t13 + (-(t275 * t476 + t236 + 0.2e1 * t462 - 0.2e1 * t516 + 0.2e1 * t731) * t390 - t213 - mrSges(4,1) * t417 - t313 * t833 - (mrSges(4,1) * t821 + 0.2e1 * t673 * t894) * qJD(3)) * pkin(2)) * t391 + t425 + (-t391 * t498 + t172) * t705 + t407 * t663 + (Ifges(3,3) + t380) * qJDD(2); (t619 + t891) * t395 + (t13 + ((-mrSges(5,1) * qJDD(2) - t394 * t689 - t504) * t390 - qJDD(2) * mrSges(4,2)) * pkin(2) + (-mrSges(5,2) * t390 + t899) * t796) * t391 + t527 * t419 + (t16 + (t46 + (t53 - t635) * t395 + (0.2e1 * t58 + 0.4e1 * t530) * t391 - 0.2e1 * t784) * t396 + (t121 + t746) * t395 + t628 + t620) * t362 + t425 + t380 * qJDD(2) - t498 * t583 + ((pkin(2) * (mrSges(5,2) * t418 + t504) + t621) * t395 + t44 * t391 + (mrSges(4,2) * t418 + t390 * t76) * pkin(2) + (t395 * t548 + t899) * t737 + t644) * t396 + (((t122 + t440) * t395 + t436) * t396 + t48 + t771 * t395 + (0.2e1 * (t395 * t59 - t44 + 0.4e1 * t606) * t361 + t396 * t782 + t44) * t392 + (t21 * t395 + t103 * t398 + (t396 * t831 + (t280 * t702 + 0.2e1 * t102) * pkin(2)) * t392 + t432) * t391) * t397 + t576 * m(6); ((-0.8e1 * t608 + (-0.4e1 * t727 - 0.2e1 * t326 + ((0.4e1 * t716 + 0.2e1 * t812) * t395 + 0.2e1 * t712) * t419) * t394 + 0.4e1 * t726 + (0.4e1 * t715 + 0.2e1 * t813) * t586 - 0.2e1 * t542 - 0.2e1 * t369) * t361 + (t493 * t824 + (mrSges(6,2) * pkin(2) * t395 + (t803 * t316 + (-0.2e1 * t812 - 0.4e1 * t463) * t390) * t391) * t703 + 0.4e1 * t483 + (t207 + (-0.4e1 * pkin(4) * t710 + pkin(2)) * t615) * t395 - 0.2e1 * t901 * t706 * t389) * t396 + (t326 + (-mrSges(6,2) * t634 - t712) * t419 + t864) * t394 + (pkin(4) - t634) * t615 + t369 + t852) * t362 + (-0.2e1 * (-0.4e1 * t493 + (t390 * t860 - t569) * t703 + t501 + 0.2e1 * (pkin(4) * t769 + Ifges(6,4)) * t584 + t901 * t711) * t713 + ((t390 * t528 + t522 * t802) * t359 + (0.4e1 * t169 * t588 + (t316 * t502 - t495) * t395 + (t384 * t537 + t390 * t540) * t389 + t90 * t390 + t616 + (t500 * t395 + ((pkin(1) - 0.2e1 * t631) * t395 + t708 * t835) * t419) * mrSges(6,2)) * t394 - 0.4e1 * t179 * t588 + ((t317 * t502 - t470) * t389 + t684 + (pkin(3) * t537 + t500) * t769) * t395 + (-pkin(3) * t687 + t390 * t92 + t391 * t509) * t389 - t390 * t562 + 0.2e1 * t522) * t396 + (t391 * t528 - 0.4e1 * t521) * t714 + ((t389 * t540 + t90) * t391 * t395 + ((-mrSges(6,2) * t786 + t93) * t391 + t860 * t704) * t390 - t569 * t704 - t454 * mrSges(6,1)) * t394 + t392 * t501 + ((t390 * t509 + t391 * t92) * t389 + 0.2e1 * t521 - t391 * t562) * t395 + (t317 * t704 + (t391 * t94 + (pkin(3) * t704 - g(1) * t707) * mrSges(6,1)) * t390 + t454 * mrSges(6,2)) * t389 - Ifges(6,3) * t390 * t665) * t397 + ((-t169 + t864) * t394 + t179 + t852) * t361 + ((t392 * t528 + t207) * t714 + (t579 * t829 + ((-0.2e1 * t384 * t585 + t392 * t540) * t389 - t735 + t460 * mrSges(6,2)) * t395 + t579 + t520 * t834 + (mrSges(6,2) * t859 + t78) * t390) * t394 - 0.2e1 * t483 + ((mrSges(6,1) * t460 + t734) * t389 + 0.2e1 * t519 - t392 * t562) * t395 + (t317 * t706 + t79 * t390 + t520 * t832 + t859 * t768) * t389 - t390 * t561) * t396 + (-t384 * t477 + (-t701 + t719) * Ifges(6,4)) * t830 + (-t727 + (t316 * t709 + t78 * t391 + (t391 * t467 + t780) * mrSges(6,2)) * t395 + (t384 * t719 + t477 * t872) * t389 - mrSges(6,1) * t785 + ((-pkin(2) * t760 + t735) * t391 + mrSges(6,2) * t660) * t390 - t616 * t708 + t514 * t834 + g(1) * t767 + t610 + t144 * mrSges(6,2) + t240) * t394 + t726 + ((mrSges(6,1) * t780 + t317 * t709) * t389 + ((-mrSges(6,1) * t629 - pkin(2) * t689 + t79) * t389 - t561) * t391) * t395 + (((-pkin(2) * t765 - t734) * t391 + mrSges(6,1) * t660) * t390 + t144 * mrSges(6,1) + t242 + (0.2e1 * t514 + (-pkin(1) + t631) * qJDD(1) + t499) * mrSges(6,2)) * t389 + t333 * t477 - t369 + t753 + t349;];
tau = t1(:);