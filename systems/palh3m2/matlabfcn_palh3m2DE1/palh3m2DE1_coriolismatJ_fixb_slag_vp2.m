% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2DE1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:03:20
% EndTime: 2020-05-07 02:03:22
% DurationCPUTime: 1.66s
% Computational Cost: add. (2434->234), mult. (2480->296), div. (0->0), fcn. (1310->80), ass. (0->152)
t854 = qJ(3) + qJ(2);
t749 = sin(t854);
t761 = pkin(17) + pkin(18);
t830 = pkin(15) + t761;
t819 = (pkin(16) + t830);
t809 = qJ(2) + t819;
t803 = qJ(4) + t809;
t800 = qJ(3) + t803;
t810 = -qJ(2) + t819;
t804 = -qJ(4) + t810;
t801 = -qJ(3) + t804;
t784 = -sin(t801) / 0.4e1 + sin(t800) / 0.4e1;
t785 = cos(t801) / 0.4e1 + cos(t800) / 0.4e1;
t822 = 0.2e1 * t854;
t796 = sin(t822) / 0.2e1;
t805 = qJ(3) + t809;
t806 = -qJ(3) + t810;
t821 = -qJ(2) + t830;
t814 = -qJ(3) + t821;
t820 = qJ(2) + t830;
t815 = qJ(3) + t820;
t877 = m(5) + m(6);
t842 = t877 * pkin(4);
t818 = t842 + mrSges(4,1) + mrSges(9,1);
t825 = -pkin(10) * m(6) + mrSges(5,2) - mrSges(6,3);
t832 = cos(t854);
t838 = pkin(8) * m(6) + mrSges(5,1);
t865 = mrSges(4,2) + mrSges(9,2);
t729 = qJ(4) + t810;
t709 = -qJ(3) + t729;
t871 = cos(t709);
t730 = -qJ(4) + t809;
t708 = qJ(3) + t730;
t872 = cos(t708);
t883 = -sin(t708) / 0.4e1 + sin(t709) / 0.4e1;
t893 = (-Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2)) * t796 + ((cos(t815) / 0.2e1 + cos(t814) / 0.2e1) * mrSges(9,2) + (sin(t815) / 0.2e1 - sin(t814) / 0.2e1) * mrSges(9,1)) * pkin(3) - (Ifges(4,4) + Ifges(9,4)) * cos(t822) + (-t818 * t749 - t865 * t832) * pkin(12) + (t796 * t842 + (-t871 / 0.4e1 - t872 / 0.4e1 + t785) * mrSges(6,2) + (t784 - t883) * mrSges(6,1) + (sin(t805) / 0.2e1 - sin(t806) / 0.2e1) * t838 + (-cos(t805) / 0.2e1 + cos(t806) / 0.2e1) * t825) * pkin(4);
t880 = 2 * qJ(2);
t892 = sin(t880) / 0.2e1;
t769 = sin(qJ(2));
t772 = cos(qJ(3));
t837 = -pkin(4) * t772 + pkin(1);
t768 = sin(qJ(3));
t773 = cos(qJ(2));
t857 = t768 * t773;
t657 = pkin(4) * t857 - t837 * t769;
t858 = t768 * t769;
t737 = pkin(4) * t858;
t658 = t837 * t773 + t737;
t770 = sin(pkin(15));
t774 = cos(pkin(15));
t633 = t657 * t774 - t658 * t770;
t764 = sin(pkin(16));
t765 = cos(pkin(16));
t827 = t657 * t770 + t658 * t774;
t891 = t633 * t764 + t765 * t827;
t890 = t633 * t765 - t764 * t827;
t856 = t772 * t773;
t691 = t856 - t858;
t692 = t769 * t772 + t857;
t860 = t692 * t774;
t639 = t691 * t770 + t860;
t861 = t692 * t770;
t826 = t691 * t774 - t861;
t889 = -t639 * t764 + t765 * t826;
t888 = t639 * t765 + t764 * t826;
t879 = 2 * qJ(4);
t758 = -qJ(2) + pkin(14) - pkin(15);
t878 = 0.2e1 * t758;
t876 = pkin(4) * mrSges(6,1);
t875 = pkin(4) * mrSges(6,2);
t874 = pkin(8) * mrSges(6,1);
t873 = pkin(8) * mrSges(6,2);
t870 = cos(t758);
t867 = pkin(8) * t770;
t866 = pkin(8) * t774;
t767 = sin(qJ(4));
t863 = t767 * mrSges(6,1);
t771 = cos(qJ(4));
t862 = t771 * mrSges(6,2);
t723 = mrSges(6,1) * t771 - mrSges(6,2) * t767;
t859 = t723 * t657;
t855 = t880 + qJ(3);
t853 = qJ(4) - qJ(2);
t852 = qJD(4) * t723;
t690 = t818 * pkin(1) * sin(t855);
t728 = pkin(1) * t865 * cos(t855);
t743 = sin(t758);
t786 = sin(t804) / 0.4e1 - sin(t803) / 0.4e1;
t787 = -cos(t804) / 0.4e1 - cos(t803) / 0.4e1;
t844 = pkin(18) + pkin(15);
t833 = qJ(2) + t844;
t834 = -qJ(2) + t844;
t835 = -cos(t729) / 0.4e1 - cos(t730) / 0.4e1;
t836 = -sin(t729) / 0.4e1 + sin(t730) / 0.4e1;
t620 = -t690 - Ifges(7,4) * cos(t878) + (-Ifges(3,1) + Ifges(3,2)) * t892 - Ifges(3,4) * cos(t880) - (-Ifges(7,1) + Ifges(7,2)) * sin(t878) / 0.2e1 - t728 + (mrSges(3,1) * t769 + mrSges(3,2) * t773) * pkin(12) + ((cos(t833) / 0.2e1 - cos(t834) / 0.2e1) * mrSges(8,2) + (sin(t834) / 0.2e1 - sin(t833) / 0.2e1) * mrSges(8,1) + (sin(t821) / 0.2e1 - sin(t820) / 0.2e1) * m(9) * pkin(3) + (t787 - t835) * mrSges(6,2) + (t786 - t836) * mrSges(6,1) + (-sin(t809) / 0.2e1 + sin(t810) / 0.2e1) * t838 + (pkin(1) * t892 + pkin(12) * t769) * (m(4) + m(8) + m(9) + t877) + (cos(t809) / 0.2e1 - cos(t810) / 0.2e1) * t825) * pkin(1) + (mrSges(7,1) * t743 - mrSges(7,2) * t870) * pkin(6) + t893;
t851 = t620 * qJD(1);
t754 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t755 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t811 = qJ(4) + t819;
t798 = 2 * t811;
t812 = -qJ(4) + t819;
t799 = 2 * t812;
t823 = 2 * t819;
t816 = qJ(4) + t823;
t817 = -qJ(4) + t823;
t621 = -(t754 + t874) * sin(t817) / 0.4e1 - (-t755 - t873) * cos(t816) / 0.4e1 + (t863 / 0.2e1 + t862 / 0.2e1) * pkin(8) + (-t755 + t873) * cos(t817) / 0.4e1 + (-t754 + t874) * sin(t816) / 0.4e1 + (cos(t879) / 0.2e1 - cos(t799) / 0.4e1 - cos(t798) / 0.4e1) * Ifges(6,4) + ((-cos(t811) / 0.2e1 - cos(t812) / 0.2e1) * mrSges(6,2) + (sin(t812) / 0.2e1 - sin(t811) / 0.2e1) * mrSges(6,1)) * pkin(12) + (t784 * mrSges(6,1) + t785 * mrSges(6,2)) * pkin(4) + ((t787 + t835) * mrSges(6,2) + (t786 + t836) * mrSges(6,1)) * pkin(1) + t883 * t876 + (t872 + t871) * t875 / 0.4e1 + (sin(t799) / 0.8e1 - sin(t798) / 0.8e1 + sin(t879) / 0.4e1) * (Ifges(6,1) - Ifges(6,2));
t850 = t621 * qJD(1);
t807 = t818 * t768;
t831 = t865 * t772;
t624 = -t690 / 0.2e1 - t728 / 0.2e1 + (-t831 / 0.2e1 - t807 / 0.2e1) * pkin(1) + t893;
t849 = t624 * qJD(1);
t762 = qJ(2) + qJ(4);
t750 = sin(t762);
t751 = -sin(t853);
t752 = cos(t762);
t753 = cos(t853);
t756 = qJ(4) + t854;
t741 = sin(t756);
t757 = qJ(3) - t853;
t742 = sin(t757);
t744 = cos(t756);
t745 = cos(t757);
t795 = (-t745 / 0.4e1 + t744 / 0.4e1) * mrSges(6,2) + (t742 / 0.4e1 + t741 / 0.4e1) * mrSges(6,1);
t792 = t795 * pkin(4);
t783 = t792 + ((t753 / 0.4e1 - t752 / 0.4e1) * mrSges(6,2) + (-t751 / 0.4e1 - t750 / 0.4e1) * mrSges(6,1)) * pkin(1);
t628 = t859 / 0.2e1 + t783;
t848 = t628 * qJD(1);
t630 = t749 * (Ifges(4,6) + Ifges(9,6)) + (pkin(4) * mrSges(5,3) - Ifges(4,5) - Ifges(9,5)) * t832 + (t741 / 0.2e1 - t742 / 0.2e1) * t876 + (t744 + t745) * t875 / 0.2e1;
t847 = t630 * qJD(3);
t802 = t723 * t692;
t631 = (t802 / 0.2e1 + t795) * pkin(4);
t846 = t631 * qJD(1);
t722 = -t862 - t863;
t747 = sin(t761);
t748 = cos(t761);
t841 = t722 * pkin(4) * (-t888 * t747 + t889 * t748);
t824 = -t841 / 0.2e1;
t808 = mrSges(6,1) * t867 + t755 * t774;
t793 = -t774 * t856 + t861;
t794 = t770 * t856 + t860;
t779 = (((-t764 * t770 + t765 * t774) * t748 + t747 * (-t764 * t774 - t765 * t770)) * t737 + ((t794 * t764 + t793 * t765) * t748 + t747 * (-t764 * t793 + t794 * t765)) * pkin(4)) * t722;
t622 = t824 - t779 / 0.2e1;
t649 = (t807 + t831) * pkin(1);
t797 = -qJD(2) * t649 + qJD(4) * t622;
t684 = mrSges(6,2) * t867 + t754 * t774;
t683 = mrSges(6,1) * t866 - t755 * t770;
t682 = mrSges(6,2) * t866 - t754 * t770;
t645 = t649 * qJD(3);
t632 = -pkin(4) * t802 / 0.2e1 + t792;
t629 = -t859 / 0.2e1 + t783;
t623 = t779 / 0.2e1 + t824;
t1 = [-qJD(2) * t620 - qJD(3) * t624 - qJD(4) * t621, t629 * qJD(4) + t847 - t851 + (Ifges(3,5) * t773 + Ifges(7,5) * t870 - Ifges(3,6) * t769 + Ifges(7,6) * t743 + t630 + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t773 + (-t753 / 0.2e1 - t752 / 0.2e1) * mrSges(6,2) + (t751 / 0.2e1 - t750 / 0.2e1) * mrSges(6,1)) * pkin(1)) * qJD(2), t630 * qJD(2) + t632 * qJD(4) + t847 - t849, -t850 + t629 * qJD(2) + t632 * qJD(3) + ((-(t683 * t765 - t764 * t808) * t767 + (-t682 * t765 + t684 * t764) * t771) * t748 + (-(-t764 * t683 - t808 * t765) * t767 + (t682 * t764 + t684 * t765) * t771) * t747 - t722 * (pkin(12) + t658)) * qJD(4); qJD(4) * t628 + t851, t645, t645 - t797, t848 - t622 * qJD(3) - (t747 * t891 - t890 * t748) * t852; qJD(4) * t631 + t849, t797, 0, t846 + t622 * qJD(2) + pkin(4) * (t889 * t747 + t888 * t748) * t852; -qJD(2) * t628 - qJD(3) * t631 + t850, -t848 + (t747 * t890 + t891 * t748) * t722 * qJD(2) + t623 * qJD(3), qJD(2) * t623 - qJD(3) * t841 - t846, 0;];
Cq = t1;
