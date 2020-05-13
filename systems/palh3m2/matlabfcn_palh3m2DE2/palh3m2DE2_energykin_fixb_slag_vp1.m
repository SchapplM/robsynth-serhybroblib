% Calculate kinetic energy for
% palh3m2DE2
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
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE2_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:15
% EndTime: 2020-05-07 02:14:23
% DurationCPUTime: 6.99s
% Computational Cost: add. (8250->368), mult. (12286->548), div. (0->0), fcn. (18403->69), ass. (0->237)
t794 = qJD(2) + qJD(3);
t818 = cos(qJ(1));
t757 = t794 * t818;
t812 = sin(qJ(1));
t758 = t794 * t812;
t800 = qJ(2) + qJ(3);
t782 = sin(t800);
t785 = cos(t800);
t847 = -Icges(4,5) * t785 + Icges(4,6) * t782;
t807 = sin(pkin(18));
t808 = cos(pkin(18));
t813 = sin(pkin(15));
t819 = cos(pkin(15));
t744 = t807 * t819 + t808 * t813;
t745 = -t807 * t813 + t808 * t819;
t811 = sin(qJ(2));
t817 = cos(qJ(2));
t671 = qJ(2) + atan2(t744 * t817 + t745 * t811, t744 * t811 - t745 * t817);
t810 = sin(qJ(3));
t816 = cos(qJ(3));
t748 = t810 * t819 + t813 * t816;
t749 = -t810 * t813 + t816 * t819;
t694 = t748 * t811 - t749 * t817;
t799 = pkin(17) + pkin(18);
t778 = sin(t799);
t779 = cos(t799);
t833 = t748 * t817 + t749 * t811;
t652 = pkin(17) - atan2(t694 * t778 - t779 * t833, t694 * t779 + t778 * t833) - t671;
t650 = sin(t652);
t651 = cos(t652);
t918 = Icges(9,5) * t650 - Icges(9,6) * t651;
t921 = -(-Icges(4,3) * t818 + t812 * t847) * t757 + (Icges(4,3) * t812 + t818 * t847) * t758 + (-Icges(4,5) * t782 - Icges(4,6) * t785 + t918) * qJD(1);
t904 = t812 ^ 2;
t903 = t818 ^ 2;
t920 = Icges(3,3) + Icges(7,3);
t814 = sin(pkin(14));
t820 = cos(pkin(14));
t747 = t813 * t820 - t814 * t819;
t751 = t813 * t814 + t819 * t820;
t692 = t747 * t817 + t751 * t811;
t866 = -t747 * t811 + t751 * t817;
t919 = Icges(3,5) * t817 + Icges(7,5) * t866 - Icges(3,6) * t811 - Icges(7,6) * t692;
t915 = Icges(7,4) * t866;
t788 = qJ(2) + pkin(15) + pkin(18);
t772 = pkin(17) + qJ(3) + t788;
t771 = pkin(16) + t772;
t874 = atan2(-sin(t771), cos(t771)) + t800;
t863 = -qJ(4) - t874;
t708 = qJ(1) - t863;
t700 = cos(t708);
t862 = -qJ(4) + t874;
t709 = qJ(1) - t862;
t701 = cos(t709);
t801 = qJ(1) + qJ(4);
t786 = cos(t801);
t795 = qJD(1) + qJD(4);
t914 = (-t701 / 0.4e1 - t700 / 0.4e1 - t786 / 0.2e1) * t795;
t710 = qJ(1) + t862;
t698 = sin(t710);
t711 = qJ(1) + t863;
t699 = sin(t711);
t802 = qJ(1) - qJ(4);
t784 = sin(t802);
t796 = qJD(1) - qJD(4);
t913 = (t699 / 0.4e1 + t698 / 0.4e1 - t784 / 0.2e1) * t796;
t912 = t818 * t812;
t911 = Icges(3,5) * t811 + Icges(7,5) * t692 + Icges(3,6) * t817 + Icges(7,6) * t866;
t910 = t812 * t919 - t818 * t920;
t909 = t812 * t920 + t818 * t919;
t890 = Icges(7,4) * t692;
t673 = Icges(7,2) * t866 + t890;
t674 = Icges(7,1) * t692 + t915;
t894 = Icges(3,4) * t811;
t760 = Icges(3,2) * t817 + t894;
t893 = Icges(3,4) * t817;
t761 = Icges(3,1) * t811 + t893;
t908 = -t673 * t692 + t674 * t866 - t760 * t811 + t761 * t817;
t850 = -Icges(7,2) * t692 + t915;
t661 = Icges(7,6) * t812 + t818 * t850;
t854 = Icges(7,1) * t866 - t890;
t663 = Icges(7,5) * t812 + t818 * t854;
t852 = -Icges(3,2) * t811 + t893;
t728 = Icges(3,6) * t812 + t818 * t852;
t856 = Icges(3,1) * t817 - t894;
t730 = Icges(3,5) * t812 + t818 * t856;
t907 = -t661 * t692 + t663 * t866 - t728 * t811 + t730 * t817;
t662 = -Icges(7,6) * t818 + t812 * t850;
t664 = -Icges(7,5) * t818 + t812 * t854;
t727 = -Icges(3,6) * t818 + t812 * t852;
t729 = -Icges(3,5) * t818 + t812 * t856;
t906 = t662 * t692 - t664 * t866 + t727 * t811 - t729 * t817;
t713 = qJ(1) - t874;
t705 = sin(t713);
t902 = t705 / 0.2e1;
t712 = qJ(1) + t874;
t706 = cos(t712);
t901 = -t706 / 0.2e1;
t791 = t812 * rSges(6,2);
t900 = t791 / 0.2e1;
t897 = pkin(1) * qJD(2);
t896 = t812 * rSges(6,1);
t895 = t818 * rSges(6,1);
t892 = Icges(4,4) * t782;
t891 = Icges(4,4) * t785;
t805 = sin(pkin(16));
t806 = cos(pkin(16));
t742 = t805 * t819 + t806 * t813;
t743 = -t805 * t813 + t806 * t819;
t690 = t742 * t816 + t743 * t810;
t691 = -t742 * t810 + t743 * t816;
t665 = -t690 * t811 + t691 * t817;
t836 = t690 * t817 + t691 * t811;
t647 = atan2(-t665 * t778 - t779 * t836, t665 * t779 - t778 * t836) + t800;
t645 = sin(t647);
t887 = t645 * t812;
t886 = t645 * t818;
t809 = sin(qJ(4));
t885 = t809 * t812;
t884 = t809 * t818;
t815 = cos(qJ(4));
t883 = t812 * t815;
t882 = t815 * t818;
t881 = qJD(1) * t812;
t880 = qJD(1) * t818;
t879 = qJD(2) * t812;
t878 = qJD(2) * t818;
t877 = qJD(4) * t645;
t737 = atan2(sin(t788), -cos(t788));
t876 = pkin(17) - t737;
t875 = t811 * t897;
t704 = sin(t712);
t873 = t704 / 0.2e1 + t902;
t872 = t902 - t704 / 0.2e1;
t707 = cos(t713);
t871 = -t707 / 0.2e1 + t901;
t870 = t707 / 0.2e1 + t901;
t689 = -atan2(-sin(t772), -cos(t772)) + t876;
t865 = t812 * t875;
t864 = t818 * t875;
t777 = t817 * t897;
t735 = -pkin(4) * t785 * t794 + t777;
t767 = rSges(3,1) * t817 - rSges(3,2) * t811;
t859 = -rSges(4,1) * t785 + rSges(4,2) * t782;
t734 = qJ(2) + t737;
t858 = rSges(8,1) * cos(t734) - rSges(8,2) * sin(t734);
t688 = -qJ(2) + t689;
t857 = -rSges(9,1) * cos(t688) - rSges(9,2) * sin(t688);
t855 = -Icges(4,1) * t785 + t892;
t851 = Icges(4,2) * t782 - t891;
t845 = -Icges(9,5) * t651 - Icges(9,6) * t650;
t646 = cos(t647);
t640 = t646 * t885 - t882;
t641 = -t646 * t883 - t884;
t642 = t646 * t884 + t883;
t643 = -t646 * t882 + t885;
t844 = (Icges(6,5) * t641 + Icges(6,6) * t640 - Icges(6,3) * t887) * t812 + (Icges(6,5) * t643 + Icges(6,6) * t642 - Icges(6,3) * t886) * t818;
t837 = -sin(t689) * (rSges(9,1) * t817 - rSges(9,2) * t811) + cos(t689) * (rSges(9,1) * t811 + rSges(9,2) * t817);
t764 = -rSges(7,1) * t819 + rSges(7,2) * t813;
t769 = rSges(7,1) * t813 + rSges(7,2) * t819;
t831 = t764 * t814 + t769 * t820;
t773 = pkin(1) * t817 + pkin(12);
t830 = t773 * t880 - t865;
t696 = sin(t708);
t697 = sin(t709);
t783 = sin(t801);
t829 = (t697 / 0.4e1 + t696 / 0.4e1 + t783 / 0.2e1) * t795;
t702 = cos(t710);
t703 = cos(t711);
t787 = cos(t802);
t828 = (-t703 / 0.4e1 - t702 / 0.4e1 + t787 / 0.2e1) * t796;
t827 = t844 * t645;
t826 = pkin(3) * cos(qJ(2) - t876) + t857;
t717 = -Icges(4,6) * t818 + t812 * t851;
t718 = Icges(4,6) * t812 + t818 * t851;
t719 = -Icges(4,5) * t818 + t812 * t855;
t720 = Icges(4,5) * t812 + t818 * t855;
t739 = -Icges(4,2) * t785 - t892;
t740 = -Icges(4,1) * t782 - t891;
t824 = (t718 * t782 - t720 * t785) * t758 - (t717 * t782 - t719 * t785) * t757 + (t739 * t782 - t740 * t785) * qJD(1);
t822 = qJD(2) ^ 2;
t821 = pkin(10) + rSges(6,3);
t804 = qJ(1) - qJ(2);
t803 = qJ(1) + qJ(2);
t798 = qJD(1) - qJD(2);
t797 = qJD(1) + qJD(2);
t792 = t818 * rSges(6,2);
t790 = qJ(1) - t800;
t789 = qJ(1) + t800;
t781 = qJD(1) - t794;
t780 = qJD(1) + t794;
t776 = pkin(12) * t880;
t775 = -t896 / 0.2e1;
t770 = t818 * t773;
t768 = rSges(2,1) * t818 - rSges(2,2) * t812;
t763 = rSges(2,1) * t812 + rSges(2,2) * t818;
t762 = rSges(3,1) * t811 + rSges(3,2) * t817;
t750 = -t810 * t811 + t816 * t817;
t746 = t810 * t817 + t811 * t816;
t741 = -rSges(4,1) * t782 - rSges(4,2) * t785;
t722 = rSges(4,3) * t812 + t818 * t859;
t721 = -rSges(4,3) * t818 + t812 * t859;
t714 = -t764 * t820 + t769 * t814;
t683 = t776 - t762 * t879 + qJD(1) * (rSges(3,3) * t812 + t767 * t818);
t682 = -t762 * t878 + (t818 * rSges(3,3) + (-pkin(12) - t767) * t812) * qJD(1);
t681 = t714 * t817 - t811 * t831;
t680 = t714 * t811 + t817 * t831;
t679 = rSges(7,3) * t812 + t681 * t818;
t678 = -rSges(7,3) * t818 + t681 * t812;
t676 = qJD(1) * t722 - t741 * t758 + t830;
t675 = -t864 - t741 * t757 + (-t812 * t773 - t721) * qJD(1);
t670 = cos(t671);
t669 = sin(t671);
t668 = -t865 + (rSges(8,3) * t812 + t818 * t858 + t770) * qJD(1);
t667 = -t864 + (t818 * rSges(8,3) + (-t773 - t858) * t812) * qJD(1);
t666 = t721 * t758 + t722 * t757 + t777;
t658 = t794 * t857 + t777;
t657 = -t680 * t879 + (-pkin(6) * t818 + t679) * qJD(1);
t656 = -t680 * t878 + (pkin(6) * t812 - t678) * qJD(1);
t655 = (t678 * t812 + t679 * t818) * qJD(2);
t654 = qJD(1) * (rSges(5,1) * t871 - rSges(5,2) * t872 + t812 * rSges(5,3)) + (t746 * t758 - t750 * t880) * pkin(4) + t830;
t653 = -t773 * t881 - t864 - qJD(1) * (-rSges(5,1) * t873 + rSges(5,2) * t870 - t818 * rSges(5,3)) + (t746 * t757 + t750 * t881) * pkin(4);
t649 = -t865 + t837 * t758 + (rSges(9,3) * t812 + t818 * t826 + t770) * qJD(1);
t648 = -t864 + t837 * t757 + (t818 * rSges(9,3) + (-t773 - t826) * t812) * qJD(1);
t644 = qJD(4) * t646 + qJD(1);
t632 = Icges(9,3) * t812 + t818 * t845;
t631 = -Icges(9,3) * t818 + t812 * t845;
t630 = t776 + (-t781 * cos(t790) / 0.2e1 - t780 * cos(t789) / 0.2e1) * pkin(4) + (t798 * cos(t804) / 0.2e1 + t797 * cos(t803) / 0.2e1) * pkin(1) + (pkin(8) * t871 + t821 * t872) * qJD(1) + (t829 - t913) * rSges(6,2) + (t828 + t914) * rSges(6,1);
t629 = (t781 * sin(t790) / 0.2e1 + t780 * sin(t789) / 0.2e1) * pkin(4) + (-t798 * sin(t804) / 0.2e1 - t797 * sin(t803) / 0.2e1) * pkin(1) + (pkin(8) * t873 - t812 * pkin(12) + t821 * t870) * qJD(1) + (t828 - t914) * rSges(6,2) + (t829 + t913) * rSges(6,1);
t628 = ((t792 + t896) * t787 + (t791 - t895) * t784 + (t792 - t896) * t786 + t783 * (t791 + t895) + (t702 + t703) * (t775 - t792 / 0.2e1) + (t701 + t700) * (t775 + t792 / 0.2e1) + (-t698 - t699) * (t900 - t895 / 0.2e1) + (t697 + t696) * (t900 + t895 / 0.2e1) + ((-t706 + t707) * t818 + (-t704 + t705) * t812) * rSges(6,3)) * t877 / 0.2e1 + t735;
t627 = Icges(6,5) * t646 + (-Icges(6,1) * t815 + Icges(6,4) * t809) * t645;
t626 = Icges(6,6) * t646 + (-Icges(6,4) * t815 + Icges(6,2) * t809) * t645;
t625 = Icges(6,3) * t646 + (-Icges(6,5) * t815 + Icges(6,6) * t809) * t645;
t624 = Icges(6,1) * t643 + Icges(6,4) * t642 - Icges(6,5) * t886;
t623 = Icges(6,1) * t641 + Icges(6,4) * t640 - Icges(6,5) * t887;
t622 = Icges(6,4) * t643 + Icges(6,2) * t642 - Icges(6,6) * t886;
t621 = Icges(6,4) * t641 + Icges(6,2) * t640 - Icges(6,6) * t887;
t1 = m(5) * (t653 ^ 2 + t654 ^ 2 + t735 ^ 2) / 0.2e1 + m(6) * (t628 ^ 2 + t629 ^ 2 + t630 ^ 2) / 0.2e1 + m(3) * (t767 ^ 2 * t822 + t682 ^ 2 + t683 ^ 2) / 0.2e1 + m(8) * (pkin(1) ^ 2 * t817 ^ 2 * t822 + t667 ^ 2 + t668 ^ 2) / 0.2e1 + t644 * ((t646 * t625 + (t626 * t809 - t627 * t815) * t645) * t644 + ((-(t622 * t809 - t624 * t815) * t818 - (t621 * t809 - t623 * t815) * t812) * t645 - t844 * t646) * t877) / 0.2e1 + m(7) * (t655 ^ 2 + t656 ^ 2 + t657 ^ 2) / 0.2e1 + m(9) * (t648 ^ 2 + t649 ^ 2 + t658 ^ 2) / 0.2e1 + m(4) * (t666 ^ 2 + t675 ^ 2 + t676 ^ 2) / 0.2e1 + ((t909 * t904 + (t906 * t818 + (t907 - t910) * t812) * t818) * qJD(2) + (t812 * t911 + t818 * t908) * qJD(1)) * t879 / 0.2e1 - ((t910 * t903 + (t907 * t812 + (t906 - t909) * t818) * t812) * qJD(2) + (t812 * t908 - t818 * t911) * qJD(1)) * t878 / 0.2e1 - (t818 * ((-t625 * t886 + t626 * t642 + t627 * t643) * t644 + (-(t621 * t642 + t623 * t643) * t812 + (-t622 * t642 - t624 * t643 + t827) * t818) * t877) + t812 * ((-t625 * t887 + t626 * t640 + t627 * t641) * t644 + (-(t622 * t640 + t624 * t641) * t818 + (-t621 * t640 - t623 * t641 + t827) * t812) * t877)) * t877 / 0.2e1 + (t646 ^ 2 * Icges(5,2) + (Icges(5,1) * t645 + 0.2e1 * Icges(5,4) * t646) * t645 + Icges(8,2) * t670 ^ 2 + (Icges(8,1) * t669 + 0.2e1 * Icges(8,4) * t670) * t669 + Icges(2,3) + m(2) * (t763 ^ 2 + t768 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t718 * t785 - t720 * t782) * t758 - (-t717 * t785 - t719 * t782) * t757 + ((-t662 * t866 - t664 * t692 - t727 * t817 - t729 * t811) * t818 + (t661 * t866 + t663 * t692 + t728 * t817 + t730 * t811) * t812) * qJD(2) + (t651 ^ 2 * Icges(9,2) + t673 * t866 + t692 * t674 - t785 * t739 - t782 * t740 + t817 * t760 + t811 * t761 + (Icges(9,1) * t650 - 0.2e1 * Icges(9,4) * t651) * t650) * qJD(1) + (t903 + t904) * t794 * t918) * qJD(1) / 0.2e1 - ((t903 * t631 - t632 * t912) * t794 + t824 * t812 - t921 * t818) * t757 / 0.2e1 + ((-t631 * t912 + t904 * t632) * t794 + t824 * t818 + t921 * t812) * t758 / 0.2e1;
T = t1;
