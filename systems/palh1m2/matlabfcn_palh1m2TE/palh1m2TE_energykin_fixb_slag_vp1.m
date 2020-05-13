% Calculate kinetic energy for
% palh1m2TE
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:10:31
% EndTime: 2020-05-01 20:10:44
% DurationCPUTime: 9.69s
% Computational Cost: add. (2847->430), mult. (4903->691), div. (0->0), fcn. (5736->24), ass. (0->265)
t907 = -Icges(3,4) - Icges(10,4);
t906 = Icges(3,1) + Icges(10,1);
t905 = Icges(3,2) + Icges(10,2);
t780 = cos(qJ(2));
t904 = t907 * t780;
t774 = sin(qJ(2));
t903 = t907 * t774;
t902 = Icges(3,5) + Icges(10,5);
t901 = -Icges(3,6) - Icges(10,6);
t900 = -t905 * t780 + t903;
t899 = -t906 * t774 + t904;
t898 = Icges(9,3) + Icges(4,3);
t775 = sin(qJ(1));
t781 = cos(qJ(1));
t897 = t900 * t775 + t901 * t781;
t896 = t901 * t775 - t900 * t781;
t895 = t899 * t775 - t902 * t781;
t894 = t902 * t775 + t899 * t781;
t893 = t905 * t774 + t904;
t892 = t906 * t780 + t903;
t768 = sin(pkin(19));
t771 = cos(pkin(19));
t714 = t768 * t774 - t771 * t780;
t718 = t768 * t780 + t771 * t774;
t773 = sin(qJ(3));
t779 = cos(qJ(3));
t653 = t714 * t779 + t718 * t773;
t850 = t779 * t780;
t721 = -t773 * t774 + t850;
t851 = t774 * t779;
t723 = t773 * t780 + t851;
t891 = -Icges(4,5) * t721 + Icges(9,5) * t653 + Icges(4,6) * t723;
t890 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t776 = sin(pkin(18));
t777 = sin(pkin(17));
t782 = cos(pkin(18));
t783 = cos(pkin(17));
t722 = t776 * t783 - t777 * t782;
t724 = t776 * t777 + t782 * t783;
t658 = t722 * t774 + t724 * t780;
t659 = t722 * t780 - t724 * t774;
t889 = Icges(7,5) * t659 - Icges(7,6) * t658 - t902 * t774 + t901 * t780;
t767 = sin(pkin(20));
t770 = cos(pkin(20));
t716 = -t782 * t767 + t770 * t776;
t720 = t767 * t776 + t770 * t782;
t762 = pkin(22) + pkin(21);
t752 = sin(t762);
t753 = cos(t762);
t645 = t716 * t752 + t720 * t753;
t887 = Icges(5,4) * t645;
t886 = Icges(7,4) * t658;
t861 = Icges(9,4) * t653;
t885 = t781 * qJD(1);
t869 = pkin(1) * t774;
t763 = qJD(2) + qJD(3);
t728 = t763 * t781;
t713 = t768 * t779 + t771 * t773;
t717 = t768 * t773 - t771 * t779;
t652 = -t713 * t780 + t717 * t774;
t655 = -t714 * t773 + t718 * t779;
t755 = qJD(2) * t775;
t727 = qJD(3) * t775 + t755;
t884 = (t898 * t781 + (-Icges(9,6) * t652 + t891) * t775) * t728 + ((-Icges(9,6) * t655 - t891) * t781 + t898 * t775) * t727 + (Icges(4,5) * t723 + Icges(9,5) * t655 + Icges(4,6) * t721 - Icges(9,6) * t653) * qJD(1);
t883 = Icges(7,5) * t658 + Icges(7,6) * t659 + t901 * t774 + t902 * t780;
t882 = t889 * t775 - t890 * t781;
t881 = t890 * t775 + t889 * t781;
t610 = Icges(7,2) * t659 + t886;
t862 = Icges(7,4) * t659;
t611 = Icges(7,1) * t658 + t862;
t880 = -t610 * t658 + t611 * t659 - t892 * t774 + t893 * t780;
t823 = -Icges(7,2) * t658 + t862;
t597 = Icges(7,6) * t775 + t781 * t823;
t828 = Icges(7,1) * t659 - t886;
t599 = Icges(7,5) * t775 + t781 * t828;
t879 = -t597 * t658 + t599 * t659 - t894 * t774 + t896 * t780;
t596 = -Icges(7,6) * t781 + t775 * t823;
t598 = -Icges(7,5) * t781 + t775 * t828;
t878 = t596 * t658 - t598 * t659 + t895 * t774 + t897 * t780;
t877 = t775 ^ 2;
t876 = t781 ^ 2;
t871 = pkin(15) * t775;
t870 = (-t713 * t774 - t717 * t780) * pkin(2);
t868 = Icges(9,1) * t653;
t865 = Icges(4,4) * t721;
t864 = Icges(4,4) * t723;
t646 = -t716 * t753 + t720 * t752;
t863 = Icges(5,4) * t646;
t860 = Icges(9,4) * t655;
t592 = -Icges(9,5) * t781 + (Icges(9,4) * t652 - t868) * t775;
t856 = t592 * t653;
t593 = Icges(9,5) * t775 + (-t860 - t868) * t781;
t855 = t593 * t653;
t604 = Icges(9,1) * t655 - t861;
t854 = t604 * t653;
t853 = t646 * t775;
t852 = t646 * t781;
t746 = t763 * t775;
t765 = sin(pkin(22));
t769 = cos(pkin(22));
t715 = -t765 * t776 - t769 * t782;
t719 = t765 * t782 - t769 * t776;
t849 = (cos(pkin(21)) * t715 + t719 * sin(pkin(21))) * pkin(4);
t706 = t869 * t775;
t707 = t869 * t781;
t844 = qJD(2) * t781;
t648 = -t706 * t755 - t707 * t844;
t750 = pkin(5) * t773 + pkin(1);
t848 = pkin(5) * t850 - t750 * t774 + t869;
t751 = pkin(15) * t885;
t847 = -qJD(1) * t707 + t751;
t764 = qJ(3) + qJ(2);
t757 = sin(t764);
t846 = Icges(11,4) * t757;
t758 = cos(t764);
t845 = Icges(11,4) * t758;
t843 = qJD(4) * t646;
t842 = pkin(1) * qJD(2) * t780;
t703 = -t746 + t727;
t839 = t706 - t871;
t725 = -rSges(9,1) * t771 + rSges(9,2) * t768;
t726 = rSges(9,1) * t768 + rSges(9,2) * t771;
t838 = -t725 * t773 + t726 * t779;
t837 = t781 * t842;
t661 = t848 * t775;
t836 = -t661 + t839;
t835 = -rSges(3,1) * t774 - rSges(3,2) * t780;
t772 = sin(qJ(4));
t778 = cos(qJ(4));
t834 = rSges(6,1) * t778 - rSges(6,2) * t772;
t833 = -rSges(10,1) * t774 - rSges(10,2) * t780;
t832 = rSges(11,1) * t758 - rSges(11,2) * t757;
t830 = Icges(4,1) * t721 - t864;
t825 = -Icges(4,2) * t723 + t865;
t628 = rSges(6,3) * t720 - t716 * t834;
t629 = rSges(6,3) * t716 + t720 * t834;
t813 = t628 * t752 - t629 * t753;
t736 = -t776 * rSges(5,1) - t782 * rSges(5,2);
t743 = rSges(5,1) * t782 - rSges(5,2) * t776;
t668 = -t736 * t767 + t743 * t770;
t801 = t736 * t770 + t743 * t767;
t812 = -t668 * t753 + t752 * t801;
t744 = -pkin(9) * t776 + pkin(11) * t782;
t745 = pkin(9) * t782 + pkin(11) * t776;
t672 = t744 * t770 + t745 * t767;
t673 = -t744 * t767 + t745 * t770;
t811 = t672 * t752 - t673 * t753;
t798 = -Icges(11,2) * t757 + t845;
t677 = -Icges(11,6) * t781 + t775 * t798;
t799 = Icges(11,1) * t758 - t846;
t679 = -Icges(11,5) * t781 + t775 * t799;
t810 = t677 * t757 - t679 * t758;
t678 = Icges(11,6) * t775 + t781 * t798;
t680 = Icges(11,5) * t775 + t781 * t799;
t809 = -t678 * t757 + t680 * t758;
t709 = Icges(11,2) * t758 + t846;
t710 = Icges(11,1) * t757 + t845;
t804 = -t709 * t757 + t710 * t758;
t741 = -t776 * rSges(7,1) + rSges(7,2) * t782;
t742 = rSges(7,1) * t782 + rSges(7,2) * t776;
t800 = t741 * t777 - t742 * t783;
t797 = Icges(11,5) * t758 - Icges(11,6) * t757;
t683 = pkin(5) * t851 + (-pkin(1) + t750) * t780;
t796 = -t728 * t683 - t837;
t662 = t848 * t781;
t795 = t727 * t661 + t662 * t728 + t648;
t794 = -(rSges(11,1) * t757 + rSges(11,2) * t758) * t763 - t842;
t793 = -t775 * t842 + t847;
t792 = t772 * t645;
t791 = t645 * t778;
t787 = qJD(1) * t662 - t683 * t727 + t793;
t640 = -Icges(4,6) * t781 + t775 * t825;
t641 = Icges(4,6) * t775 + t781 * t825;
t642 = -Icges(4,5) * t781 + t775 * t830;
t643 = Icges(4,5) * t775 + t781 * t830;
t664 = Icges(4,2) * t721 + t864;
t665 = Icges(4,1) * t723 + t865;
t785 = (-t641 * t723 + t643 * t721) * t727 - (-t640 * t723 + t642 * t721) * t728 + (-t664 * t723 + t665 * t721) * qJD(1);
t740 = rSges(2,1) * t781 - rSges(2,2) * t775;
t739 = rSges(3,1) * t780 - rSges(3,2) * t774;
t738 = rSges(10,1) * t780 - rSges(10,2) * t774;
t737 = rSges(6,1) * t772 + rSges(6,2) * t778;
t735 = rSges(2,1) * t775 + rSges(2,2) * t781;
t708 = Icges(11,5) * t757 + Icges(11,6) * t758;
t699 = rSges(3,3) * t775 + t781 * t835;
t698 = rSges(10,3) * t775 + t781 * t833;
t697 = -rSges(3,3) * t781 + t775 * t835;
t696 = -rSges(10,3) * t781 + t775 * t833;
t682 = rSges(11,3) * t775 + t781 * t832;
t681 = -rSges(11,3) * t781 + t775 * t832;
t676 = Icges(11,3) * t775 + t781 * t797;
t675 = -Icges(11,3) * t781 + t775 * t797;
t674 = (rSges(4,1) * t774 + rSges(4,2) * t780) * t779 + t773 * (rSges(4,1) * t780 - rSges(4,2) * t774);
t671 = -t741 * t783 - t742 * t777;
t670 = (rSges(4,1) * t779 - rSges(4,2) * t773) * t780 - t774 * (rSges(4,1) * t773 + rSges(4,2) * t779);
t669 = (-rSges(8,1) * t782 + rSges(8,2) * t776) * t769 - t765 * (rSges(8,1) * t776 + rSges(8,2) * t782);
t666 = -t725 * t779 - t726 * t773;
t651 = rSges(4,3) * t775 + t670 * t781;
t650 = -rSges(4,3) * t781 + t670 * t775;
t649 = t652 * pkin(2);
t636 = qJD(1) * t699 - t739 * t755 + t751;
t635 = -t739 * t844 + (-t697 - t871) * qJD(1);
t634 = qJD(4) * t645 + qJD(1);
t632 = t870 * t781;
t631 = t870 * t775;
t630 = (t697 * t775 + t699 * t781) * qJD(2);
t627 = t772 * t775 - t781 * t791;
t626 = t775 * t778 + t781 * t792;
t625 = -t772 * t781 - t775 * t791;
t624 = t775 * t792 - t778 * t781;
t623 = t671 * t774 - t780 * t800;
t622 = t671 * t780 + t774 * t800;
t621 = t775 * t843;
t620 = t781 * t843 + t703;
t619 = t774 * t666 + t780 * t838;
t618 = t666 * t780 - t774 * t838;
t615 = rSges(7,3) * t775 + t622 * t781;
t614 = -rSges(7,3) * t781 + t622 * t775;
t613 = qJD(1) * (rSges(8,3) * t775 + t669 * t781) + t793;
t612 = -t837 + (t781 * rSges(8,3) + t706 + (-pkin(15) - t669) * t775) * qJD(1);
t607 = t811 * t775;
t606 = rSges(9,3) * t775 + t618 * t781;
t605 = -rSges(9,3) * t781 + t618 * t775;
t603 = -Icges(9,2) * t653 + t860;
t601 = -rSges(5,3) * t781 + t775 * t812;
t591 = Icges(9,6) * t775 + (-Icges(9,2) * t655 - t861) * t781;
t590 = -Icges(9,6) * t781 + (Icges(9,2) * t652 - t861) * t775;
t587 = Icges(5,1) * t646 - t887;
t586 = -Icges(5,2) * t645 + t863;
t584 = Icges(5,5) * t775 + t781 * (-Icges(5,1) * t645 - t863);
t582 = Icges(5,6) * t775 + (-Icges(5,2) * t646 - t887) * t781;
t578 = (t681 * t775 + t682 * t781) * t763 + t648;
t577 = qJD(1) * t651 - t674 * t727 + t793;
t576 = -t837 - t674 * t728 + (-t650 + t839) * qJD(1);
t575 = t628 * t753 + t629 * t752;
t574 = Icges(6,5) * t645 + (Icges(6,1) * t778 - Icges(6,4) * t772) * t646;
t573 = Icges(6,6) * t645 + (Icges(6,4) * t778 - Icges(6,2) * t772) * t646;
t572 = Icges(6,3) * t645 + (Icges(6,5) * t778 - Icges(6,6) * t772) * t646;
t571 = t650 * t727 + t651 * t728 + t648;
t570 = t737 * t775 + t781 * t813;
t569 = -t737 * t781 + t775 * t813;
t568 = t794 * t775 + (t781 * t849 + t682) * qJD(1) + t847;
t567 = t794 * t781 + (-t681 + t706 + (-pkin(15) - t849) * t775) * qJD(1);
t566 = -t738 * t755 + t649 * t727 + t751 + (t632 + t698) * qJD(1);
t565 = -t738 * t844 + t649 * t728 + (-t631 - t696 - t871) * qJD(1);
t564 = -t623 * t755 + (-pkin(14) * t781 + t615) * qJD(1);
t563 = -t623 * t844 + (pkin(14) * t775 - t614) * qJD(1);
t562 = (t614 * t775 + t615 * t781) * qJD(2);
t561 = qJD(1) * t606 - t619 * t727 + t751;
t560 = -t619 * t728 + (-t605 - t871) * qJD(1);
t559 = Icges(6,1) * t627 + Icges(6,4) * t626 + Icges(6,5) * t852;
t558 = Icges(6,1) * t625 + Icges(6,4) * t624 + Icges(6,5) * t853;
t557 = Icges(6,4) * t627 + Icges(6,2) * t626 + Icges(6,6) * t852;
t556 = Icges(6,4) * t625 + Icges(6,2) * t624 + Icges(6,6) * t853;
t555 = Icges(6,5) * t627 + Icges(6,6) * t626 + Icges(6,3) * t852;
t554 = Icges(6,5) * t625 + Icges(6,6) * t624 + Icges(6,3) * t853;
t553 = t605 * t727 + t606 * t728;
t552 = t631 * t727 + t632 * t728 + (t696 * t775 + t698 * t781) * qJD(2);
t551 = qJD(1) * (rSges(5,3) * t775 + t781 * t812) - (t752 * t668 + t753 * t801) * t703 + t787;
t550 = (-t601 + t836) * qJD(1) + t796;
t549 = t601 * t703 + t795;
t548 = t811 * t885 + t570 * t634 - t575 * t620 - (t672 * t753 + t673 * t752) * t703 + t787;
t547 = -t569 * t634 + t575 * t621 + (-t607 + t836) * qJD(1) + t796;
t546 = t569 * t620 - t570 * t621 + t607 * t703 + t795;
t1 = t634 * ((t554 * t621 + t555 * t620 + t572 * t634) * t645 + ((-t557 * t772 + t559 * t778) * t620 + (-t556 * t772 + t558 * t778) * t621 + (-t573 * t772 + t574 * t778) * t634) * t646) / 0.2e1 + m(9) * (t553 ^ 2 + t560 ^ 2 + t561 ^ 2) / 0.2e1 + m(7) * (t562 ^ 2 + t563 ^ 2 + t564 ^ 2) / 0.2e1 + m(10) * (t552 ^ 2 + t565 ^ 2 + t566 ^ 2) / 0.2e1 + m(4) * (t571 ^ 2 + t576 ^ 2 + t577 ^ 2) / 0.2e1 + m(11) * (t567 ^ 2 + t568 ^ 2 + t578 ^ 2) / 0.2e1 + ((t775 * t708 + t781 * t804) * qJD(1) + (t877 * t676 + (t810 * t781 + (-t675 + t809) * t775) * t781) * t763) * t746 / 0.2e1 + t620 * ((t555 * t852 + t626 * t557 + t627 * t559) * t620 + (t554 * t852 + t556 * t626 + t558 * t627) * t621 + (t572 * t852 + t573 * t626 + t574 * t627) * t634) / 0.2e1 + t621 * ((t555 * t853 + t557 * t624 + t559 * t625) * t620 + (t554 * t853 + t624 * t556 + t625 * t558) * t621 + (t572 * t853 + t573 * t624 + t574 * t625) * t634) / 0.2e1 + t703 * (((Icges(5,5) * t646 - Icges(5,6) * t645) * qJD(1) + (Icges(5,3) * t775 + (-Icges(5,5) * t645 - Icges(5,6) * t646) * t781) * t703) * t775 + ((-t582 * t646 - t584 * t645) * t703 + (-t586 * t646 - t587 * t645) * qJD(1)) * t781) / 0.2e1 + m(6) * (t546 ^ 2 + t547 ^ 2 + t548 ^ 2) / 0.2e1 + m(5) * (t549 ^ 2 + t550 ^ 2 + t551 ^ 2) / 0.2e1 + m(3) * (t630 ^ 2 + t635 ^ 2 + t636 ^ 2) / 0.2e1 + m(8) * (t612 ^ 2 + t613 ^ 2 + t648 ^ 2) / 0.2e1 + (((-t591 * t655 - t855) * t727 - (-t590 * t655 - t856) * t728 + (-t603 * t655 - t854) * qJD(1) + t785) * t781 + t884 * t775) * t727 / 0.2e1 + (Icges(8,1) * t719 ^ 2 + (0.2e1 * Icges(8,4) * t719 + Icges(8,2) * t715) * t715 + Icges(2,3) + m(2) * (t735 ^ 2 + t740 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t881 * t877 + (t878 * t781 + (t879 - t882) * t775) * t781) * qJD(2) + (t883 * t775 + t880 * t781) * qJD(1)) * t755 / 0.2e1 - ((t882 * t876 + (t879 * t775 + (t878 - t881) * t781) * t775) * qJD(2) + (t880 * t775 - t883 * t781) * qJD(1)) * t844 / 0.2e1 + (((t678 * t758 + t680 * t757) * t775 - (t677 * t758 + t679 * t757) * t781) * t763 + (-t582 * t645 + t584 * t646) * t703 - (-t590 * t653 + t592 * t655 + t640 * t721 + t642 * t723) * t728 + (-t591 * t653 + t593 * t655 + t641 * t721 + t643 * t723) * t727 + ((-t596 * t659 - t598 * t658 + t897 * t774 - t895 * t780) * t781 + (t597 * t659 + t599 * t658 + t896 * t774 + t894 * t780) * t775) * qJD(2) + (-t645 * t586 + t646 * t587 - t653 * t603 + t655 * t604 + t659 * t610 + t658 * t611 + t721 * t664 + t723 * t665 + t758 * t709 + t757 * t710 + t893 * t774 + t892 * t780) * qJD(1)) * qJD(1) / 0.2e1 - ((-t781 * t708 + t775 * t804) * qJD(1) + (t876 * t675 + (t809 * t775 + (-t676 + t810) * t781) * t775) * t763 - t884 * t781 + ((t591 * t652 - t855) * t727 - (t590 * t652 - t856) * t728 + (t603 * t652 - t854) * qJD(1) + t785) * t775) * t728 / 0.2e1;
T = t1;
