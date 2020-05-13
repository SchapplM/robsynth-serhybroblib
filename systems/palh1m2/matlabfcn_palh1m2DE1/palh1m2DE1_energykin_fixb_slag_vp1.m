% Calculate kinetic energy for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE1_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:56:03
% EndTime: 2020-05-01 20:56:12
% DurationCPUTime: 9.63s
% Computational Cost: add. (2198->376), mult. (4014->630), div. (0->0), fcn. (4725->24), ass. (0->223)
t896 = -Icges(3,4) - Icges(10,4);
t895 = Icges(3,1) + Icges(10,1);
t894 = Icges(3,2) + Icges(10,2);
t769 = cos(qJ(2));
t893 = t896 * t769;
t763 = sin(qJ(2));
t892 = t896 * t763;
t891 = Icges(3,5) + Icges(10,5);
t890 = -Icges(3,6) - Icges(10,6);
t889 = -t894 * t769 + t892;
t888 = -t895 * t763 + t893;
t757 = sin(pkin(19));
t760 = cos(pkin(19));
t762 = sin(qJ(3));
t768 = cos(qJ(3));
t703 = t757 * t768 + t760 * t762;
t707 = -t757 * t762 + t760 * t768;
t655 = -t703 * t769 - t707 * t763;
t704 = t757 * t763 - t760 * t769;
t708 = t757 * t769 + t760 * t763;
t656 = t704 * t768 + t708 * t762;
t658 = -t704 * t762 + t708 * t768;
t837 = t768 * t769;
t711 = -t762 * t763 + t837;
t838 = t763 * t768;
t713 = t762 * t769 + t838;
t752 = qJD(2) + qJD(3);
t764 = sin(qJ(1));
t718 = t752 * t764;
t770 = cos(qJ(1));
t719 = t752 * t770;
t753 = qJ(3) + qJ(2);
t747 = sin(t753);
t748 = cos(t753);
t876 = Icges(11,5) * t747 + Icges(11,6) * t748;
t879 = -Icges(4,5) * t711 + Icges(9,5) * t656 + Icges(4,6) * t713;
t886 = Icges(9,3) + Icges(4,3);
t887 = (t886 * t770 + (-Icges(9,6) * t655 + t879) * t764) * t719 + ((-Icges(9,6) * t658 - t879) * t770 + t886 * t764) * t718 + (Icges(4,5) * t713 + Icges(9,5) * t658 + Icges(4,6) * t711 - Icges(9,6) * t656 + t876) * qJD(1);
t863 = t764 ^ 2;
t862 = t770 ^ 2;
t849 = Icges(9,4) * t656;
t885 = t889 * t764 + t890 * t770;
t884 = t890 * t764 - t889 * t770;
t883 = t888 * t764 - t891 * t770;
t882 = t891 * t764 + t888 * t770;
t881 = t894 * t763 + t893;
t880 = t895 * t769 + t892;
t878 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t765 = sin(pkin(18));
t766 = sin(pkin(17));
t771 = cos(pkin(18));
t772 = cos(pkin(17));
t712 = t765 * t772 - t766 * t771;
t714 = t765 * t766 + t771 * t772;
t661 = t712 * t763 + t714 * t769;
t662 = t712 * t769 - t714 * t763;
t877 = Icges(7,5) * t662 - Icges(7,6) * t661 - t891 * t763 + t890 * t769;
t873 = Icges(7,4) * t661;
t872 = t764 * t770;
t869 = t877 * t764 - t878 * t770;
t868 = Icges(7,5) * t661 + Icges(7,6) * t662 + t890 * t763 + t891 * t769;
t867 = t878 * t764 + t877 * t770;
t621 = Icges(7,2) * t662 + t873;
t850 = Icges(7,4) * t662;
t622 = Icges(7,1) * t661 + t850;
t866 = -t621 * t661 + t622 * t662 - t880 * t763 + t881 * t769;
t811 = -Icges(7,2) * t661 + t850;
t610 = -Icges(7,6) * t770 + t764 * t811;
t815 = Icges(7,1) * t662 - t873;
t612 = -Icges(7,5) * t770 + t764 * t815;
t865 = t610 * t661 - t612 * t662 + t883 * t763 + t885 * t769;
t611 = Icges(7,6) * t764 + t770 * t811;
t613 = Icges(7,5) * t764 + t770 * t815;
t864 = -t611 * t661 + t613 * t662 - t882 * t763 + t884 * t769;
t857 = pkin(1) * qJD(2);
t856 = pkin(5) * qJD(3);
t741 = pkin(1) * t763 - pkin(15);
t855 = Icges(9,1) * t656;
t852 = Icges(4,4) * t711;
t851 = Icges(4,4) * t713;
t848 = Icges(9,4) * t658;
t606 = Icges(9,5) * t764 + (-t848 - t855) * t770;
t844 = t606 * t656;
t607 = -Icges(9,5) * t770 + (Icges(9,4) * t655 - t855) * t764;
t843 = t607 * t656;
t617 = Icges(9,1) * t658 - t849;
t842 = t617 * t656;
t756 = sin(pkin(20));
t759 = cos(pkin(20));
t706 = -t771 * t756 + t759 * t765;
t710 = t756 * t765 + t759 * t771;
t751 = pkin(22) + pkin(21);
t744 = sin(t751);
t745 = cos(t751);
t651 = -t706 * t745 + t710 * t744;
t841 = t651 * t764;
t840 = t651 * t770;
t742 = pkin(5) * t762 + pkin(1);
t784 = pkin(5) * t837 - t742 * t763;
t652 = t784 * qJD(2) + t711 * t856;
t834 = qJD(1) * t764;
t833 = qJD(1) * t770;
t832 = qJD(2) * t764;
t831 = qJD(2) * t770;
t830 = qJD(4) * t651;
t829 = t763 * t857;
t828 = t769 * t857;
t761 = sin(qJ(4));
t767 = cos(qJ(4));
t820 = rSges(6,1) * t767 - rSges(6,2) * t761;
t637 = rSges(6,3) * t710 - t706 * t820;
t638 = rSges(6,3) * t706 + t710 * t820;
t827 = (t637 * t745 + t638 * t744) * t830;
t716 = -rSges(9,1) * t760 + rSges(9,2) * t757;
t717 = rSges(9,1) * t757 + rSges(9,2) * t760;
t822 = -t716 * t762 + t717 * t768;
t821 = t764 * t828;
t730 = rSges(3,1) * t763 + rSges(3,2) * t769;
t819 = -rSges(10,1) * t763 - rSges(10,2) * t769;
t818 = rSges(11,1) * t748 - rSges(11,2) * t747;
t816 = Icges(4,1) * t711 - t851;
t812 = -Icges(4,2) * t713 + t852;
t650 = t706 * t744 + t710 * t745;
t781 = t761 * t650;
t633 = t764 * t781 - t767 * t770;
t779 = t650 * t767;
t634 = -t761 * t770 - t764 * t779;
t635 = t764 * t767 + t770 * t781;
t636 = t761 * t764 - t770 * t779;
t805 = (Icges(6,5) * t634 + Icges(6,6) * t633 + Icges(6,3) * t841) * t764 + (Icges(6,5) * t636 + Icges(6,6) * t635 + Icges(6,3) * t840) * t770;
t801 = t637 * t744 - t638 * t745;
t720 = rSges(5,1) * t756 - rSges(5,2) * t759;
t721 = rSges(5,1) * t759 + rSges(5,2) * t756;
t800 = -(t720 * t771 - t721 * t765) * t744 + (t720 * t765 + t721 * t771) * t745;
t737 = -t765 * rSges(7,1) + rSges(7,2) * t771;
t738 = rSges(7,1) * t771 + rSges(7,2) * t765;
t789 = t737 * t766 - t738 * t772;
t754 = sin(pkin(22));
t758 = cos(pkin(22));
t709 = t754 * t771 - t758 * t765;
t786 = Icges(11,5) * t748 - Icges(11,6) * t747;
t785 = t741 * t834 - t770 * t828;
t783 = -(rSges(11,1) * t747 + rSges(11,2) * t748) * t752 - t828;
t782 = -qJD(2) * (pkin(5) * t838 + t742 * t769) - t713 * t856;
t780 = t805 * t651;
t778 = t782 * t764;
t645 = -Icges(4,6) * t770 + t764 * t812;
t646 = Icges(4,6) * t764 + t770 * t812;
t647 = -Icges(4,5) * t770 + t764 * t816;
t648 = Icges(4,5) * t764 + t770 * t816;
t664 = Icges(4,2) * t711 + t851;
t665 = Icges(4,1) * t713 + t852;
t775 = (-t646 * t713 + t648 * t711) * t718 - (-t645 * t713 + t647 * t711) * t719 + (-t664 * t713 + t665 * t711) * qJD(1);
t773 = qJD(2) ^ 2;
t743 = pkin(15) * t833;
t736 = rSges(2,1) * t770 - rSges(2,2) * t764;
t735 = rSges(3,1) * t769 - rSges(3,2) * t763;
t734 = rSges(10,1) * t769 - rSges(10,2) * t763;
t732 = rSges(6,1) * t761 + rSges(6,2) * t767;
t731 = rSges(2,1) * t764 + rSges(2,2) * t770;
t729 = pkin(9) * t759 - pkin(11) * t756;
t728 = pkin(9) * t756 + pkin(11) * t759;
t705 = -t754 * t765 - t758 * t771;
t694 = -pkin(15) - t784;
t680 = t694 * t834;
t679 = rSges(11,3) * t764 + t770 * t818;
t678 = -rSges(11,3) * t770 + t764 * t818;
t673 = Icges(11,3) * t764 + t770 * t786;
t672 = -Icges(11,3) * t770 + t764 * t786;
t671 = (rSges(4,1) * t763 + rSges(4,2) * t769) * t768 + t762 * (rSges(4,1) * t769 - rSges(4,2) * t763);
t670 = -t737 * t772 - t738 * t766;
t669 = (rSges(4,1) * t768 - rSges(4,2) * t762) * t769 - t763 * (rSges(4,1) * t762 + rSges(4,2) * t768);
t666 = -t716 * t768 - t717 * t762;
t654 = rSges(4,3) * t764 + t669 * t770;
t653 = -rSges(4,3) * t770 + t669 * t764;
t642 = -t735 * t831 + (t770 * rSges(3,3) + (-pkin(15) + t730) * t764) * qJD(1);
t641 = t743 - t735 * t832 + qJD(1) * (rSges(3,3) * t764 - t730 * t770);
t640 = qJD(4) * t650 + qJD(1);
t639 = (-cos(pkin(21)) * t705 - t709 * sin(pkin(21))) * pkin(4) + t741;
t632 = -t821 + (t764 * rSges(8,3) + (-t741 + (-rSges(8,1) * t771 + rSges(8,2) * t765) * t758 - t754 * (rSges(8,1) * t765 + rSges(8,2) * t771)) * t770) * qJD(1);
t631 = qJD(1) * (t770 * rSges(8,3) + t764 * ((rSges(8,1) * t758 + rSges(8,2) * t754) * t771 + t765 * (rSges(8,1) * t754 - rSges(8,2) * t758))) + t785;
t630 = t763 * t670 - t769 * t789;
t629 = t670 * t769 + t763 * t789;
t628 = t763 * t666 + t769 * t822;
t627 = t666 * t769 - t763 * t822;
t626 = (t728 * t765 + t729 * t771) * t745 - t744 * (t728 * t771 - t729 * t765);
t625 = -t829 + (t678 * t764 + t679 * t770) * t752;
t624 = rSges(7,3) * t764 + t629 * t770;
t623 = -rSges(7,3) * t770 + t629 * t764;
t619 = rSges(9,3) * t764 + t627 * t770;
t618 = -rSges(9,3) * t770 + t627 * t764;
t616 = -Icges(9,2) * t656 + t848;
t614 = t819 * qJD(2) + t752 * pkin(2) * (-t703 * t763 + t707 * t769);
t605 = -Icges(9,6) * t770 + (Icges(9,2) * t655 - t849) * t764;
t604 = Icges(9,6) * t764 + (-Icges(9,2) * t658 - t849) * t770;
t601 = -t821 - t671 * t718 + (-t770 * t741 + t654) * qJD(1);
t600 = -qJD(1) * t653 - t671 * t719 + t785;
t599 = t653 * t718 + t654 * t719 - t829;
t598 = t783 * t770 + (t639 * t764 - t678) * qJD(1);
t597 = t783 * t764 + (-t639 * t770 + t679) * qJD(1);
t595 = -pkin(2) * t719 * t658 - t734 * t831 + (t770 * rSges(10,3) + (pkin(2) * t656 - pkin(15) - t819) * t764) * qJD(1);
t594 = t743 - t734 * t832 + qJD(1) * (rSges(10,3) * t764 + t770 * t819) + (-t656 * t833 - t658 * t718) * pkin(2);
t593 = Icges(6,5) * t650 + (Icges(6,1) * t767 - Icges(6,4) * t761) * t651;
t592 = Icges(6,6) * t650 + (Icges(6,4) * t767 - Icges(6,2) * t761) * t651;
t591 = Icges(6,3) * t650 + (Icges(6,5) * t767 - Icges(6,6) * t761) * t651;
t590 = -t732 * t770 + t764 * t801;
t589 = t732 * t764 + t770 * t801;
t588 = -t630 * t832 + (-pkin(14) * t770 + t624) * qJD(1);
t587 = -t630 * t831 + (pkin(14) * t764 - t623) * qJD(1);
t586 = t778 + (t764 * rSges(5,3) + (-t694 - t800) * t770) * qJD(1);
t585 = t680 + t800 * t834 + (rSges(5,3) * qJD(1) + t782) * t770;
t584 = (t623 * t764 + t624 * t770) * qJD(2);
t583 = qJD(1) * t619 - t628 * t718 + t743;
t582 = -t628 * t719 + (-t764 * pkin(15) - t618) * qJD(1);
t581 = Icges(6,1) * t636 + Icges(6,4) * t635 + Icges(6,5) * t840;
t580 = Icges(6,1) * t634 + Icges(6,4) * t633 + Icges(6,5) * t841;
t579 = Icges(6,4) * t636 + Icges(6,2) * t635 + Icges(6,6) * t840;
t578 = Icges(6,4) * t634 + Icges(6,2) * t633 + Icges(6,6) * t841;
t575 = t618 * t718 + t619 * t719;
t574 = (-t589 * t764 + t590 * t770) * t830 + t652;
t573 = t589 * t640 + t778 + (-t827 + (-t626 - t694) * qJD(1)) * t770;
t572 = -t590 * t640 + t680 + t782 * t770 + (qJD(1) * t626 + t827) * t764;
t1 = m(9) * (t575 ^ 2 + t582 ^ 2 + t583 ^ 2) / 0.2e1 + m(7) * (t584 ^ 2 + t587 ^ 2 + t588 ^ 2) / 0.2e1 + m(4) * (t599 ^ 2 + t600 ^ 2 + t601 ^ 2) / 0.2e1 + m(6) * (t572 ^ 2 + t573 ^ 2 + t574 ^ 2) / 0.2e1 + m(5) * (t585 ^ 2 + t586 ^ 2 + t652 ^ 2) / 0.2e1 + m(3) * (t730 ^ 2 * t773 + t641 ^ 2 + t642 ^ 2) / 0.2e1 + m(8) * (pkin(1) ^ 2 * t763 ^ 2 * t773 + t631 ^ 2 + t632 ^ 2) / 0.2e1 + m(11) * (t597 ^ 2 + t598 ^ 2 + t625 ^ 2) / 0.2e1 + t640 * ((t650 * t591 + (-t592 * t761 + t593 * t767) * t651) * t640 + (((-t579 * t761 + t581 * t767) * t770 + (-t578 * t761 + t580 * t767) * t764) * t651 + t805 * t650) * t830) / 0.2e1 + m(10) * (t594 ^ 2 + t595 ^ 2 + t614 ^ 2) / 0.2e1 + (t770 * ((t591 * t840 + t592 * t635 + t593 * t636) * t640 + ((t578 * t635 + t580 * t636) * t764 + (t635 * t579 + t636 * t581 + t780) * t770) * t830) + t764 * ((t591 * t841 + t592 * t633 + t593 * t634) * t640 + ((t579 * t633 + t581 * t634) * t770 + (t633 * t578 + t634 * t580 + t780) * t764) * t830)) * t830 / 0.2e1 + ((t867 * t863 + (t865 * t770 + (t864 - t869) * t764) * t770) * qJD(2) + (t868 * t764 + t866 * t770) * qJD(1)) * t832 / 0.2e1 - ((t869 * t862 + (t864 * t764 + (t865 - t867) * t770) * t764) * qJD(2) + (t866 * t764 - t868 * t770) * qJD(1)) * t831 / 0.2e1 + (m(2) * (t731 ^ 2 + t736 ^ 2) + Icges(8,1) * t709 ^ 2 + (0.2e1 * Icges(8,4) * t709 + Icges(8,2) * t705) * t705 + Icges(5,1) * t651 ^ 2 - (0.2e1 * Icges(5,4) * t651 - Icges(5,2) * t650) * t650 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t605 * t656 - t607 * t658 - t645 * t711 - t647 * t713) * t719 + (-t604 * t656 + t606 * t658 + t646 * t711 + t648 * t713) * t718 + ((-t610 * t662 - t612 * t661 + t763 * t885 - t769 * t883) * t770 + (t611 * t662 + t613 * t661 + t763 * t884 + t769 * t882) * t764) * qJD(2) + (Icges(11,2) * t748 ^ 2 - t656 * t616 + t658 * t617 + t662 * t621 + t661 * t622 + t711 * t664 + t713 * t665 + t881 * t763 + t880 * t769 + (Icges(11,1) * t747 + 0.2e1 * Icges(11,4) * t748) * t747) * qJD(1) + (t862 + t863) * t752 * t876) * qJD(1) / 0.2e1 + ((-t672 * t872 + t863 * t673) * t752 + ((-t604 * t658 - t844) * t718 - (-t605 * t658 - t843) * t719 + (-t616 * t658 - t842) * qJD(1) + t775) * t770 + t887 * t764) * t718 / 0.2e1 - ((t862 * t672 - t673 * t872) * t752 + ((t604 * t655 - t844) * t718 - (t605 * t655 - t843) * t719 + (t616 * t655 - t842) * qJD(1) + t775) * t764 - t887 * t770) * t719 / 0.2e1;
T = t1;
