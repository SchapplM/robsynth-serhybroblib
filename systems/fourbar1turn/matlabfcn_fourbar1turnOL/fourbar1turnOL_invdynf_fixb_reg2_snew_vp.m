% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% fourbar1turnOL
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = fourbar1turnOL_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:18
% EndTime: 2020-04-12 19:41:19
% DurationCPUTime: 0.80s
% Computational Cost: add. (1064->129), mult. (2473->180), div. (0->0), fcn. (1807->8), ass. (0->112)
t684 = sin(qJ(3));
t685 = sin(qJ(2));
t688 = cos(qJ(3));
t689 = cos(qJ(2));
t640 = (t684 * t685 - t688 * t689) * qJD(1);
t712 = t640 ^ 2;
t641 = (-t684 * t689 - t685 * t688) * qJD(1);
t711 = t641 ^ 2;
t678 = qJD(2) + qJD(3);
t710 = t678 ^ 2;
t709 = t640 * t641;
t708 = t678 * t640;
t707 = t678 * t641;
t686 = sin(qJ(1));
t690 = cos(qJ(1));
t664 = g(1) * t686 - t690 * g(2);
t706 = t686 * t664;
t683 = sin(qJ(4));
t679 = t683 ^ 2;
t687 = cos(qJ(4));
t681 = t687 ^ 2;
t705 = t679 + t681;
t680 = t685 ^ 2;
t682 = t689 ^ 2;
t704 = t680 + t682;
t703 = qJD(1) * qJD(2);
t702 = qJD(1) * qJD(4);
t701 = t683 * qJDD(1);
t700 = t685 * qJDD(1);
t699 = t687 * qJDD(1);
t693 = qJD(1) ^ 2;
t672 = t689 * t693 * t685;
t662 = qJDD(2) + t672;
t698 = -qJDD(2) - qJDD(3);
t697 = t685 * t703;
t696 = t689 * t703;
t665 = -g(1) * t690 - g(2) * t686;
t639 = -t685 * g(3) + t689 * t665;
t692 = qJD(2) ^ 2;
t669 = -t682 * t693 - t692;
t650 = t696 + t700;
t676 = t689 * qJDD(1);
t652 = t676 - t697;
t695 = -t641 * qJD(3) + t650 * t684 - t688 * t652;
t638 = -t689 * g(3) - t685 * t665;
t694 = -t640 * qJD(3) + t650 * t688 + t652 * t684;
t691 = qJD(4) ^ 2;
t670 = t687 * t693 * t683;
t668 = -t681 * t693 - t691;
t667 = -t680 * t693 - t692;
t666 = -t679 * t693 - t691;
t663 = -qJDD(2) + t672;
t661 = -qJDD(4) + t670;
t660 = qJDD(4) + t670;
t659 = t704 * t693;
t658 = t705 * t693;
t657 = -qJDD(1) * t686 - t690 * t693;
t656 = qJDD(1) * t690 - t686 * t693;
t655 = t704 * qJDD(1);
t654 = t705 * qJDD(1);
t653 = t676 - 0.2e1 * t697;
t651 = -0.2e1 * t683 * t702 + t699;
t649 = 0.2e1 * t696 + t700;
t648 = 0.2e1 * t687 * t702 + t701;
t647 = -pkin(1) * t693 + t665;
t646 = qJDD(1) * pkin(1) + t664;
t645 = t690 * t664;
t634 = -g(3) * t683 + t647 * t687;
t633 = -g(3) * t687 - t647 * t683;
t632 = -t710 - t711;
t631 = t663 * t689 - t667 * t685;
t630 = -t662 * t685 + t669 * t689;
t629 = t661 * t687 - t666 * t683;
t628 = -t660 * t683 + t668 * t687;
t627 = t663 * t685 + t667 * t689;
t626 = t662 * t689 + t669 * t685;
t625 = t661 * t683 + t666 * t687;
t624 = t660 * t687 + t668 * t683;
t623 = pkin(2) * t669 + t639;
t622 = pkin(2) * t662 + t638;
t621 = t664 + (t652 - t697) * pkin(2);
t620 = t698 + t709;
t619 = -t698 + t709;
t618 = -t710 - t712;
t617 = -t638 * t685 + t639 * t689;
t616 = t638 * t689 + t639 * t685;
t615 = -t633 * t683 + t634 * t687;
t614 = t633 * t687 + t634 * t683;
t613 = -t711 - t712;
t612 = -t622 * t684 - t623 * t688;
t611 = -t622 * t688 + t623 * t684;
t610 = -t620 * t688 + t632 * t684;
t609 = -t620 * t684 - t632 * t688;
t608 = t694 + t708;
t607 = t694 - t708;
t606 = t695 - t707;
t605 = t695 + t707;
t604 = -t618 * t688 + t619 * t684;
t603 = -t618 * t684 - t619 * t688;
t602 = t611 * t684 - t612 * t688;
t601 = -t611 * t688 - t612 * t684;
t600 = -t609 * t685 + t610 * t689;
t599 = t609 * t689 + t610 * t685;
t598 = -t605 * t688 + t608 * t684;
t597 = -t605 * t684 - t608 * t688;
t596 = -t603 * t685 + t604 * t689;
t595 = t603 * t689 + t604 * t685;
t594 = -t601 * t685 + t602 * t689;
t593 = t601 * t689 + t602 * t685;
t592 = -t597 * t685 + t598 * t689;
t591 = t597 * t689 + t598 * t685;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t657, -t656, 0, t665 * t690 - t706, 0, 0, 0, 0, 0, 0, t630 * t690 - t653 * t686, t631 * t690 + t649 * t686, t655 * t690 - t659 * t686, t617 * t690 - t706, 0, 0, 0, 0, 0, 0, t596 * t690 - t606 * t686, t600 * t690 - t607 * t686, t592 * t690 + t613 * t686, t594 * t690 - t621 * t686, 0, 0, 0, 0, 0, 0, t628 * t690 - t651 * t686, t629 * t690 + t648 * t686, t654 * t690 - t658 * t686, t615 * t690 - t646 * t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t656, t657, 0, t665 * t686 + t645, 0, 0, 0, 0, 0, 0, t630 * t686 + t653 * t690, t631 * t686 - t649 * t690, t655 * t686 + t659 * t690, t617 * t686 + t645, 0, 0, 0, 0, 0, 0, t596 * t686 + t606 * t690, t600 * t686 + t607 * t690, t592 * t686 - t613 * t690, t594 * t686 + t621 * t690, 0, 0, 0, 0, 0, 0, t628 * t686 + t651 * t690, t629 * t686 - t648 * t690, t654 * t686 + t658 * t690, t615 * t686 + t646 * t690; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t626, t627, 0, t616, 0, 0, 0, 0, 0, 0, t595, t599, t591, t593, 0, 0, 0, 0, 0, 0, t624, t625, 0, t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t693, -qJDD(1), 0, t665, 0, 0, 0, 0, 0, 0, t630, t631, t655, t617, 0, 0, 0, 0, 0, 0, t596, t600, t592, t594, 0, 0, 0, 0, 0, 0, t628, t629, t654, t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t693, 0, t664, 0, 0, 0, 0, 0, 0, t653, -t649, t659, t664, 0, 0, 0, 0, 0, 0, t606, t607, -t613, t621, 0, 0, 0, 0, 0, 0, t651, -t648, t658, t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t626, t627, 0, t616, 0, 0, 0, 0, 0, 0, t595, t599, t591, t593, 0, 0, 0, 0, 0, 0, t624, t625, 0, t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t669, t663, t676, t639, 0, 0, 0, 0, 0, 0, t604, t610, t598, t602, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t662, t667, -t700, t638, 0, 0, 0, 0, 0, 0, t603, t609, t597, t601, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t653, t649, -t659, -t664, 0, 0, 0, 0, 0, 0, -t606, -t607, t613, -t621, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t618, t620, t605, t612, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t619, t632, t608, t611, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t606, -t607, t613, -t621, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t668, t661, t699, t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t660, t666, -t701, t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t651, t648, -t658, -t646;];
f_new_reg = t1;
