% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m2DE2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh3m2DE2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m2DE2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobiaD_transl_sym_varpar: pkin has to be [18x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(12);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->25), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
	t46 = qJD(2) + qJD(3);
	t47 = qJ(2) + qJ(3);
	t45 = cos(t47);
	t65 = r_i_i_C(2) * t45;
	t44 = sin(t47);
	t66 = r_i_i_C(1) * t44;
	t54 = t65 + t66;
	t48 = sin(qJ(2));
	t67 = pkin(1) * t48;
	t55 = qJD(2) * t67;
	t68 = t54 * t46 - t55;
	t64 = t44 * t46;
	t63 = t45 * t46;
	t62 = r_i_i_C(1) * t64 + r_i_i_C(2) * t63;
	t49 = sin(qJ(1));
	t61 = qJD(1) * t49;
	t51 = cos(qJ(1));
	t60 = qJD(1) * t51;
	t50 = cos(qJ(2));
	t59 = qJD(2) * t50;
	t58 = r_i_i_C(1) * t63;
	t57 = r_i_i_C(2) * t64;
	t56 = qJD(1) * t66;
	t53 = -t50 * pkin(1) + r_i_i_C(1) * t45 - r_i_i_C(2) * t44 - pkin(12);
	t52 = t51 * t56 + t60 * t65 + (-t57 + t58) * t49;
	t38 = t51 * t58;
	t1 = [-t68 * t49 + (-r_i_i_C(3) * t49 + t53 * t51) * qJD(1), t38 + (-pkin(1) * t59 - t57) * t51 + (-t54 + t67) * t61, -t49 * t56 + t38 + (-t45 * t61 - t51 * t64) * r_i_i_C(2), 0; t68 * t51 + (r_i_i_C(3) * t51 + t53 * t49) * qJD(1), (-t48 * t60 - t49 * t59) * pkin(1) + t52, t52, 0; 0, -t55 + t62, t62, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:25
	% EndTime: 2020-05-07 02:13:27
	% DurationCPUTime: 2.06s
	% Computational Cost: add. (18993->73), mult. (35262->124), div. (516->4), fcn. (50293->17), ass. (0->77)
	t175 = pkin(17) + pkin(18);
	t170 = sin(t175);
	t171 = cos(t175);
	t177 = sin(pkin(16));
	t178 = cos(pkin(16));
	t182 = sin(pkin(15));
	t186 = cos(pkin(15));
	t163 = t177 * t186 + t178 * t182;
	t164 = -t177 * t182 + t178 * t186;
	t179 = sin(qJ(3));
	t183 = cos(qJ(3));
	t157 = t179 * t163 - t164 * t183;
	t180 = sin(qJ(2));
	t184 = cos(qJ(2));
	t197 = t163 * t183 + t179 * t164;
	t247 = t157 * t180 - t197 * t184;
	t248 = t157 * t184 + t180 * t197;
	t250 = t248 * t170 + t171 * t247;
	t135 = t250 ^ 2;
	t255 = -t247 * t170 + t171 * t248;
	t137 = 0.1e1 / t255 ^ 2;
	t134 = t135 * t137 + 0.1e1;
	t132 = 0.1e1 / t134;
	t136 = 0.1e1 / t255;
	t174 = qJD(2) + qJD(3);
	t224 = t137 * t250;
	t155 = t197 * qJD(3);
	t246 = t157 * qJD(3);
	t190 = qJD(2) * t248 + t180 * t155 + t184 * t246;
	t254 = -qJD(2) * t247 + t155 * t184 - t180 * t246;
	t251 = t170 * t254 + t171 * t190;
	t259 = -t190 * t170 + t254 * t171;
	t123 = (-t136 * t251 + t224 * t259) * t132 + t174;
	t193 = t136 * t255 + t250 * t224;
	t262 = -t193 * t132 + 0.1e1;
	t265 = t123 * t262;
	t176 = qJ(2) + qJ(3);
	t129 = atan2(t250, -t255) + t176;
	t128 = cos(t129);
	t264 = t128 * t265;
	t127 = sin(t129);
	t263 = t127 * t265;
	t245 = t259 * t137;
	t227 = t136 * t245;
	t256 = t251 * t224;
	t260 = 0.2e1 * t193 * (-t135 * t227 + t256) / t134 ^ 2;
	t258 = -0.2e1 * t250;
	t257 = t259 * t136;
	t172 = sin(t176);
	t234 = pkin(4) * t174;
	t168 = t172 * t234;
	t233 = pkin(1) * qJD(2);
	t192 = -t180 * t233 + t168;
	t200 = r_i_i_C(1) * t127 + r_i_i_C(2) * t128;
	t240 = t200 * t123 + t192;
	t239 = -pkin(1) * t180 + t200 * t262;
	t235 = pkin(4) * t172;
	t181 = sin(qJ(1));
	t230 = t123 * t181;
	t185 = cos(qJ(1));
	t229 = t123 * t185;
	t205 = qJD(1) * t235;
	t173 = cos(t176);
	t206 = t173 * t234;
	t210 = t181 * t206 + t185 * t205;
	t209 = qJD(1) * t181;
	t208 = qJD(1) * t185;
	t202 = t227 * t258;
	t122 = ((-t137 * t251 - t202) * t250 - t257 + t255 * t245 - t256) * t132 + t260;
	t196 = t122 * t200;
	t121 = t260 + (-t257 - t250 * t202 + (t251 * t258 + t255 * t259) * t137) * t132;
	t195 = t121 * t127 + t264;
	t194 = -t121 * t128 + t263;
	t191 = -t184 * pkin(1) + pkin(4) * t173 + r_i_i_C(1) * t128 - r_i_i_C(2) * t127 - pkin(12);
	t189 = t195 * r_i_i_C(1) - t194 * r_i_i_C(2) - t184 * t233;
	t166 = t185 * t206;
	t1 = [-t240 * t181 + (-r_i_i_C(3) * t181 + t191 * t185) * qJD(1), t166 + (-t235 - t239) * t209 + t189 * t185, -t181 * t205 + t166 + t185 * t196 + ((-t127 * t229 - t128 * t209) * r_i_i_C(2) + (-t127 * t209 + t128 * t229) * r_i_i_C(1)) * t262, 0; t240 * t185 + (r_i_i_C(3) * t185 + t191 * t181) * qJD(1), t189 * t181 + t239 * t208 + t210, t181 * t196 + ((-t127 * t230 + t128 * t208) * r_i_i_C(2) + (t127 * t208 + t128 * t230) * r_i_i_C(1)) * t262 + t210, 0; 0, t194 * r_i_i_C(1) + t195 * r_i_i_C(2) + t192, t168 + (t122 * t127 + t264) * r_i_i_C(2) + (-t122 * t128 + t263) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:44
	% EndTime: 2020-05-07 02:13:48
	% DurationCPUTime: 3.79s
	% Computational Cost: add. (37686->129), mult. (61626->201), div. (1422->6), fcn. (87433->35), ass. (0->126)
	t787 = qJ(2) + qJ(3);
	t716 = pkin(17) + pkin(18);
	t705 = sin(t716);
	t706 = cos(t716);
	t719 = sin(pkin(16));
	t720 = cos(pkin(16));
	t725 = sin(pkin(15));
	t730 = cos(pkin(15));
	t690 = t719 * t730 + t720 * t725;
	t691 = -t719 * t725 + t720 * t730;
	t722 = sin(qJ(3));
	t727 = cos(qJ(3));
	t667 = t722 * t690 - t691 * t727;
	t723 = sin(qJ(2));
	t728 = cos(qJ(2));
	t753 = t690 * t727 + t722 * t691;
	t831 = t667 * t723 - t753 * t728;
	t833 = t667 * t728 + t723 * t753;
	t835 = t833 * t705 + t706 * t831;
	t840 = -t831 * t705 + t706 * t833;
	t616 = atan2(t835, -t840) + t787;
	t614 = sin(t616);
	t615 = cos(t616);
	t711 = qJ(1) + t787;
	t713 = qJD(2) + qJD(3);
	t778 = pkin(4) * (qJD(1) + t713) / 0.2e1;
	t686 = sin(t711) * t778;
	t712 = qJ(1) - t787;
	t703 = sin(t712);
	t724 = sin(qJ(1));
	t640 = t835 ^ 2;
	t642 = 0.1e1 / t840 ^ 2;
	t621 = t640 * t642 + 0.1e1;
	t619 = 0.1e1 / t621;
	t641 = 0.1e1 / t840;
	t807 = t642 * t835;
	t665 = t753 * qJD(3);
	t830 = t667 * qJD(3);
	t740 = qJD(2) * t833 + t723 * t665 + t728 * t830;
	t839 = -qJD(2) * t831 + t665 * t728 - t723 * t830;
	t836 = t705 * t839 + t706 * t740;
	t844 = -t740 * t705 + t839 * t706;
	t610 = (-t641 * t836 + t807 * t844) * t619 + t713;
	t721 = sin(qJ(4));
	t726 = cos(qJ(4));
	t756 = r_i_i_C(1) * t721 + r_i_i_C(2) * t726;
	t745 = r_i_i_C(3) * t610 - t756 * qJD(4);
	t729 = cos(qJ(1));
	t757 = r_i_i_C(1) * t726 - r_i_i_C(2) * t721;
	t749 = t757 * t729;
	t702 = pkin(16) + pkin(15) + t716 + t787;
	t699 = sin(t702);
	t700 = cos(t702);
	t793 = t699 ^ 2 / t700 ^ 2;
	t682 = 0.1e1 + t793;
	t679 = 0.1e1 / t682;
	t659 = -t679 * t793 - t679 + 0.1e1;
	t655 = (-t679 * t682 + 0.1e1) * t713;
	t654 = -qJD(1) + t655;
	t678 = atan2(-t699, t700) + t787;
	t677 = -qJ(1) + t678;
	t802 = t654 * cos(t677);
	t758 = t659 * t802 / 0.2e1;
	t803 = t654 * sin(t677);
	t766 = -t803 / 0.2e1;
	t653 = qJD(1) + t655;
	t676 = qJ(1) + t678;
	t804 = t653 * cos(t676);
	t768 = -t804 / 0.2e1;
	t818 = pkin(4) * (qJD(1) - t713);
	t777 = -t818 / 0.2e1;
	t785 = qJD(1) * t729;
	t814 = t610 * t724;
	t815 = pkin(10) * t659;
	t816 = pkin(8) * t659;
	t769 = t653 * sin(t676) / 0.2e1;
	t832 = t659 * t769;
	t748 = t641 * t840 + t835 * t807;
	t847 = -t748 * t619 + 0.1e1;
	t850 = ((-r_i_i_C(3) * t785 + t757 * t814) * t615 + (qJD(1) * t749 + t745 * t724) * t614) * t847 + t832 * pkin(8) + pkin(10) * t758 + t703 * t777 + t766 * t816 + t768 * t815 + t686;
	t687 = cos(t711) * t778;
	t704 = cos(t712);
	t765 = t803 / 0.2e1;
	t767 = t804 / 0.2e1;
	t786 = qJD(1) * t724;
	t849 = ((r_i_i_C(3) * t786 + t610 * t749) * t615 + (t745 * t729 - t757 * t786) * t614) * t847 + pkin(8) * t758 + t832 * pkin(10) + t704 * t777 + t765 * t815 + t767 * t816 + t687;
	t746 = -r_i_i_C(3) * t615 + t757 * t614;
	t783 = qJD(4) * t615;
	t848 = (t746 * t610 + t756 * t783) * t847 + pkin(4) * t713 * sin(t787) + (sin(t678) * pkin(8) - cos(t678) * pkin(10)) * t655 * t659;
	t829 = t844 * t642;
	t810 = t641 * t829;
	t841 = t836 * t807;
	t845 = 0.2e1 * t748 * (-t640 * t810 + t841) / t621 ^ 2;
	t843 = -0.2e1 * t835;
	t842 = t844 * t641;
	t760 = qJD(1) * t615 + qJD(4);
	t812 = t610 * t729;
	t821 = t614 * t812 + t760 * t724;
	t819 = pkin(1) * (qJD(1) - qJD(2));
	t813 = t610 * t726;
	t784 = qJD(4) * t614;
	t781 = -pkin(1) * (qJD(1) + qJD(2)) / 0.2e1;
	t780 = -t819 / 0.2e1;
	t779 = t819 / 0.2e1;
	t776 = t818 / 0.2e1;
	t764 = -t802 / 0.2e1;
	t763 = t810 * t843;
	t761 = qJD(1) + t783;
	t751 = t761 * t729;
	t750 = t760 * t729;
	t747 = -r_i_i_C(3) * t614 - t757 * t615;
	t742 = t729 * t746;
	t741 = t746 * t724;
	t718 = qJ(1) - qJ(2);
	t717 = qJ(1) + qJ(2);
	t710 = cos(t718);
	t709 = sin(t718);
	t693 = cos(t717) * t781;
	t692 = sin(t717) * t781;
	t609 = t726 * t750 + (-t614 * t813 - t761 * t721) * t724;
	t608 = t761 * t726 * t724 + (-t614 * t814 + t750) * t721;
	t607 = t721 * t751 + t821 * t726;
	t606 = -t821 * t721 + t726 * t751;
	t605 = ((-t642 * t836 - t763) * t835 - t842 + t840 * t829 - t841) * t619 + t845;
	t604 = t845 + (-t842 - t835 * t763 + (t836 * t843 + t840 * t844) * t642) * t619;
	t1 = [t609 * r_i_i_C(1) - t608 * r_i_i_C(2) + t710 * t780 + t693 + t704 * t776 + t687 - pkin(12) * t785 + (t614 * t785 + t615 * t814) * r_i_i_C(3) + (t769 + t766) * pkin(10) + (t764 + t767) * pkin(8), t604 * t742 + t710 * t779 + t693 + t849, t605 * t742 + t849, t606 * r_i_i_C(1) - t607 * r_i_i_C(2); t607 * r_i_i_C(1) + t606 * r_i_i_C(2) + t692 + t709 * t780 + t686 + t703 * t776 - pkin(12) * t786 + (t614 * t786 - t615 * t812) * r_i_i_C(3) + (t764 + t768) * pkin(10) + (t769 + t765) * pkin(8), t604 * t741 + t709 * t779 + t692 + t850, t605 * t741 + t850, t608 * r_i_i_C(1) + t609 * r_i_i_C(2); 0, -qJD(2) * t723 * pkin(1) + t747 * t604 + t848, t747 * t605 + t848, (t615 * t813 - t721 * t784) * r_i_i_C(2) + (t610 * t615 * t721 + t726 * t784) * r_i_i_C(1);];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (71->17), mult. (178->36), div. (0->0), fcn. (200->8), ass. (0->19)
	t49 = sin(pkin(15));
	t50 = sin(pkin(14));
	t53 = cos(pkin(15));
	t54 = cos(pkin(14));
	t45 = t49 * t54 - t53 * t50;
	t46 = t49 * t50 + t53 * t54;
	t47 = sin(qJ(2));
	t51 = cos(qJ(2));
	t40 = -t51 * t45 - t47 * t46;
	t38 = t40 * qJD(2);
	t63 = -t45 * t47 + t46 * t51;
	t39 = t63 * qJD(2);
	t57 = t38 * r_i_i_C(1) - r_i_i_C(2) * t39;
	t48 = sin(qJ(1));
	t59 = qJD(1) * t48;
	t52 = cos(qJ(1));
	t58 = qJD(1) * t52;
	t55 = -r_i_i_C(1) * t63 - r_i_i_C(2) * t40 + pkin(6);
	t1 = [-t57 * t48 + (-r_i_i_C(3) * t48 + t55 * t52) * qJD(1), (-t38 * t52 + t59 * t63) * r_i_i_C(2) + (-t39 * t52 - t40 * t59) * r_i_i_C(1), 0, 0; t57 * t52 + (r_i_i_C(3) * t52 + t55 * t48) * qJD(1), (-t38 * t48 - t58 * t63) * r_i_i_C(2) + (-t39 * t48 + t40 * t58) * r_i_i_C(1), 0, 0; 0, t57, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:53
	% EndTime: 2020-05-07 02:13:53
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (1570->33), mult. (3340->62), div. (252->4), fcn. (4607->11), ass. (0->40)
	t72 = sin(pkin(18));
	t73 = cos(pkin(18));
	t76 = sin(pkin(15));
	t79 = cos(pkin(15));
	t69 = t79 * t72 + t76 * t73;
	t70 = -t76 * t72 + t79 * t73;
	t74 = sin(qJ(2));
	t77 = cos(qJ(2));
	t66 = t74 * t69 - t77 * t70;
	t63 = 0.1e1 / t66 ^ 2;
	t67 = t77 * t69 + t74 * t70;
	t91 = t67 ^ 2 * t63;
	t59 = 0.1e1 + t91;
	t57 = 0.1e1 / t59;
	t62 = 0.1e1 / t66;
	t86 = t62 * t66 + t91;
	t52 = -t86 * t57 + 0.1e1;
	t56 = qJ(2) + atan2(t67, t66);
	t54 = sin(t56);
	t55 = cos(t56);
	t87 = r_i_i_C(1) * t54 + r_i_i_C(2) * t55;
	t98 = qJD(1) * (pkin(1) * t74 + t87 * t52);
	t60 = t66 * qJD(2);
	t61 = t67 * qJD(2);
	t93 = t63 * t67;
	t51 = qJD(2) + (-t60 * t62 - t61 * t93) * t57;
	t90 = pkin(1) * qJD(2);
	t88 = t74 * t90;
	t97 = t87 * t51 + t88;
	t94 = t51 * t52;
	t92 = t62 * t91;
	t89 = t60 * t93;
	t85 = -t77 * pkin(1) - r_i_i_C(1) * t55 + r_i_i_C(2) * t54 - pkin(12);
	t50 = 0.2e1 * t86 / t59 ^ 2 * (-t61 * t92 - t89) + (0.2e1 * t89 + (t63 * t66 - t62 + 0.2e1 * t92) * t61) * t57;
	t84 = -t50 * t54 - t55 * t94;
	t83 = -t50 * t55 + t54 * t94;
	t82 = t84 * r_i_i_C(1) + t83 * r_i_i_C(2) - t77 * t90;
	t78 = cos(qJ(1));
	t75 = sin(qJ(1));
	t1 = [t97 * t75 + (-r_i_i_C(3) * t75 + t85 * t78) * qJD(1), t75 * t98 + t82 * t78, 0, 0; -t97 * t78 + (r_i_i_C(3) * t78 + t85 * t75) * qJD(1), t82 * t75 - t78 * t98, 0, 0; 0, -t83 * r_i_i_C(1) + t84 * r_i_i_C(2) - t88, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:54
	% EndTime: 2020-05-07 02:13:56
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (12798->77), mult. (23149->145), div. (930->8), fcn. (31677->18), ass. (0->80)
	t203 = pkin(17) + pkin(18);
	t201 = sin(t203);
	t202 = cos(t203);
	t208 = sin(pkin(15));
	t209 = cos(qJ(3));
	t212 = cos(pkin(15));
	t259 = sin(qJ(3));
	t196 = t209 * t208 + t259 * t212;
	t206 = sin(qJ(2));
	t210 = cos(qJ(2));
	t265 = t259 * t208 - t209 * t212;
	t266 = t206 * t196 + t210 * t265;
	t267 = t196 * t210 - t206 * t265;
	t274 = t201 * t267 + t202 * t266;
	t156 = 0.1e1 / t274 ^ 2;
	t275 = -t201 * t266 + t202 * t267;
	t282 = t275 ^ 2 * t156;
	t193 = t265 * qJD(3);
	t194 = t196 * qJD(3);
	t169 = qJD(2) * t267 - t193 * t206 + t194 * t210;
	t261 = qJD(2) * t266 + t193 * t210 + t194 * t206;
	t148 = t169 * t202 - t201 * t261;
	t155 = 0.1e1 / t274;
	t271 = t148 * t155;
	t281 = t271 * t282;
	t151 = 0.1e1 + t282;
	t149 = 0.1e1 / t151;
	t250 = t156 * t275;
	t280 = -t155 * t274 - t250 * t275;
	t134 = t280 * t149;
	t147 = t169 * t201 + t261 * t202;
	t130 = -0.2e1 * t280 * (-t147 * t250 - t281) / t151 ^ 2 + (-t271 + 0.2e1 * t281 + (0.2e1 * t147 * t275 + t148 * t274) * t156) * t149;
	t204 = sin(pkin(18));
	t205 = cos(pkin(18));
	t195 = t212 * t204 + t208 * t205;
	t227 = -t208 * t204 + t212 * t205;
	t182 = t206 * t195 - t210 * t227;
	t178 = 0.1e1 / t182;
	t179 = 0.1e1 / t182 ^ 2;
	t183 = t210 * t195 + t206 * t227;
	t247 = t183 ^ 2 * t179;
	t276 = t178 * t182 + t247;
	t176 = t182 * qJD(2);
	t173 = 0.1e1 + t247;
	t171 = 0.1e1 / t173;
	t177 = t183 * qJD(2);
	t237 = t171 * t177 * t179;
	t249 = t171 * t178;
	t137 = t176 * t249 + t183 * t237 - qJD(2);
	t132 = (-t147 * t155 - t148 * t250) * t149 + t137;
	t170 = pkin(17) - atan2(t183, t182) - qJ(2);
	t164 = sin(t170);
	t141 = -atan2(-t275, t274) + t170;
	t139 = sin(t141);
	t140 = cos(t141);
	t230 = r_i_i_C(1) * t139 - r_i_i_C(2) * t140;
	t258 = pkin(1) * qJD(2);
	t238 = t206 * t258;
	t263 = -pkin(3) * t137 * t164 + t230 * t132 - t238;
	t138 = t276 * t171 - 0.1e1;
	t133 = t138 + t134;
	t253 = t138 * t164;
	t262 = -pkin(1) * t206 - pkin(3) * t253 + t230 * t133;
	t257 = t132 * t139;
	t256 = t132 * t140;
	t207 = sin(qJ(1));
	t255 = t132 * t207;
	t211 = cos(qJ(1));
	t254 = t132 * t211;
	t240 = qJD(1) * t207;
	t239 = qJD(1) * t211;
	t225 = t130 * t230;
	t131 = t177 * t249 - t182 * t237 + (-0.2e1 * t171 + 0.2e1 * t276 / t173 ^ 2) * (t183 * t179 * t176 + t178 * t177 * t247);
	t129 = t131 + t130;
	t224 = t129 * t139 + t133 * t256;
	t223 = -t129 * t140 + t133 * t257;
	t165 = cos(t170);
	t219 = -t210 * pkin(1) - pkin(3) * t165 + r_i_i_C(1) * t140 + r_i_i_C(2) * t139 - pkin(12);
	t217 = -t210 * t258 + t223 * r_i_i_C(2) + t224 * r_i_i_C(1) + (-t137 * t138 * t165 - t131 * t164) * pkin(3);
	t1 = [-t263 * t207 + (-r_i_i_C(3) * t207 + t219 * t211) * qJD(1), t217 * t211 - t262 * t240, t211 * t225 + ((t139 * t254 + t140 * t240) * r_i_i_C(2) + (-t139 * t240 + t140 * t254) * r_i_i_C(1)) * t134, 0; t263 * t211 + (r_i_i_C(3) * t211 + t219 * t207) * qJD(1), t217 * t207 + t262 * t239, t207 * t225 + ((t139 * t255 - t140 * t239) * r_i_i_C(2) + (t139 * t239 + t140 * t255) * r_i_i_C(1)) * t134, 0; 0, -t238 + t224 * r_i_i_C(2) - t223 * r_i_i_C(1) + (-t131 * t165 + t137 * t253) * pkin(3), (t130 * t139 + t134 * t256) * r_i_i_C(2) + (t130 * t140 - t134 * t257) * r_i_i_C(1), 0;];
	JaD_transl = t1;
end