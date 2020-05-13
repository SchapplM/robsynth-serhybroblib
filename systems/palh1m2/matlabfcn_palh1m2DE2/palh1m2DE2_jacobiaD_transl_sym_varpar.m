% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m2DE2
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
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh1m2DE2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m2DE2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobiaD_transl_sym_varpar: pkin has to be [22x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
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
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16;
	t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 - pkin(15);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (-t16 * t26 + t18 * t23) * r_i_i_C(2) + (t16 * t23 + t18 * t26) * r_i_i_C(1), 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t25 + t18 * t24) * r_i_i_C(2) + (t16 * t24 - t18 * t25) * r_i_i_C(1), 0, 0; 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t37 = qJD(2) + qJD(3);
	t38 = qJ(2) + qJ(3);
	t36 = cos(t38);
	t56 = r_i_i_C(2) * t36;
	t35 = sin(t38);
	t58 = r_i_i_C(1) * t35;
	t47 = t56 + t58;
	t45 = t47 * t37;
	t41 = cos(qJ(2));
	t59 = pkin(1) * t41;
	t60 = qJD(2) * t59 + t45;
	t57 = r_i_i_C(2) * t35;
	t55 = t36 * t37;
	t40 = sin(qJ(1));
	t54 = qJD(1) * t40;
	t42 = cos(qJ(1));
	t53 = qJD(1) * t42;
	t39 = sin(qJ(2));
	t52 = qJD(2) * t39;
	t51 = r_i_i_C(1) * t55;
	t50 = t37 * t57;
	t49 = qJD(1) * t56;
	t46 = pkin(1) * t39 - r_i_i_C(1) * t36 - pkin(15) + t57;
	t44 = t40 * t49 + t54 * t58 + (t50 - t51) * t42;
	t30 = t40 * t50;
	t1 = [t60 * t40 + (-r_i_i_C(3) * t40 + t46 * t42) * qJD(1), (t41 * t54 + t42 * t52) * pkin(1) + t44, t44, 0; -t60 * t42 + (r_i_i_C(3) * t42 + t46 * t40) * qJD(1), t30 + (pkin(1) * t52 - t51) * t40 + (-t47 - t59) * t53, -t42 * t49 + t30 + (-t35 * t53 - t40 * t55) * r_i_i_C(1), 0; 0, -t60, -t45, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:31
	% EndTime: 2020-05-02 21:08:34
	% DurationCPUTime: 2.16s
	% Computational Cost: add. (18993->66), mult. (35250->109), div. (516->4), fcn. (50285->17), ass. (0->76)
	t179 = pkin(22) + pkin(21);
	t175 = sin(t179);
	t176 = cos(t179);
	t182 = cos(pkin(20));
	t186 = sin(pkin(18));
	t232 = sin(pkin(20));
	t235 = cos(pkin(18));
	t170 = t186 * t182 - t235 * t232;
	t171 = t235 * t182 + t186 * t232;
	t183 = sin(qJ(3));
	t187 = cos(qJ(3));
	t160 = t170 * t187 - t183 * t171;
	t188 = cos(qJ(2));
	t209 = qJD(2) * t188;
	t184 = sin(qJ(2));
	t210 = qJD(2) * t184;
	t241 = t183 * t170 + t171 * t187;
	t239 = t241 * qJD(3);
	t251 = t160 * qJD(3);
	t264 = -t184 * t239 + t251 * t188;
	t250 = -t160 * t209 + t210 * t241 - t264;
	t257 = t184 * t160 + t188 * t241;
	t260 = qJD(2) * t257 + t184 * t251 + t239 * t188;
	t129 = -t250 * t175 + t260 * t176;
	t223 = t260 * t175;
	t130 = t250 * t176 + t223;
	t202 = t160 * t188 - t184 * t241;
	t262 = t257 * t175 - t202 * t176;
	t134 = t262 ^ 2;
	t261 = t202 * t175 + t176 * t257;
	t136 = 0.1e1 / t261 ^ 2;
	t133 = t134 * t136 + 0.1e1;
	t131 = 0.1e1 / t133;
	t135 = 0.1e1 / t261;
	t180 = qJD(2) + qJD(3);
	t259 = t136 * t262;
	t122 = (-t129 * t135 - t130 * t259) * t131 + t180;
	t194 = t135 * t261 + t262 * t259;
	t269 = -t194 * t131 + 0.1e1;
	t272 = t122 * t269;
	t181 = qJ(2) + qJ(3);
	t128 = atan2(t262, -t261) + t181;
	t127 = cos(t128);
	t271 = t127 * t272;
	t126 = sin(t128);
	t270 = t126 * t272;
	t177 = sin(t181);
	t203 = r_i_i_C(1) * t126 + r_i_i_C(2) * t127;
	t237 = pkin(5) * t177 + t203 * t269;
	t245 = t130 * t136;
	t228 = t135 * t245;
	t255 = t129 * t259;
	t267 = 0.2e1 * t194 * (t134 * t228 + t255) / t133 ^ 2;
	t266 = -0.2e1 * t262;
	t265 = t262 * t228 * t266;
	t233 = pkin(5) * t180;
	t207 = t177 * t233;
	t168 = -pkin(1) * t209 - t207;
	t238 = t203 * t122 - t168;
	t236 = t188 * pkin(1) + t237;
	t185 = sin(qJ(1));
	t212 = qJD(1) * t185;
	t189 = cos(qJ(1));
	t211 = qJD(1) * t189;
	t178 = cos(t181);
	t206 = t178 * t233;
	t200 = t184 * pkin(1) - pkin(5) * t178 - r_i_i_C(1) * t127 + r_i_i_C(2) * t126 - pkin(15);
	t120 = t267 + (-((t202 * qJD(2) + t264) * t176 - t223) * t135 + t265 + (t129 * t266 - t261 * t130) * t136) * t131;
	t198 = -t120 * t126 - t271;
	t197 = -t120 * t127 + t270;
	t121 = (t130 * t135 - t261 * t245 - 0.2e1 * t255 + t265) * t131 + t267;
	t196 = -t121 * t126 - t271;
	t195 = -t121 * t127 + t270;
	t193 = pkin(1) * t210 + t198 * r_i_i_C(1) + t197 * r_i_i_C(2) - t206;
	t192 = t196 * r_i_i_C(1) + t195 * r_i_i_C(2) - t206;
	t1 = [t238 * t185 + (-r_i_i_C(3) * t185 + t200 * t189) * qJD(1), t193 * t189 + t236 * t212, t192 * t189 + t237 * t212, 0; -t238 * t189 + (r_i_i_C(3) * t189 + t200 * t185) * qJD(1), t193 * t185 - t236 * t211, t192 * t185 - t237 * t211, 0; 0, -t197 * r_i_i_C(1) + t198 * r_i_i_C(2) + t168, -t195 * r_i_i_C(1) + t196 * r_i_i_C(2) - t207, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:46
	% EndTime: 2020-05-02 21:08:50
	% DurationCPUTime: 3.67s
	% Computational Cost: add. (32774->112), mult. (60438->178), div. (882->4), fcn. (86521->30), ass. (0->112)
	t737 = qJ(2) + qJ(3);
	t680 = pkin(22) + pkin(21);
	t671 = sin(t680);
	t672 = cos(t680);
	t686 = cos(pkin(20));
	t691 = sin(pkin(18));
	t758 = sin(pkin(20));
	t763 = cos(pkin(18));
	t658 = t691 * t686 - t763 * t758;
	t659 = t763 * t686 + t691 * t758;
	t688 = sin(qJ(3));
	t693 = cos(qJ(3));
	t646 = t658 * t693 - t688 * t659;
	t694 = cos(qJ(2));
	t689 = sin(qJ(2));
	t768 = t688 * t658 + t659 * t693;
	t739 = t689 * t768;
	t712 = t646 * t694 - t739;
	t784 = t689 * t646 + t694 * t768;
	t788 = t712 * t671 + t672 * t784;
	t789 = t784 * t671 - t712 * t672;
	t614 = atan2(t789, -t788) + t737;
	t612 = sin(t614);
	t613 = cos(t614);
	t679 = qJ(1) - t737;
	t669 = sin(t679);
	t690 = sin(qJ(1));
	t695 = cos(qJ(1));
	t731 = qJD(2) * t694;
	t766 = t768 * qJD(3);
	t778 = t646 * qJD(3);
	t791 = -t689 * t766 + t778 * t694;
	t777 = qJD(2) * t739 - t646 * t731 - t791;
	t787 = qJD(2) * t784 + t689 * t778 + t766 * t694;
	t615 = -t777 * t671 + t787 * t672;
	t748 = t787 * t671;
	t616 = t777 * t672 + t748;
	t620 = t789 ^ 2;
	t622 = 0.1e1 / t788 ^ 2;
	t619 = t620 * t622 + 0.1e1;
	t617 = 0.1e1 / t619;
	t621 = 0.1e1 / t788;
	t681 = qJD(2) + qJD(3);
	t786 = t622 * t789;
	t608 = (-t615 * t621 - t616 * t786) * t617 + t681;
	t687 = sin(qJ(4));
	t692 = cos(qJ(4));
	t714 = r_i_i_C(1) * t687 + r_i_i_C(2) * t692;
	t704 = -r_i_i_C(3) * t608 + t714 * qJD(4);
	t715 = r_i_i_C(1) * t692 - r_i_i_C(2) * t687;
	t708 = t715 * t695;
	t760 = pkin(5) * (qJD(1) - t681);
	t721 = t760 / 0.2e1;
	t756 = t608 * t690;
	t707 = t621 * t788 + t789 * t786;
	t796 = -t707 * t617 + 0.1e1;
	t799 = t669 * t721 + ((r_i_i_C(3) * qJD(1) * t695 - t715 * t756) * t613 + (-qJD(1) * t708 + t704 * t690) * t612) * t796;
	t670 = cos(t679);
	t733 = qJD(1) * t690;
	t798 = t670 * t721 + ((-r_i_i_C(3) * t733 - t608 * t708) * t613 + (t704 * t695 + t715 * t733) * t612) * t796;
	t706 = r_i_i_C(3) * t613 - t715 * t612;
	t729 = qJD(4) * t613;
	t797 = (t706 * t608 - t714 * t729) * t796 - pkin(5) * t681 * sin(t737);
	t772 = t616 * t622;
	t753 = t621 * t772;
	t782 = t615 * t786;
	t794 = 0.2e1 * t707 * (t620 * t753 + t782) / t619 ^ 2;
	t793 = -0.2e1 * t789;
	t792 = t789 * t753 * t793;
	t713 = pkin(18) - pkin(20) - t680;
	t666 = qJ(1) + t713;
	t765 = sin(t666) / 0.2e1;
	t667 = -qJ(1) + t713;
	t764 = cos(t667) / 0.2e1;
	t762 = pkin(1) * (qJD(1) + qJD(2));
	t761 = pkin(1) * (qJD(1) - qJD(2));
	t759 = t612 * r_i_i_C(3);
	t757 = t608 * t613;
	t755 = t608 * t692;
	t738 = t692 * t695;
	t678 = qJ(1) + t737;
	t723 = -pkin(5) * (qJD(1) + t681) / 0.2e1;
	t654 = sin(t678) * t723;
	t684 = qJ(1) + qJ(2);
	t735 = t654 - cos(t684) * t762 / 0.2e1;
	t655 = cos(t678) * t723;
	t734 = t655 + sin(t684) * t762 / 0.2e1;
	t730 = qJD(4) * t612;
	t726 = r_i_i_C(3) * t757;
	t725 = -t761 / 0.2e1;
	t724 = t761 / 0.2e1;
	t722 = -t760 / 0.2e1;
	t719 = -pkin(15) - t759;
	t717 = -qJD(1) + t729;
	t716 = qJD(1) * t613 - qJD(4);
	t710 = t717 * t687;
	t705 = t715 * t613 + t759;
	t703 = t608 * t612 * t695 + t716 * t690;
	t702 = t706 * t690;
	t701 = t706 * t695;
	t685 = qJ(1) - qJ(2);
	t677 = cos(t685);
	t676 = sin(t685);
	t664 = cos(t666);
	t663 = sin(t667);
	t607 = -t716 * t738 + (t612 * t755 + t710) * t690;
	t606 = t717 * t692 * t690 + (-t612 * t756 + t716 * t695) * t687;
	t605 = t703 * t692 + t695 * t710;
	t604 = t703 * t687 - t717 * t738;
	t603 = (t616 * t621 - t788 * t772 - 0.2e1 * t782 + t792) * t617 + t794;
	t602 = t794 + (-((t712 * qJD(2) + t791) * t672 - t748) * t621 + t792 + (t615 * t793 - t788 * t616) * t622) * t617;
	t1 = [t607 * r_i_i_C(1) + t606 * r_i_i_C(2) - t690 * t726 + t676 * t725 + t670 * t722 + (t719 * t695 + (t663 / 0.2e1 + t765) * pkin(11) + (t664 / 0.2e1 + t764) * pkin(9)) * qJD(1) + t734, t602 * t701 + t676 * t724 + t734 + t798, t603 * t701 + t655 + t798, t604 * r_i_i_C(1) + t605 * r_i_i_C(2); -t605 * r_i_i_C(1) + t604 * r_i_i_C(2) + t695 * t726 + t677 * t724 + t669 * t722 + (t719 * t690 + (-t664 / 0.2e1 + t764) * pkin(11) + (-t663 / 0.2e1 + t765) * pkin(9)) * qJD(1) + t735, t602 * t702 + t677 * t725 + t735 + t799, t603 * t702 + t654 + t799, -t606 * r_i_i_C(1) + t607 * r_i_i_C(2); 0, -pkin(1) * t731 + t705 * t602 + t797, t705 * t603 + t797, (-t613 * t755 + t687 * t730) * r_i_i_C(2) + (-t687 * t757 - t692 * t730) * r_i_i_C(1);];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (71->17), mult. (178->36), div. (0->0), fcn. (200->8), ass. (0->19)
	t49 = sin(pkin(18));
	t50 = sin(pkin(17));
	t53 = cos(pkin(18));
	t54 = cos(pkin(17));
	t45 = t49 * t54 - t53 * t50;
	t46 = t50 * t49 + t53 * t54;
	t47 = sin(qJ(2));
	t51 = cos(qJ(2));
	t40 = -t47 * t45 - t46 * t51;
	t38 = t40 * qJD(2);
	t41 = t45 * t51 - t46 * t47;
	t39 = t41 * qJD(2);
	t56 = t38 * r_i_i_C(1) - r_i_i_C(2) * t39;
	t48 = sin(qJ(1));
	t58 = qJD(1) * t48;
	t52 = cos(qJ(1));
	t57 = qJD(1) * t52;
	t55 = -r_i_i_C(1) * t41 - r_i_i_C(2) * t40 + pkin(14);
	t1 = [-t56 * t48 + (-r_i_i_C(3) * t48 + t55 * t52) * qJD(1), (-t38 * t52 + t41 * t58) * r_i_i_C(2) + (-t39 * t52 - t40 * t58) * r_i_i_C(1), 0, 0; t56 * t52 + (r_i_i_C(3) * t52 + t55 * t48) * qJD(1), (-t38 * t48 - t41 * t57) * r_i_i_C(2) + (-t39 * t48 + t40 * t57) * r_i_i_C(1), 0, 0; 0, t56, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->19), mult. (50->32), div. (0->0), fcn. (27->8), ass. (0->18)
	t38 = r_i_i_C(1) / 0.2e1;
	t37 = -r_i_i_C(2) / 0.2e1;
	t31 = cos(qJ(2));
	t36 = qJD(1) * t31;
	t29 = sin(qJ(2));
	t35 = qJD(2) * t29;
	t34 = pkin(18) - pkin(22);
	t33 = qJD(2) * t31 * pkin(1);
	t32 = cos(qJ(1));
	t30 = sin(qJ(1));
	t28 = -qJ(1) + t34;
	t27 = qJ(1) + t34;
	t26 = -t29 * pkin(1) + pkin(15);
	t25 = cos(t28);
	t24 = sin(t27);
	t23 = qJD(1) * cos(t27) / 0.2e1;
	t22 = -qJD(1) * sin(t28) / 0.2e1;
	t1 = [t30 * t33 + t23 * r_i_i_C(1) + t22 * r_i_i_C(2) + (-t30 * r_i_i_C(3) + t24 * t37 + t25 * t38 - t32 * t26) * qJD(1), (t30 * t36 + t32 * t35) * pkin(1), 0, 0; -t32 * t33 + t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + (t32 * r_i_i_C(3) + t24 * t38 + t25 * t37 - t30 * t26) * qJD(1), (t30 * t35 - t32 * t36) * pkin(1), 0, 0; 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (716->30), mult. (1642->68), div. (254->4), fcn. (2068->9), ass. (0->35)
	t65 = cos(pkin(19));
	t67 = cos(qJ(3));
	t78 = sin(pkin(19));
	t85 = sin(qJ(3));
	t61 = t85 * t65 + t67 * t78;
	t58 = 0.1e1 / t61 ^ 2;
	t62 = t67 * t65 - t85 * t78;
	t60 = t62 ^ 2;
	t55 = t60 * t58 + 0.1e1;
	t86 = 0.1e1 / t61;
	t53 = 0.1e1 / t55;
	t56 = t61 * qJD(3);
	t57 = t62 * qJD(3);
	t79 = t58 * t62;
	t47 = qJD(2) + (t56 * t86 + t57 * t79) * t53;
	t48 = t55 * t53;
	t84 = t47 * t48;
	t66 = sin(qJ(1));
	t83 = t47 * t66;
	t68 = cos(qJ(1));
	t82 = t47 * t68;
	t80 = t57 * t86 * t58;
	t81 = (-t56 * t79 - t60 * t80) / t55 ^ 2;
	t77 = qJD(1) * t66;
	t76 = qJD(1) * t68;
	t52 = qJ(2) + atan2(-t62, t61);
	t50 = sin(t52);
	t51 = cos(t52);
	t75 = r_i_i_C(1) * t51 - r_i_i_C(2) * t50;
	t74 = r_i_i_C(1) * t50 + r_i_i_C(2) * t51 - pkin(15);
	t73 = t75 * t68;
	t72 = (t50 * t76 + t51 * t83) * r_i_i_C(2) + (t50 * t83 - t51 * t76) * r_i_i_C(1);
	t71 = (-t50 * t77 + t51 * t82) * r_i_i_C(2) + (t50 * t82 + t51 * t77) * r_i_i_C(1);
	t46 = -0.2e1 * t81 + 0.2e1 * (-t53 * t56 * t58 + (-t53 * t80 - t58 * t81) * t62) * t62;
	t1 = [t75 * t83 + (-r_i_i_C(3) * t66 + t74 * t68) * qJD(1), t71, -t46 * t73 + t71 * t48, 0; -t47 * t73 + (r_i_i_C(3) * t68 + t74 * t66) * qJD(1), t72, -t75 * t66 * t46 + t72 * t48, 0; 0, -t75 * t47, (-t46 * t51 + t50 * t84) * r_i_i_C(2) + (-t46 * t50 - t51 * t84) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (1824->39), mult. (4087->67), div. (669->4), fcn. (5236->12), ass. (0->48)
	t111 = cos(pkin(19));
	t113 = cos(qJ(3));
	t138 = sin(pkin(19));
	t146 = sin(qJ(3));
	t105 = t111 * t146 + t113 * t138;
	t107 = t113 * t111 - t146 * t138;
	t89 = qJ(2) + atan2(-t107, t105);
	t165 = pkin(2) * sin(t89);
	t147 = pkin(2) * cos(t89);
	t149 = 0.1e1 / t105 ^ 2;
	t152 = t107 ^ 2 * t149;
	t158 = 0.1e1 + t152;
	t161 = 0.1e1 / t158;
	t78 = t158 * t161;
	t164 = t78 * t147;
	t101 = 0.1e1 / t105;
	t97 = t105 * qJD(3);
	t162 = t161 * t101 * t97;
	t134 = qJD(2) + t162;
	t140 = t107 * t149;
	t98 = t107 * qJD(3);
	t77 = t140 * t161 * t98 + t134;
	t163 = t77 * t165;
	t137 = t101 * t152;
	t154 = t97 * t140;
	t85 = atan2(-t107, -t105) + t89;
	t82 = sin(t85);
	t83 = cos(t85);
	t127 = r_i_i_C(1) * t83 - r_i_i_C(2) * t82;
	t124 = t101 * t105 + t152;
	t75 = -t124 * t161 + t78;
	t151 = t127 * t75 - t164;
	t74 = t134 - t162;
	t118 = t127 * t74 - t77 * t147;
	t145 = t74 * t75;
	t112 = sin(qJ(1));
	t133 = qJD(1) * t112;
	t114 = cos(qJ(1));
	t132 = qJD(1) * t114;
	t126 = -r_i_i_C(1) * t82 - r_i_i_C(2) * t83;
	t123 = -pkin(15) + t126 + t165;
	t72 = -0.2e1 * t124 / t158 ^ 2 * (t98 * t137 + t154) + (0.2e1 * t154 + (t105 * t149 - t101 + 0.2e1 * t137) * t98) * t161;
	t122 = t145 * t83 + t72 * t82;
	t121 = -t145 * t82 + t72 * t83;
	t120 = -t127 + t147;
	t119 = t126 * t74 + t163;
	t117 = t121 * r_i_i_C(1) - t122 * r_i_i_C(2) + t78 * t163;
	t1 = [-t118 * t112 + (-r_i_i_C(3) * t112 + t114 * t123) * qJD(1), t114 * t119 + t120 * t133, t117 * t114 - t151 * t133, 0; t118 * t114 + (r_i_i_C(3) * t114 + t112 * t123) * qJD(1), t112 * t119 - t120 * t132, t117 * t112 + t151 * t132, 0; 0, t118, t122 * r_i_i_C(1) + t121 * r_i_i_C(2) - t77 * t164, 0;];
	JaD_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:33
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (13050->79), mult. (23662->153), div. (957->8), fcn. (32457->18), ass. (0->84)
	t222 = pkin(22) + pkin(21);
	t220 = sin(t222);
	t221 = cos(t222);
	t277 = sin(qJ(3));
	t278 = sin(pkin(18));
	t279 = cos(qJ(3));
	t280 = cos(pkin(18));
	t213 = -t277 * t278 - t279 * t280;
	t214 = -t277 * t280 + t279 * t278;
	t224 = sin(qJ(2));
	t226 = cos(qJ(2));
	t244 = t213 * t226 - t224 * t214;
	t284 = t224 * t213 + t226 * t214;
	t299 = t220 * t284 - t244 * t221;
	t167 = 0.1e1 / t299 ^ 2;
	t293 = t220 * t244 + t221 * t284;
	t302 = t293 ^ 2 * t167;
	t209 = t213 * qJD(3);
	t210 = t214 * qJD(3);
	t233 = t244 * qJD(2) + t209 * t226 - t224 * t210;
	t234 = qJD(2) * t284 + t209 * t224 + t210 * t226;
	t160 = t220 * t233 + t234 * t221;
	t166 = 0.1e1 / t299;
	t300 = t160 * t166;
	t301 = t300 * t302;
	t163 = 0.1e1 + t302;
	t161 = 0.1e1 / t163;
	t267 = t167 * t293;
	t297 = -t166 * t299 - t267 * t293;
	t144 = t297 * t161;
	t298 = -t220 * t234 + t221 * t233;
	t140 = -0.2e1 * (t267 * t298 - t301) / t163 ^ 2 * t297 + (-t300 + 0.2e1 * t301 + (t160 * t299 - 0.2e1 * t293 * t298) * t167) * t161;
	t223 = sin(pkin(22));
	t275 = cos(pkin(22));
	t211 = t278 * t223 + t275 * t280;
	t212 = t280 * t223 - t278 * t275;
	t196 = t224 * t211 + t212 * t226;
	t191 = 0.1e1 / t196;
	t192 = 0.1e1 / t196 ^ 2;
	t285 = t211 * t226 - t212 * t224;
	t263 = t192 * t285;
	t294 = t191 * t196 + t285 * t263;
	t188 = t285 * qJD(2);
	t190 = t196 * qJD(2);
	t194 = t285 ^ 2;
	t185 = t194 * t192 + 0.1e1;
	t183 = 0.1e1 / t185;
	t253 = t183 * t188 * t285;
	t266 = t183 * t191;
	t147 = t190 * t266 + t192 * t253 - qJD(2);
	t142 = (-t160 * t267 + t166 * t298) * t161 + t147;
	t177 = pkin(21) - atan2(t285, t196) - qJ(2);
	t176 = cos(t177);
	t151 = -atan2(-t293, t299) + t177;
	t149 = sin(t151);
	t150 = cos(t151);
	t246 = r_i_i_C(1) * t150 + r_i_i_C(2) * t149;
	t276 = pkin(1) * qJD(2);
	t254 = t226 * t276;
	t282 = -pkin(4) * t147 * t176 + t246 * t142 + t254;
	t148 = t294 * t183 - 0.1e1;
	t143 = t148 + t144;
	t270 = t148 * t176;
	t281 = pkin(1) * t226 - pkin(4) * t270 + t246 * t143;
	t274 = t142 * t149;
	t273 = t142 * t150;
	t225 = sin(qJ(1));
	t272 = t142 * t225;
	t227 = cos(qJ(1));
	t271 = t142 * t227;
	t265 = t183 * t192;
	t256 = qJD(1) * t225;
	t255 = qJD(1) * t227;
	t251 = t190 * t263;
	t241 = t140 * t246;
	t193 = t191 * t192;
	t141 = -t183 * t251 - 0.2e1 * t294 * (-t194 * t193 * t188 - t251) / t185 ^ 2 + (-t190 * t265 - 0.2e1 * t193 * t253) * t285 + (-t196 * t265 + t266) * t188;
	t139 = t140 + t141;
	t240 = t139 * t149 + t143 * t273;
	t239 = -t139 * t150 + t143 * t274;
	t175 = sin(t177);
	t235 = t224 * pkin(1) - pkin(4) * t175 + r_i_i_C(1) * t149 - r_i_i_C(2) * t150 - pkin(15);
	t232 = t224 * t276 - t240 * r_i_i_C(2) + t239 * r_i_i_C(1) + (-t147 * t148 * t175 + t141 * t176) * pkin(4);
	t1 = [t282 * t225 + (-r_i_i_C(3) * t225 + t235 * t227) * qJD(1), t232 * t227 + t281 * t256, -t227 * t241 + ((t149 * t256 - t150 * t271) * r_i_i_C(2) + (t149 * t271 + t150 * t256) * r_i_i_C(1)) * t144, 0; -t282 * t227 + (r_i_i_C(3) * t227 + t235 * t225) * qJD(1), t232 * t225 - t281 * t255, -t225 * t241 + ((-t149 * t255 - t150 * t272) * r_i_i_C(2) + (t149 * t272 - t150 * t255) * r_i_i_C(1)) * t144, 0; 0, -t254 + t239 * r_i_i_C(2) + t240 * r_i_i_C(1) + (-t141 * t175 - t147 * t270) * pkin(4), (-t140 * t150 + t144 * t274) * r_i_i_C(2) + (t140 * t149 + t144 * t273) * r_i_i_C(1), 0;];
	JaD_transl = t1;
end