% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m2DE2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m2DE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m2DE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobia_transl_sym_varpar: pkin has to be [22x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = -pkin(15) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) + t5 * t2, t7 * t4, 0, 0; t2 * r_i_i_C(3) - t5 * t4, t7 * t2, 0, 0; 0, t6, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t5 = qJ(2) + qJ(3);
	t3 = sin(t5);
	t4 = cos(t5);
	t13 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t16 = t13 - sin(qJ(2)) * pkin(1);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(15) + t16;
	t10 = -cos(qJ(2)) * pkin(1) + t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t9 * r_i_i_C(3) - t11 * t7, t10 * t9, t12 * t9, 0; t7 * r_i_i_C(3) + t11 * t9, t10 * t7, t12 * t7, 0; 0, t16, t13, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (2493->20), mult. (4333->30), div. (72->0), fcn. (6505->17), ass. (0->27)
	t32 = qJ(2) + qJ(3);
	t41 = pkin(5) * sin(t32);
	t33 = cos(pkin(20));
	t37 = sin(pkin(18));
	t47 = sin(pkin(20));
	t51 = cos(pkin(18));
	t22 = t37 * t33 - t51 * t47;
	t23 = t51 * t33 + t37 * t47;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t17 = t22 * t38 - t34 * t23;
	t35 = sin(qJ(2));
	t39 = cos(qJ(2));
	t55 = t34 * t22 + t23 * t38;
	t60 = t35 * t17 + t55 * t39;
	t54 = -t17 * t39 + t35 * t55;
	t27 = pkin(5) * cos(t32);
	t46 = -t35 * pkin(1) + t27;
	t31 = pkin(22) + pkin(21);
	t28 = sin(t31);
	t29 = cos(t31);
	t5 = atan2(t60 * t28 + t54 * t29, t54 * t28 - t29 * t60) + t32;
	t43 = pkin(15) + t46 + r_i_i_C(1) * cos(t5) - r_i_i_C(2) * sin(t5);
	t42 = -t39 * pkin(1) - t41;
	t40 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t40 * r_i_i_C(3) - t43 * t36, t42 * t40, -t41 * t40, 0; t36 * r_i_i_C(3) + t43 * t40, t42 * t36, -t41 * t36, 0; 0, t46, t27, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (4032->53), mult. (6929->69), div. (108->0), fcn. (10365->30), ass. (0->58)
	t53 = cos(pkin(20));
	t58 = sin(pkin(18));
	t71 = sin(pkin(20));
	t82 = cos(pkin(18));
	t26 = t58 * t53 - t82 * t71;
	t27 = t82 * t53 + t58 * t71;
	t55 = sin(qJ(3));
	t60 = cos(qJ(3));
	t21 = t26 * t60 - t55 * t27;
	t56 = sin(qJ(2));
	t61 = cos(qJ(2));
	t92 = t55 * t26 + t27 * t60;
	t97 = t56 * t21 + t92 * t61;
	t91 = -t21 * t61 + t56 * t92;
	t90 = -pkin(1) / 0.2e1;
	t50 = pkin(22) + pkin(21);
	t67 = pkin(18) - pkin(20) - t50;
	t37 = -qJ(1) + t67;
	t88 = -sin(t37) / 0.2e1;
	t36 = qJ(1) + t67;
	t87 = -cos(t36) / 0.2e1;
	t52 = qJ(1) - qJ(2);
	t86 = pkin(1) * sin(t52);
	t85 = pkin(1) * cos(t52);
	t70 = qJ(2) + qJ(3);
	t49 = qJ(1) - t70;
	t84 = pkin(5) * sin(t49);
	t83 = pkin(5) * cos(t49);
	t54 = sin(qJ(4));
	t57 = sin(qJ(1));
	t77 = t57 * t54;
	t59 = cos(qJ(4));
	t76 = t57 * t59;
	t62 = cos(qJ(1));
	t75 = t62 * t54;
	t74 = t62 * t59;
	t48 = qJ(1) + t70;
	t32 = -pkin(5) * sin(t48) / 0.2e1;
	t51 = qJ(1) + qJ(2);
	t73 = t32 + cos(t51) * t90;
	t34 = pkin(5) * cos(t48) / 0.2e1;
	t72 = t34 + sin(t51) * t90;
	t44 = sin(t50);
	t45 = cos(t50);
	t9 = atan2(t97 * t44 + t91 * t45, t91 * t44 - t45 * t97) + t70;
	t7 = sin(t9);
	t69 = t7 * r_i_i_C(3) + pkin(15);
	t41 = pkin(5) * cos(t70);
	t35 = -t83 / 0.2e1;
	t33 = t84 / 0.2e1;
	t31 = cos(t37);
	t28 = sin(t36);
	t8 = cos(t9);
	t6 = t8 * t74 + t77;
	t5 = -t8 * t75 + t76;
	t4 = -t8 * t76 + t75;
	t3 = t8 * t77 + t74;
	t1 = [t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t85 / 0.2e1 - t84 / 0.2e1 - t69 * t57 + (t31 / 0.2e1 + t87) * pkin(11) + (t28 / 0.2e1 + t88) * pkin(9) + t73, -t85 / 0.2e1 + t33 + t73, t32 + t33, t5 * r_i_i_C(1) - t6 * r_i_i_C(2); t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t86 / 0.2e1 + t83 / 0.2e1 + t69 * t62 + (-t28 / 0.2e1 + t88) * pkin(11) + (-t31 / 0.2e1 + t87) * pkin(9) + t72, -t86 / 0.2e1 + t35 + t72, t34 + t35, -t3 * r_i_i_C(1) + t4 * r_i_i_C(2); 0, -t56 * pkin(1) + t41, t41, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t59) * t7;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (39->9), mult. (82->18), div. (0->0), fcn. (112->8), ass. (0->16)
	t11 = cos(qJ(2));
	t10 = sin(pkin(17));
	t13 = cos(pkin(18));
	t14 = cos(pkin(17));
	t9 = sin(pkin(18));
	t5 = -t13 * t10 + t9 * t14;
	t6 = t10 * t9 + t14 * t13;
	t7 = sin(qJ(2));
	t1 = -t6 * t11 - t7 * t5;
	t2 = t5 * t11 - t7 * t6;
	t19 = t2 * r_i_i_C(1) + r_i_i_C(2) * t1;
	t16 = r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t15 = -pkin(14) + t19;
	t12 = cos(qJ(1));
	t8 = sin(qJ(1));
	t3 = [t12 * r_i_i_C(3) - t15 * t8, t16 * t12, 0, 0; t8 * r_i_i_C(3) + t15 * t12, t16 * t8, 0, 0; 0, t19, 0, 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (31->17), mult. (23->18), div. (0->0), fcn. (19->8), ass. (0->13)
	t14 = sin(qJ(2)) * pkin(1);
	t13 = cos(qJ(2)) * pkin(1);
	t12 = pkin(18) - pkin(22);
	t11 = cos(qJ(1));
	t9 = sin(qJ(1));
	t7 = -qJ(1) + t12;
	t6 = qJ(1) + t12;
	t5 = pkin(15) - t14;
	t4 = cos(t6);
	t3 = sin(t7);
	t2 = -cos(t7) / 0.2e1;
	t1 = sin(t6) / 0.2e1;
	t8 = [(-t3 / 0.2e1 + t1) * r_i_i_C(1) + (t4 / 0.2e1 + t2) * r_i_i_C(2) + t11 * r_i_i_C(3) - t9 * t5, -t11 * t13, 0, 0; (t2 - t4 / 0.2e1) * r_i_i_C(1) + (t1 + t3 / 0.2e1) * r_i_i_C(2) + t9 * r_i_i_C(3) + t11 * t5, -t9 * t13, 0, 0; 0, -t14, 0, 0;];
	Ja_transl = t8;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (120->8), mult. (222->14), div. (30->0), fcn. (344->9), ass. (0->15)
	t21 = sin(qJ(3));
	t20 = sin(pkin(19));
	t11 = cos(pkin(19));
	t13 = cos(qJ(3));
	t4 = qJ(2) + atan2(-t11 * t13 + t20 * t21, t11 * t21 + t13 * t20);
	t2 = sin(t4);
	t3 = cos(t4);
	t19 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t2;
	t18 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t17 = -pkin(15) - t18;
	t12 = sin(qJ(1));
	t16 = t19 * t12;
	t14 = cos(qJ(1));
	t15 = t19 * t14;
	t1 = [t14 * r_i_i_C(3) + t12 * t17, t15, t15, 0; t12 * r_i_i_C(3) - t14 * t17, t16, t16, 0; 0, t18, t18, 0;];
	Ja_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (302->17), mult. (525->18), div. (81->0), fcn. (848->12), ass. (0->18)
	t21 = cos(pkin(19));
	t23 = cos(qJ(3));
	t31 = sin(pkin(19));
	t34 = sin(qJ(3));
	t16 = t34 * t21 + t23 * t31;
	t18 = t23 * t21 - t34 * t31;
	t9 = qJ(2) + atan2(-t18, t16);
	t5 = atan2(-t18, -t16) + t9;
	t3 = sin(t5);
	t36 = pkin(2) * sin(t9);
	t4 = cos(t5);
	t38 = t3 * r_i_i_C(1) + t4 * r_i_i_C(2) - t36;
	t35 = pkin(2) * cos(t9);
	t28 = pkin(15) + t38;
	t27 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3 - t35;
	t24 = cos(qJ(1));
	t22 = sin(qJ(1));
	t1 = [t24 * r_i_i_C(3) - t28 * t22, t27 * t24, -t35 * t24, 0; t22 * r_i_i_C(3) + t28 * t24, t27 * t22, -t35 * t22, 0; 0, t38, -t36, 0;];
	Ja_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (1787->30), mult. (2887->43), div. (126->2), fcn. (4395->18), ass. (0->36)
	t61 = sin(qJ(3));
	t62 = sin(pkin(18));
	t63 = cos(qJ(3));
	t64 = cos(pkin(18));
	t34 = -t61 * t62 - t63 * t64;
	t35 = -t61 * t64 + t63 * t62;
	t45 = sin(qJ(2));
	t47 = cos(qJ(2));
	t55 = t45 * t34 + t47 * t35;
	t44 = sin(pkin(22));
	t56 = cos(pkin(22));
	t32 = t62 * t44 + t56 * t64;
	t33 = t64 * t44 - t62 * t56;
	t21 = t45 * t32 + t33 * t47;
	t22 = t32 * t47 - t33 * t45;
	t15 = pkin(21) - atan2(t22, t21) - qJ(2);
	t66 = pkin(4) * sin(t15);
	t65 = t45 * pkin(1);
	t59 = t22 ^ 2 / t21 ^ 2;
	t54 = -t34 * t47 + t45 * t35;
	t43 = pkin(22) + pkin(21);
	t41 = sin(t43);
	t42 = cos(t43);
	t6 = -atan2(t41 * t54 - t42 * t55, t41 * t55 + t54 * t42) + t15;
	t4 = sin(t6);
	t5 = cos(t6);
	t53 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t4;
	t52 = r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t16 = 0.1e1 / (0.1e1 + t59);
	t3 = t16 * t59 + t16 - 0.1e1;
	t50 = pkin(15) - t65 - t52 + t66;
	t1 = -0.1e1 + t3;
	t49 = pkin(4) * cos(t15) * t3 - pkin(1) * t47 + t53 * t1;
	t48 = cos(qJ(1));
	t46 = sin(qJ(1));
	t2 = [t48 * r_i_i_C(3) - t50 * t46, t49 * t48, -t48 * t53, 0; t46 * r_i_i_C(3) + t50 * t48, t49 * t46, -t46 * t53, 0; 0, t52 * t1 - t3 * t66 - t65, -t52, 0;];
	Ja_transl = t2;
end