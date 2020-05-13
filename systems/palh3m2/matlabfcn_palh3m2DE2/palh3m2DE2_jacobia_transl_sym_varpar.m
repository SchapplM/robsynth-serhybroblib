% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m2DE2
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
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m2DE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m2DE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobia_transl_sym_varpar: pkin has to be [18x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
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
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(12) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, t6 * t4, 0, 0; t2 * r_i_i_C(3) + t5 * t4, t6 * t2, 0, 0; 0, t7, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->9), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->13)
	t10 = qJ(2) + qJ(3);
	t7 = sin(t10);
	t8 = cos(t10);
	t23 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t15 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t22 = t15 + cos(qJ(2)) * pkin(1);
	t12 = sin(qJ(1));
	t18 = t23 * t12;
	t13 = cos(qJ(1));
	t17 = t23 * t13;
	t16 = sin(qJ(2)) * pkin(1);
	t14 = -pkin(12) - t22;
	t1 = [t13 * r_i_i_C(3) + t14 * t12, -t13 * t16 + t17, t17, 0; t12 * r_i_i_C(3) - t14 * t13, -t12 * t16 + t18, t18, 0; 0, t22, t15, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (2493->21), mult. (4337->30), div. (72->0), fcn. (6509->17), ass. (0->29)
	t35 = sin(qJ(2));
	t39 = cos(qJ(2));
	t32 = sin(pkin(16));
	t33 = cos(pkin(16));
	t37 = sin(pkin(15));
	t41 = cos(pkin(15));
	t20 = t32 * t41 + t33 * t37;
	t21 = -t32 * t37 + t33 * t41;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t45 = t34 * t20 - t21 * t38;
	t49 = -t20 * t38 - t34 * t21;
	t59 = -t35 * t49 + t45 * t39;
	t31 = qJ(2) + qJ(3);
	t54 = pkin(4) * cos(t31);
	t58 = -t39 * pkin(1) + t54;
	t46 = t35 * t45 + t39 * t49;
	t55 = pkin(4) * sin(t31);
	t30 = pkin(17) + pkin(18);
	t25 = sin(t30);
	t26 = cos(t30);
	t5 = atan2(t25 * t59 + t26 * t46, t46 * t25 - t59 * t26) + t31;
	t43 = -pkin(12) + r_i_i_C(1) * cos(t5) - r_i_i_C(2) * sin(t5) + t58;
	t42 = pkin(1) * t35;
	t40 = cos(qJ(1));
	t36 = sin(qJ(1));
	t23 = t40 * t55;
	t22 = t36 * t55;
	t1 = [t40 * r_i_i_C(3) + t43 * t36, -t42 * t40 + t23, t23, 0; t36 * r_i_i_C(3) - t43 * t40, -t42 * t36 + t22, t22, 0; 0, -t58, -t54, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:29
	% EndTime: 2020-05-07 02:13:30
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (5088->68), mult. (7105->89), div. (208->2), fcn. (10589->35), ass. (0->70)
	t97 = qJ(2) + qJ(3);
	t65 = qJ(1) - t97;
	t109 = pkin(4) * sin(t65);
	t66 = pkin(17) + pkin(18);
	t50 = pkin(16) + pkin(15) + t66 + t97;
	t48 = sin(t50);
	t49 = cos(t50);
	t41 = atan2(-t48, t49);
	t64 = qJ(1) + t97;
	t37 = t41 + t64;
	t31 = sin(t37);
	t38 = -t41 + t65;
	t34 = cos(t38);
	t51 = pkin(4) * sin(t64) / 0.2e1;
	t112 = -cos(t37) / 0.2e1;
	t104 = t48 ^ 2 / t49 ^ 2;
	t40 = 0.1e1 / (0.1e1 + t104);
	t26 = -t40 * t104 - t40 + 0.1e1;
	t93 = t26 * t112;
	t114 = -t26 / 0.2e1;
	t32 = sin(t38);
	t94 = t32 * t114;
	t95 = pkin(10) * t114;
	t96 = pkin(8) * t26 / 0.2e1;
	t123 = pkin(8) * t94 + pkin(10) * t93 + t31 * t96 + t34 * t95 - t109 / 0.2e1 + t51;
	t108 = pkin(4) * cos(t65);
	t53 = -pkin(4) * cos(t64) / 0.2e1;
	t122 = pkin(8) * t93 + pkin(10) * t94 + t31 * t95 + t34 * t96 + t108 / 0.2e1 + t53;
	t39 = t41 + t97;
	t121 = -pkin(4) * cos(t97) + (-pkin(8) * cos(t39) - pkin(10) * sin(t39)) * t26;
	t73 = sin(qJ(2));
	t78 = cos(qJ(2));
	t69 = sin(pkin(16));
	t70 = cos(pkin(16));
	t75 = sin(pkin(15));
	t80 = cos(pkin(15));
	t44 = t69 * t80 + t70 * t75;
	t45 = -t69 * t75 + t70 * t80;
	t72 = sin(qJ(3));
	t77 = cos(qJ(3));
	t88 = t72 * t44 - t45 * t77;
	t91 = -t44 * t77 - t72 * t45;
	t118 = -t73 * t91 + t88 * t78;
	t89 = t73 * t88 + t78 * t91;
	t113 = t32 / 0.2e1;
	t68 = qJ(1) - qJ(2);
	t111 = pkin(1) * sin(t68);
	t110 = pkin(1) * cos(t68);
	t71 = sin(qJ(4));
	t74 = sin(qJ(1));
	t101 = t74 * t71;
	t76 = cos(qJ(4));
	t100 = t74 * t76;
	t79 = cos(qJ(1));
	t99 = t79 * t71;
	t98 = t79 * t76;
	t59 = sin(t66);
	t60 = cos(t66);
	t9 = atan2(t118 * t59 + t60 * t89, -t118 * t60 + t89 * t59) + t97;
	t7 = sin(t9);
	t92 = t7 * r_i_i_C(3) - pkin(12);
	t67 = qJ(1) + qJ(2);
	t56 = pkin(1) * cos(t67) / 0.2e1;
	t55 = -pkin(1) * sin(t67) / 0.2e1;
	t8 = cos(t9);
	t6 = t8 * t98 - t101;
	t5 = t8 * t99 + t100;
	t4 = t8 * t100 + t99;
	t3 = t8 * t101 - t98;
	t1 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t111 / 0.2e1 + t55 + t109 / 0.2e1 + t51 + t92 * t74 + (t112 + t34 / 0.2e1) * pkin(10) + (t113 + t31 / 0.2e1) * pkin(8), t111 / 0.2e1 + t55 + t123, t123, t5 * r_i_i_C(1) + t6 * r_i_i_C(2); -t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t56 + t110 / 0.2e1 + t53 - t108 / 0.2e1 - t92 * t79 + (t113 - t31 / 0.2e1) * pkin(10) + (t112 - t34 / 0.2e1) * pkin(8), t56 - t110 / 0.2e1 + t122, t122, t3 * r_i_i_C(1) + t4 * r_i_i_C(2); 0, t78 * pkin(1) + t121, t121, (r_i_i_C(1) * t71 + r_i_i_C(2) * t76) * t7;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (39->9), mult. (82->18), div. (0->0), fcn. (112->8), ass. (0->16)
	t11 = cos(qJ(2));
	t10 = sin(pkin(14));
	t13 = cos(pkin(15));
	t14 = cos(pkin(14));
	t9 = sin(pkin(15));
	t5 = -t13 * t10 + t9 * t14;
	t6 = t9 * t10 + t13 * t14;
	t7 = sin(qJ(2));
	t1 = -t11 * t5 - t7 * t6;
	t2 = t6 * t11 - t7 * t5;
	t19 = t2 * r_i_i_C(1) + r_i_i_C(2) * t1;
	t16 = r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t15 = -pkin(6) + t19;
	t12 = cos(qJ(1));
	t8 = sin(qJ(1));
	t3 = [t12 * r_i_i_C(3) - t15 * t8, t16 * t12, 0, 0; t8 * r_i_i_C(3) + t15 * t12, t16 * t8, 0, 0; 0, t19, 0, 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (252->15), mult. (479->25), div. (36->2), fcn. (719->11), ass. (0->23)
	t14 = sin(pkin(18));
	t15 = cos(pkin(18));
	t18 = sin(pkin(15));
	t21 = cos(pkin(15));
	t10 = t21 * t14 + t18 * t15;
	t11 = -t18 * t14 + t21 * t15;
	t16 = sin(qJ(2));
	t19 = cos(qJ(2));
	t8 = t16 * t10 - t19 * t11;
	t9 = t19 * t10 + t16 * t11;
	t25 = t9 ^ 2 / t8 ^ 2;
	t4 = qJ(2) + atan2(t9, t8);
	t2 = sin(t4);
	t3 = cos(t4);
	t24 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t13 = t19 * pkin(1);
	t23 = t13 + pkin(12) + t24;
	t5 = 0.1e1 / (0.1e1 + t25);
	t1 = -t5 * t25 - t5 + 0.1e1;
	t22 = -pkin(1) * t16 + (-r_i_i_C(1) * t2 - r_i_i_C(2) * t3) * t1;
	t20 = cos(qJ(1));
	t17 = sin(qJ(1));
	t6 = [t20 * r_i_i_C(3) - t23 * t17, t22 * t20, 0, 0; t17 * r_i_i_C(3) + t23 * t20, t22 * t17, 0, 0; 0, t24 * t1 + t13, 0, 0;];
	Ja_transl = t6;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:53
	% EndTime: 2020-05-07 02:13:54
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (1760->30), mult. (2833->43), div. (126->2), fcn. (4305->18), ass. (0->36)
	t44 = sin(pkin(15));
	t45 = cos(qJ(3));
	t48 = cos(pkin(15));
	t61 = sin(qJ(3));
	t32 = t45 * t44 + t61 * t48;
	t33 = -t61 * t44 + t45 * t48;
	t42 = sin(qJ(2));
	t46 = cos(qJ(2));
	t24 = -t42 * t32 + t46 * t33;
	t64 = -t32 * t46 - t42 * t33;
	t40 = sin(pkin(18));
	t41 = cos(pkin(18));
	t31 = t48 * t40 + t44 * t41;
	t52 = -t44 * t40 + t48 * t41;
	t22 = t42 * t31 - t46 * t52;
	t23 = t46 * t31 + t42 * t52;
	t15 = pkin(17) - atan2(t23, t22) - qJ(2);
	t62 = pkin(3) * cos(t15);
	t59 = t23 ^ 2 / t22 ^ 2;
	t39 = pkin(17) + pkin(18);
	t36 = sin(t39);
	t37 = cos(t39);
	t6 = -atan2(-t36 * t24 + t64 * t37, -t24 * t37 - t36 * t64) + t15;
	t4 = sin(t6);
	t5 = cos(t6);
	t54 = r_i_i_C(1) * t5 + r_i_i_C(2) * t4;
	t53 = r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t16 = 0.1e1 / (0.1e1 + t59);
	t3 = t16 * t59 + t16 - 0.1e1;
	t38 = t46 * pkin(1);
	t50 = t38 + pkin(12) - t54 + t62;
	t1 = t3 - 0.1e1;
	t49 = -pkin(3) * sin(t15) * t3 - pkin(1) * t42 + t53 * t1;
	t47 = cos(qJ(1));
	t43 = sin(qJ(1));
	t2 = [t47 * r_i_i_C(3) - t50 * t43, t49 * t47, -t47 * t53, 0; t43 * r_i_i_C(3) + t50 * t47, t49 * t43, -t43 * t53, 0; 0, t54 * t1 - t3 * t62 + t38, -t54, 0;];
	Ja_transl = t2;
end