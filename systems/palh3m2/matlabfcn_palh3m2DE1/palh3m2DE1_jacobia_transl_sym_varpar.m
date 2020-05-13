% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m2DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m2DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_jacobia_transl_sym_varpar: pkin has to be [18x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
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
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
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
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->10), mult. (71->18), div. (0->0), fcn. (87->6), ass. (0->16)
	t15 = cos(qJ(2));
	t12 = sin(qJ(3));
	t13 = sin(qJ(2));
	t18 = cos(qJ(3));
	t7 = t12 * t13 - t18 * t15;
	t8 = t12 * t15 + t18 * t13;
	t20 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t26 = pkin(1) * t15 + t20;
	t25 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7;
	t14 = sin(qJ(1));
	t22 = t25 * t14;
	t16 = cos(qJ(1));
	t21 = t25 * t16;
	t19 = pkin(1) * t13;
	t17 = pkin(12) + t26;
	t1 = [t16 * r_i_i_C(3) - t17 * t14, -t16 * t19 + t21, t21, 0; t14 * r_i_i_C(3) + t17 * t16, -t14 * t19 + t22, t22, 0; 0, t26, t20, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (41->15), mult. (72->26), div. (0->0), fcn. (84->12), ass. (0->25)
	t21 = cos(qJ(3));
	t27 = pkin(4) * t21;
	t11 = pkin(1) - t27;
	t22 = cos(qJ(2));
	t17 = sin(qJ(3));
	t18 = sin(qJ(2));
	t9 = t18 * t17 * pkin(4);
	t28 = t11 * t22 + t9;
	t26 = pkin(4) * t22;
	t14 = pkin(17) + pkin(18);
	t12 = sin(t14);
	t13 = cos(t14);
	t15 = sin(pkin(16));
	t16 = cos(pkin(16));
	t20 = sin(pkin(15));
	t24 = cos(pkin(15));
	t6 = t15 * t24 + t16 * t20;
	t7 = -t15 * t20 + t16 * t24;
	t25 = r_i_i_C(1) * (t12 * t6 - t7 * t13) - r_i_i_C(2) * (t7 * t12 + t6 * t13) + pkin(12) + t28;
	t23 = cos(qJ(1));
	t19 = sin(qJ(1));
	t10 = t17 * t26;
	t5 = t18 * t27 + t10;
	t4 = -t11 * t18 + t10;
	t1 = [t23 * r_i_i_C(3) - t25 * t19, t4 * t23, t5 * t23, 0; t19 * r_i_i_C(3) + t25 * t23, t4 * t19, t5 * t19, 0; 0, t28, -t21 * t26 + t9, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:36
	% EndTime: 2020-05-07 01:57:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (112->29), mult. (198->54), div. (0->0), fcn. (242->14), ass. (0->37)
	t31 = cos(qJ(3));
	t41 = pkin(4) * t31;
	t32 = cos(qJ(2));
	t40 = pkin(4) * t32;
	t26 = sin(qJ(3));
	t27 = sin(qJ(2));
	t17 = t27 * t26 * pkin(4);
	t19 = pkin(1) - t41;
	t39 = t19 * t32 + t17;
	t23 = sin(pkin(16));
	t24 = cos(pkin(16));
	t29 = sin(pkin(15));
	t34 = cos(pkin(15));
	t12 = t23 * t34 + t24 * t29;
	t13 = -t23 * t29 + t24 * t34;
	t22 = pkin(17) + pkin(18);
	t20 = sin(t22);
	t21 = cos(t22);
	t5 = t12 * t21 + t13 * t20;
	t38 = -t12 * t20 + t13 * t21;
	t28 = sin(qJ(1));
	t37 = t38 * t28;
	t33 = cos(qJ(1));
	t36 = t38 * t33;
	t14 = t29 * pkin(8) + t34 * pkin(10);
	t15 = t34 * pkin(8) - t29 * pkin(10);
	t35 = r_i_i_C(3) * t5 + t20 * (t14 * t24 + t23 * t15) - t21 * (-t23 * t14 + t15 * t24) + pkin(12) + t39;
	t30 = cos(qJ(4));
	t25 = sin(qJ(4));
	t18 = t26 * t40;
	t11 = t27 * t41 + t18;
	t10 = -t19 * t27 + t18;
	t4 = -t28 * t25 + t30 * t36;
	t3 = t25 * t36 + t28 * t30;
	t2 = t33 * t25 + t30 * t37;
	t1 = t25 * t37 - t33 * t30;
	t6 = [t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t35 * t28, t10 * t33, t11 * t33, t3 * r_i_i_C(1) + t4 * r_i_i_C(2); -t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t35 * t33, t10 * t28, t11 * t28, t1 * r_i_i_C(1) + t2 * r_i_i_C(2); 0, t39, -t31 * t40 + t17, (-r_i_i_C(1) * t25 - r_i_i_C(2) * t30) * t5;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:37
	% EndTime: 2020-05-07 01:57:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (39->9), mult. (82->18), div. (0->0), fcn. (112->8), ass. (0->16)
	t12 = cos(qJ(2));
	t11 = sin(pkin(15));
	t14 = cos(pkin(14));
	t18 = sin(pkin(14));
	t19 = cos(pkin(15));
	t6 = t14 * t11 - t18 * t19;
	t7 = t18 * t11 + t14 * t19;
	t9 = sin(qJ(2));
	t23 = t7 * t12 - t9 * t6;
	t16 = -t12 * t6 - t9 * t7;
	t22 = r_i_i_C(1) * t23 + r_i_i_C(2) * t16;
	t17 = r_i_i_C(1) * t16 - r_i_i_C(2) * t23;
	t15 = -pkin(6) + t22;
	t13 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [t13 * r_i_i_C(3) - t15 * t10, t17 * t13, 0, 0; t10 * r_i_i_C(3) + t15 * t13, t17 * t10, 0, 0; 0, t22, 0, 0;];
	Ja_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:38
	% EndTime: 2020-05-07 01:57:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->9), mult. (27->14), div. (0->0), fcn. (31->8), ass. (0->10)
	t13 = sin(qJ(2)) * pkin(1);
	t11 = cos(pkin(15));
	t4 = pkin(1) * cos(qJ(2));
	t5 = sin(pkin(18));
	t6 = cos(pkin(18));
	t9 = sin(pkin(15));
	t12 = r_i_i_C(1) * (t6 * t11 - t5 * t9) + r_i_i_C(2) * (t5 * t11 + t6 * t9) - t4 - pkin(12);
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t10 * r_i_i_C(3) + t12 * t8, -t10 * t13, 0, 0; t8 * r_i_i_C(3) - t12 * t10, -t8 * t13, 0, 0; 0, t4, 0, 0;];
	Ja_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:38
	% EndTime: 2020-05-07 01:57:39
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (41->13), mult. (59->21), div. (0->0), fcn. (61->12), ass. (0->17)
	t10 = qJ(3) + qJ(2);
	t7 = sin(t10);
	t8 = cos(t10);
	t28 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t20 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t27 = t20 + pkin(1) * cos(qJ(2));
	t15 = sin(qJ(1));
	t23 = t28 * t15;
	t17 = cos(qJ(1));
	t22 = t28 * t17;
	t21 = sin(qJ(2)) * pkin(1);
	t11 = sin(pkin(18));
	t13 = cos(pkin(18));
	t16 = sin(pkin(15));
	t18 = cos(pkin(15));
	t19 = -pkin(12) - (-(-t11 * t16 + t13 * t18) * cos(pkin(17)) + (t18 * t11 + t16 * t13) * sin(pkin(17))) * pkin(3) - t27;
	t1 = [t17 * r_i_i_C(3) + t19 * t15, -t17 * t21 + t22, t22, 0; t15 * r_i_i_C(3) - t19 * t17, -t15 * t21 + t23, t23, 0; 0, t27, t20, 0;];
	Ja_transl = t1;
end