% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m2DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m2DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_jacobia_transl_sym_varpar: pkin has to be [22x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
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
	% StartTime: 2020-05-01 21:04:21
	% EndTime: 2020-05-01 21:04:21
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
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->12), mult. (71->20), div. (0->0), fcn. (87->6), ass. (0->17)
	t11 = cos(qJ(3));
	t12 = cos(qJ(2));
	t8 = sin(qJ(3));
	t9 = sin(qJ(2));
	t5 = t11 * t12 - t9 * t8;
	t6 = t11 * t9 + t12 * t8;
	t16 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t20 = -pkin(1) * t9 + t16;
	t18 = r_i_i_C(1) * t6;
	t17 = r_i_i_C(2) * t5;
	t15 = -pkin(1) * t12 - t17;
	t14 = pkin(15) + t20;
	t13 = cos(qJ(1));
	t10 = sin(qJ(1));
	t2 = t13 * t18;
	t1 = t10 * t18;
	t3 = [t13 * r_i_i_C(3) - t14 * t10, t15 * t13 - t2, -t13 * t17 - t2, 0; t10 * r_i_i_C(3) + t14 * t13, t15 * t10 - t1, -t10 * t17 - t1, 0; 0, t20, t16, 0;];
	Ja_transl = t3;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (41->15), mult. (72->27), div. (0->0), fcn. (84->12), ass. (0->24)
	t15 = sin(qJ(3));
	t26 = t15 * pkin(5);
	t16 = sin(qJ(2));
	t19 = cos(qJ(3));
	t25 = t16 * t19;
	t20 = cos(qJ(2));
	t8 = pkin(5) * t20 * t19;
	t9 = pkin(1) + t26;
	t24 = -t9 * t16 + t8;
	t12 = pkin(22) + pkin(21);
	t10 = sin(t12);
	t11 = cos(t12);
	t13 = sin(pkin(20));
	t14 = cos(pkin(20));
	t18 = sin(pkin(18));
	t22 = cos(pkin(18));
	t6 = -t22 * t13 + t18 * t14;
	t7 = t18 * t13 + t22 * t14;
	t23 = r_i_i_C(1) * (t6 * t10 + t11 * t7) + r_i_i_C(2) * (t10 * t7 - t6 * t11) - pkin(15) - t24;
	t21 = cos(qJ(1));
	t17 = sin(qJ(1));
	t5 = (-t15 * t20 - t25) * pkin(5);
	t4 = -pkin(5) * t25 - t9 * t20;
	t1 = [t21 * r_i_i_C(3) + t23 * t17, t4 * t21, t5 * t21, 0; t17 * r_i_i_C(3) - t23 * t21, t4 * t17, t5 * t17, 0; 0, t24, -t16 * t26 + t8, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:21
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (112->29), mult. (198->55), div. (0->0), fcn. (242->14), ass. (0->36)
	t23 = sin(qJ(3));
	t38 = t23 * pkin(5);
	t24 = sin(qJ(2));
	t28 = cos(qJ(3));
	t37 = t24 * t28;
	t29 = cos(qJ(2));
	t15 = pkin(5) * t29 * t28;
	t16 = pkin(1) + t38;
	t36 = -t16 * t24 + t15;
	t20 = sin(pkin(20));
	t21 = cos(pkin(20));
	t26 = sin(pkin(18));
	t31 = cos(pkin(18));
	t11 = -t31 * t20 + t26 * t21;
	t12 = t26 * t20 + t31 * t21;
	t19 = pkin(22) + pkin(21);
	t17 = sin(t19);
	t18 = cos(t19);
	t35 = t11 * t17 + t12 * t18;
	t25 = sin(qJ(1));
	t34 = t35 * t25;
	t30 = cos(qJ(1));
	t33 = t35 * t30;
	t13 = -t26 * pkin(9) + t31 * pkin(11);
	t14 = t31 * pkin(9) + t26 * pkin(11);
	t5 = -t11 * t18 + t17 * t12;
	t32 = r_i_i_C(3) * t5 + t17 * (t13 * t21 + t20 * t14) - t18 * (-t20 * t13 + t14 * t21) + pkin(15) + t36;
	t27 = cos(qJ(4));
	t22 = sin(qJ(4));
	t10 = (-t23 * t29 - t37) * pkin(5);
	t9 = -pkin(5) * t37 - t16 * t29;
	t4 = -t25 * t22 + t27 * t33;
	t3 = t22 * t33 + t25 * t27;
	t2 = t30 * t22 + t27 * t34;
	t1 = t22 * t34 - t30 * t27;
	t6 = [t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t32 * t25, t9 * t30, t10 * t30, t3 * r_i_i_C(1) + t4 * r_i_i_C(2); -t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t32 * t30, t9 * t25, t10 * t25, t1 * r_i_i_C(1) + t2 * r_i_i_C(2); 0, t36, -t24 * t38 + t15, (-r_i_i_C(1) * t22 - r_i_i_C(2) * t27) * t5;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (39->9), mult. (82->18), div. (0->0), fcn. (112->8), ass. (0->16)
	t12 = cos(qJ(2));
	t11 = sin(pkin(18));
	t14 = cos(pkin(17));
	t18 = sin(pkin(17));
	t19 = cos(pkin(18));
	t6 = t11 * t14 - t19 * t18;
	t7 = t18 * t11 + t14 * t19;
	t9 = sin(qJ(2));
	t23 = t6 * t12 - t7 * t9;
	t16 = -t7 * t12 - t9 * t6;
	t22 = r_i_i_C(1) * t23 + r_i_i_C(2) * t16;
	t17 = r_i_i_C(1) * t16 - r_i_i_C(2) * t23;
	t15 = -pkin(14) + t22;
	t13 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [t13 * r_i_i_C(3) - t15 * t10, t17 * t13, 0, 0; t10 * r_i_i_C(3) + t15 * t13, t17 * t10, 0, 0; 0, t22, 0, 0;];
	Ja_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->10), mult. (27->14), div. (0->0), fcn. (31->8), ass. (0->10)
	t14 = pkin(1) * sin(qJ(2));
	t13 = cos(qJ(2)) * pkin(1);
	t11 = cos(pkin(18));
	t4 = sin(pkin(22));
	t5 = cos(pkin(22));
	t8 = sin(pkin(18));
	t12 = r_i_i_C(1) * (t5 * t11 + t8 * t4) + r_i_i_C(2) * (t11 * t4 - t8 * t5) - pkin(15) + t14;
	t10 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t10 * r_i_i_C(3) + t12 * t7, -t10 * t13, 0, 0; t7 * r_i_i_C(3) - t12 * t10, -t7 * t13, 0, 0; 0, -t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (60->16), mult. (128->28), div. (0->0), fcn. (176->8), ass. (0->22)
	t15 = sin(pkin(19));
	t16 = cos(pkin(19));
	t17 = sin(qJ(3));
	t20 = cos(qJ(3));
	t10 = t20 * t15 + t17 * t16;
	t12 = t17 * t15 - t20 * t16;
	t18 = sin(qJ(2));
	t19 = sin(qJ(1));
	t21 = cos(qJ(2));
	t11 = t15 * t18 - t16 * t21;
	t13 = t15 * t21 + t16 * t18;
	t26 = -t11 * t17 + t13 * t20;
	t34 = t26 * r_i_i_C(1);
	t39 = ((t18 * t10 + t21 * t12) * r_i_i_C(2) - t34) * t19;
	t22 = cos(qJ(1));
	t25 = -t11 * t20 - t17 * t13;
	t38 = (-r_i_i_C(2) * t25 - t34) * t22;
	t37 = r_i_i_C(1) * t25;
	t35 = r_i_i_C(2) * t26;
	t33 = -t35 + t37;
	t27 = pkin(15) + t37;
	t1 = [t22 * r_i_i_C(3) + (-t27 + t35) * t19, t38, t38, 0; t19 * r_i_i_C(3) + ((-t10 * t21 + t18 * t12) * r_i_i_C(2) + t27) * t22, t39, t39, 0; 0, t33, t33, 0;];
	Ja_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (38->11), mult. (106->22), div. (0->0), fcn. (114->8), ass. (0->16)
	t10 = sin(qJ(2));
	t13 = cos(qJ(2));
	t12 = cos(qJ(3));
	t7 = sin(pkin(19));
	t8 = cos(pkin(19));
	t9 = sin(qJ(3));
	t4 = t12 * t7 + t9 * t8;
	t5 = t12 * t8 - t9 * t7;
	t15 = (-t10 * t4 + t13 * t5) * pkin(2);
	t23 = -t10 * r_i_i_C(1) - t13 * r_i_i_C(2) + t15;
	t2 = (-t10 * t5 - t13 * t4) * pkin(2);
	t17 = -pkin(15) - t23;
	t16 = -r_i_i_C(1) * t13 + r_i_i_C(2) * t10 + t2;
	t14 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [t14 * r_i_i_C(3) + t17 * t11, t16 * t14, t2 * t14, 0; t11 * r_i_i_C(3) - t17 * t14, t16 * t11, t2 * t11, 0; 0, t23, t15, 0;];
	Ja_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 21:04:22
	% EndTime: 2020-05-01 21:04:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (41->12), mult. (59->21), div. (0->0), fcn. (61->12), ass. (0->15)
	t5 = qJ(3) + qJ(2);
	t3 = sin(t5);
	t4 = cos(t5);
	t18 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t21 = t18 - pkin(1) * sin(qJ(2));
	t17 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = sin(pkin(18));
	t14 = cos(pkin(18));
	t6 = sin(pkin(22));
	t8 = cos(pkin(22));
	t16 = pkin(15) + (-(t11 * t6 + t8 * t14) * cos(pkin(21)) + (-t11 * t8 + t14 * t6) * sin(pkin(21))) * pkin(4) + t21;
	t15 = -cos(qJ(2)) * pkin(1) + t17;
	t13 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [t13 * r_i_i_C(3) - t16 * t10, t15 * t13, t17 * t13, 0; t10 * r_i_i_C(3) + t16 * t13, t15 * t10, t17 * t10, 0; 0, t21, t18, 0;];
	Ja_transl = t1;
end