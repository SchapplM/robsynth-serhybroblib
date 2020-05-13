% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_transl [3x13]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m1OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(3,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_jacobia_transl_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m1OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_jacobia_transl_sym_varpar: pkin has to be [20x1] (double)');
Ja_transl=NaN(3,13);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = -pkin(15) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) + t5 * t2, t7 * t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t2 * r_i_i_C(3) - t5 * t4, t7 * t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.08s
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
	t1 = [t9 * r_i_i_C(3) - t11 * t7, t10 * t9, t12 * t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t7 * r_i_i_C(3) + t11 * t9, t10 * t7, t12 * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t16, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (80->11), mult. (59->18), div. (0->0), fcn. (59->8), ass. (0->14)
	t9 = qJ(2) + qJ(3);
	t8 = qJ(4) + t9;
	t5 = sin(t8);
	t6 = cos(t8);
	t18 = t6 * r_i_i_C(1) - r_i_i_C(2) * t5;
	t17 = t18 + pkin(5) * cos(t9);
	t22 = t17 - sin(qJ(2)) * pkin(1);
	t16 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t13 = t16 - pkin(5) * sin(t9);
	t15 = pkin(15) + t22;
	t14 = -cos(qJ(2)) * pkin(1) + t13;
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [r_i_i_C(3) * t12 - t15 * t11, t14 * t12, t13 * t12, t16 * t12, 0, 0, 0, 0, 0, 0, 0, 0, 0; t11 * r_i_i_C(3) + t15 * t12, t14 * t11, t13 * t11, t16 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t22, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (200->30), mult. (160->41), div. (0->0), fcn. (168->10), ass. (0->32)
	t22 = qJ(2) + qJ(3);
	t21 = qJ(4) + t22;
	t18 = sin(t21);
	t19 = cos(t21);
	t23 = sin(qJ(5));
	t44 = r_i_i_C(2) * t23;
	t51 = pkin(11) + r_i_i_C(3);
	t52 = t18 * t44 + t19 * t51;
	t49 = t19 * pkin(9) + t51 * t18;
	t26 = cos(qJ(5));
	t45 = r_i_i_C(1) * t26;
	t32 = (-pkin(9) - t45) * t18;
	t29 = t32 - pkin(5) * sin(t22);
	t17 = pkin(5) * cos(t22);
	t43 = sin(qJ(2)) * pkin(1);
	t48 = pkin(15) + t17 - t43 + t49;
	t25 = sin(qJ(1));
	t40 = t25 * t23;
	t39 = t25 * t26;
	t27 = cos(qJ(1));
	t38 = t27 * t23;
	t37 = t27 * t26;
	t36 = t52 * t25;
	t34 = t52 * t27;
	t31 = -cos(qJ(2)) * pkin(1) + t29;
	t30 = (-t44 + t45) * t19 + t49;
	t28 = t17 + t30;
	t4 = t19 * t37 + t40;
	t3 = -t19 * t38 + t39;
	t2 = -t19 * t39 + t38;
	t1 = t19 * t40 + t37;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t48 * t25, t31 * t27 + t34, t29 * t27 + t34, t27 * t32 + t34, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t48 * t27, t31 * t25 + t36, t29 * t25 + t36, t25 * t32 + t36, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0; 0, t28 - t43, t28, t30, (-r_i_i_C(1) * t23 - r_i_i_C(2) * t26) * t18, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(6));
	t3 = cos(qJ(6));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = -pkin(14) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, 0, 0, 0, 0, t6 * t4, 0, 0, 0, 0, 0, 0, 0; t2 * r_i_i_C(3) + t5 * t4, 0, 0, 0, 0, t6 * t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->12), mult. (39->16), div. (0->0), fcn. (39->6), ass. (0->14)
	t6 = qJ(2) + qJ(7);
	t4 = sin(t6);
	t5 = cos(t6);
	t12 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t17 = -t12 - sin(qJ(2)) * pkin(1);
	t16 = r_i_i_C(1) * t5;
	t15 = r_i_i_C(2) * t4;
	t13 = -cos(qJ(2)) * pkin(1) - t16;
	t11 = -pkin(15) - t17;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t2 = t10 * t15;
	t1 = t8 * t15;
	t3 = [t10 * r_i_i_C(3) + t11 * t8, t13 * t10 + t2, 0, 0, 0, 0, -t10 * t16 + t2, 0, 0, 0, 0, 0, 0; t8 * r_i_i_C(3) - t11 * t10, t13 * t8 + t1, 0, 0, 0, 0, -t8 * t16 + t1, 0, 0, 0, 0, 0, 0; 0, t17, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (28->6), mult. (32->10), div. (0->0), fcn. (32->4), ass. (0->11)
	t6 = qJ(2) + qJ(8);
	t4 = sin(t6);
	t5 = cos(t6);
	t10 = -r_i_i_C(1) * t5 + r_i_i_C(2) * t4;
	t3 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t9 = -pkin(15) - t3;
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t2 = t10 * t8;
	t1 = t10 * t7;
	t11 = [t8 * r_i_i_C(3) + t9 * t7, t2, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0; t7 * r_i_i_C(3) - t9 * t8, t1, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0; 0, t3, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0;];
	Ja_transl = t11;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (75->12), mult. (54->16), div. (0->0), fcn. (54->6), ass. (0->17)
	t14 = qJ(2) + qJ(8);
	t13 = qJ(9) + t14;
	t10 = cos(t13);
	t9 = sin(t13);
	t21 = t9 * r_i_i_C(1) + t10 * r_i_i_C(2);
	t3 = -pkin(2) * sin(t14) + t21;
	t22 = r_i_i_C(2) * t9;
	t19 = r_i_i_C(1) * t10;
	t18 = -pkin(2) * cos(t14) - t22;
	t17 = pkin(15) + t3;
	t16 = cos(qJ(1));
	t15 = sin(qJ(1));
	t5 = t16 * t19;
	t4 = t15 * t19;
	t2 = t18 * t16 + t5;
	t1 = t18 * t15 + t4;
	t6 = [t16 * r_i_i_C(3) - t17 * t15, t2, 0, 0, 0, 0, 0, t2, -t16 * t22 + t5, 0, 0, 0, 0; t15 * r_i_i_C(3) + t17 * t16, t1, 0, 0, 0, 0, 0, t1, -t15 * t22 + t4, 0, 0, 0, 0; 0, t3, 0, 0, 0, 0, 0, t3, t21, 0, 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:17
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (110->14), mult. (59->18), div. (0->0), fcn. (59->8), ass. (0->16)
	t13 = pkin(19) - qJ(7) - qJ(2);
	t12 = -qJ(10) + t13;
	t10 = cos(t12);
	t9 = sin(t12);
	t19 = -r_i_i_C(1) * t9 + t10 * r_i_i_C(2);
	t18 = t19 + pkin(4) * sin(t13);
	t28 = t18 - sin(qJ(2)) * pkin(1);
	t27 = r_i_i_C(1) * t10 + r_i_i_C(2) * t9;
	t15 = sin(qJ(1));
	t24 = t27 * t15;
	t16 = cos(qJ(1));
	t23 = t27 * t16;
	t21 = pkin(4) * cos(t13);
	t17 = -pkin(15) - t28;
	t2 = -t21 - cos(qJ(2)) * pkin(1);
	t1 = [r_i_i_C(3) * t16 + t17 * t15, t16 * t2 + t23, 0, 0, 0, 0, -t16 * t21 + t23, 0, 0, t23, 0, 0, 0; t15 * r_i_i_C(3) - t17 * t16, t15 * t2 + t24, 0, 0, 0, 0, -t15 * t21 + t24, 0, 0, t24, 0, 0, 0; 0, t28, 0, 0, 0, 0, t18, 0, 0, t19, 0, 0, 0;];
	Ja_transl = t1;
end