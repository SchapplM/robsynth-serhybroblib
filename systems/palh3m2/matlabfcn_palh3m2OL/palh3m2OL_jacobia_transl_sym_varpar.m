% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m2OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% Ja_transl [3x10]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m2OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(10,1),uint8(0),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_jacobia_transl_sym_varpar: qJ has to be [10x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m2OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_jacobia_transl_sym_varpar: pkin has to be [16x1] (double)');
Ja_transl=NaN(3,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(12) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, t6 * t4, 0, 0, 0, 0, 0, 0, 0, 0; t2 * r_i_i_C(3) + t5 * t4, t6 * t2, 0, 0, 0, 0, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:50
	% DurationCPUTime: 0.08s
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
	t1 = [t13 * r_i_i_C(3) + t12 * t14, -t13 * t16 + t17, t17, 0, 0, 0, 0, 0, 0, 0; t12 * r_i_i_C(3) - t13 * t14, -t12 * t16 + t18, t18, 0, 0, 0, 0, 0, 0, 0; 0, t22, t15, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (80->14), mult. (59->18), div. (0->0), fcn. (59->8), ass. (0->16)
	t14 = qJ(2) + qJ(3);
	t12 = qJ(4) + t14;
	t8 = sin(t12);
	t9 = cos(t12);
	t19 = -t9 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t17 = t19 - pkin(4) * cos(t14);
	t28 = cos(qJ(2)) * pkin(1) + t17;
	t27 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9;
	t15 = sin(qJ(1));
	t23 = t27 * t15;
	t16 = cos(qJ(1));
	t22 = t27 * t16;
	t21 = pkin(4) * sin(t14);
	t18 = -pkin(12) - t28;
	t2 = t21 - sin(qJ(2)) * pkin(1);
	t1 = [t16 * r_i_i_C(3) + t18 * t15, t16 * t2 + t22, t16 * t21 + t22, t22, 0, 0, 0, 0, 0, 0; t15 * r_i_i_C(3) - t18 * t16, t15 * t2 + t23, t15 * t21 + t23, t23, 0, 0, 0, 0, 0, 0; 0, t28, t17, t19, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (200->32), mult. (160->41), div. (0->0), fcn. (168->10), ass. (0->32)
	t21 = cos(qJ(5));
	t45 = r_i_i_C(1) * t21 + pkin(8);
	t18 = qJ(2) + qJ(3);
	t16 = qJ(4) + t18;
	t13 = cos(t16);
	t17 = cos(qJ(2)) * pkin(1);
	t12 = sin(t16);
	t42 = pkin(10) + r_i_i_C(3);
	t29 = t42 * t12;
	t39 = pkin(4) * cos(t18);
	t44 = t13 * pkin(8) - pkin(12) - t17 + t29 + t39;
	t43 = t45 * t12;
	t19 = sin(qJ(5));
	t36 = r_i_i_C(2) * t19;
	t27 = -t12 * t36 - t42 * t13;
	t25 = t27 + pkin(4) * sin(t18);
	t20 = sin(qJ(1));
	t41 = t43 * t20;
	t22 = cos(qJ(1));
	t35 = t43 * t22;
	t34 = t20 * t19;
	t33 = t20 * t21;
	t32 = t22 * t19;
	t31 = t22 * t21;
	t26 = -sin(qJ(2)) * pkin(1) + t25;
	t24 = -t29 + (t36 - t45) * t13;
	t23 = t24 - t39;
	t4 = t13 * t31 - t34;
	t3 = t13 * t32 + t33;
	t2 = t13 * t33 + t32;
	t1 = t13 * t34 - t31;
	t5 = [t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t44 * t20, t26 * t22 + t35, t25 * t22 + t35, t27 * t22 + t35, t3 * r_i_i_C(1) + t4 * r_i_i_C(2), 0, 0, 0, 0, 0; -t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t44 * t22, t26 * t20 + t41, t25 * t20 + t41, t27 * t20 + t41, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, t17 + t23, t23, t24, (r_i_i_C(1) * t19 + r_i_i_C(2) * t21) * t12, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(6));
	t3 = cos(qJ(6));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = -pkin(6) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, 0, 0, 0, 0, t6 * t4, 0, 0, 0, 0; t2 * r_i_i_C(3) + t5 * t4, 0, 0, 0, 0, t6 * t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, t7, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t6 = qJ(2) + qJ(7);
	t3 = sin(t6);
	t4 = cos(t6);
	t13 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t15 = t13 + cos(qJ(2)) * pkin(1);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(12) + t15;
	t10 = -sin(qJ(2)) * pkin(1) + t12;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t9 * r_i_i_C(3) - t11 * t8, t10 * t9, 0, 0, 0, 0, t12 * t9, 0, 0, 0; t8 * r_i_i_C(3) + t11 * t9, t10 * t8, 0, 0, 0, 0, t12 * t8, 0, 0, 0; 0, t15, 0, 0, 0, 0, t13, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (110->17), mult. (59->20), div. (0->0), fcn. (59->8), ass. (0->17)
	t10 = pkin(15) - qJ(7) - qJ(2);
	t9 = -qJ(8) + t10;
	t6 = sin(t9);
	t7 = cos(t9);
	t16 = r_i_i_C(1) * t7 + t6 * r_i_i_C(2);
	t14 = -t16 + pkin(3) * cos(t10);
	t22 = cos(qJ(2)) * pkin(1) + t14;
	t20 = r_i_i_C(1) * t6;
	t17 = -t20 + pkin(3) * sin(t10);
	t19 = r_i_i_C(2) * t7;
	t18 = -sin(qJ(2)) * pkin(1) + t17;
	t15 = -pkin(12) - t22;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t4 = t13 * t19;
	t3 = t12 * t19;
	t1 = [t13 * r_i_i_C(3) + t15 * t12, t18 * t13 + t4, 0, 0, 0, 0, t17 * t13 + t4, -t13 * t20 + t4, 0, 0; t12 * r_i_i_C(3) - t15 * t13, t18 * t12 + t3, 0, 0, 0, 0, t17 * t12 + t3, -t12 * t20 + t3, 0, 0; 0, t22, 0, 0, 0, 0, t14, -t16, 0, 0;];
	Ja_transl = t1;
end