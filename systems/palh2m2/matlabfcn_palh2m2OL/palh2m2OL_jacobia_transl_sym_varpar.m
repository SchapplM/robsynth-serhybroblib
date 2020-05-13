% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m2OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh2m2OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m2OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
Ja_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, t6 * t4, 0, 0, 0, 0; t2 * r_i_i_C(3) + t5 * t4, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t13 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t15 = t13 + cos(qJ(2)) * pkin(4);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(1) + t15;
	t10 = -sin(qJ(2)) * pkin(4) + t12;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t9 * r_i_i_C(3) - t11 * t8, t10 * t9, t12 * t9, 0, 0, 0; t8 * r_i_i_C(3) + t11 * t9, t10 * t8, t12 * t8, 0, 0, 0; 0, t15, t13, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (80->11), mult. (59->18), div. (0->0), fcn. (59->8), ass. (0->14)
	t10 = qJ(2) + qJ(3);
	t8 = qJ(4) + t10;
	t5 = sin(t8);
	t6 = cos(t8);
	t18 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t17 = t18 + pkin(2) * cos(t10);
	t21 = t17 + cos(qJ(2)) * pkin(4);
	t16 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t13 = t16 - pkin(2) * sin(t10);
	t15 = pkin(1) + t21;
	t14 = -sin(qJ(2)) * pkin(4) + t13;
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [t12 * r_i_i_C(3) - t15 * t11, t14 * t12, t13 * t12, t16 * t12, 0, 0; t11 * r_i_i_C(3) + t15 * t12, t14 * t11, t13 * t11, t16 * t11, 0, 0; 0, t21, t17, t18, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (156->14), mult. (82->22), div. (0->0), fcn. (82->10), ass. (0->17)
	t14 = qJ(2) + qJ(3);
	t12 = qJ(4) + t14;
	t10 = qJ(5) + t12;
	t6 = sin(t10);
	t7 = cos(t10);
	t24 = t7 * r_i_i_C(1) - r_i_i_C(2) * t6;
	t23 = t24 + pkin(5) * cos(t12);
	t21 = t23 + pkin(2) * cos(t14);
	t27 = cos(qJ(2)) * pkin(4) + t21;
	t22 = -r_i_i_C(1) * t6 - r_i_i_C(2) * t7;
	t17 = t22 - pkin(5) * sin(t12);
	t18 = -pkin(2) * sin(t14) + t17;
	t20 = pkin(1) + t27;
	t19 = -sin(qJ(2)) * pkin(4) + t18;
	t16 = cos(qJ(1));
	t15 = sin(qJ(1));
	t1 = [r_i_i_C(3) * t16 - t20 * t15, t19 * t16, t18 * t16, t17 * t16, t22 * t16, 0; t15 * r_i_i_C(3) + t20 * t16, t19 * t15, t18 * t15, t17 * t15, t22 * t15, 0; 0, t27, t21, t23, t24, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (295->34), mult. (180->45), div. (0->0), fcn. (188->12), ass. (0->35)
	t21 = qJ(2) + qJ(3);
	t19 = qJ(4) + t21;
	t17 = qJ(5) + t19;
	t13 = sin(t17);
	t14 = cos(t17);
	t24 = cos(qJ(6));
	t41 = r_i_i_C(1) * t24;
	t31 = -r_i_i_C(3) * t14 + (-pkin(3) - t41) * t13;
	t26 = t31 - pkin(5) * sin(t19);
	t28 = -pkin(2) * sin(t21) + t26;
	t43 = t14 * pkin(3) - t13 * r_i_i_C(3);
	t22 = sin(qJ(6));
	t40 = r_i_i_C(2) * t22;
	t23 = sin(qJ(1));
	t38 = t23 * t22;
	t37 = t23 * t24;
	t25 = cos(qJ(1));
	t36 = t25 * t22;
	t35 = t25 * t24;
	t34 = t13 * t40;
	t12 = pkin(5) * cos(t19);
	t15 = pkin(2) * cos(t21);
	t20 = cos(qJ(2)) * pkin(4);
	t33 = t15 + t12 + t20 + pkin(1) + t43;
	t32 = t43 + (-t40 + t41) * t14;
	t30 = t12 + t32;
	t29 = -sin(qJ(2)) * pkin(4) + t28;
	t27 = t15 + t30;
	t9 = t25 * t34;
	t8 = t23 * t34;
	t5 = t14 * t35 - t38;
	t4 = -t14 * t36 - t37;
	t3 = -t14 * t37 - t36;
	t2 = t14 * t38 - t35;
	t1 = [t3 * r_i_i_C(1) + t2 * r_i_i_C(2) - t33 * t23, t29 * t25 + t9, t28 * t25 + t9, t26 * t25 + t9, t31 * t25 + t9, t4 * r_i_i_C(1) - t5 * r_i_i_C(2); t5 * r_i_i_C(1) + t4 * r_i_i_C(2) + t33 * t25, t29 * t23 + t8, t28 * t23 + t8, t26 * t23 + t8, t31 * t23 + t8, -t2 * r_i_i_C(1) + t3 * r_i_i_C(2); 0, t20 + t27, t27, t30, t32, (-r_i_i_C(1) * t22 - r_i_i_C(2) * t24) * t13;];
	Ja_transl = t1;
end