% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh2m1OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m1OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) - t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [-t4 * r_i_i_C(3) - t2 * t5, t6 * t4, 0, 0, 0; -t2 * r_i_i_C(3) + t4 * t5, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t5 = qJ(2) + qJ(3);
	t3 = sin(t5);
	t4 = cos(t5);
	t13 = -t4 * r_i_i_C(1) + t3 * r_i_i_C(2);
	t16 = t13 - cos(qJ(2)) * pkin(2);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(1) - t16;
	t10 = -sin(qJ(2)) * pkin(2) + t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [-t9 * r_i_i_C(3) - t11 * t7, t10 * t9, t12 * t9, 0, 0; -t7 * r_i_i_C(3) + t11 * t9, t10 * t7, t12 * t7, 0, 0; 0, t16, t13, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (80->11), mult. (59->18), div. (0->0), fcn. (59->8), ass. (0->14)
	t9 = qJ(2) + qJ(3);
	t8 = qJ(4) + t9;
	t4 = sin(t8);
	t5 = cos(t8);
	t18 = -t5 * r_i_i_C(1) + t4 * r_i_i_C(2);
	t16 = t18 - pkin(3) * cos(t9);
	t23 = t16 - cos(qJ(2)) * pkin(2);
	t17 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t13 = t17 - pkin(3) * sin(t9);
	t15 = pkin(1) - t23;
	t14 = -sin(qJ(2)) * pkin(2) + t13;
	t12 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t12 * r_i_i_C(3) - t15 * t10, t14 * t12, t13 * t12, t17 * t12, 0; -t10 * r_i_i_C(3) + t15 * t12, t14 * t10, t13 * t10, t17 * t10, 0; 0, t23, t16, t18, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (200->32), mult. (160->41), div. (0->0), fcn. (168->10), ass. (0->32)
	t19 = qJ(2) + qJ(3);
	t18 = qJ(4) + t19;
	t14 = sin(t18);
	t15 = cos(t18);
	t20 = sin(qJ(5));
	t43 = r_i_i_C(2) * t20;
	t49 = pkin(6) + r_i_i_C(3);
	t50 = t14 * t43 + t15 * t49;
	t32 = t49 * t14;
	t42 = cos(qJ(2)) * pkin(2);
	t44 = pkin(3) * cos(t19);
	t48 = t15 * pkin(4) + pkin(1) + t32 + t42 + t44;
	t22 = cos(qJ(5));
	t31 = -r_i_i_C(1) * t22 - pkin(4);
	t29 = t31 * t14;
	t27 = t29 - pkin(3) * sin(t19);
	t21 = sin(qJ(1));
	t39 = t21 * t20;
	t38 = t21 * t22;
	t24 = cos(qJ(1));
	t37 = t24 * t20;
	t36 = t24 * t22;
	t35 = t50 * t21;
	t33 = t50 * t24;
	t28 = -sin(qJ(2)) * pkin(2) + t27;
	t26 = -t32 + (t31 + t43) * t15;
	t25 = t26 - t44;
	t4 = t15 * t36 - t39;
	t3 = -t15 * t37 - t38;
	t2 = -t15 * t38 - t37;
	t1 = t15 * t39 - t36;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t48 * t21, t28 * t24 + t33, t27 * t24 + t33, t24 * t29 + t33, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t48 * t24, t28 * t21 + t35, t27 * t21 + t35, t21 * t29 + t35, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t25 - t42, t25, t26, (r_i_i_C(1) * t20 + r_i_i_C(2) * t22) * t14;];
	Ja_transl = t5;
end