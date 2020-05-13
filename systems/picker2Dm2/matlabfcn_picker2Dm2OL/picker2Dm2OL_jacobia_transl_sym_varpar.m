% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% picker2Dm2OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% Ja_transl [3x12]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = picker2Dm2OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(12,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_jacobia_transl_sym_varpar: qJ has to be [12x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'picker2Dm2OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm2OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
Ja_transl=NaN(3,12);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t6 = qJ(1) + qJ(2);
	t4 = sin(t6);
	t5 = cos(t6);
	t8 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t7 = -t5 * r_i_i_C(1) + t4 * r_i_i_C(2);
	t1 = [sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t7, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t10 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t9 = -pkin(2) * cos(t7) + t11;
	t8 = t10 + pkin(2) * sin(t7);
	t1 = [sin(qJ(1)) * pkin(1) + t8, t8, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t9, t9, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t9 = qJ(1) + qJ(2);
	t8 = qJ(4) + t9;
	t5 = sin(t8);
	t6 = cos(t8);
	t13 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t12 = pkin(3) * sin(t9) + t13;
	t11 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t10 = -pkin(3) * cos(t9) + t11;
	t1 = [sin(qJ(1)) * pkin(1) + t12, t12, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t10, t10, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(8) + qJ(5);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 0, 0, 0, -t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (32->5), mult. (14->6), div. (0->0), fcn. (14->4), ass. (0->6)
	t6 = qJ(1) + qJ(2) + qJ(6);
	t4 = sin(t6);
	t5 = cos(t6);
	t2 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t3 = [sin(qJ(1)) * pkin(1) + t1, t1, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t2, t2, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(7));
	t1 = sin(qJ(7));
	t3 = [0, 0, 0, 0, 0, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(8);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [sin(qJ(1)) * pkin(1) + t5, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t6, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (84->11), mult. (28->10), div. (0->0), fcn. (28->8), ass. (0->12)
	t12 = qJ(1) + qJ(2);
	t11 = qJ(3) + t12;
	t9 = qJ(9) + t11;
	t5 = sin(t9);
	t6 = cos(t9);
	t18 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t17 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t16 = -pkin(6) * sin(t11) + t18;
	t15 = t17 + pkin(6) * cos(t11);
	t14 = t16 + pkin(2) * sin(t12);
	t13 = -pkin(2) * cos(t12) + t15;
	t1 = [sin(qJ(1)) * pkin(1) + t14, t14, t16, 0, 0, 0, 0, 0, t18, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t13, t13, t15, 0, 0, 0, 0, 0, t17, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (84->11), mult. (28->10), div. (0->0), fcn. (28->8), ass. (0->12)
	t10 = qJ(1) + qJ(2);
	t9 = qJ(4) + t10;
	t7 = qJ(10) + t9;
	t3 = sin(t7);
	t4 = cos(t7);
	t16 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t15 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t14 = -pkin(4) * cos(t9) + t16;
	t13 = t15 + pkin(4) * sin(t9);
	t12 = t13 + pkin(3) * sin(t10);
	t11 = -pkin(3) * cos(t10) + t14;
	t1 = [sin(qJ(1)) * pkin(1) + t12, t12, 0, t13, 0, 0, 0, 0, 0, t15, 0, 0; -cos(qJ(1)) * pkin(1) + t11, t11, 0, t14, 0, 0, 0, 0, 0, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
end