% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m1DE
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh2m1DE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m1DE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1DE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) - t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [-t4 * r_i_i_C(3) - t2 * t5, t6 * t4, 0, 0; -t2 * r_i_i_C(3) + t4 * t5, t6 * t2, 0, 0; 0, t7, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.08s
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
	t1 = [-t9 * r_i_i_C(3) - t11 * t7, t10 * t9, t12 * t9, 0; -t7 * r_i_i_C(3) + t11 * t9, t10 * t7, t12 * t7, 0; 0, t16, t13, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->9), mult. (23->12), div. (0->0), fcn. (23->6), ass. (0->9)
	t5 = qJ(2) + qJ(3);
	t12 = pkin(3) * sin(t5);
	t11 = pkin(3) * cos(t5);
	t9 = -cos(qJ(2)) * pkin(2) - t11;
	t10 = r_i_i_C(1) + pkin(1) - t9;
	t8 = cos(qJ(1));
	t6 = sin(qJ(1));
	t2 = -t12 - sin(qJ(2)) * pkin(2);
	t1 = [-t8 * r_i_i_C(3) - t10 * t6, t8 * t2, -t8 * t12, 0; -t6 * r_i_i_C(3) + t10 * t8, t6 * t2, -t6 * t12, 0; 0, t9, -t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:40
	% EndTime: 2020-05-02 23:52:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->13), mult. (47->18), div. (0->0), fcn. (55->8), ass. (0->15)
	t10 = qJ(2) + qJ(3);
	t22 = pkin(3) * cos(t10);
	t24 = -cos(qJ(2)) * pkin(2) - t22;
	t23 = pkin(3) * sin(t10);
	t11 = sin(qJ(4));
	t13 = sin(qJ(1));
	t14 = cos(qJ(4));
	t16 = cos(qJ(1));
	t5 = t13 * t11 - t16 * t14;
	t6 = -t16 * t11 - t13 * t14;
	t21 = t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t20 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t18 = pkin(1) + pkin(4) - t24;
	t17 = -sin(qJ(2)) * pkin(2) - t23;
	t1 = [-t13 * t18 + t21, t17 * t16, -t16 * t23, t21; t16 * t18 + t20, t17 * t13, -t13 * t23, t20; 0, t24, -t22, 0;];
	Ja_transl = t1;
end