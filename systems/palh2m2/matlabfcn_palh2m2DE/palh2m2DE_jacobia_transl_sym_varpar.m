% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m2DE
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh2m2DE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m2DE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2DE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t2 * t5, t6 * t4, 0, 0; t2 * r_i_i_C(3) + t4 * t5, t6 * t2, 0, 0; 0, t7, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->6)
	t7 = sin(qJ(2)) * pkin(4);
	t2 = cos(qJ(2)) * pkin(4);
	t6 = r_i_i_C(1) + t2 + pkin(1);
	t5 = cos(qJ(1));
	t4 = sin(qJ(1));
	t1 = [t5 * r_i_i_C(3) - t4 * t6, -t5 * t7, 0, 0; t4 * r_i_i_C(3) + t5 * t6, -t4 * t7, 0, 0; 0, t2, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->7), mult. (29->14), div. (0->0), fcn. (29->6), ass. (0->10)
	t11 = sin(qJ(2)) * pkin(4);
	t3 = sin(qJ(3));
	t6 = cos(qJ(3));
	t10 = t6 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t6;
	t2 = cos(qJ(2)) * pkin(4);
	t8 = pkin(2) + t2 + pkin(1) + t10;
	t7 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [t7 * r_i_i_C(3) - t8 * t5, -t7 * t11, t9 * t7, 0; t5 * r_i_i_C(3) + t8 * t7, -t5 * t11, t9 * t5, 0; 0, t2, t10, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->7), mult. (20->12), div. (0->0), fcn. (20->6), ass. (0->8)
	t10 = sin(qJ(2)) * pkin(4);
	t9 = sin(qJ(3)) * pkin(5);
	t2 = cos(qJ(3)) * pkin(5);
	t3 = cos(qJ(2)) * pkin(4);
	t8 = r_i_i_C(1) + pkin(2) + t2 + t3 + pkin(1);
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t7 * r_i_i_C(3) - t8 * t6, -t7 * t10, -t7 * t9, 0; t6 * r_i_i_C(3) + t8 * t7, -t6 * t10, -t6 * t9, 0; 0, t3, t2, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->11), mult. (40->18), div. (0->0), fcn. (48->8), ass. (0->14)
	t10 = sin(qJ(4));
	t13 = sin(qJ(1));
	t14 = cos(qJ(4));
	t15 = cos(qJ(1));
	t6 = t10 * t13 - t14 * t15;
	t7 = -t10 * t15 - t13 * t14;
	t19 = t7 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t18 = -t6 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t17 = sin(qJ(2)) * pkin(4);
	t16 = sin(qJ(3)) * pkin(5);
	t9 = cos(qJ(2)) * pkin(4);
	t8 = cos(qJ(3)) * pkin(5);
	t5 = t9 + t8 + pkin(1) + pkin(2) + pkin(3);
	t1 = [-t13 * t5 + t19, -t15 * t17, -t15 * t16, t19; t15 * t5 + t18, -t13 * t17, -t13 * t16, t18; 0, t9, t8, 0;];
	Ja_transl = t1;
end