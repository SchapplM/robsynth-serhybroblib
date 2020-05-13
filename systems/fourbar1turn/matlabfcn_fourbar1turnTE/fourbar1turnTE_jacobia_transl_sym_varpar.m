% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% fourbar1turnTE
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% Ja_transl [3x2]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = fourbar1turnTE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_jacobia_transl_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbar1turnTE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnTE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
Ja_transl=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0; 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (7->4), mult. (20->10), div. (0->0), fcn. (20->4), ass. (0->7)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t7 = [t4 * r_i_i_C(3) - t2 * t6, t5 * t4; t2 * r_i_i_C(3) + t4 * t6, t5 * t2; 0, t6;];
	Ja_transl = t7;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (1306->56), mult. (1877->101), div. (88->4), fcn. (571->6), ass. (0->47)
	t28 = pkin(1) ^ 2;
	t24 = cos(qJ(2));
	t46 = t24 * pkin(2);
	t42 = -0.2e1 * pkin(1) * t46 + t28;
	t18 = pkin(2) ^ 2 + t42;
	t15 = pkin(3) ^ 2 - pkin(4) ^ 2 + t18;
	t19 = pkin(1) * t24 - pkin(2);
	t22 = sin(qJ(2));
	t52 = -pkin(3) - pkin(4);
	t13 = (pkin(2) - t52) * (pkin(2) + t52) + t42;
	t51 = -pkin(3) + pkin(4);
	t14 = (pkin(2) - t51) * (pkin(2) + t51) + t42;
	t29 = sqrt(-t13 * t14);
	t43 = t22 * t29;
	t7 = -pkin(1) * t43 - t19 * t15;
	t56 = -t7 / 0.2e1;
	t55 = t22 / 0.2e1;
	t54 = -t24 / 0.2e1;
	t53 = t24 / 0.2e1;
	t23 = sin(qJ(1));
	t16 = 0.1e1 / t18;
	t26 = 0.1e1 / pkin(3);
	t45 = t16 * t26;
	t37 = t45 * t53;
	t40 = t22 * t45;
	t12 = pkin(1) * t22 * t15;
	t8 = -t19 * t29 + t12;
	t50 = (t7 * t37 - t8 * t40 / 0.2e1) * t23;
	t25 = cos(qJ(1));
	t49 = (t8 * t37 + t7 * t40 / 0.2e1) * t25;
	t47 = t22 * pkin(2);
	t48 = 0.1e1 / t29 * (-t13 - t14) * pkin(1) * t47;
	t44 = t22 * t24;
	t41 = pkin(1) * pkin(2) / t18 ^ 2;
	t1 = t12 + (-t29 * t24 + (-0.2e1 * t19 * pkin(2) - t48) * t22) * pkin(1);
	t39 = -t1 / 0.2e1 + t8 / 0.2e1;
	t21 = t22 ^ 2;
	t2 = -t19 * t48 + 0.2e1 * t28 * t21 * pkin(2) + (t24 * t15 + t43) * pkin(1);
	t38 = t56 - t2 / 0.2e1;
	t36 = t21 * t7 + t8 * t44;
	t35 = -t21 * t8 + t7 * t44;
	t34 = t1 * t54 + t2 * t55;
	t33 = t1 * t55 + t2 * t53;
	t32 = t7 * t54 + t8 * t55;
	t31 = t22 * t56 + t8 * t54;
	t30 = (t35 * r_i_i_C(1) - t36 * r_i_i_C(2)) * t41;
	t3 = [t50 * r_i_i_C(1) + t25 * r_i_i_C(3) + (r_i_i_C(2) * t31 * t45 - t46) * t23, t49 * r_i_i_C(1) + (-t47 + (t30 + (t34 * r_i_i_C(1) + (-t32 + t33) * r_i_i_C(2)) * t16) * t26) * t25; t49 * r_i_i_C(2) + t23 * r_i_i_C(3) + (r_i_i_C(1) * t32 * t45 + t46) * t25, t50 * r_i_i_C(2) + (-t47 + (t30 + ((-t31 + t34) * r_i_i_C(1) + t33 * r_i_i_C(2)) * t16) * t26) * t23; 0, t46 + ((t36 * r_i_i_C(1) + t35 * r_i_i_C(2)) * t41 + ((t38 * r_i_i_C(1) + t39 * r_i_i_C(2)) * t24 + (t39 * r_i_i_C(1) - t38 * r_i_i_C(2)) * t22) * t16) * t26;];
	Ja_transl = t3;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (519->33), mult. (748->58), div. (32->4), fcn. (213->6), ass. (0->33)
	t22 = pkin(2) ^ 2;
	t19 = cos(qJ(2));
	t32 = pkin(2) * t19;
	t29 = (-0.2e1 * t32 + pkin(1)) * pkin(1);
	t14 = t22 + t29;
	t12 = 0.1e1 / t14;
	t21 = 0.1e1 / pkin(4);
	t15 = pkin(1) - t32;
	t34 = -pkin(3) + pkin(4);
	t10 = (pkin(2) - t34) * (pkin(2) + t34) + t29;
	t35 = -pkin(3) - pkin(4);
	t9 = (pkin(2) - t35) * (pkin(2) + t35) + t29;
	t24 = sqrt(-t9 * t10);
	t11 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t14;
	t17 = sin(qJ(2));
	t33 = pkin(2) * t17;
	t8 = t11 * t33;
	t4 = t15 * t24 + t8;
	t37 = r_i_i_C(2) * t4;
	t30 = t17 * t24;
	t3 = -pkin(2) * t30 + t15 * t11;
	t38 = r_i_i_C(1) * t3;
	t41 = (t38 / 0.2e1 + t37 / 0.2e1) * t12 * t21 - pkin(1);
	t40 = 0.2e1 * pkin(1);
	t28 = pkin(1) * t33;
	t36 = 0.1e1 / t24 * (-t10 - t9) * t28;
	t39 = -t8 / 0.2e1 - (-t19 * t24 + (t15 * t40 - t36) * t17) * pkin(2) / 0.2e1;
	t27 = 0.1e1 / t14 ^ 2 * t28;
	t2 = t15 * t36 + t22 * t17 ^ 2 * t40 + (t19 * t11 + t30) * pkin(2);
	t25 = ((r_i_i_C(1) * t39 - t2 * r_i_i_C(2) / 0.2e1) * t12 + (t37 + t38) * t27) * t21;
	t20 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [t20 * r_i_i_C(3) + t41 * t18, t20 * t25; t18 * r_i_i_C(3) - t41 * t20, t18 * t25; 0, ((t2 * r_i_i_C(1) / 0.2e1 + r_i_i_C(2) * t39) * t12 + (-r_i_i_C(1) * t4 + r_i_i_C(2) * t3) * t27) * t21;];
	Ja_transl = t1;
end