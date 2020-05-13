% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = fourbar1turnDE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_jacobia_transl_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbar1turnDE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
Ja_transl=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
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
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
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
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:04
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (6150->66), mult. (8941->120), div. (432->6), fcn. (2523->8), ass. (0->56)
	t35 = pkin(1) ^ 2;
	t30 = cos(qJ(2));
	t65 = t30 * pkin(2);
	t54 = -0.2e1 * pkin(1) * t65 + t35;
	t24 = pkin(2) ^ 2 + t54;
	t22 = 0.1e1 / t24 ^ 2;
	t38 = pkin(3) ^ 2;
	t72 = t22 / t38;
	t20 = -pkin(4) ^ 2 + t24 + t38;
	t28 = sin(qJ(2));
	t17 = pkin(1) * t28 * t20;
	t25 = pkin(1) * t30 - pkin(2);
	t69 = -pkin(3) - pkin(4);
	t18 = (pkin(2) - t69) * (pkin(2) + t69) + t54;
	t68 = pkin(4) - pkin(3);
	t19 = (pkin(2) - t68) * (pkin(2) + t68) + t54;
	t37 = sqrt(-t18 * t19);
	t13 = -t25 * t37 + t17;
	t66 = t28 * pkin(2);
	t53 = pkin(1) * t66;
	t58 = 0.1e1 / t37 * (-t18 - t19) * t53;
	t5 = t17 + (-t37 * t30 + (-0.2e1 * t25 * pkin(2) - t58) * t28) * pkin(1);
	t63 = t13 - t5;
	t56 = t28 * t37;
	t12 = -pkin(1) * t56 - t20 * t25;
	t27 = t28 ^ 2;
	t6 = -t25 * t58 + 0.2e1 * t35 * t27 * pkin(2) + (t20 * t30 + t56) * pkin(1);
	t64 = t12 + t6;
	t71 = t64 * r_i_i_C(1) - t63 * r_i_i_C(2);
	t21 = 0.1e1 / t24;
	t55 = t12 ^ 2 + t13 ^ 2;
	t9 = t55 * t72;
	t7 = t9 ^ (-0.1e1 / 0.2e1);
	t70 = 0.2e1 * (-0.2e1 * t21 * t53 * t55 + t12 * t5 + t13 * t6) * t7 / t9 * t72;
	t31 = cos(qJ(1));
	t32 = 0.1e1 / pkin(3);
	t52 = t21 * t32 * t7;
	t49 = t30 * t52;
	t62 = t12 * t28;
	t67 = (t13 * t49 + t52 * t62) * t31;
	t61 = t12 * t30;
	t60 = t13 * t28;
	t59 = t13 * t30;
	t57 = t28 * t30;
	t51 = t63 * r_i_i_C(1);
	t50 = pkin(1) * pkin(2) * t22 * t7;
	t48 = t60 - t61;
	t47 = t12 * t27 + t13 * t57;
	t46 = t61 / 0.2e1 - t60 / 0.2e1;
	t45 = -t59 / 0.2e1 - t62 / 0.2e1;
	t44 = 0.2e1 * t12 * t57 - 0.2e1 * t13 * t27;
	t43 = (t46 * r_i_i_C(1) + t45 * r_i_i_C(2)) * t70;
	t42 = (r_i_i_C(1) * t44 - 0.2e1 * r_i_i_C(2) * t47) * t50;
	t29 = sin(qJ(1));
	t2 = t29 * t12 * t49;
	t1 = [t2 * r_i_i_C(1) + t31 * r_i_i_C(3) + (-t65 + (-r_i_i_C(1) * t60 + (-t59 - t62) * r_i_i_C(2)) * t52) * t29, t67 * r_i_i_C(1) + (-t66 + (t42 + (t43 + ((t28 * t6 - t30 * t5) * r_i_i_C(1) + (t28 * t5 + t30 * t6 - t48) * r_i_i_C(2)) * t7) * t21) * t32) * t31; t67 * r_i_i_C(2) + t29 * r_i_i_C(3) + (r_i_i_C(1) * t48 * t52 + t65) * t31, t2 * r_i_i_C(2) + (-t66 + (t42 + (t43 + ((t6 * r_i_i_C(2) + t51) * t30 + t71 * t28) * t7) * t21) * t32) * t29; 0, t65 + ((0.2e1 * r_i_i_C(1) * t47 + r_i_i_C(2) * t44) * t50 + ((-r_i_i_C(1) * t45 + r_i_i_C(2) * t46) * t70 + (-t71 * t30 + (t64 * r_i_i_C(2) + t51) * t28) * t7) * t21) * t32;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:04
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (2683->37), mult. (3902->74), div. (180->6), fcn. (1075->8), ass. (0->42)
	t30 = pkin(2) ^ 2;
	t26 = cos(qJ(2));
	t49 = pkin(2) * t26;
	t46 = (-0.2e1 * t49 + pkin(1)) * pkin(1);
	t21 = t30 + t46;
	t34 = pkin(4) ^ 2;
	t17 = -pkin(3) ^ 2 + t21 + t34;
	t24 = sin(qJ(2));
	t50 = pkin(2) * t24;
	t14 = t17 * t50;
	t22 = pkin(1) - t49;
	t53 = -pkin(3) - pkin(4);
	t15 = (pkin(2) - t53) * (pkin(2) + t53) + t46;
	t52 = -pkin(3) + pkin(4);
	t16 = (pkin(2) - t52) * (pkin(2) + t52) + t46;
	t33 = sqrt(-t15 * t16);
	t10 = t22 * t33 + t14;
	t19 = 0.1e1 / t21 ^ 2;
	t47 = t24 * t33;
	t9 = -pkin(2) * t47 + t17 * t22;
	t57 = t19 * (t10 ^ 2 + t9 ^ 2);
	t18 = 0.1e1 / t21;
	t28 = 0.1e1 / pkin(4);
	t29 = 0.1e1 / t34;
	t6 = t29 * t57;
	t4 = t6 ^ (-0.1e1 / 0.2e1);
	t42 = r_i_i_C(1) * t9 + r_i_i_C(2) * t10;
	t56 = t42 * t4 * t28 * t18 - pkin(1);
	t55 = 0.2e1 * pkin(1);
	t45 = pkin(1) * t50;
	t48 = 0.1e1 / t33 * (-t15 - t16) * t45;
	t2 = t14 + (-t26 * t33 + (t22 * t55 - t48) * t24) * pkin(2);
	t3 = t22 * t48 + t30 * t24 ^ 2 * t55 + (t17 * t26 + t47) * pkin(2);
	t44 = 0.2e1 * t19;
	t54 = ((t10 * t3 + t2 * t9) * t44 - 0.4e1 * t18 * t45 * t57) * t29 * t4 / t6;
	t41 = -t2 * t4 + t9 * t54 / 0.2e1;
	t40 = t4 * t44 * t45;
	t39 = t3 * t4 - t10 * t54 / 0.2e1;
	t38 = t28 * (t42 * t40 + (t41 * r_i_i_C(1) - t39 * r_i_i_C(2)) * t18);
	t27 = cos(qJ(1));
	t25 = sin(qJ(1));
	t1 = [r_i_i_C(3) * t27 + t56 * t25, t27 * t38; r_i_i_C(3) * t25 - t56 * t27, t25 * t38; 0, ((-r_i_i_C(1) * t10 + r_i_i_C(2) * t9) * t40 + (t39 * r_i_i_C(1) + t41 * r_i_i_C(2)) * t18) * t28;];
	Ja_transl = t1;
end