% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% fourbar1DE1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% Ja_transl [3x1]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = fourbar1DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_jacobia_transl_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbar1DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_jacobia_transl_sym_varpar: pkin has to be [4x1] (double)');
Ja_transl=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:33
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:33
	% EndTime: 2020-04-24 19:57:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2; r_i_i_C(1) * t2 - r_i_i_C(2) * t1; 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:33
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (288->32), mult. (410->48), div. (16->4), fcn. (112->4), ass. (0->27)
	t34 = -2 * pkin(1);
	t18 = cos(qJ(1));
	t29 = pkin(2) * t18;
	t14 = -pkin(1) + t29;
	t17 = sin(qJ(1));
	t27 = (pkin(1) ^ 2) + t29 * t34;
	t31 = -pkin(3) - pkin(4);
	t8 = (pkin(2) - t31) * (pkin(2) + t31) + t27;
	t30 = -pkin(3) + pkin(4);
	t9 = (pkin(2) - t30) * (pkin(2) + t30) + t27;
	t22 = sqrt(-t8 * t9);
	t28 = t17 * pkin(2);
	t32 = 0.1e1 / t22 * (-t8 - t9) * pkin(1) * t28;
	t20 = pkin(2) ^ 2;
	t13 = t20 + t27;
	t10 = pkin(3) ^ 2 - pkin(4) ^ 2 + t13;
	t7 = t10 * t28;
	t33 = t7 / 0.2e1 + (t18 * t22 + (t14 * t34 + t32) * t17) * pkin(2) / 0.2e1;
	t26 = 0.1e1 / t13 ^ 2 * t28;
	t11 = 0.1e1 / t13;
	t25 = t20 * t17 ^ 2 * t11;
	t24 = t14 * t22 + t7;
	t4 = t22 * t28;
	t23 = t10 * t29 + t14 * t32 - t4;
	t19 = 0.1e1 / pkin(3);
	t2 = -t14 * t10 + t4;
	t1 = [-t28 + ((r_i_i_C(1) * t33 + t23 * r_i_i_C(2) / 0.2e1) * t11 + (r_i_i_C(2) * t25 + (-t2 * r_i_i_C(1) - t24 * r_i_i_C(2)) * t26) * pkin(1)) * t19; t29 + ((-t23 * r_i_i_C(1) / 0.2e1 + r_i_i_C(2) * t33) * t11 + (-r_i_i_C(1) * t25 + (t24 * r_i_i_C(1) - t2 * r_i_i_C(2)) * t26) * pkin(1)) * t19; 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:33
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (286->29), mult. (408->47), div. (16->4), fcn. (110->4), ass. (0->24)
	t19 = pkin(2) ^ 2;
	t17 = cos(qJ(1));
	t28 = pkin(2) * t17;
	t27 = (-0.2e1 * t28 + pkin(1)) * pkin(1);
	t12 = t19 + t27;
	t10 = 0.1e1 / t12;
	t13 = -pkin(1) + t28;
	t16 = sin(qJ(1));
	t31 = -pkin(3) - pkin(4);
	t7 = (pkin(2) - t31) * (pkin(2) + t31) + t27;
	t30 = -pkin(3) + pkin(4);
	t8 = (pkin(2) - t30) * (pkin(2) + t30) + t27;
	t21 = sqrt(-t7 * t8);
	t29 = pkin(2) * t16;
	t26 = 0.1e1 / t12 ^ 2 * t29;
	t9 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t12;
	t34 = (-t13 * t21 + t9 * t29) * t26 - t19 * t16 ^ 2 * t10;
	t32 = 0.1e1 / t21 * (-t7 - t8) * pkin(1) * t29;
	t33 = (t17 * t21 + (0.2e1 * t13 * pkin(1) + t32 - t9) * t16) * pkin(2) / 0.2e1;
	t4 = t21 * t29;
	t24 = (t13 * t9 + t4) * t26;
	t22 = -t13 * t32 + t9 * t28 + t4;
	t18 = 0.1e1 / pkin(4);
	t1 = [((r_i_i_C(1) * t33 - t22 * r_i_i_C(2) / 0.2e1) * t10 + (-r_i_i_C(1) * t24 + t34 * r_i_i_C(2)) * pkin(1)) * t18; ((t22 * r_i_i_C(1) / 0.2e1 + r_i_i_C(2) * t33) * t10 + (-t34 * r_i_i_C(1) - r_i_i_C(2) * t24) * pkin(1)) * t18; 0;];
	Ja_transl = t1;
end