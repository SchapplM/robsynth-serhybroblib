% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% fourbarprisDE2
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% Ja_transl [3x1]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = fourbarprisDE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_jacobia_transl_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbarprisDE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisDE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_jacobia_transl_sym_varpar: pkin has to be [3x1] (double)');
Ja_transl=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:45:26
	% EndTime: 2020-05-07 09:45:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:45:27
	% EndTime: 2020-05-07 09:45:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (154->19), mult. (70->25), div. (16->4), fcn. (8->2), ass. (0->16)
	t11 = qJ(1) + pkin(3);
	t16 = -pkin(2) + t11;
	t7 = (pkin(1) + t16);
	t17 = (-pkin(2) - t11);
	t8 = (pkin(1) + t17);
	t19 = (t7 * t8);
	t6 = (pkin(1) - t16);
	t18 = (t6 * t19);
	t10 = 1 / (t11 ^ 2);
	t9 = 0.1e1 / t11;
	t15 = -t11 * t9 - ((-pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) * t10) / 0.2e1;
	t5 = pkin(1) - t17;
	t13 = sqrt(-(t5 * t18));
	t14 = t10 * t13 / 0.2e1 - t9 / t13 * (-t18 + (t19 + (t7 - t8) * t6) * t5) / 0.4e1;
	t12 = 1 / pkin(1);
	t1 = [(t15 * r_i_i_C(1) - t14 * r_i_i_C(2)) * t12; (t14 * r_i_i_C(1) + t15 * r_i_i_C(2)) * t12; 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:45:27
	% EndTime: 2020-05-07 09:45:27
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (251->32), mult. (119->40), div. (26->4), fcn. (14->2), ass. (0->22)
	t30 = r_i_i_C(3) / 0.2e1;
	t11 = qJ(1) + pkin(3);
	t21 = -pkin(2) + t11;
	t7 = pkin(1) + t21;
	t22 = -pkin(2) - t11;
	t8 = pkin(1) + t22;
	t28 = t7 * t8;
	t6 = pkin(1) - t21;
	t23 = t6 * t28;
	t5 = pkin(1) - t22;
	t17 = sqrt(-t5 * t23);
	t29 = (-t23 + (t28 + (t7 - t8) * t6) * t5) / t17;
	t26 = (pkin(3) * qJ(1));
	t24 = (pkin(2) ^ 2 - pkin(3) ^ 2);
	t15 = (pkin(1) ^ 2);
	t20 = t15 + t24;
	t16 = 0.1e1 / pkin(1);
	t12 = (qJ(1) ^ 2);
	t10 = 0.1e1 / t11 ^ 2;
	t9 = 0.1e1 / t11;
	t4 = t12 + t15 - t24 + 2 * t26;
	t1 = [((-r_i_i_C(1) * t29 / 0.4e1 - t11 * r_i_i_C(3) - 0.3e1 / 0.2e1 * t12 + t20 / 0.2e1 - (2 * t26)) * t9 + (-((-t12 + t15) * pkin(3)) - (-t12 + t20) * qJ(1) / 0.2e1 + t17 * r_i_i_C(1) / 0.2e1 + t4 * t30) * t10) * t16; (r_i_i_C(1) + (-t17 / 0.2e1 + (-r_i_i_C(3) / 0.4e1 - qJ(1) / 0.4e1) * t29) * t9 + (-t4 * r_i_i_C(1) / 0.2e1 + (t30 + qJ(1) / 0.2e1) * t17) * t10) * t16; 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:45:27
	% EndTime: 2020-05-07 09:45:27
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (108->13), mult. (38->14), div. (8->2), fcn. (4->2), ass. (0->12)
	t7 = -pkin(3) - qJ(1);
	t10 = -pkin(2) - t7;
	t5 = pkin(1) + t10;
	t11 = -pkin(2) + t7;
	t6 = pkin(1) + t11;
	t14 = t5 * t6;
	t4 = pkin(1) - t10;
	t12 = t4 * t14;
	t3 = pkin(1) - t11;
	t15 = (-t12 + (t14 + (t5 - t6) * t4) * t3) * (-t3 * t12) ^ (-0.1e1 / 0.2e1);
	t13 = 0.1e1 / pkin(2) / pkin(1);
	t1 = [(t7 * r_i_i_C(1) + r_i_i_C(2) * t15 / 0.4e1) * t13; (-r_i_i_C(1) * t15 / 0.4e1 + t7 * r_i_i_C(2)) * t13; 0;];
	Ja_transl = t1;
end