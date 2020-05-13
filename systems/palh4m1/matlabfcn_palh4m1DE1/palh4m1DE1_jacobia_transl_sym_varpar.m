% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh4m1DE1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:26
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh4m1DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh4m1DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Ja_transl = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->17)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t3 = cos(qJ(1));
	unknown(1,1) = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(2,1) = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (6178->85), mult. (6932->210), div. (336->12), fcn. (1849->8), ass. (0->94)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = pkin(1) * t3;
	t6 = 0.2e1 * pkin(2) * t4;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t15 * t11);
	t19 = pkin(1) * t17 * t2;
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t25 * t20 - t19;
	t28 = t27 * t1;
	t29 = 0.1e1 / t22;
	t30 = t29 * t28;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t34 = t2 * pkin(1);
	t35 = t25 * t34;
	t36 = -t17 * t20 + t35;
	t37 = t36 ^ 2;
	t38 = 0.1e1 / t23;
	t39 = t38 * t37;
	t40 = t31 ^ 2;
	t41 = 0.1e1 / t40;
	t43 = t27 ^ 2;
	t44 = t38 * t43;
	t46 = t41 * t39 + t41 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = t48 * t32;
	t50 = r_i_i_C(1) * t49;
	t52 = t36 * t1;
	t53 = t29 * t52;
	t54 = r_i_i_C(2) * t49;
	t56 = cos(qJ(1));
	t60 = 0.1e1 / t17;
	t61 = t60 * t2;
	t66 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t71 = -t66 * pkin(1) * t61 / 0.2e1 - 0.2e1 * t22 * t20;
	t75 = t27 * t56;
	t78 = t29 * t75;
	t80 = 0.1e1 / t47 / t46;
	t81 = t80 * t32;
	t82 = t38 * t36;
	t83 = -t60 * t20;
	t87 = t66 * t83 / 0.2e1 + 0.2e1 * t22 * t34;
	t91 = 0.1e1 / t23 / t22;
	t94 = t38 * t27;
	t99 = -t41 * t91 * t37 - t41 * t91 * t43 + t71 * t41 * t94 + t87 * t41 * t82;
	t101 = 0.2e1 * t99 * r_i_i_C(1) * t81;
	t107 = t36 * t56;
	t110 = t29 * t107;
	t112 = 0.2e1 * t99 * r_i_i_C(2) * t81;
	t121 = pkin(1) * pkin(2);
	t123 = t15 * pkin(2) * t34 + t121 * t2 * t11;
	t130 = -t123 * pkin(1) * t61 + pkin(1) * t17 * t3 + 0.2e1 * t121 * t2 * t20 - t35;
	t134 = t41 * t29;
	t137 = pkin(2) * t34;
	t138 = t137 * r_i_i_C(1) * t48;
	t144 = t2 ^ 2;
	t148 = -0.2e1 * pkin(2) * t144 * t7 + t123 * t83 - t25 * t4 - t19;
	t153 = 0.1e1 / t40 / t31;
	t163 = 0.2e1 * t130 * t41 * t94 + 0.4e1 * t137 * t153 * t39 + 0.4e1 * t137 * t153 * t44 + 0.2e1 * t148 * t41 * t82;
	t165 = t163 * r_i_i_C(1) * t81;
	t173 = t137 * r_i_i_C(2) * t48;
	t177 = t163 * r_i_i_C(2) * t81;
	t221 = t29 * t36;
	t222 = t32 * t221;
	t223 = r_i_i_C(1) * t80;
	t230 = t29 * t27;
	t231 = t32 * t230;
	t232 = r_i_i_C(2) * t80;
	t239 = t48 * t41;
	unknown(1,1) = pkin(9) * t1 + r_i_i_C(3) * t56 - t50 * t30 + t54 * t53;
	unknown(1,2) = t50 * t29 * t71 * t56 - t50 * t38 * t75 - t101 * t78 / 0.2e1 - t54 * t29 * t87 * t56 + t54 * t38 * t107 + t112 * t110 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t50 * t29 * t130 * t56 + 0.2e1 * t138 * t134 * t75 - t165 * t78 / 0.2e1 - t54 * t29 * t148 * t56 - 0.2e1 * t173 * t134 * t107 + t177 * t110 / 0.2e1;
	unknown(2,1) = -pkin(9) * t56 + r_i_i_C(3) * t1 - t54 * t110 + t50 * t78;
	unknown(2,2) = t50 * t29 * t71 * t1 - t50 * t38 * t28 - t101 * t30 / 0.2e1 - t54 * t29 * t87 * t1 + t54 * t38 * t52 + t112 * t53 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t50 * t29 * t130 * t1 + 0.2e1 * t138 * t134 * t28 - t165 * t30 / 0.2e1 - t54 * t29 * t148 * t1 - 0.2e1 * t173 * t134 * t52 + t177 * t53 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t99 * t223 * t222 - t99 * t232 * t231 + t50 * t29 * t87 + t54 * t29 * t71 - t50 * t82 - t54 * t94;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t50 * t29 * t148 + 0.2e1 * t121 * t2 * r_i_i_C(1) * t239 * t221 - t163 * t223 * t222 / 0.2e1 + t54 * t29 * t130 + 0.2e1 * t121 * t2 * r_i_i_C(2) * t239 * t230 - t163 * t232 * t231 / 0.2e1;
	Ja_transl = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (9003->101), mult. (10095->249), div. (469->12), fcn. (2692->8), ass. (0->102)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t15 * t11);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t28 = t27 * t1;
	t29 = 0.1e1 / t22;
	t30 = t29 * t28;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = t27 ^ 2;
	t34 = 0.1e1 / t23;
	t35 = t34 * t33;
	t36 = t31 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t34 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = t48 * t32;
	t50 = r_i_i_C(1) * t49;
	t52 = cos(qJ(1));
	t54 = t42 * t1;
	t55 = t29 * t54;
	t56 = r_i_i_C(3) * t49;
	t61 = 0.1e1 / t17;
	t62 = t61 * t4;
	t67 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t71 = t67 * t62 / 0.2e1 + 0.2e1 * t22 * t20;
	t75 = t27 * t52;
	t78 = t29 * t75;
	t80 = 0.1e1 / t47 / t46;
	t81 = t80 * t32;
	t82 = t34 * t27;
	t86 = 0.1e1 / t23 / t22;
	t89 = t34 * t42;
	t90 = t61 * t19;
	t95 = -t67 * pkin(1) * t90 / 0.2e1 + 0.2e1 * t22 * t4;
	t100 = -t37 * t86 * t33 - t37 * t86 * t43 + t71 * t37 * t82 + t95 * t37 * t89;
	t102 = 0.2e1 * t100 * r_i_i_C(1) * t81;
	t105 = t95 * t52;
	t108 = t42 * t52;
	t111 = t29 * t108;
	t113 = 0.2e1 * t100 * r_i_i_C(3) * t81;
	t117 = 0.2e1 * t100 * t81;
	t124 = pkin(1) * pkin(2);
	t126 = t15 * pkin(2) * t20 + t124 * t19 * t11;
	t130 = t19 ^ 2;
	t134 = -0.2e1 * pkin(2) * t130 * t7 + t126 * t62 - t25 * t3 - t40;
	t138 = t37 * t29;
	t141 = pkin(2) * t20;
	t142 = t141 * r_i_i_C(1) * t48;
	t149 = 0.1e1 / t36 / t31;
	t161 = -t126 * pkin(1) * t90 + pkin(1) * t17 * t2 - 0.2e1 * t124 * t19 * t4 - t26;
	t168 = 0.2e1 * t134 * t37 * t82 + 0.4e1 * t141 * t149 * t35 + 0.4e1 * t141 * t149 * t44 + 0.2e1 * t161 * t37 * t89;
	t170 = t168 * r_i_i_C(1) * t81;
	t173 = t161 * t52;
	t178 = t141 * r_i_i_C(3) * t48;
	t182 = t168 * r_i_i_C(3) * t81;
	t188 = t124 * t19 * t48;
	t191 = t168 * t81;
	t208 = t95 * t1;
	t227 = t161 * t1;
	t245 = t29 * t42;
	t246 = t32 * t245;
	t247 = r_i_i_C(1) * t80;
	t254 = t29 * t27;
	t255 = t32 * t254;
	t256 = r_i_i_C(3) * t80;
	t262 = t32 * t27;
	t269 = t48 * t37;
	unknown(1,1) = pkin(9) * t1 - r_i_i_C(2) * t52 - t50 * t30 - t49 * t54 - t56 * t55;
	unknown(1,2) = t50 * t29 * t71 * t52 - t50 * t34 * t75 - t102 * t78 / 0.2e1 + t56 * t29 * t105 - t56 * t34 * t108 - t113 * t111 / 0.2e1 + t49 * t105 - t117 * t108 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t50 * t29 * t134 * t52 + 0.2e1 * t142 * t138 * t75 - t170 * t78 / 0.2e1 + t56 * t29 * t173 + 0.2e1 * t178 * t138 * t108 - t182 * t111 / 0.2e1 + t49 * t173 + 0.2e1 * t188 * t37 * t108 - t191 * t108 / 0.2e1;
	unknown(2,1) = -pkin(9) * t52 - r_i_i_C(2) * t1 + t49 * t108 + t56 * t111 + t50 * t78;
	unknown(2,2) = t50 * t29 * t71 * t1 - t50 * t34 * t28 - t102 * t30 / 0.2e1 + t56 * t29 * t208 - t56 * t34 * t54 - t113 * t55 / 0.2e1 + t49 * t208 - t117 * t54 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t50 * t29 * t134 * t1 + 0.2e1 * t142 * t138 * t28 - t170 * t30 / 0.2e1 + t56 * t29 * t227 + 0.2e1 * t178 * t138 * t54 - t182 * t55 / 0.2e1 + t49 * t227 + 0.2e1 * t188 * t37 * t54 - t191 * t54 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t100 * t247 * t246 - t100 * t256 * t255 - t100 * t80 * t262 - t50 * t29 * t95 + t56 * t29 * t71 + t48 * t32 * t71 + t50 * t89 - t56 * t82;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t50 * t29 * t161 - 0.2e1 * t124 * t19 * r_i_i_C(1) * t269 * t245 + t168 * t247 * t246 / 0.2e1 + t56 * t29 * t134 + 0.2e1 * t124 * t19 * r_i_i_C(3) * t269 * t254 - t168 * t256 * t255 / 0.2e1 + t48 * t32 * t134 + 0.2e1 * t141 * t48 * t37 * t27 - t168 * t80 * t262 / 0.2e1;
	Ja_transl = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (24807->200), mult. (27489->516), div. (1925->16), fcn. (6951->10), ass. (0->176)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t28 = t27 * t1;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t28;
	t34 = t27 ^ 2;
	t35 = t29 * t34;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t29 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t50 = t49 * t48;
	t51 = 0.1e1 / t24;
	t52 = t51 * t29;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t58 = t51 * t29 * t55 - t16 * t52;
	t59 = sqrt(t58);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t50;
	t64 = t42 * t1;
	t65 = t32 * t64;
	t66 = t54 * t48;
	t67 = t60 * t49;
	t68 = t67 * t66;
	t76 = cos(qJ(1));
	t78 = t48 * t31;
	t82 = 0.1e1 / t17;
	t83 = t82 * t4;
	t85 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t87 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t88 = -t85 - t87;
	t92 = t88 * t83 / 0.2e1 + 0.2e1 * t22 * t20;
	t94 = t32 * t92 * t76;
	t96 = t27 * t76;
	t98 = 0.1e1 / t23 / t22;
	t99 = t31 * t98;
	t100 = t99 * t96;
	t103 = t32 * t96;
	t105 = 0.1e1 / t47 / t46;
	t106 = t49 * t105;
	t107 = t29 * t27;
	t112 = t29 * t42;
	t113 = t82 * t19;
	t118 = -t88 * pkin(1) * t113 / 0.2e1 + 0.2e1 * t22 * t4;
	t123 = t92 * t37 * t107 + t118 * t37 * t112 - t37 * t98 * t34 - t37 * t98 * t43;
	t125 = 0.2e1 * t123 * t61 * t106;
	t128 = t60 * t82;
	t130 = t88 * t128 * t50;
	t134 = 0.1e1 / t59 / t58;
	t135 = t134 * t17;
	t141 = t29 * t54;
	t148 = 0.4e1 * t22 * t51 * t141 + 0.2e1 * t16 * t51 * t98 - 0.2e1 * t51 * t98 * t55 - t85 * t52 - t87 * t52;
	t150 = t148 * t135 * t50;
	t153 = t118 * t76;
	t154 = t32 * t153;
	t156 = t42 * t76;
	t157 = t99 * t156;
	t160 = t32 * t156;
	t161 = t54 * t105;
	t163 = 0.2e1 * t123 * t67 * t161;
	t167 = 0.2e1 * t67 * t22 * t48;
	t169 = t134 * t49;
	t171 = t148 * t169 * t66;
	t196 = t105 * t31;
	t197 = 0.2e1 * t123 * t196;
	t204 = pkin(1) * pkin(2);
	t206 = t15 * pkin(2) * t20 + t204 * t19 * t11;
	t210 = t19 ^ 2;
	t214 = -0.2e1 * pkin(2) * t210 * t7 + t206 * t83 - t25 * t3 - t40;
	t216 = t32 * t214 * t76;
	t219 = t48 * t37 * t29;
	t220 = t219 * t96;
	t221 = t17 * t49;
	t223 = pkin(2) * t20;
	t224 = t223 * t60 * t221;
	t231 = 0.1e1 / t36 / t30;
	t243 = -t206 * pkin(1) * t113 + pkin(1) * t17 * t2 - 0.2e1 * t204 * t19 * t4 - t26;
	t250 = 0.2e1 * t214 * t37 * t107 + 0.2e1 * t243 * t37 * t112 + 0.4e1 * t223 * t231 * t35 + 0.4e1 * t223 * t231 * t44;
	t252 = t250 * t61 * t106;
	t256 = 0.2e1 * t206 * t128 * t50;
	t269 = 0.2e1 * t15 * t204 * t19 * t52 + 0.2e1 * t223 * t11 * t52 + 0.4e1 * t223 * t51 * t141;
	t271 = t269 * t135 * t50;
	t274 = t243 * t76;
	t275 = t32 * t274;
	t277 = t219 * t156;
	t278 = t49 * t54;
	t280 = t223 * t60 * t278;
	t284 = t250 * t67 * t161;
	t287 = t48 * t32;
	t291 = t60 * t49 * pkin(2) * t20;
	t295 = t269 * t169 * t66;
	t324 = t204 * t19 * t48;
	t327 = t250 * t196;
	t344 = t32 * t92 * t1;
	t346 = t99 * t28;
	t355 = t118 * t1;
	t356 = t32 * t355;
	t358 = t99 * t64;
	t392 = t32 * t214 * t1;
	t394 = t219 * t28;
	t403 = t243 * t1;
	t404 = t32 * t403;
	t406 = t219 * t64;
	t447 = t31 * t29 * t118;
	t450 = t31 * t98 * t42;
	t453 = t196 * t112;
	t454 = 0.2e1 * t123 * t60;
	t455 = t454 * t221;
	t458 = t78 * t112;
	t459 = t82 * t49;
	t461 = t88 * t60 * t459;
	t464 = t148 * t134;
	t465 = t464 * t221;
	t469 = t31 * t29 * t92;
	t472 = t31 * t98 * t27;
	t475 = t196 * t107;
	t476 = t454 * t278;
	t481 = t78 * t107;
	t482 = t464 * t278;
	t509 = t31 * t27;
	t515 = t31 * t29 * t243;
	t517 = t48 * t37;
	t518 = t49 * t517;
	t520 = t223 * t61;
	t523 = t250 * t60;
	t524 = t523 * t221;
	t528 = 0.2e1 * t206 * t60 * t459;
	t531 = t269 * t134;
	t532 = t531 * t221;
	t536 = t31 * t29 * t214;
	t538 = t54 * t517;
	t542 = t523 * t278;
	t547 = t531 * t278;
	unknown(1,1) = r_i_i_C(1) * (t62 * t33 - t68 * t65) + r_i_i_C(2) * (t68 * t33 + t62 * t65) + r_i_i_C(3) * t76 - t78 * t64 + pkin(9) * t1;
	unknown(1,2) = r_i_i_C(1) * (-t62 * t94 + 0.2e1 * t62 * t100 + t125 * t103 / 0.2e1 - t130 * t103 / 0.2e1 + t150 * t103 / 0.2e1 + t68 * t154 - 0.2e1 * t68 * t157 - t163 * t160 / 0.2e1 + t167 * t160 - t171 * t160 / 0.2e1) + r_i_i_C(2) * (-t68 * t94 + 0.2e1 * t68 * t100 + t163 * t103 / 0.2e1 - t167 * t103 + t171 * t103 / 0.2e1 - t62 * t154 + 0.2e1 * t62 * t157 + t125 * t160 / 0.2e1 - t130 * t160 / 0.2e1 + t150 * t160 / 0.2e1) + t78 * t153 - t197 * t156 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = r_i_i_C(1) * (-t62 * t216 - 0.2e1 * t224 * t220 + t252 * t103 / 0.2e1 - t256 * t103 / 0.2e1 + t271 * t103 / 0.2e1 + t68 * t275 + 0.2e1 * t280 * t277 - t284 * t160 / 0.2e1 + 0.2e1 * t291 * t287 * t156 - t295 * t160 / 0.2e1) + r_i_i_C(2) * (-t68 * t216 - 0.2e1 * t280 * t220 + t284 * t103 / 0.2e1 - 0.2e1 * t291 * t287 * t96 + t295 * t103 / 0.2e1 - t62 * t275 - 0.2e1 * t224 * t277 + t252 * t160 / 0.2e1 - t256 * t160 / 0.2e1 + t271 * t160 / 0.2e1) + t78 * t274 + 0.2e1 * t324 * t37 * t156 - t327 * t156 / 0.2e1;
	unknown(2,1) = r_i_i_C(1) * (-t62 * t103 + t68 * t160) + r_i_i_C(2) * (-t68 * t103 - t62 * t160) + r_i_i_C(3) * t1 + t78 * t156 - pkin(9) * t76;
	unknown(2,2) = r_i_i_C(1) * (-t62 * t344 + 0.2e1 * t62 * t346 + t125 * t33 / 0.2e1 - t130 * t33 / 0.2e1 + t150 * t33 / 0.2e1 + t68 * t356 - 0.2e1 * t68 * t358 - t163 * t65 / 0.2e1 + t167 * t65 - t171 * t65 / 0.2e1) + r_i_i_C(2) * (-t68 * t344 + 0.2e1 * t68 * t346 + t163 * t33 / 0.2e1 - t167 * t33 + t171 * t33 / 0.2e1 - t62 * t356 + 0.2e1 * t62 * t358 + t125 * t65 / 0.2e1 - t130 * t65 / 0.2e1 + t150 * t65 / 0.2e1) + t78 * t355 - t197 * t64 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = r_i_i_C(1) * (-t62 * t392 - 0.2e1 * t224 * t394 + t252 * t33 / 0.2e1 - t256 * t33 / 0.2e1 + t271 * t33 / 0.2e1 + t68 * t404 + 0.2e1 * t280 * t406 - t284 * t65 / 0.2e1 + 0.2e1 * t291 * t287 * t64 - t295 * t65 / 0.2e1) + r_i_i_C(2) * (-t68 * t392 - 0.2e1 * t280 * t394 + t284 * t33 / 0.2e1 - 0.2e1 * t291 * t287 * t28 + t295 * t33 / 0.2e1 - t62 * t404 - 0.2e1 * t224 * t406 + t252 * t65 / 0.2e1 - t256 * t65 / 0.2e1 + t271 * t65 / 0.2e1) + t78 * t403 + 0.2e1 * t324 * t37 * t64 - t327 * t64 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = r_i_i_C(1) * (t62 * t447 - 0.2e1 * t62 * t450 - t455 * t453 / 0.2e1 + t461 * t458 / 0.2e1 - t465 * t458 / 0.2e1 + t68 * t469 - 0.2e1 * t68 * t472 - t476 * t475 / 0.2e1 + t167 * t31 * t107 - t482 * t481 / 0.2e1) + r_i_i_C(2) * (t68 * t447 - 0.2e1 * t68 * t450 - t476 * t453 / 0.2e1 + t167 * t31 * t112 - t482 * t458 / 0.2e1 - t62 * t469 + 0.2e1 * t62 * t472 + t455 * t475 / 0.2e1 - t461 * t481 / 0.2e1 + t465 * t481 / 0.2e1) + t48 * t31 * t92 - t123 * t105 * t509;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = r_i_i_C(1) * (t62 * t515 + 0.2e1 * t520 * t518 * t112 - t524 * t453 / 0.2e1 + t528 * t458 / 0.2e1 - t532 * t458 / 0.2e1 + t68 * t536 + 0.2e1 * t291 * t538 * t107 - t542 * t475 / 0.2e1 + 0.2e1 * t291 * t481 - t547 * t481 / 0.2e1) + r_i_i_C(2) * (t68 * t515 + 0.2e1 * t291 * t538 * t112 - t542 * t453 / 0.2e1 + 0.2e1 * t291 * t458 - t547 * t458 / 0.2e1 - t62 * t536 - 0.2e1 * t520 * t518 * t107 + t524 * t475 / 0.2e1 - t528 * t481 / 0.2e1 + t532 * t481 / 0.2e1) + t48 * t31 * t214 + 0.2e1 * t223 * t48 * t37 * t27 - t250 * t105 * t509 / 0.2e1;
	Ja_transl = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:51
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (60615->233), mult. (67106->574), div. (4877->16), fcn. (16971->12), ass. (0->200)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t28 = t27 * t1;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t28;
	t34 = t27 ^ 2;
	t35 = t29 * t34;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t29 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t50 = t49 * t48;
	t51 = 0.1e1 / t24;
	t52 = t51 * t29;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t58 = t51 * t29 * t55 - t16 * t52;
	t59 = sqrt(t58);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t50;
	t64 = t42 * t1;
	t65 = t32 * t64;
	t66 = t54 * t48;
	t67 = t60 * t49;
	t68 = t67 * t66;
	t70 = t62 * t33 - t68 * t65;
	t71 = cos(qJ(3));
	t75 = t68 * t33 + t62 * t65;
	t76 = sin(qJ(3));
	t84 = cos(qJ(1));
	t87 = t48 * t31;
	t91 = 0.1e1 / t17;
	t92 = t91 * t4;
	t94 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t96 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t97 = -t94 - t96;
	t101 = t97 * t92 / 0.2e1 + 0.2e1 * t22 * t20;
	t103 = t32 * t101 * t84;
	t105 = t27 * t84;
	t107 = 0.1e1 / t23 / t22;
	t108 = t31 * t107;
	t109 = t108 * t105;
	t112 = t32 * t105;
	t114 = 0.1e1 / t47 / t46;
	t115 = t49 * t114;
	t116 = t29 * t27;
	t121 = t29 * t42;
	t122 = t91 * t19;
	t127 = -t97 * pkin(1) * t122 / 0.2e1 + 0.2e1 * t22 * t4;
	t132 = t101 * t37 * t116 - t37 * t107 * t34 - t37 * t107 * t43 + t127 * t37 * t121;
	t134 = 0.2e1 * t132 * t61 * t115;
	t137 = t60 * t91;
	t139 = t97 * t137 * t50;
	t143 = 0.1e1 / t59 / t58;
	t144 = t143 * t17;
	t150 = t29 * t54;
	t157 = 0.2e1 * t16 * t51 * t107 - 0.2e1 * t51 * t107 * t55 + 0.4e1 * t22 * t51 * t150 - t94 * t52 - t96 * t52;
	t159 = t157 * t144 * t50;
	t162 = t127 * t84;
	t163 = t32 * t162;
	t165 = t42 * t84;
	t166 = t108 * t165;
	t169 = t32 * t165;
	t170 = t54 * t114;
	t172 = 0.2e1 * t132 * t67 * t170;
	t176 = 0.2e1 * t67 * t22 * t48;
	t178 = t143 * t49;
	t180 = t157 * t178 * t66;
	t183 = -t62 * t103 + 0.2e1 * t62 * t109 + t134 * t112 / 0.2e1 - t139 * t112 / 0.2e1 + t159 * t112 / 0.2e1 + t68 * t163 - 0.2e1 * t68 * t166 - t172 * t169 / 0.2e1 + t176 * t169 - t180 * t169 / 0.2e1;
	t202 = -t68 * t103 + 0.2e1 * t68 * t109 + t172 * t112 / 0.2e1 - t176 * t112 + t180 * t112 / 0.2e1 - t62 * t163 + 0.2e1 * t62 * t166 + t134 * t169 / 0.2e1 - t139 * t169 / 0.2e1 + t159 * t169 / 0.2e1;
	t212 = t114 * t31;
	t213 = 0.2e1 * t132 * t212;
	t219 = -t62 * t112 + t68 * t169;
	t223 = -t68 * t112 - t62 * t169;
	t225 = -t76 * t219 + t71 * t223;
	t229 = -t71 * t219 - t76 * t223;
	t235 = pkin(1) * pkin(2);
	t237 = t15 * pkin(2) * t20 + t235 * t19 * t11;
	t241 = t19 ^ 2;
	t245 = -0.2e1 * pkin(2) * t241 * t7 + t237 * t92 - t25 * t3 - t40;
	t247 = t32 * t245 * t84;
	t250 = t48 * t37 * t29;
	t251 = t250 * t105;
	t252 = t17 * t49;
	t254 = pkin(2) * t20;
	t255 = t254 * t60 * t252;
	t262 = 0.1e1 / t36 / t30;
	t274 = -t237 * pkin(1) * t122 + pkin(1) * t17 * t2 - 0.2e1 * t235 * t19 * t4 - t26;
	t281 = 0.2e1 * t245 * t37 * t116 + 0.2e1 * t274 * t37 * t121 + 0.4e1 * t254 * t262 * t35 + 0.4e1 * t254 * t262 * t44;
	t283 = t281 * t61 * t115;
	t287 = 0.2e1 * t237 * t137 * t50;
	t300 = 0.2e1 * t15 * t235 * t19 * t52 + 0.2e1 * t254 * t11 * t52 + 0.4e1 * t254 * t51 * t150;
	t302 = t300 * t144 * t50;
	t305 = t274 * t84;
	t306 = t32 * t305;
	t308 = t250 * t165;
	t309 = t49 * t54;
	t311 = t254 * t60 * t309;
	t315 = t281 * t67 * t170;
	t318 = t48 * t32;
	t322 = t60 * t49 * pkin(2) * t20;
	t326 = t300 * t178 * t66;
	t329 = -t62 * t247 - 0.2e1 * t255 * t251 + t283 * t112 / 0.2e1 - t287 * t112 / 0.2e1 + t302 * t112 / 0.2e1 + t68 * t306 + 0.2e1 * t311 * t308 - t315 * t169 / 0.2e1 + 0.2e1 * t322 * t318 * t165 - t326 * t169 / 0.2e1;
	t350 = -t68 * t247 - 0.2e1 * t311 * t251 + t315 * t112 / 0.2e1 - 0.2e1 * t322 * t318 * t105 + t326 * t112 / 0.2e1 - t62 * t306 - 0.2e1 * t255 * t308 + t283 * t169 / 0.2e1 - t287 * t169 / 0.2e1 + t302 * t169 / 0.2e1;
	t362 = t235 * t19 * t48;
	t365 = t281 * t212;
	t377 = t32 * t101 * t1;
	t379 = t108 * t28;
	t388 = t127 * t1;
	t389 = t32 * t388;
	t391 = t108 * t64;
	t399 = -t62 * t377 + 0.2e1 * t62 * t379 + t134 * t33 / 0.2e1 - t139 * t33 / 0.2e1 + t159 * t33 / 0.2e1 + t68 * t389 - 0.2e1 * t68 * t391 - t172 * t65 / 0.2e1 + t176 * t65 - t180 * t65 / 0.2e1;
	t418 = -t68 * t377 + 0.2e1 * t68 * t379 + t172 * t33 / 0.2e1 - t176 * t33 + t180 * t33 / 0.2e1 - t62 * t389 + 0.2e1 * t62 * t391 + t134 * t65 / 0.2e1 - t139 * t65 / 0.2e1 + t159 * t65 / 0.2e1;
	t441 = t32 * t245 * t1;
	t443 = t250 * t28;
	t452 = t274 * t1;
	t453 = t32 * t452;
	t455 = t250 * t64;
	t465 = -t62 * t441 - 0.2e1 * t255 * t443 + t283 * t33 / 0.2e1 - t287 * t33 / 0.2e1 + t302 * t33 / 0.2e1 + t68 * t453 + 0.2e1 * t311 * t455 - t315 * t65 / 0.2e1 + 0.2e1 * t322 * t318 * t64 - t326 * t65 / 0.2e1;
	t486 = -t68 * t441 - 0.2e1 * t311 * t443 + t315 * t33 / 0.2e1 - 0.2e1 * t322 * t318 * t28 + t326 * t33 / 0.2e1 - t62 * t453 - 0.2e1 * t255 * t455 + t283 * t65 / 0.2e1 - t287 * t65 / 0.2e1 + t302 * t65 / 0.2e1;
	t503 = t31 * t29 * t127;
	t506 = t31 * t107 * t42;
	t509 = t212 * t121;
	t510 = 0.2e1 * t132 * t60;
	t511 = t510 * t252;
	t514 = t87 * t121;
	t515 = t91 * t49;
	t517 = t97 * t60 * t515;
	t520 = t157 * t143;
	t521 = t520 * t252;
	t525 = t31 * t29 * t101;
	t528 = t31 * t107 * t27;
	t531 = t212 * t116;
	t532 = t510 * t309;
	t535 = t31 * t116;
	t537 = t87 * t116;
	t538 = t520 * t309;
	t541 = t62 * t503 - 0.2e1 * t62 * t506 - t511 * t509 / 0.2e1 + t517 * t514 / 0.2e1 - t521 * t514 / 0.2e1 + t68 * t525 - 0.2e1 * t68 * t528 - t532 * t531 / 0.2e1 + t176 * t535 - t538 * t537 / 0.2e1;
	t548 = t31 * t121;
	t561 = t68 * t503 - 0.2e1 * t68 * t506 - t532 * t509 / 0.2e1 + t176 * t548 - t538 * t514 / 0.2e1 - t62 * t525 + 0.2e1 * t62 * t528 + t511 * t531 / 0.2e1 - t517 * t537 / 0.2e1 + t521 * t537 / 0.2e1;
	t572 = t31 * t27;
	t579 = t68 * t535 + t62 * t548;
	t583 = -t62 * t535 + t68 * t548;
	t593 = t31 * t29 * t274;
	t595 = t48 * t37;
	t596 = t49 * t595;
	t598 = t254 * t61;
	t601 = t281 * t60;
	t602 = t601 * t252;
	t606 = 0.2e1 * t237 * t60 * t515;
	t609 = t300 * t143;
	t610 = t609 * t252;
	t614 = t31 * t29 * t245;
	t616 = t54 * t595;
	t620 = t601 * t309;
	t625 = t609 * t309;
	t628 = t62 * t593 + 0.2e1 * t598 * t596 * t121 - t602 * t509 / 0.2e1 + t606 * t514 / 0.2e1 - t610 * t514 / 0.2e1 + t68 * t614 + 0.2e1 * t322 * t616 * t116 - t620 * t531 / 0.2e1 + 0.2e1 * t322 * t537 - t625 * t537 / 0.2e1;
	t650 = t68 * t593 + 0.2e1 * t322 * t616 * t121 - t620 * t509 / 0.2e1 + 0.2e1 * t322 * t514 - t625 * t514 / 0.2e1 - t62 * t614 - 0.2e1 * t598 * t596 * t116 + t602 * t531 / 0.2e1 - t606 * t537 / 0.2e1 + t610 * t537 / 0.2e1;
	unknown(1,1) = r_i_i_C(1) * (t71 * t70 + t76 * t75) + r_i_i_C(2) * (-t76 * t70 + t71 * t75) + r_i_i_C(3) * t84 + pkin(4) * t70 - t87 * t64 + pkin(9) * t1;
	unknown(1,2) = r_i_i_C(1) * (t71 * t183 + t76 * t202) + r_i_i_C(2) * (-t76 * t183 + t71 * t202) + pkin(4) * t183 + t87 * t162 - t213 * t165 / 0.2e1;
	unknown(1,3) = r_i_i_C(1) * t225 + r_i_i_C(2) * t229;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = r_i_i_C(1) * (t71 * t329 + t76 * t350) + r_i_i_C(2) * (-t76 * t329 + t71 * t350) + pkin(4) * t329 + t87 * t305 + 0.2e1 * t362 * t37 * t165 - t365 * t165 / 0.2e1;
	unknown(2,1) = pkin(4) * t219 - pkin(9) * t84 - r_i_i_C(1) * t229 + r_i_i_C(2) * t225 + r_i_i_C(3) * t1 + t87 * t165;
	unknown(2,2) = r_i_i_C(1) * (t71 * t399 + t76 * t418) + r_i_i_C(2) * (-t76 * t399 + t71 * t418) + pkin(4) * t399 + t87 * t388 - t213 * t64 / 0.2e1;
	unknown(2,3) = r_i_i_C(1) * (t76 * t70 - t71 * t75) + r_i_i_C(2) * (t71 * t70 + t76 * t75);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = r_i_i_C(1) * (t71 * t465 + t76 * t486) + r_i_i_C(2) * (-t76 * t465 + t71 * t486) + pkin(4) * t465 + t87 * t452 + 0.2e1 * t362 * t37 * t64 - t365 * t64 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = r_i_i_C(1) * (t71 * t541 + t76 * t561) + r_i_i_C(2) * (-t76 * t541 + t71 * t561) + pkin(4) * t541 + t48 * t31 * t101 - t132 * t114 * t572;
	unknown(3,3) = r_i_i_C(1) * (-t76 * t579 + t71 * t583) + r_i_i_C(2) * (-t71 * t579 - t76 * t583);
	unknown(3,4) = 0.0e0;
	unknown(3,5) = r_i_i_C(1) * (t71 * t628 + t76 * t650) + r_i_i_C(2) * (-t76 * t628 + t71 * t650) + pkin(4) * t628 + t48 * t31 * t245 + 0.2e1 * t254 * t48 * t37 * t27 - t281 * t114 * t572 / 0.2e1;
	Ja_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:52
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (110254->262), mult. (122036->630), div. (8989->16), fcn. (30917->14), ass. (0->224)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t28 = t27 * t1;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t28;
	t34 = t27 ^ 2;
	t35 = t29 * t34;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t29 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t50 = t49 * t48;
	t51 = 0.1e1 / t24;
	t52 = t51 * t29;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t58 = t51 * t29 * t55 - t16 * t52;
	t59 = sqrt(t58);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t50;
	t64 = t42 * t1;
	t65 = t32 * t64;
	t66 = t54 * t48;
	t67 = t60 * t49;
	t68 = t67 * t66;
	t70 = t62 * t33 - t68 * t65;
	t71 = sin(qJ(3));
	t75 = t68 * t33 + t62 * t65;
	t76 = cos(qJ(3));
	t78 = -t71 * t70 + t76 * t75;
	t79 = cos(qJ(4));
	t81 = cos(qJ(1));
	t82 = sin(qJ(4));
	t83 = t82 * t81;
	t87 = t79 * t81;
	t92 = t76 * t70 + t71 * t75;
	t96 = t48 * t31;
	t100 = 0.1e1 / t17;
	t101 = t100 * t4;
	t103 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t105 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t106 = -t103 - t105;
	t110 = t106 * t101 / 0.2e1 + 0.2e1 * t22 * t20;
	t112 = t32 * t110 * t81;
	t114 = t27 * t81;
	t116 = 0.1e1 / t23 / t22;
	t117 = t31 * t116;
	t118 = t117 * t114;
	t121 = t32 * t114;
	t123 = 0.1e1 / t47 / t46;
	t124 = t49 * t123;
	t125 = t29 * t27;
	t130 = t29 * t42;
	t131 = t100 * t19;
	t136 = -t106 * pkin(1) * t131 / 0.2e1 + 0.2e1 * t22 * t4;
	t141 = t110 * t37 * t125 - t37 * t116 * t34 - t37 * t116 * t43 + t136 * t37 * t130;
	t143 = 0.2e1 * t141 * t61 * t124;
	t146 = t60 * t100;
	t148 = t106 * t146 * t50;
	t152 = 0.1e1 / t59 / t58;
	t153 = t152 * t17;
	t159 = t29 * t54;
	t166 = 0.2e1 * t16 * t51 * t116 - 0.2e1 * t51 * t116 * t55 + 0.4e1 * t22 * t51 * t159 - t103 * t52 - t105 * t52;
	t168 = t166 * t153 * t50;
	t171 = t136 * t81;
	t172 = t32 * t171;
	t174 = t42 * t81;
	t175 = t117 * t174;
	t178 = t32 * t174;
	t179 = t54 * t123;
	t181 = 0.2e1 * t141 * t67 * t179;
	t185 = 0.2e1 * t67 * t22 * t48;
	t187 = t152 * t49;
	t189 = t166 * t187 * t66;
	t192 = -t62 * t112 + 0.2e1 * t62 * t118 + t143 * t121 / 0.2e1 - t148 * t121 / 0.2e1 + t168 * t121 / 0.2e1 + t68 * t172 - 0.2e1 * t68 * t175 - t181 * t178 / 0.2e1 + t185 * t178 - t189 * t178 / 0.2e1;
	t211 = -t68 * t112 + 0.2e1 * t68 * t118 + t181 * t121 / 0.2e1 - t185 * t121 + t189 * t121 / 0.2e1 - t62 * t172 + 0.2e1 * t62 * t175 + t143 * t178 / 0.2e1 - t148 * t178 / 0.2e1 + t168 * t178 / 0.2e1;
	t213 = -t71 * t192 + t76 * t211;
	t220 = t76 * t192 + t71 * t211;
	t225 = t123 * t31;
	t226 = 0.2e1 * t141 * t225;
	t232 = -t62 * t121 + t68 * t178;
	t236 = -t68 * t121 - t62 * t178;
	t238 = -t76 * t232 - t71 * t236;
	t245 = -t71 * t232 + t76 * t236;
	t251 = -t79 * t1 + t82 * t245;
	t255 = t82 * t1 + t79 * t245;
	t261 = pkin(1) * pkin(2);
	t263 = t15 * pkin(2) * t20 + t261 * t19 * t11;
	t267 = t19 ^ 2;
	t271 = -0.2e1 * pkin(2) * t267 * t7 + t263 * t101 - t25 * t3 - t40;
	t273 = t32 * t271 * t81;
	t276 = t48 * t37 * t29;
	t277 = t276 * t114;
	t278 = t17 * t49;
	t280 = pkin(2) * t20;
	t281 = t280 * t60 * t278;
	t288 = 0.1e1 / t36 / t30;
	t300 = -t263 * pkin(1) * t131 + pkin(1) * t17 * t2 - 0.2e1 * t261 * t19 * t4 - t26;
	t307 = 0.2e1 * t271 * t37 * t125 + 0.2e1 * t300 * t37 * t130 + 0.4e1 * t280 * t288 * t35 + 0.4e1 * t280 * t288 * t44;
	t309 = t307 * t61 * t124;
	t313 = 0.2e1 * t263 * t146 * t50;
	t326 = 0.2e1 * t15 * t261 * t19 * t52 + 0.2e1 * t280 * t11 * t52 + 0.4e1 * t280 * t51 * t159;
	t328 = t326 * t153 * t50;
	t331 = t300 * t81;
	t332 = t32 * t331;
	t334 = t276 * t174;
	t335 = t49 * t54;
	t337 = t280 * t60 * t335;
	t341 = t307 * t67 * t179;
	t344 = t48 * t32;
	t348 = t60 * t49 * pkin(2) * t20;
	t352 = t326 * t187 * t66;
	t355 = -t62 * t273 - 0.2e1 * t281 * t277 + t309 * t121 / 0.2e1 - t313 * t121 / 0.2e1 + t328 * t121 / 0.2e1 + t68 * t332 + 0.2e1 * t337 * t334 - t341 * t178 / 0.2e1 + 0.2e1 * t348 * t344 * t174 - t352 * t178 / 0.2e1;
	t376 = -t68 * t273 - 0.2e1 * t337 * t277 + t341 * t121 / 0.2e1 - 0.2e1 * t348 * t344 * t114 + t352 * t121 / 0.2e1 - t62 * t332 - 0.2e1 * t281 * t334 + t309 * t178 / 0.2e1 - t313 * t178 / 0.2e1 + t328 * t178 / 0.2e1;
	t378 = -t71 * t355 + t76 * t376;
	t385 = t76 * t355 + t71 * t376;
	t392 = t261 * t19 * t48;
	t395 = t307 * t225;
	t408 = t32 * t110 * t1;
	t410 = t117 * t28;
	t419 = t136 * t1;
	t420 = t32 * t419;
	t422 = t117 * t64;
	t430 = -t62 * t408 + 0.2e1 * t62 * t410 + t143 * t33 / 0.2e1 - t148 * t33 / 0.2e1 + t168 * t33 / 0.2e1 + t68 * t420 - 0.2e1 * t68 * t422 - t181 * t65 / 0.2e1 + t185 * t65 - t189 * t65 / 0.2e1;
	t449 = -t68 * t408 + 0.2e1 * t68 * t410 + t181 * t33 / 0.2e1 - t185 * t33 + t189 * t33 / 0.2e1 - t62 * t420 + 0.2e1 * t62 * t422 + t143 * t65 / 0.2e1 - t148 * t65 / 0.2e1 + t168 * t65 / 0.2e1;
	t451 = -t71 * t430 + t76 * t449;
	t458 = t76 * t430 + t71 * t449;
	t468 = t76 * t70 + t71 * t75;
	t475 = t71 * t70 - t76 * t75;
	t487 = t32 * t271 * t1;
	t489 = t276 * t28;
	t498 = t300 * t1;
	t499 = t32 * t498;
	t501 = t276 * t64;
	t511 = -t62 * t487 - 0.2e1 * t281 * t489 + t309 * t33 / 0.2e1 - t313 * t33 / 0.2e1 + t328 * t33 / 0.2e1 + t68 * t499 + 0.2e1 * t337 * t501 - t341 * t65 / 0.2e1 + 0.2e1 * t348 * t344 * t64 - t352 * t65 / 0.2e1;
	t532 = -t68 * t487 - 0.2e1 * t337 * t489 + t341 * t33 / 0.2e1 - 0.2e1 * t348 * t344 * t28 + t352 * t33 / 0.2e1 - t62 * t499 - 0.2e1 * t281 * t501 + t309 * t65 / 0.2e1 - t313 * t65 / 0.2e1 + t328 * t65 / 0.2e1;
	t534 = -t71 * t511 + t76 * t532;
	t541 = t76 * t511 + t71 * t532;
	t553 = t31 * t29 * t136;
	t556 = t31 * t116 * t42;
	t559 = t225 * t130;
	t560 = 0.2e1 * t141 * t60;
	t561 = t560 * t278;
	t564 = t96 * t130;
	t565 = t100 * t49;
	t567 = t106 * t60 * t565;
	t570 = t166 * t152;
	t571 = t570 * t278;
	t575 = t31 * t29 * t110;
	t578 = t31 * t116 * t27;
	t581 = t225 * t125;
	t582 = t560 * t335;
	t585 = t31 * t125;
	t587 = t96 * t125;
	t588 = t570 * t335;
	t591 = t62 * t553 - 0.2e1 * t62 * t556 - t561 * t559 / 0.2e1 + t567 * t564 / 0.2e1 - t571 * t564 / 0.2e1 + t68 * t575 - 0.2e1 * t68 * t578 - t582 * t581 / 0.2e1 + t185 * t585 - t588 * t587 / 0.2e1;
	t598 = t31 * t130;
	t611 = t68 * t553 - 0.2e1 * t68 * t556 - t582 * t559 / 0.2e1 + t185 * t598 - t588 * t564 / 0.2e1 - t62 * t575 + 0.2e1 * t62 * t578 + t561 * t581 / 0.2e1 - t567 * t587 / 0.2e1 + t571 * t587 / 0.2e1;
	t613 = -t71 * t591 + t76 * t611;
	t620 = t76 * t591 + t71 * t611;
	t626 = t31 * t27;
	t633 = t68 * t585 + t62 * t598;
	t637 = -t62 * t585 + t68 * t598;
	t639 = -t76 * t633 - t71 * t637;
	t646 = -t71 * t633 + t76 * t637;
	t656 = t31 * t29 * t300;
	t658 = t48 * t37;
	t659 = t49 * t658;
	t661 = t280 * t61;
	t664 = t307 * t60;
	t665 = t664 * t278;
	t669 = 0.2e1 * t263 * t60 * t565;
	t672 = t326 * t152;
	t673 = t672 * t278;
	t677 = t31 * t29 * t271;
	t679 = t54 * t658;
	t683 = t664 * t335;
	t688 = t672 * t335;
	t691 = t62 * t656 + 0.2e1 * t661 * t659 * t130 - t665 * t559 / 0.2e1 + t669 * t564 / 0.2e1 - t673 * t564 / 0.2e1 + t68 * t677 + 0.2e1 * t348 * t679 * t125 - t683 * t581 / 0.2e1 + 0.2e1 * t348 * t587 - t688 * t587 / 0.2e1;
	t713 = t68 * t656 + 0.2e1 * t348 * t679 * t130 - t683 * t559 / 0.2e1 + 0.2e1 * t348 * t564 - t688 * t564 / 0.2e1 - t62 * t677 - 0.2e1 * t661 * t659 * t125 + t665 * t581 / 0.2e1 - t669 * t587 / 0.2e1 + t673 * t587 / 0.2e1;
	t715 = -t71 * t691 + t76 * t713;
	t722 = t76 * t691 + t71 * t713;
	unknown(1,1) = r_i_i_C(1) * (-t79 * t78 - t83) + r_i_i_C(2) * (t82 * t78 - t87) + r_i_i_C(3) * t92 + pkin(5) * t92 + pkin(4) * t70 - t96 * t64 + pkin(9) * t1;
	unknown(1,2) = -r_i_i_C(1) * t79 * t213 + r_i_i_C(2) * t82 * t213 + r_i_i_C(3) * t220 + pkin(5) * t220 + pkin(4) * t192 + t96 * t171 - t226 * t174 / 0.2e1;
	unknown(1,3) = -r_i_i_C(1) * t79 * t238 + r_i_i_C(2) * t82 * t238 + pkin(5) * t245 + r_i_i_C(3) * t245;
	unknown(1,4) = r_i_i_C(1) * t251 + r_i_i_C(2) * t255;
	unknown(1,5) = -r_i_i_C(1) * t79 * t378 + r_i_i_C(2) * t82 * t378 + r_i_i_C(3) * t385 + pkin(5) * t385 + pkin(4) * t355 + t96 * t331 + 0.2e1 * t392 * t37 * t174 - t395 * t174 / 0.2e1;
	unknown(2,1) = pkin(4) * t232 - pkin(5) * t238 - pkin(9) * t81 - r_i_i_C(1) * t255 + r_i_i_C(2) * t251 - r_i_i_C(3) * t238 + t96 * t174;
	unknown(2,2) = -r_i_i_C(1) * t79 * t451 + r_i_i_C(2) * t82 * t451 + r_i_i_C(3) * t458 + pkin(5) * t458 + pkin(4) * t430 + t96 * t419 - t226 * t64 / 0.2e1;
	unknown(2,3) = -r_i_i_C(1) * t79 * t468 + r_i_i_C(2) * t82 * t468 + pkin(5) * t475 + r_i_i_C(3) * t475;
	unknown(2,4) = r_i_i_C(1) * (t82 * t475 + t87) + r_i_i_C(2) * (t79 * t475 - t83);
	unknown(2,5) = -r_i_i_C(1) * t79 * t534 + r_i_i_C(2) * t82 * t534 + r_i_i_C(3) * t541 + pkin(5) * t541 + pkin(4) * t511 + t96 * t498 + 0.2e1 * t392 * t37 * t64 - t395 * t64 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -r_i_i_C(1) * t79 * t613 + r_i_i_C(2) * t82 * t613 + t48 * t31 * t110 - t141 * t123 * t626 + pkin(4) * t591 + pkin(5) * t620 + r_i_i_C(3) * t620;
	unknown(3,3) = -r_i_i_C(1) * t79 * t639 + r_i_i_C(2) * t82 * t639 + pkin(5) * t646 + r_i_i_C(3) * t646;
	unknown(3,4) = r_i_i_C(1) * t82 * t646 + r_i_i_C(2) * t79 * t646;
	unknown(3,5) = -r_i_i_C(1) * t79 * t715 + r_i_i_C(2) * t82 * t715 + r_i_i_C(3) * t722 + pkin(5) * t722 + pkin(4) * t691 + t48 * t31 * t271 + 0.2e1 * t280 * t48 * t37 * t27 - t307 * t123 * t626 / 0.2e1;
	Ja_transl = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (22->18), div. (0->0), fcn. (22->4), ass. (0->23)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t1 * t2;
	t5 = cos(qJ(5));
	t6 = t1 * t5;
	t8 = cos(qJ(1));
	t12 = t8 * t5;
	t14 = t8 * t2;
	unknown(1,1) = -t1 * pkin(8) + t3 * r_i_i_C(1) + t6 * r_i_i_C(2) + t8 * r_i_i_C(3);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t12 * r_i_i_C(1) + t14 * r_i_i_C(2);
	unknown(2,1) = t8 * pkin(8) - t14 * r_i_i_C(1) - t12 * r_i_i_C(2) + t1 * r_i_i_C(3);
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t6 * r_i_i_C(1) + t3 * r_i_i_C(2);
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t2 * r_i_i_C(1) - t5 * r_i_i_C(2);
	Ja_transl = unknown;
else
	Ja_transl=NaN(3,5);
end