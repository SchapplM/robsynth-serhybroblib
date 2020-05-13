% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh4m1DE2
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
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh4m1DE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh4m1DE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:31
	% EndTime: 2020-04-11 22:53:32
	% DurationCPUTime: 0.04s
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
	% StartTime: 2020-04-11 22:53:31
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
	% DurationCPUTime: 0.25s
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.45s
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (35027->251), mult. (37574->635), div. (2621->20), fcn. (9889->13), ass. (0->225)
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
	t50 = 0.1e1 / pkin(3);
	t51 = t50 * t29;
	t53 = t6 - t7 - t21 + t23 + t24;
	t56 = atan2(t17 * t51, t50 * t29 * t53);
	t57 = t56 + qJ(3);
	t58 = sin(t57);
	t59 = t58 * t49;
	t61 = t42 * t1;
	t62 = t29 * t61;
	t63 = cos(t57);
	t64 = t63 * t49;
	t66 = t59 * t30 - t64 * t62;
	t70 = t64 * t30 + t59 * t62;
	t72 = cos(qJ(1));
	t74 = t32 * t34;
	t77 = 0.1e1 / t24;
	t78 = t77 * t34;
	t80 = t53 ^ 2;
	t83 = t77 * t34 * t80 - t16 * t78;
	t84 = sqrt(t83);
	t85 = 0.1e1 / t84;
	t87 = pkin(4) * t85 * t17;
	t88 = t87 * t50 * t48;
	t90 = t74 * t61;
	t93 = pkin(4) * t85 * t50;
	t94 = t93 * t53 * t48;
	t99 = 0.1e1 / t17;
	t100 = t99 * t4;
	t102 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t104 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t105 = -t102 - t104;
	t109 = t105 * t100 / 0.2e1 + 0.2e1 * t22 * t20;
	t110 = t109 * t72;
	t111 = t29 * t110;
	t113 = t27 * t72;
	t114 = t34 * t113;
	t116 = t29 * t113;
	t118 = 0.1e1 / t47 / t46;
	t119 = t118 * t32;
	t120 = t34 * t27;
	t124 = 0.1e1 / t23 / t22;
	t127 = t34 * t42;
	t128 = t99 * t19;
	t133 = -t105 * pkin(1) * t128 / 0.2e1 + 0.2e1 * t22 * t4;
	t138 = t109 * t37 * t120 - t37 * t124 * t33 - t37 * t124 * t43 + t133 * t37 * t127;
	t140 = 0.2e1 * t138 * t58 * t119;
	t149 = 0.1e1 / t53;
	t152 = 0.1e1 / t80;
	t155 = 0.1e1 / (-t152 * t16 + 0.1e1);
	t160 = t34 * t53;
	t166 = t155 * t152 * t17;
	t168 = t155 * pkin(3) * t22 * t149 * (-t17 * t50 * t34 + t105 * t99 * t51 / 0.2e1) - t166 * pkin(3) * t22 * (0.2e1 * t50 * t29 * t22 - t50 * t160);
	t170 = t63 * t168 * t49;
	t172 = t133 * t72;
	t173 = t29 * t172;
	t175 = t42 * t72;
	t176 = t34 * t175;
	t178 = t29 * t175;
	t180 = 0.2e1 * t138 * t63 * t119;
	t184 = t58 * t168 * t49;
	t202 = t32 * t124;
	t206 = t118 * t74;
	t207 = t206 * t113;
	t208 = t17 * t50;
	t209 = pkin(4) * t85;
	t210 = 0.2e1 * t138 * t209;
	t211 = t210 * t208;
	t214 = t48 * t74;
	t215 = t214 * t113;
	t216 = t99 * t50;
	t218 = t105 * t209 * t216;
	t223 = pkin(4) / t84 / t83;
	t236 = (0.2e1 * t16 * t77 * t124 - 0.2e1 * t77 * t124 * t80 + 0.4e1 * t22 * t77 * t160 - t102 * t78 - t104 * t78) * t223;
	t237 = t236 * t208;
	t245 = t206 * t175;
	t246 = t50 * t53;
	t247 = t210 * t246;
	t250 = t74 * t175;
	t252 = 0.2e1 * t93 * t22 * t48;
	t254 = t214 * t175;
	t255 = t236 * t246;
	t259 = 0.2e1 * t138 * t119;
	t262 = r_i_i_C(1) * (-t59 * t111 + t59 * t114 + t140 * t116 / 0.2e1 - t170 * t116 + t64 * t173 - t64 * t176 - t180 * t178 / 0.2e1 - t184 * t178) + r_i_i_C(2) * (-t64 * t111 + t64 * t114 + t180 * t116 / 0.2e1 + t184 * t116 - t59 * t173 + t59 * t176 + t140 * t178 / 0.2e1 - t170 * t178) - t88 * t74 * t110 + 0.2e1 * t88 * t202 * t113 + t211 * t207 / 0.2e1 - t218 * t215 / 0.2e1 + t237 * t215 / 0.2e1 + t94 * t74 * t172 - 0.2e1 * t94 * t202 * t175 - t247 * t245 / 0.2e1 + t252 * t250 - t255 * t254 / 0.2e1 + t49 * t172 - t259 * t175 / 0.2e1;
	t265 = -t64 * t116 - t59 * t178;
	t269 = t59 * t116 - t64 * t178;
	t275 = pkin(1) * pkin(2);
	t277 = t15 * pkin(2) * t20 + t275 * t19 * t11;
	t281 = t19 ^ 2;
	t285 = -0.2e1 * pkin(2) * t281 * t7 + t277 * t100 - t25 * t3 - t40;
	t286 = t285 * t72;
	t287 = t29 * t286;
	t289 = t37 * t29;
	t290 = t289 * t113;
	t292 = pkin(2) * t20;
	t293 = t292 * t58 * t48;
	t300 = 0.1e1 / t36 / t31;
	t312 = -t277 * pkin(1) * t128 + pkin(1) * t17 * t2 - 0.2e1 * t275 * t19 * t4 - t26;
	t319 = 0.2e1 * t285 * t37 * t120 + 0.2e1 * t312 * t37 * t127 + 0.4e1 * t292 * t300 * t35 + 0.4e1 * t292 * t300 * t44;
	t321 = t319 * t58 * t119;
	t330 = t155 * t149 * t277 * t99 - 0.2e1 * t166 * t292;
	t332 = t63 * t330 * t49;
	t334 = t312 * t72;
	t335 = t29 * t334;
	t337 = t289 * t175;
	t339 = t292 * t63 * t48;
	t343 = t319 * t63 * t119;
	t347 = t58 * t330 * t49;
	t367 = t48 * t37;
	t368 = t50 * t367;
	t370 = t292 * t87;
	t373 = t319 * t209;
	t374 = t373 * t208;
	t378 = 0.2e1 * t277 * t209 * t216;
	t392 = (0.2e1 * t15 * t275 * t19 * t78 + 0.2e1 * t292 * t11 * t78 + 0.4e1 * t292 * t77 * t160) * t223;
	t393 = t392 * t208;
	t398 = t53 * t367;
	t400 = t292 * t93;
	t403 = t373 * t246;
	t408 = t392 * t246;
	t414 = t275 * t19 * t48;
	t417 = t319 * t119;
	t420 = r_i_i_C(1) * (-t59 * t287 - 0.2e1 * t293 * t290 + t321 * t116 / 0.2e1 - t332 * t116 + t64 * t335 + 0.2e1 * t339 * t337 - t343 * t178 / 0.2e1 - t347 * t178) + r_i_i_C(2) * (-t64 * t287 - 0.2e1 * t339 * t290 + t343 * t116 / 0.2e1 + t347 * t116 - t59 * t335 - 0.2e1 * t293 * t337 + t321 * t178 / 0.2e1 - t332 * t178) - t88 * t74 * t286 - 0.2e1 * t370 * t368 * t114 + t374 * t207 / 0.2e1 - t378 * t215 / 0.2e1 + t393 * t215 / 0.2e1 + t94 * t74 * t334 + 0.2e1 * t400 * t398 * t176 - t403 * t245 / 0.2e1 + 0.2e1 * t400 * t254 - t408 * t254 / 0.2e1 + t49 * t334 + 0.2e1 * t414 * t37 * t175 - t417 * t175 / 0.2e1;
	t430 = t109 * t1;
	t431 = t29 * t430;
	t433 = t34 * t28;
	t438 = t133 * t1;
	t439 = t29 * t438;
	t441 = t34 * t61;
	t465 = t206 * t28;
	t468 = t214 * t28;
	t478 = t206 * t61;
	t482 = t214 * t61;
	t488 = r_i_i_C(1) * (-t59 * t431 + t59 * t433 + t140 * t30 / 0.2e1 - t170 * t30 + t64 * t439 - t64 * t441 - t180 * t62 / 0.2e1 - t184 * t62) + r_i_i_C(2) * (-t64 * t431 + t64 * t433 + t180 * t30 / 0.2e1 + t184 * t30 - t59 * t439 + t59 * t441 + t140 * t62 / 0.2e1 - t170 * t62) - t88 * t74 * t430 + 0.2e1 * t88 * t202 * t28 + t211 * t465 / 0.2e1 - t218 * t468 / 0.2e1 + t237 * t468 / 0.2e1 + t94 * t74 * t438 - 0.2e1 * t94 * t202 * t61 - t247 * t478 / 0.2e1 + t252 * t90 - t255 * t482 / 0.2e1 + t49 * t438 - t259 * t61 / 0.2e1;
	t492 = t285 * t1;
	t493 = t29 * t492;
	t495 = t289 * t28;
	t501 = t312 * t1;
	t502 = t29 * t501;
	t504 = t289 * t61;
	t554 = r_i_i_C(1) * (-t59 * t493 - 0.2e1 * t293 * t495 + t321 * t30 / 0.2e1 - t332 * t30 + t64 * t502 + 0.2e1 * t339 * t504 - t343 * t62 / 0.2e1 - t347 * t62) + r_i_i_C(2) * (-t64 * t493 - 0.2e1 * t339 * t495 + t343 * t30 / 0.2e1 + t347 * t30 - t59 * t502 - 0.2e1 * t293 * t504 + t321 * t62 / 0.2e1 - t332 * t62) - t88 * t74 * t492 - 0.2e1 * t370 * t368 * t433 + t374 * t465 / 0.2e1 - t378 * t468 / 0.2e1 + t393 * t468 / 0.2e1 + t94 * t74 * t501 + 0.2e1 * t400 * t398 * t441 - t403 * t478 / 0.2e1 + 0.2e1 * t400 * t482 - t408 * t482 / 0.2e1 + t49 * t501 + 0.2e1 * t414 * t37 * t61 - t417 * t61 / 0.2e1;
	t555 = t29 * t133;
	t558 = t29 * t42;
	t559 = t32 * t558;
	t560 = t58 * t118;
	t561 = 0.2e1 * t138 * t560;
	t564 = t168 * t48;
	t565 = t63 * t564;
	t567 = t29 * t109;
	t570 = t29 * t27;
	t571 = t32 * t570;
	t572 = t63 * t118;
	t573 = 0.2e1 * t138 * t572;
	t576 = t58 * t564;
	t594 = t209 * t208;
	t600 = t119 * t127;
	t603 = t49 * t127;
	t610 = t209 * t246;
	t616 = t119 * t120;
	t619 = t49 * t120;
	t627 = t32 * t27;
	t631 = r_i_i_C(1) * (t59 * t555 - t59 * t127 - t561 * t559 / 0.2e1 + t565 * t559 + t64 * t567 - t64 * t120 - t573 * t571 / 0.2e1 - t576 * t571) + r_i_i_C(2) * (t64 * t555 - t64 * t127 - t573 * t559 / 0.2e1 - t576 * t559 - t59 * t567 + t59 * t120 + t561 * t571 / 0.2e1 - t565 * t571) + t594 * t49 * t34 * t133 - 0.2e1 * t594 * t49 * t124 * t42 - t211 * t600 / 0.2e1 + t218 * t603 / 0.2e1 - t237 * t603 / 0.2e1 + t610 * t49 * t34 * t109 - 0.2e1 * t610 * t49 * t124 * t27 - t247 * t616 / 0.2e1 + 0.2e1 * t209 * t50 * t22 * t619 - t255 * t619 / 0.2e1 + t48 * t32 * t109 - t138 * t118 * t627;
	t641 = t29 * t312;
	t643 = t367 * t558;
	t645 = t275 * t19 * t58;
	t648 = t319 * t560;
	t651 = t330 * t48;
	t652 = t63 * t651;
	t654 = t29 * t285;
	t656 = t367 * t570;
	t658 = t275 * t19 * t63;
	t661 = t319 * t572;
	t664 = t58 * t651;
	t718 = r_i_i_C(1) * (t59 * t641 + 0.2e1 * t645 * t643 - t648 * t559 / 0.2e1 + t652 * t559 + t64 * t654 + 0.2e1 * t658 * t656 - t661 * t571 / 0.2e1 - t664 * t571) + r_i_i_C(2) * (t64 * t641 + 0.2e1 * t658 * t643 - t661 * t559 / 0.2e1 - t664 * t559 - t59 * t654 - 0.2e1 * t645 * t656 + t648 * t571 / 0.2e1 - t652 * t571) + t594 * t49 * t34 * t312 + 0.2e1 * t370 * t368 * t127 - t374 * t600 / 0.2e1 + t378 * t603 / 0.2e1 - t393 * t603 / 0.2e1 + t610 * t49 * t34 * t285 + 0.2e1 * t400 * t398 * t120 - t403 * t616 / 0.2e1 + 0.2e1 * t93 * t275 * t19 * t49 * t120 - t408 * t619 / 0.2e1 + t48 * t32 * t285 + 0.2e1 * t292 * t48 * t37 * t27 - t319 * t118 * t627 / 0.2e1;
	unknown(1,1) = t88 * t74 * t28 + pkin(9) * t1 + r_i_i_C(1) * t66 + r_i_i_C(2) * t70 + r_i_i_C(3) * t72 - t49 * t61 - t94 * t90;
	unknown(1,2) = t262;
	unknown(1,3) = r_i_i_C(1) * t265 + r_i_i_C(2) * t269;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t420;
	unknown(2,1) = -t88 * t74 * t113 - pkin(9) * t72 - r_i_i_C(1) * t269 + r_i_i_C(2) * t265 + r_i_i_C(3) * t1 + t49 * t175 + t94 * t250;
	unknown(2,2) = t488;
	unknown(2,3) = -r_i_i_C(1) * t70 + r_i_i_C(2) * t66;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t554;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t631;
	unknown(3,3) = r_i_i_C(1) * (t64 * t558 - t59 * t570) + r_i_i_C(2) * (-t59 * t558 - t64 * t570);
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t718;
	Ja_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (57542->280), mult. (61156->691), div. (4333->20), fcn. (16325->15), ass. (0->245)
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
	t50 = 0.1e1 / pkin(3);
	t51 = t50 * t29;
	t53 = t6 - t7 - t21 + t23 + t24;
	t56 = atan2(t17 * t51, t50 * t29 * t53);
	t57 = t56 + qJ(3);
	t58 = cos(t57);
	t59 = t58 * t49;
	t61 = t42 * t1;
	t62 = t29 * t61;
	t63 = sin(t57);
	t64 = t63 * t49;
	t66 = t59 * t30 + t64 * t62;
	t67 = cos(qJ(4));
	t69 = cos(qJ(1));
	t70 = sin(qJ(4));
	t71 = t70 * t69;
	t75 = t67 * t69;
	t80 = t64 * t30 - t59 * t62;
	t83 = t32 * t34;
	t86 = 0.1e1 / t24;
	t87 = t86 * t34;
	t89 = t53 ^ 2;
	t92 = t86 * t34 * t89 - t16 * t87;
	t93 = sqrt(t92);
	t94 = 0.1e1 / t93;
	t96 = pkin(4) * t94 * t17;
	t97 = t96 * t50 * t48;
	t99 = t83 * t61;
	t102 = pkin(4) * t94 * t50;
	t103 = t102 * t53 * t48;
	t108 = 0.1e1 / t17;
	t109 = t108 * t4;
	t111 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t113 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t114 = -t111 - t113;
	t118 = t114 * t109 / 0.2e1 + 0.2e1 * t22 * t20;
	t119 = t118 * t69;
	t120 = t29 * t119;
	t122 = t27 * t69;
	t123 = t34 * t122;
	t125 = t29 * t122;
	t127 = 0.1e1 / t47 / t46;
	t128 = t127 * t32;
	t129 = t34 * t27;
	t133 = 0.1e1 / t23 / t22;
	t136 = t34 * t42;
	t137 = t108 * t19;
	t142 = -t114 * pkin(1) * t137 / 0.2e1 + 0.2e1 * t22 * t4;
	t147 = t118 * t37 * t129 - t37 * t133 * t33 - t37 * t133 * t43 + t142 * t37 * t136;
	t149 = 0.2e1 * t147 * t58 * t128;
	t158 = 0.1e1 / t53;
	t161 = 0.1e1 / t89;
	t164 = 0.1e1 / (-t161 * t16 + 0.1e1);
	t169 = t34 * t53;
	t175 = t164 * t161 * t17;
	t177 = t164 * pkin(3) * t22 * t158 * (-t17 * t50 * t34 + t114 * t108 * t51 / 0.2e1) - t175 * pkin(3) * t22 * (0.2e1 * t50 * t29 * t22 - t50 * t169);
	t179 = t63 * t177 * t49;
	t181 = t142 * t69;
	t182 = t29 * t181;
	t184 = t42 * t69;
	t185 = t34 * t184;
	t187 = t29 * t184;
	t189 = 0.2e1 * t147 * t63 * t128;
	t193 = t58 * t177 * t49;
	t195 = -t59 * t120 + t59 * t123 + t149 * t125 / 0.2e1 + t179 * t125 - t64 * t182 + t64 * t185 + t189 * t187 / 0.2e1 - t193 * t187;
	t210 = -t64 * t120 + t64 * t123 + t189 * t125 / 0.2e1 - t193 * t125 + t59 * t182 - t59 * t185 - t149 * t187 / 0.2e1 - t179 * t187;
	t215 = t32 * t133;
	t219 = t127 * t83;
	t220 = t219 * t122;
	t221 = t17 * t50;
	t222 = pkin(4) * t94;
	t223 = 0.2e1 * t147 * t222;
	t224 = t223 * t221;
	t227 = t48 * t83;
	t228 = t227 * t122;
	t229 = t108 * t50;
	t231 = t114 * t222 * t229;
	t236 = pkin(4) / t93 / t92;
	t249 = (0.2e1 * t16 * t86 * t133 - 0.2e1 * t86 * t133 * t89 + 0.4e1 * t22 * t86 * t169 - t111 * t87 - t113 * t87) * t236;
	t250 = t249 * t221;
	t258 = t219 * t184;
	t259 = t50 * t53;
	t260 = t223 * t259;
	t263 = t83 * t184;
	t265 = 0.2e1 * t102 * t22 * t48;
	t267 = t227 * t184;
	t268 = t249 * t259;
	t272 = 0.2e1 * t147 * t128;
	t275 = -r_i_i_C(1) * t67 * t195 + r_i_i_C(2) * t70 * t195 + r_i_i_C(3) * t210 + pkin(5) * t210 - t97 * t83 * t119 + 0.2e1 * t97 * t215 * t122 + t224 * t220 / 0.2e1 - t231 * t228 / 0.2e1 + t250 * t228 / 0.2e1 + t103 * t83 * t181 - 0.2e1 * t103 * t215 * t184 - t260 * t258 / 0.2e1 + t265 * t263 - t268 * t267 / 0.2e1 + t49 * t181 - t272 * t184 / 0.2e1;
	t278 = t64 * t125 - t59 * t187;
	t285 = -t59 * t125 - t64 * t187;
	t291 = -t67 * t1 + t70 * t285;
	t295 = t70 * t1 + t67 * t285;
	t301 = pkin(1) * pkin(2);
	t303 = t15 * pkin(2) * t20 + t301 * t19 * t11;
	t307 = t19 ^ 2;
	t311 = -0.2e1 * pkin(2) * t307 * t7 + t303 * t109 - t25 * t3 - t40;
	t312 = t311 * t69;
	t313 = t29 * t312;
	t315 = t37 * t29;
	t316 = t315 * t122;
	t318 = pkin(2) * t20;
	t319 = t318 * t58 * t48;
	t326 = 0.1e1 / t36 / t31;
	t338 = -t303 * pkin(1) * t137 + pkin(1) * t17 * t2 - 0.2e1 * t301 * t19 * t4 - t26;
	t345 = 0.2e1 * t311 * t37 * t129 + 0.2e1 * t338 * t37 * t136 + 0.4e1 * t318 * t326 * t35 + 0.4e1 * t318 * t326 * t44;
	t347 = t345 * t58 * t128;
	t356 = t164 * t158 * t303 * t108 - 0.2e1 * t175 * t318;
	t358 = t63 * t356 * t49;
	t360 = t338 * t69;
	t361 = t29 * t360;
	t363 = t315 * t184;
	t365 = t318 * t63 * t48;
	t369 = t345 * t63 * t128;
	t373 = t58 * t356 * t49;
	t375 = -t59 * t313 - 0.2e1 * t319 * t316 + t347 * t125 / 0.2e1 + t358 * t125 - t64 * t361 - 0.2e1 * t365 * t363 + t369 * t187 / 0.2e1 - t373 * t187;
	t392 = -t64 * t313 - 0.2e1 * t365 * t316 + t369 * t125 / 0.2e1 - t373 * t125 + t59 * t361 + 0.2e1 * t319 * t363 - t347 * t187 / 0.2e1 - t358 * t187;
	t397 = t48 * t37;
	t398 = t50 * t397;
	t400 = t318 * t96;
	t403 = t345 * t222;
	t404 = t403 * t221;
	t408 = 0.2e1 * t303 * t222 * t229;
	t422 = (0.2e1 * t15 * t301 * t19 * t87 + 0.2e1 * t318 * t11 * t87 + 0.4e1 * t318 * t86 * t169) * t236;
	t423 = t422 * t221;
	t428 = t53 * t397;
	t430 = t318 * t102;
	t433 = t403 * t259;
	t438 = t422 * t259;
	t444 = t301 * t19 * t48;
	t447 = t345 * t128;
	t450 = -r_i_i_C(1) * t67 * t375 + r_i_i_C(2) * t70 * t375 + r_i_i_C(3) * t392 + pkin(5) * t392 - t97 * t83 * t312 - 0.2e1 * t400 * t398 * t123 + t404 * t220 / 0.2e1 - t408 * t228 / 0.2e1 + t423 * t228 / 0.2e1 + t103 * t83 * t360 + 0.2e1 * t430 * t428 * t185 - t433 * t258 / 0.2e1 + 0.2e1 * t430 * t267 - t438 * t267 / 0.2e1 + t49 * t360 + 0.2e1 * t444 * t37 * t184 - t447 * t184 / 0.2e1;
	t461 = t118 * t1;
	t462 = t29 * t461;
	t464 = t34 * t28;
	t469 = t142 * t1;
	t470 = t29 * t469;
	t472 = t34 * t61;
	t477 = -t59 * t462 + t59 * t464 + t149 * t30 / 0.2e1 + t179 * t30 - t64 * t470 + t64 * t472 + t189 * t62 / 0.2e1 - t193 * t62;
	t492 = -t64 * t462 + t64 * t464 + t189 * t30 / 0.2e1 - t193 * t30 + t59 * t470 - t59 * t472 - t149 * t62 / 0.2e1 - t179 * t62;
	t500 = t219 * t28;
	t503 = t227 * t28;
	t513 = t219 * t61;
	t517 = t227 * t61;
	t523 = -r_i_i_C(1) * t67 * t477 + r_i_i_C(2) * t70 * t477 + r_i_i_C(3) * t492 + pkin(5) * t492 - t97 * t83 * t461 + 0.2e1 * t97 * t215 * t28 + t224 * t500 / 0.2e1 - t231 * t503 / 0.2e1 + t250 * t503 / 0.2e1 + t103 * t83 * t469 - 0.2e1 * t103 * t215 * t61 - t260 * t513 / 0.2e1 + t265 * t99 - t268 * t517 / 0.2e1 + t49 * t469 - t272 * t61 / 0.2e1;
	t538 = t311 * t1;
	t539 = t29 * t538;
	t541 = t315 * t28;
	t547 = t338 * t1;
	t548 = t29 * t547;
	t550 = t315 * t61;
	t556 = -t59 * t539 - 0.2e1 * t319 * t541 + t347 * t30 / 0.2e1 + t358 * t30 - t64 * t548 - 0.2e1 * t365 * t550 + t369 * t62 / 0.2e1 - t373 * t62;
	t573 = -t64 * t539 - 0.2e1 * t365 * t541 + t369 * t30 / 0.2e1 - t373 * t30 + t59 * t548 + 0.2e1 * t319 * t550 - t347 * t62 / 0.2e1 - t358 * t62;
	t604 = -r_i_i_C(1) * t67 * t556 + r_i_i_C(2) * t70 * t556 + r_i_i_C(3) * t573 + pkin(5) * t573 - t97 * t83 * t538 - 0.2e1 * t400 * t398 * t464 + t404 * t500 / 0.2e1 - t408 * t503 / 0.2e1 + t423 * t503 / 0.2e1 + t103 * t83 * t547 + 0.2e1 * t430 * t428 * t472 - t433 * t513 / 0.2e1 + 0.2e1 * t430 * t517 - t438 * t517 / 0.2e1 + t49 * t547 + 0.2e1 * t444 * t37 * t61 - t447 * t61 / 0.2e1;
	t605 = t29 * t142;
	t608 = t29 * t42;
	t609 = t32 * t608;
	t610 = t58 * t127;
	t611 = 0.2e1 * t147 * t610;
	t614 = t177 * t48;
	t615 = t63 * t614;
	t617 = t29 * t118;
	t620 = t29 * t27;
	t621 = t32 * t620;
	t622 = t63 * t127;
	t623 = 0.2e1 * t147 * t622;
	t626 = t58 * t614;
	t628 = t59 * t605 - t59 * t136 - t611 * t609 / 0.2e1 - t615 * t609 - t64 * t617 + t64 * t129 + t623 * t621 / 0.2e1 - t626 * t621;
	t643 = t64 * t605 - t64 * t136 - t623 * t609 / 0.2e1 + t626 * t609 + t59 * t617 - t59 * t129 - t611 * t621 / 0.2e1 - t615 * t621;
	t648 = t222 * t221;
	t654 = t128 * t136;
	t657 = t49 * t136;
	t664 = t222 * t259;
	t670 = t128 * t129;
	t673 = t49 * t129;
	t681 = t32 * t27;
	t685 = -r_i_i_C(1) * t67 * t628 + r_i_i_C(2) * t70 * t628 + r_i_i_C(3) * t643 + pkin(5) * t643 + t648 * t49 * t34 * t142 - 0.2e1 * t648 * t49 * t133 * t42 - t224 * t654 / 0.2e1 + t231 * t657 / 0.2e1 - t250 * t657 / 0.2e1 + t664 * t49 * t34 * t118 - 0.2e1 * t664 * t49 * t133 * t27 - t260 * t670 / 0.2e1 + 0.2e1 * t222 * t50 * t22 * t673 - t268 * t673 / 0.2e1 + t48 * t32 * t118 - t147 * t127 * t681;
	t688 = -t59 * t620 - t64 * t608;
	t695 = t59 * t608 - t64 * t620;
	t704 = t29 * t338;
	t706 = t397 * t608;
	t708 = t301 * t19 * t58;
	t711 = t345 * t610;
	t714 = t356 * t48;
	t715 = t63 * t714;
	t717 = t29 * t311;
	t719 = t397 * t620;
	t721 = t301 * t19 * t63;
	t724 = t345 * t622;
	t727 = t58 * t714;
	t729 = t59 * t704 + 0.2e1 * t708 * t706 - t711 * t609 / 0.2e1 - t715 * t609 - t64 * t717 - 0.2e1 * t721 * t719 + t724 * t621 / 0.2e1 - t727 * t621;
	t746 = t64 * t704 + 0.2e1 * t721 * t706 - t724 * t609 / 0.2e1 + t727 * t609 + t59 * t717 + 0.2e1 * t708 * t719 - t711 * t621 / 0.2e1 - t715 * t621;
	t785 = -r_i_i_C(1) * t67 * t729 + r_i_i_C(2) * t70 * t729 + r_i_i_C(3) * t746 + pkin(5) * t746 + t648 * t49 * t34 * t338 + 0.2e1 * t400 * t398 * t136 - t404 * t654 / 0.2e1 + t408 * t657 / 0.2e1 - t423 * t657 / 0.2e1 + t664 * t49 * t34 * t311 + 0.2e1 * t430 * t428 * t129 - t433 * t670 / 0.2e1 + 0.2e1 * t102 * t301 * t19 * t49 * t129 - t438 * t673 / 0.2e1 + t48 * t32 * t311 + 0.2e1 * t318 * t48 * t37 * t27 - t345 * t127 * t681 / 0.2e1;
	unknown(1,1) = r_i_i_C(1) * (-t67 * t66 - t71) + r_i_i_C(2) * (t70 * t66 - t75) + r_i_i_C(3) * t80 + pkin(5) * t80 + t97 * t83 * t28 - t103 * t99 - t49 * t61 + pkin(9) * t1;
	unknown(1,2) = t275;
	unknown(1,3) = -r_i_i_C(1) * t67 * t278 + r_i_i_C(2) * t70 * t278 + pkin(5) * t285 + r_i_i_C(3) * t285;
	unknown(1,4) = r_i_i_C(1) * t291 + r_i_i_C(2) * t295;
	unknown(1,5) = t450;
	unknown(2,1) = -t97 * t83 * t122 - pkin(5) * t278 - pkin(9) * t69 - r_i_i_C(1) * t295 + r_i_i_C(2) * t291 - r_i_i_C(3) * t278 + t103 * t263 + t49 * t184;
	unknown(2,2) = t523;
	unknown(2,3) = -r_i_i_C(1) * t67 * t80 + r_i_i_C(2) * t70 * t80 - pkin(5) * t66 - r_i_i_C(3) * t66;
	unknown(2,4) = r_i_i_C(1) * (-t70 * t66 + t75) + r_i_i_C(2) * (-t67 * t66 - t71);
	unknown(2,5) = t604;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t685;
	unknown(3,3) = -r_i_i_C(1) * t67 * t688 + r_i_i_C(2) * t70 * t688 + pkin(5) * t695 + r_i_i_C(3) * t695;
	unknown(3,4) = r_i_i_C(1) * t70 * t695 + r_i_i_C(2) * t67 * t695;
	unknown(3,5) = t785;
	Ja_transl = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:31
	% EndTime: 2020-04-11 22:53:32
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