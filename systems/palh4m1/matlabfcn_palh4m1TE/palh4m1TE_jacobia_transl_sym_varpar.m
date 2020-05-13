% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh4m1TE
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
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh4m1TE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh4m1TE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1TE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:17
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (1098->66), mult. (1212->144), div. (56->5), fcn. (323->6), ass. (0->68)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = t3 * pkin(1);
	t6 = 0.2e1 * t4 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = t2 * t17 * pkin(1);
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t20 * t25 - t19;
	t28 = t1 * t27;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t29 * t31;
	t33 = t32 * r_i_i_C(1);
	t37 = pkin(1) * t2;
	t38 = t37 * t25;
	t39 = -t20 * t17 + t38;
	t40 = t1 * t39;
	t41 = t32 * r_i_i_C(2);
	t44 = cos(qJ(1));
	t48 = 0.1e1 / t17;
	t49 = t2 * t48;
	t54 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t59 = -t49 * pkin(1) * t54 / 0.2e1 - 0.2e1 * t20 * t22;
	t62 = t44 * t27;
	t63 = 0.1e1 / t23;
	t64 = t63 * t31;
	t65 = t64 * r_i_i_C(1);
	t67 = -t20 * t48;
	t71 = t67 * t54 / 0.2e1 + 0.2e1 * t37 * t22;
	t74 = t44 * t39;
	t75 = t64 * r_i_i_C(2);
	t83 = pkin(1) * pkin(2);
	t85 = t37 * pkin(2) * t15 + t11 * t2 * t83;
	t92 = t3 * t17 * pkin(1) - t49 * pkin(1) * t85 + 0.2e1 * t20 * t2 * t83 - t38;
	t96 = t30 ^ 2;
	t97 = 0.1e1 / t96;
	t98 = t29 * t97;
	t101 = r_i_i_C(1) * t2 * t83;
	t106 = t2 ^ 2;
	t110 = -0.2e1 * t7 * t106 * pkin(2) - t4 * t25 + t67 * t85 - t19;
	t116 = r_i_i_C(2) * t2 * t83;
	t145 = t31 * r_i_i_C(1);
	t150 = t31 * r_i_i_C(2);
	unknown(1,1) = -t28 * t33 / 0.2e1 + t40 * t41 / 0.2e1 + t44 * r_i_i_C(3) + t1 * pkin(9);
	unknown(1,2) = t44 * t59 * t33 / 0.2e1 - t44 * t71 * t41 / 0.2e1 - t62 * t65 / 0.2e1 + t74 * t75 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t44 * t92 * t33 / 0.2e1 + t62 * t98 * t101 - t44 * t110 * t41 / 0.2e1 - t74 * t98 * t116;
	unknown(2,1) = t62 * t33 / 0.2e1 - t74 * t41 / 0.2e1 + t1 * r_i_i_C(3) - t44 * pkin(9);
	unknown(2,2) = t1 * t59 * t33 / 0.2e1 - t1 * t71 * t41 / 0.2e1 - t28 * t65 / 0.2e1 + t40 * t75 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t1 * t92 * t33 / 0.2e1 + t28 * t98 * t101 - t1 * t110 * t41 / 0.2e1 - t40 * t98 * t116;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t71 * t29 * t145 / 0.2e1 - t39 * t63 * t145 / 0.2e1 - t27 * t63 * t150 / 0.2e1 + t59 * t29 * t150 / 0.2e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t110 * t29 * t145 / 0.2e1 + t39 * t29 * t97 * t101 + t92 * t29 * t150 / 0.2e1 + t27 * t29 * t97 * t116;
	Ja_transl = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (1554->76), mult. (1710->167), div. (67->5), fcn. (457->6), ass. (0->73)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t29 * t31;
	t33 = t32 * r_i_i_C(1);
	t36 = cos(qJ(1));
	t39 = t19 * t17 * pkin(1);
	t41 = t4 * t25 - t39;
	t42 = t1 * t41;
	t43 = t32 * r_i_i_C(3);
	t50 = 0.1e1 / t17;
	t51 = t4 * t50;
	t56 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t60 = t51 * t56 / 0.2e1 + 0.2e1 * t20 * t22;
	t63 = t36 * t27;
	t64 = 0.1e1 / t23;
	t65 = t64 * t31;
	t66 = t65 * r_i_i_C(1);
	t68 = t19 * t50;
	t73 = -t68 * pkin(1) * t56 / 0.2e1 + 0.2e1 * t4 * t22;
	t74 = t36 * t73;
	t76 = t36 * t41;
	t77 = t65 * r_i_i_C(3);
	t84 = pkin(1) * pkin(2);
	t86 = t20 * pkin(2) * t15 + t11 * t19 * t84;
	t90 = t19 ^ 2;
	t94 = -0.2e1 * t7 * t90 * pkin(2) - t3 * t25 + t51 * t86 - t39;
	t98 = t30 ^ 2;
	t99 = 0.1e1 / t98;
	t100 = t29 * t99;
	t103 = r_i_i_C(1) * t19 * t84;
	t113 = t2 * t17 * pkin(1) - t68 * pkin(1) * t86 - 0.2e1 * t4 * t19 * t84 - t26;
	t114 = t36 * t113;
	t119 = r_i_i_C(3) * t19 * t84;
	t124 = t20 * pkin(2);
	t139 = t1 * t73;
	t149 = t1 * t113;
	t160 = t31 * r_i_i_C(1);
	t165 = t31 * r_i_i_C(3);
	unknown(1,1) = -t28 * t33 / 0.2e1 - t36 * r_i_i_C(2) - t42 * t43 / 0.2e1 - t42 * t31 / 0.2e1 + t1 * pkin(9);
	unknown(1,2) = t36 * t60 * t33 / 0.2e1 + t74 * t31 / 0.2e1 + t74 * t43 / 0.2e1 - t63 * t66 / 0.2e1 - t76 * t77 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t36 * t94 * t33 / 0.2e1 + t63 * t100 * t103 + t114 * t43 / 0.2e1 + t76 * t100 * t119 + t114 * t31 / 0.2e1 + t76 * t99 * t124;
	unknown(2,1) = t63 * t33 / 0.2e1 - t1 * r_i_i_C(2) + t76 * t43 / 0.2e1 + t76 * t31 / 0.2e1 - t36 * pkin(9);
	unknown(2,2) = t1 * t60 * t33 / 0.2e1 + t139 * t31 / 0.2e1 + t139 * t43 / 0.2e1 - t28 * t66 / 0.2e1 - t42 * t77 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t1 * t94 * t33 / 0.2e1 + t28 * t100 * t103 + t149 * t43 / 0.2e1 + t42 * t100 * t119 + t149 * t31 / 0.2e1 + t42 * t99 * t124;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t73 * t29 * t160 / 0.2e1 + t41 * t64 * t160 / 0.2e1 - t27 * t64 * t165 / 0.2e1 + t60 * t29 * t165 / 0.2e1 + t60 * t31 / 0.2e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t113 * t29 * t160 / 0.2e1 - t41 * t29 * t99 * t103 + t94 * t29 * t165 / 0.2e1 + t27 * t29 * t99 * t119 + t94 * t31 / 0.2e1 + t27 * t99 * t124;
	Ja_transl = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:17
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (4382->132), mult. (4690->326), div. (251->7), fcn. (1214->6), ass. (0->112)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t53 = cos(qJ(1));
	t59 = 0.1e1 / t17;
	t60 = t4 * t59;
	t65 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t69 = t60 * t65 / 0.2e1 + 0.2e1 * t20 * t22;
	t71 = t53 * t69 * t29;
	t74 = t53 * t27;
	t76 = 0.1e1 / t23 / t22;
	t77 = t74 * t76;
	t80 = t74 * t29;
	t82 = t34 * t59 * t65;
	t85 = t19 * t59;
	t90 = -t85 * pkin(1) * t65 / 0.2e1 + 0.2e1 * t4 * t22;
	t91 = t53 * t90;
	t92 = t91 * t29;
	t95 = t53 * t40;
	t96 = t95 * t76;
	t99 = t95 * t29;
	t101 = 0.2e1 * t32 * t22 * t33;
	t126 = pkin(1) * pkin(2);
	t128 = t20 * pkin(2) * t15 + t11 * t19 * t126;
	t132 = t19 ^ 2;
	t136 = -0.2e1 * t7 * t132 * pkin(2) + t60 * t128 - t3 * t25 - t38;
	t138 = t53 * t136 * t29;
	t141 = t31 ^ 2;
	t142 = 0.1e1 / t141;
	t143 = t29 * t142;
	t144 = t74 * t143;
	t146 = t20 * pkin(2);
	t147 = t33 * t17 * t146;
	t151 = 0.2e1 * t34 * t59 * t128;
	t162 = -t85 * pkin(1) * t128 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t126 - t26;
	t163 = t53 * t162;
	t164 = t163 * t29;
	t167 = t95 * t143;
	t169 = t43 * t33 * t146;
	t172 = t29 * t32;
	t175 = t20 * pkin(2) * t33;
	t214 = t1 * t69 * t29;
	t217 = t28 * t76;
	t222 = t1 * t90;
	t223 = t222 * t29;
	t226 = t41 * t76;
	t251 = t1 * t136 * t29;
	t254 = t28 * t143;
	t259 = t1 * t162;
	t260 = t259 * t29;
	t263 = t41 * t143;
	t291 = t90 * t29;
	t294 = t40 * t76;
	t297 = t40 * t29;
	t298 = t297 * t32;
	t299 = t33 * t59;
	t300 = t299 * t65;
	t303 = t69 * t29;
	t306 = t27 * t76;
	t309 = t27 * t29;
	t324 = t309 * t32;
	t332 = t162 * t29;
	t335 = t142 * t33;
	t337 = t37 * t126;
	t340 = 0.2e1 * t299 * t128;
	t343 = t136 * t29;
	t346 = t142 * t43;
	unknown(1,1) = (t30 * t35 - t42 * t45) * r_i_i_C(1) / 0.4e1 + (t30 * t45 + t42 * t35) * r_i_i_C(2) / 0.4e1 + t53 * r_i_i_C(3) - t41 * t32 / 0.2e1 + t1 * pkin(9);
	unknown(1,2) = (-t71 * t35 / 0.4e1 + t77 * t35 / 0.2e1 - t80 * t82 / 0.8e1 + t92 * t45 / 0.4e1 - t96 * t45 / 0.2e1 + t99 * t101 / 0.4e1) * r_i_i_C(1) + (-t71 * t45 / 0.4e1 + t77 * t45 / 0.2e1 - t80 * t101 / 0.4e1 - t92 * t35 / 0.4e1 + t96 * t35 / 0.2e1 - t99 * t82 / 0.8e1) * r_i_i_C(2) + t91 * t32 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = (-t138 * t35 / 0.4e1 - t144 * t147 / 0.2e1 - t80 * t151 / 0.8e1 + t164 * t45 / 0.4e1 + t167 * t169 / 0.2e1 + t95 * t172 * t175 / 0.2e1) * r_i_i_C(1) + (-t138 * t45 / 0.4e1 - t144 * t169 / 0.2e1 - t74 * t172 * t175 / 0.2e1 - t164 * t35 / 0.4e1 - t167 * t147 / 0.2e1 - t99 * t151 / 0.8e1) * r_i_i_C(2) + t163 * t32 / 0.2e1 + t95 * t142 * t146;
	unknown(2,1) = (-t80 * t35 + t99 * t45) * r_i_i_C(1) / 0.4e1 + (-t99 * t35 - t80 * t45) * r_i_i_C(2) / 0.4e1 + t1 * r_i_i_C(3) + t95 * t32 / 0.2e1 - t53 * pkin(9);
	unknown(2,2) = (-t214 * t35 / 0.4e1 + t217 * t35 / 0.2e1 - t30 * t82 / 0.8e1 + t223 * t45 / 0.4e1 - t226 * t45 / 0.2e1 + t42 * t101 / 0.4e1) * r_i_i_C(1) + (-t214 * t45 / 0.4e1 + t217 * t45 / 0.2e1 - t30 * t101 / 0.4e1 - t223 * t35 / 0.4e1 + t226 * t35 / 0.2e1 - t42 * t82 / 0.8e1) * r_i_i_C(2) + t222 * t32 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = (-t251 * t35 / 0.4e1 - t254 * t147 / 0.2e1 - t30 * t151 / 0.8e1 + t260 * t45 / 0.4e1 + t263 * t169 / 0.2e1 + t41 * t172 * t175 / 0.2e1) * r_i_i_C(1) + (-t251 * t45 / 0.4e1 - t254 * t169 / 0.2e1 - t28 * t172 * t175 / 0.2e1 - t260 * t35 / 0.4e1 - t263 * t147 / 0.2e1 - t42 * t151 / 0.8e1) * r_i_i_C(2) + t259 * t32 / 0.2e1 + t41 * t142 * t146;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = (t291 * t35 / 0.4e1 - t294 * t35 / 0.2e1 + t298 * t300 / 0.8e1 + t303 * t45 / 0.4e1 - t306 * t45 / 0.2e1 + t309 * t101 / 0.4e1) * r_i_i_C(1) + (t291 * t45 / 0.4e1 - t294 * t45 / 0.2e1 + t297 * t101 / 0.4e1 - t303 * t35 / 0.4e1 + t306 * t35 / 0.2e1 - t324 * t300 / 0.8e1) * r_i_i_C(2) + t69 * t32 / 0.2e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = (t332 * t35 / 0.4e1 + t297 * t335 * t337 / 0.2e1 + t298 * t340 / 0.8e1 + t343 * t45 / 0.4e1 + t309 * t346 * t175 / 0.2e1 + t324 * t175 / 0.2e1) * r_i_i_C(1) + (t332 * t45 / 0.4e1 + t297 * t346 * t175 / 0.2e1 + t298 * t175 / 0.2e1 - t343 * t35 / 0.4e1 - t309 * t335 * t337 / 0.2e1 - t324 * t340 / 0.8e1) * r_i_i_C(2) + t136 * t32 / 0.2e1 + t27 * t142 * t146;
	Ja_transl = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:17
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (11210->165), mult. (11994->398), div. (683->7), fcn. (3140->8), ass. (0->134)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t47 = t30 * t35 - t42 * t45;
	t48 = cos(qJ(3));
	t52 = t30 * t45 + t42 * t35;
	t53 = sin(qJ(3));
	t61 = cos(qJ(1));
	t68 = 0.1e1 / t17;
	t69 = t4 * t68;
	t74 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t78 = t69 * t74 / 0.2e1 + 0.2e1 * t20 * t22;
	t80 = t61 * t78 * t29;
	t83 = t61 * t27;
	t85 = 0.1e1 / t23 / t22;
	t86 = t83 * t85;
	t89 = t83 * t29;
	t91 = t34 * t68 * t74;
	t94 = t19 * t68;
	t99 = -t94 * pkin(1) * t74 / 0.2e1 + 0.2e1 * t4 * t22;
	t100 = t61 * t99;
	t101 = t100 * t29;
	t104 = t61 * t40;
	t105 = t104 * t85;
	t108 = t104 * t29;
	t110 = 0.2e1 * t32 * t22 * t33;
	t113 = -t80 * t35 / 0.4e1 + t86 * t35 / 0.2e1 - t89 * t91 / 0.8e1 + t101 * t45 / 0.4e1 - t105 * t45 / 0.2e1 + t108 * t110 / 0.4e1;
	t127 = -t80 * t45 / 0.4e1 + t86 * t45 / 0.2e1 - t89 * t110 / 0.4e1 - t101 * t35 / 0.4e1 + t105 * t35 / 0.2e1 - t108 * t91 / 0.8e1;
	t141 = t108 * t45 - t89 * t35;
	t145 = -t108 * t35 - t89 * t45;
	t147 = -t141 * t53 / 0.4e1 + t145 * t48 / 0.4e1;
	t151 = -t141 * t48 / 0.4e1 - t145 * t53 / 0.4e1;
	t157 = pkin(1) * pkin(2);
	t159 = t20 * pkin(2) * t15 + t11 * t19 * t157;
	t163 = t19 ^ 2;
	t167 = -0.2e1 * t7 * t163 * pkin(2) + t69 * t159 - t3 * t25 - t38;
	t169 = t61 * t167 * t29;
	t172 = t31 ^ 2;
	t173 = 0.1e1 / t172;
	t174 = t29 * t173;
	t175 = t83 * t174;
	t177 = t20 * pkin(2);
	t178 = t33 * t17 * t177;
	t182 = 0.2e1 * t34 * t68 * t159;
	t193 = -t94 * pkin(1) * t159 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t157 - t26;
	t194 = t61 * t193;
	t195 = t194 * t29;
	t198 = t104 * t174;
	t200 = t43 * t33 * t177;
	t203 = t29 * t32;
	t206 = t20 * pkin(2) * t33;
	t209 = -t169 * t35 / 0.4e1 - t175 * t178 / 0.2e1 - t89 * t182 / 0.8e1 + t195 * t45 / 0.4e1 + t198 * t200 / 0.2e1 + t104 * t203 * t206 / 0.2e1;
	t224 = -t169 * t45 / 0.4e1 - t175 * t200 / 0.2e1 - t83 * t203 * t206 / 0.2e1 - t195 * t35 / 0.4e1 - t198 * t178 / 0.2e1 - t108 * t182 / 0.8e1;
	t247 = t1 * t78 * t29;
	t250 = t28 * t85;
	t255 = t1 * t99;
	t256 = t255 * t29;
	t259 = t41 * t85;
	t264 = -t247 * t35 / 0.4e1 + t250 * t35 / 0.2e1 - t30 * t91 / 0.8e1 + t256 * t45 / 0.4e1 - t259 * t45 / 0.2e1 + t42 * t110 / 0.4e1;
	t278 = -t247 * t45 / 0.4e1 + t250 * t45 / 0.2e1 - t30 * t110 / 0.4e1 - t256 * t35 / 0.4e1 + t259 * t35 / 0.2e1 - t42 * t91 / 0.8e1;
	t300 = t1 * t167 * t29;
	t303 = t28 * t174;
	t308 = t1 * t193;
	t309 = t308 * t29;
	t312 = t41 * t174;
	t318 = -t300 * t35 / 0.4e1 - t303 * t178 / 0.2e1 - t30 * t182 / 0.8e1 + t309 * t45 / 0.4e1 + t312 * t200 / 0.2e1 + t41 * t203 * t206 / 0.2e1;
	t333 = -t300 * t45 / 0.4e1 - t303 * t200 / 0.2e1 - t28 * t203 * t206 / 0.2e1 - t309 * t35 / 0.4e1 - t312 * t178 / 0.2e1 - t42 * t182 / 0.8e1;
	t347 = t99 * t29;
	t350 = t40 * t85;
	t353 = t40 * t29;
	t354 = t353 * t32;
	t355 = t33 * t68;
	t356 = t355 * t74;
	t359 = t78 * t29;
	t362 = t27 * t85;
	t365 = t27 * t29;
	t368 = t347 * t35 / 0.4e1 - t350 * t35 / 0.2e1 + t354 * t356 / 0.8e1 + t359 * t45 / 0.4e1 - t362 * t45 / 0.2e1 + t365 * t110 / 0.4e1;
	t380 = t365 * t32;
	t383 = t347 * t45 / 0.4e1 - t350 * t45 / 0.2e1 + t353 * t110 / 0.4e1 - t359 * t35 / 0.4e1 + t362 * t35 / 0.2e1 - t380 * t356 / 0.8e1;
	t397 = t353 * t35 + t365 * t45;
	t401 = -t365 * t35 + t353 * t45;
	t410 = t193 * t29;
	t413 = t173 * t33;
	t415 = t37 * t157;
	t418 = 0.2e1 * t355 * t159;
	t421 = t167 * t29;
	t424 = t173 * t43;
	t430 = t410 * t35 / 0.4e1 + t353 * t413 * t415 / 0.2e1 + t354 * t418 / 0.8e1 + t421 * t45 / 0.4e1 + t365 * t424 * t206 / 0.2e1 + t380 * t206 / 0.2e1;
	t446 = t410 * t45 / 0.4e1 + t353 * t424 * t206 / 0.2e1 + t354 * t206 / 0.2e1 - t421 * t35 / 0.4e1 - t365 * t413 * t415 / 0.2e1 - t380 * t418 / 0.8e1;
	unknown(1,1) = (t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1) * r_i_i_C(1) + (-t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1) * r_i_i_C(2) + t61 * r_i_i_C(3) + t47 * pkin(4) / 0.4e1 - t41 * t32 / 0.2e1 + t1 * pkin(9);
	unknown(1,2) = (t113 * t48 + t127 * t53) * r_i_i_C(1) + (-t113 * t53 + t127 * t48) * r_i_i_C(2) + t113 * pkin(4) + t100 * t32 / 0.2e1;
	unknown(1,3) = t147 * r_i_i_C(1) + t151 * r_i_i_C(2);
	unknown(1,4) = 0.0e0;
	unknown(1,5) = (t209 * t48 + t224 * t53) * r_i_i_C(1) + (-t209 * t53 + t224 * t48) * r_i_i_C(2) + t209 * pkin(4) + t194 * t32 / 0.2e1 + t104 * t173 * t177;
	unknown(2,1) = -t151 * r_i_i_C(1) + t147 * r_i_i_C(2) + t1 * r_i_i_C(3) + t141 * pkin(4) / 0.4e1 + t104 * t32 / 0.2e1 - t61 * pkin(9);
	unknown(2,2) = (t264 * t48 + t278 * t53) * r_i_i_C(1) + (-t264 * t53 + t278 * t48) * r_i_i_C(2) + t264 * pkin(4) + t255 * t32 / 0.2e1;
	unknown(2,3) = (t47 * t53 / 0.4e1 - t52 * t48 / 0.4e1) * r_i_i_C(1) + (t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1) * r_i_i_C(2);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = (t318 * t48 + t333 * t53) * r_i_i_C(1) + (-t318 * t53 + t333 * t48) * r_i_i_C(2) + t318 * pkin(4) + t308 * t32 / 0.2e1 + t41 * t173 * t177;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = (t368 * t48 + t383 * t53) * r_i_i_C(1) + (-t368 * t53 + t383 * t48) * r_i_i_C(2) + t368 * pkin(4) + t78 * t32 / 0.2e1;
	unknown(3,3) = (-t397 * t53 / 0.4e1 + t401 * t48 / 0.4e1) * r_i_i_C(1) + (-t397 * t48 / 0.4e1 - t401 * t53 / 0.4e1) * r_i_i_C(2);
	unknown(3,4) = 0.0e0;
	unknown(3,5) = (t430 * t48 + t446 * t53) * r_i_i_C(1) + (-t430 * t53 + t446 * t48) * r_i_i_C(2) + t430 * pkin(4) + t167 * t32 / 0.2e1 + t27 * t173 * t177;
	Ja_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:17
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (20945->194), mult. (22432->454), div. (1307->7), fcn. (5926->10), ass. (0->158)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t47 = t30 * t35 - t42 * t45;
	t48 = sin(qJ(3));
	t52 = t30 * t45 + t42 * t35;
	t53 = cos(qJ(3));
	t55 = -t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1;
	t56 = cos(qJ(4));
	t58 = cos(qJ(1));
	t59 = sin(qJ(4));
	t60 = t58 * t59;
	t64 = t58 * t56;
	t69 = t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	t77 = 0.1e1 / t17;
	t78 = t4 * t77;
	t83 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t87 = t78 * t83 / 0.2e1 + 0.2e1 * t20 * t22;
	t89 = t58 * t87 * t29;
	t92 = t58 * t27;
	t94 = 0.1e1 / t23 / t22;
	t95 = t92 * t94;
	t98 = t92 * t29;
	t100 = t34 * t77 * t83;
	t103 = t19 * t77;
	t108 = -t103 * pkin(1) * t83 / 0.2e1 + 0.2e1 * t4 * t22;
	t109 = t58 * t108;
	t110 = t109 * t29;
	t113 = t58 * t40;
	t114 = t113 * t94;
	t117 = t113 * t29;
	t119 = 0.2e1 * t32 * t22 * t33;
	t122 = -t89 * t35 / 0.4e1 + t95 * t35 / 0.2e1 - t98 * t100 / 0.8e1 + t110 * t45 / 0.4e1 - t114 * t45 / 0.2e1 + t117 * t119 / 0.4e1;
	t136 = -t89 * t45 / 0.4e1 + t95 * t45 / 0.2e1 - t98 * t119 / 0.4e1 - t110 * t35 / 0.4e1 + t114 * t35 / 0.2e1 - t117 * t100 / 0.8e1;
	t138 = -t122 * t48 + t136 * t53;
	t145 = t122 * t53 + t136 * t48;
	t154 = t117 * t45 - t98 * t35;
	t158 = -t117 * t35 - t98 * t45;
	t160 = -t154 * t53 / 0.4e1 - t158 * t48 / 0.4e1;
	t167 = -t154 * t48 / 0.4e1 + t158 * t53 / 0.4e1;
	t173 = -t1 * t56 + t167 * t59;
	t177 = t1 * t59 + t167 * t56;
	t183 = pkin(1) * pkin(2);
	t185 = t20 * pkin(2) * t15 + t11 * t19 * t183;
	t189 = t19 ^ 2;
	t193 = -0.2e1 * t7 * t189 * pkin(2) + t78 * t185 - t3 * t25 - t38;
	t195 = t58 * t193 * t29;
	t198 = t31 ^ 2;
	t199 = 0.1e1 / t198;
	t200 = t29 * t199;
	t201 = t92 * t200;
	t203 = t20 * pkin(2);
	t204 = t33 * t17 * t203;
	t208 = 0.2e1 * t34 * t77 * t185;
	t219 = -t103 * pkin(1) * t185 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t183 - t26;
	t220 = t58 * t219;
	t221 = t220 * t29;
	t224 = t113 * t200;
	t226 = t43 * t33 * t203;
	t229 = t29 * t32;
	t232 = t20 * pkin(2) * t33;
	t235 = -t195 * t35 / 0.4e1 - t201 * t204 / 0.2e1 - t98 * t208 / 0.8e1 + t221 * t45 / 0.4e1 + t224 * t226 / 0.2e1 + t113 * t229 * t232 / 0.2e1;
	t250 = -t195 * t45 / 0.4e1 - t201 * t226 / 0.2e1 - t92 * t229 * t232 / 0.2e1 - t221 * t35 / 0.4e1 - t224 * t204 / 0.2e1 - t117 * t208 / 0.8e1;
	t252 = -t235 * t48 + t250 * t53;
	t259 = t235 * t53 + t250 * t48;
	t278 = t1 * t87 * t29;
	t281 = t28 * t94;
	t286 = t1 * t108;
	t287 = t286 * t29;
	t290 = t41 * t94;
	t295 = -t278 * t35 / 0.4e1 + t281 * t35 / 0.2e1 - t30 * t100 / 0.8e1 + t287 * t45 / 0.4e1 - t290 * t45 / 0.2e1 + t42 * t119 / 0.4e1;
	t309 = -t278 * t45 / 0.4e1 + t281 * t45 / 0.2e1 - t30 * t119 / 0.4e1 - t287 * t35 / 0.4e1 + t290 * t35 / 0.2e1 - t42 * t100 / 0.8e1;
	t311 = -t295 * t48 + t309 * t53;
	t318 = t295 * t53 + t309 * t48;
	t327 = t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	t334 = t47 * t48 / 0.4e1 - t52 * t53 / 0.4e1;
	t346 = t1 * t193 * t29;
	t349 = t28 * t200;
	t354 = t1 * t219;
	t355 = t354 * t29;
	t358 = t41 * t200;
	t364 = -t346 * t35 / 0.4e1 - t349 * t204 / 0.2e1 - t30 * t208 / 0.8e1 + t355 * t45 / 0.4e1 + t358 * t226 / 0.2e1 + t41 * t229 * t232 / 0.2e1;
	t379 = -t346 * t45 / 0.4e1 - t349 * t226 / 0.2e1 - t28 * t229 * t232 / 0.2e1 - t355 * t35 / 0.4e1 - t358 * t204 / 0.2e1 - t42 * t208 / 0.8e1;
	t381 = -t364 * t48 + t379 * t53;
	t388 = t364 * t53 + t379 * t48;
	t397 = t108 * t29;
	t400 = t40 * t94;
	t403 = t40 * t29;
	t404 = t403 * t32;
	t405 = t33 * t77;
	t406 = t405 * t83;
	t409 = t87 * t29;
	t412 = t27 * t94;
	t415 = t27 * t29;
	t418 = t397 * t35 / 0.4e1 - t400 * t35 / 0.2e1 + t404 * t406 / 0.8e1 + t409 * t45 / 0.4e1 - t412 * t45 / 0.2e1 + t415 * t119 / 0.4e1;
	t430 = t415 * t32;
	t433 = t397 * t45 / 0.4e1 - t400 * t45 / 0.2e1 + t403 * t119 / 0.4e1 - t409 * t35 / 0.4e1 + t412 * t35 / 0.2e1 - t430 * t406 / 0.8e1;
	t435 = -t418 * t48 + t433 * t53;
	t442 = t418 * t53 + t433 * t48;
	t451 = t403 * t35 + t415 * t45;
	t455 = -t415 * t35 + t403 * t45;
	t457 = -t451 * t53 / 0.4e1 - t455 * t48 / 0.4e1;
	t464 = -t451 * t48 / 0.4e1 + t455 * t53 / 0.4e1;
	t473 = t219 * t29;
	t476 = t199 * t33;
	t478 = t37 * t183;
	t481 = 0.2e1 * t405 * t185;
	t484 = t193 * t29;
	t487 = t199 * t43;
	t493 = t473 * t35 / 0.4e1 + t403 * t476 * t478 / 0.2e1 + t404 * t481 / 0.8e1 + t484 * t45 / 0.4e1 + t415 * t487 * t232 / 0.2e1 + t430 * t232 / 0.2e1;
	t509 = t473 * t45 / 0.4e1 + t403 * t487 * t232 / 0.2e1 + t404 * t232 / 0.2e1 - t484 * t35 / 0.4e1 - t415 * t476 * t478 / 0.2e1 - t430 * t481 / 0.8e1;
	t511 = -t493 * t48 + t509 * t53;
	t518 = t509 * t48 + t493 * t53;
	unknown(1,1) = (-t55 * t56 - t60) * r_i_i_C(1) + (t55 * t59 - t64) * r_i_i_C(2) + t69 * r_i_i_C(3) + t69 * pkin(5) + t47 * pkin(4) / 0.4e1 - t41 * t32 / 0.2e1 + t1 * pkin(9);
	unknown(1,2) = -t138 * t56 * r_i_i_C(1) + t138 * t59 * r_i_i_C(2) + t145 * r_i_i_C(3) + t145 * pkin(5) + t122 * pkin(4) + t109 * t32 / 0.2e1;
	unknown(1,3) = -t160 * t56 * r_i_i_C(1) + t160 * t59 * r_i_i_C(2) + t167 * pkin(5) + t167 * r_i_i_C(3);
	unknown(1,4) = t173 * r_i_i_C(1) + t177 * r_i_i_C(2);
	unknown(1,5) = -t252 * t56 * r_i_i_C(1) + t252 * t59 * r_i_i_C(2) + t259 * r_i_i_C(3) + t259 * pkin(5) + t235 * pkin(4) + t220 * t32 / 0.2e1 + t113 * t199 * t203;
	unknown(2,1) = -t177 * r_i_i_C(1) + t173 * r_i_i_C(2) - t160 * r_i_i_C(3) - t160 * pkin(5) + t154 * pkin(4) / 0.4e1 + t113 * t32 / 0.2e1 - t58 * pkin(9);
	unknown(2,2) = -t311 * t56 * r_i_i_C(1) + t311 * t59 * r_i_i_C(2) + t318 * r_i_i_C(3) + t318 * pkin(5) + t295 * pkin(4) + t286 * t32 / 0.2e1;
	unknown(2,3) = -t327 * t56 * r_i_i_C(1) + t327 * t59 * r_i_i_C(2) + t334 * pkin(5) + t334 * r_i_i_C(3);
	unknown(2,4) = (t334 * t59 + t64) * r_i_i_C(1) + (t334 * t56 - t60) * r_i_i_C(2);
	unknown(2,5) = -t381 * t56 * r_i_i_C(1) + t381 * t59 * r_i_i_C(2) + t388 * r_i_i_C(3) + t388 * pkin(5) + t364 * pkin(4) + t354 * t32 / 0.2e1 + t41 * t199 * t203;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t435 * t56 * r_i_i_C(1) + t435 * t59 * r_i_i_C(2) + t442 * r_i_i_C(3) + t442 * pkin(5) + t418 * pkin(4) + t87 * t32 / 0.2e1;
	unknown(3,3) = -t457 * t56 * r_i_i_C(1) + t457 * t59 * r_i_i_C(2) + t464 * pkin(5) + t464 * r_i_i_C(3);
	unknown(3,4) = t464 * t59 * r_i_i_C(1) + t464 * t56 * r_i_i_C(2);
	unknown(3,5) = -t511 * t56 * r_i_i_C(1) + t511 * t59 * r_i_i_C(2) + t518 * r_i_i_C(3) + t518 * pkin(5) + t493 * pkin(4) + t193 * t32 / 0.2e1 + t27 * t199 * t203;
	Ja_transl = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:17
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.05s
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