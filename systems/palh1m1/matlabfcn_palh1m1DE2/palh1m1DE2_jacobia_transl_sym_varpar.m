% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m1DE2
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
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m1DE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m1DE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_jacobia_transl_sym_varpar: pkin has to be [23x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = -pkin(16) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) + t5 * t2, t7 * t4, 0, 0; t2 * r_i_i_C(3) - t5 * t4, t7 * t2, 0, 0; 0, t6, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->11)
	t5 = qJ(2) + qJ(3);
	t3 = sin(t5);
	t4 = cos(t5);
	t13 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t16 = t13 - sin(qJ(2)) * pkin(1);
	t12 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t11 = pkin(16) + t16;
	t10 = -cos(qJ(2)) * pkin(1) + t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t9 * r_i_i_C(3) - t11 * t7, t10 * t9, t12 * t9, 0; t7 * r_i_i_C(3) + t11 * t9, t10 * t7, t12 * t7, 0; 0, t16, t13, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:50:09
	% DurationCPUTime: 11.08s
	% Computational Cost: add. (322885->98), mult. (487477->179), div. (21464->11), fcn. (309017->21), ass. (0->116)
	t136 = -2 * pkin(1);
	t135 = -2 * pkin(5);
	t134 = pkin(4) * pkin(5);
	t68 = sin(pkin(21));
	t133 = t68 / 0.2e1;
	t73 = cos(qJ(3));
	t132 = -t73 / 0.2e1;
	t131 = -pkin(8) - pkin(3);
	t130 = -pkin(8) + pkin(3);
	t129 = (-pkin(9) - pkin(11));
	t128 = (-pkin(9) + pkin(11));
	t121 = sin(pkin(19));
	t71 = sin(qJ(2));
	t74 = cos(qJ(2));
	t76 = cos(pkin(19));
	t58 = t121 * t71 + t74 * t76;
	t123 = pkin(7) * t58;
	t57 = -t121 * t74 + t71 * t76;
	t124 = pkin(7) * t57;
	t109 = (pkin(1) ^ 2) + t124 * t136;
	t78 = pkin(7) ^ 2;
	t54 = t78 + t109;
	t51 = pkin(3) ^ 2 - pkin(8) ^ 2 + t54;
	t55 = pkin(1) - t124;
	t49 = (pkin(7) - t131) * (pkin(7) + t131) + t109;
	t50 = (pkin(7) - t130) * (pkin(7) + t130) + t109;
	t84 = sqrt(-t50 * t49);
	t45 = t123 * t51 + t55 * t84;
	t111 = t73 * t45;
	t115 = t58 * t84;
	t44 = -pkin(7) * t115 + t51 * t55;
	t70 = sin(qJ(3));
	t114 = t70 * t44;
	t52 = 0.1e1 / t54;
	t81 = 0.1e1 / pkin(3);
	t116 = t52 * t81;
	t40 = (t114 / 0.2e1 + t111 / 0.2e1) * t116;
	t112 = t73 * t44;
	t113 = t70 * t45;
	t41 = (-t112 / 0.2e1 + t113 / 0.2e1) * t116;
	t66 = pkin(23) + pkin(22);
	t63 = sin(t66);
	t64 = cos(t66);
	t33 = t40 * t64 - t41 * t63;
	t127 = pkin(4) * t33;
	t95 = t40 * t63 + t64 * t41;
	t126 = pkin(4) * t95;
	t67 = qJ(2) + qJ(3);
	t125 = pkin(5) * sin(t67);
	t110 = (pkin(5) ^ 2) - t127 * t135;
	t80 = pkin(4) ^ 2;
	t30 = t80 + t110;
	t28 = 0.1e1 / t30;
	t77 = 0.1e1 / pkin(11);
	t118 = t28 * t77;
	t27 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t30;
	t31 = pkin(5) + t127;
	t25 = ((pkin(4) - t129) * (pkin(4) + t129)) + t110;
	t26 = ((pkin(4) - t128) * (pkin(4) + t128)) + t110;
	t83 = sqrt(-t26 * t25);
	t16 = -t126 * t83 + t27 * t31;
	t17 = t126 * t27 + t31 * t83;
	t69 = cos(pkin(21));
	t14 = (-t16 * t69 / 0.2e1 + t17 * t133) * t118;
	t12 = 0.1e1 / t14 ^ 2;
	t13 = (t16 * t133 + t17 * t69 / 0.2e1) * t118;
	t122 = t77 / (t12 * t13 ^ 2 + 0.1e1);
	t120 = t12 * t13;
	t119 = t28 * t69;
	t108 = pkin(1) * t123;
	t117 = 0.2e1 / t84 * (t49 + t50) * t108;
	t107 = 0.1e1 / t30 ^ 2 * t134;
	t22 = 0.1e1 / t83;
	t106 = t22 * t31 / 0.2e1;
	t105 = -t22 * t95 / 0.2e1;
	t104 = t28 * t133;
	t103 = -t119 / 0.2e1;
	t102 = t119 / 0.2e1;
	t101 = t80 * t95 * t135;
	t61 = pkin(5) * cos(t67);
	t100 = -t71 * pkin(1) + t61;
	t99 = t135 * t31 - t27;
	t98 = 0.1e1 / t54 ^ 2 * t108;
	t97 = t68 * t107;
	t96 = t69 * t107;
	t37 = (t57 * t84 + (-t117 / 0.2e1 - t51 + t55 * t136) * t58) * pkin(7);
	t38 = t55 * t117 / 0.2e1 + t78 * t58 ^ 2 * t136 + (-t51 * t57 - t115) * pkin(7);
	t23 = ((-t70 * t37 / 0.2e1 + t38 * t132) * t52 + (-t111 - t114) * t98) * t81;
	t24 = ((t37 * t132 + t70 * t38 / 0.2e1) * t52 + (-t112 + t113) * t98) * t81;
	t19 = t23 * t64 + t24 * t63;
	t94 = t19 * t97;
	t93 = t95 * t97;
	t92 = t19 * t96;
	t91 = t95 * t96;
	t7 = atan2(t13, t14) + t67;
	t5 = sin(t7);
	t6 = cos(t7);
	t90 = r_i_i_C(1) * t6 - r_i_i_C(2) * t5;
	t89 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t88 = pkin(16) + t100 + t90;
	t87 = 0.2e1 * (t25 + t26) * t134;
	t11 = 0.1e1 / t14;
	t15 = t19 * t87;
	t20 = -t23 * t63 + t24 * t64;
	t3 = (t105 * t15 + t19 * t99 - t20 * t83) * pkin(4);
	t4 = t15 * t106 + t19 * t101 + (-t19 * t83 + t20 * t27) * pkin(4);
	t1 = 0.1e1 + ((t102 * t4 + t104 * t3 + t16 * t94 + t17 * t92) * t11 - (t103 * t3 + t104 * t4 - t16 * t92 + t17 * t94) * t120) * t122;
	t86 = -pkin(1) * t74 + t1 * t89 - t125;
	t18 = t95 * t87;
	t10 = t18 * t106 + t95 * t101 + (t27 * t33 - t83 * t95) * pkin(4);
	t9 = (t105 * t18 - t33 * t83 + t95 * t99) * pkin(4);
	t2 = 0.1e1 + ((t10 * t102 + t104 * t9 + t16 * t93 + t17 * t91) * t11 - (t10 * t104 + t103 * t9 - t16 * t91 + t17 * t93) * t120) * t122;
	t85 = t2 * t89 - t125;
	t75 = cos(qJ(1));
	t72 = sin(qJ(1));
	t8 = [t75 * r_i_i_C(3) - t72 * t88, t86 * t75, t85 * t75, 0; t72 * r_i_i_C(3) + t75 * t88, t86 * t72, t85 * t72, 0; 0, t1 * t90 + t100, t2 * t90 + t61, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:50:32
	% DurationCPUTime: 26.26s
	% Computational Cost: add. (823996->112), mult. (1243899->199), div. (54860->11), fcn. (788475->23), ass. (0->129)
	t134 = sin(pkin(19));
	t76 = sin(qJ(2));
	t80 = cos(qJ(2));
	t82 = cos(pkin(19));
	t62 = t76 * t134 + t80 * t82;
	t135 = pkin(7) * t62;
	t61 = -t80 * t134 + t76 * t82;
	t136 = pkin(7) * t61;
	t149 = -2 * pkin(1);
	t117 = (pkin(1) ^ 2) + t136 * t149;
	t84 = pkin(7) ^ 2;
	t58 = t84 + t117;
	t55 = pkin(3) ^ 2 - pkin(8) ^ 2 + t58;
	t59 = pkin(1) - t136;
	t144 = -pkin(8) - pkin(3);
	t53 = (pkin(7) - t144) * (pkin(7) + t144) + t117;
	t143 = -pkin(8) + pkin(3);
	t54 = (pkin(7) - t143) * (pkin(7) + t143) + t117;
	t90 = sqrt(-t54 * t53);
	t49 = t55 * t135 + t59 * t90;
	t79 = cos(qJ(3));
	t119 = t79 * t49;
	t127 = t62 * t90;
	t48 = -pkin(7) * t127 + t55 * t59;
	t75 = sin(qJ(3));
	t124 = t75 * t48;
	t56 = 0.1e1 / t58;
	t87 = 0.1e1 / pkin(3);
	t128 = t56 * t87;
	t44 = (t124 / 0.2e1 + t119 / 0.2e1) * t128;
	t120 = t79 * t48;
	t123 = t75 * t49;
	t45 = (-t120 / 0.2e1 + t123 / 0.2e1) * t128;
	t70 = pkin(23) + pkin(22);
	t67 = sin(t70);
	t68 = cos(t70);
	t37 = t44 * t68 - t45 * t67;
	t139 = pkin(4) * t37;
	t148 = -2 * pkin(5);
	t118 = (pkin(5) ^ 2) - t139 * t148;
	t86 = pkin(4) ^ 2;
	t34 = t86 + t118;
	t32 = 0.1e1 / t34;
	t83 = 0.1e1 / pkin(11);
	t130 = t32 * t83;
	t72 = sin(pkin(21));
	t146 = t72 / 0.2e1;
	t101 = t44 * t67 + t68 * t45;
	t138 = pkin(4) * t101;
	t31 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t34;
	t35 = pkin(5) + t139;
	t142 = -pkin(9) - pkin(11);
	t29 = (pkin(4) - t142) * (pkin(4) + t142) + t118;
	t141 = -pkin(9) + pkin(11);
	t30 = (pkin(4) - t141) * (pkin(4) + t141) + t118;
	t89 = sqrt(-t30 * t29);
	t20 = -t89 * t138 + t31 * t35;
	t21 = t31 * t138 + t35 * t89;
	t73 = cos(pkin(21));
	t17 = (t20 * t146 + t21 * t73 / 0.2e1) * t130;
	t18 = (-t20 * t73 / 0.2e1 + t21 * t146) * t130;
	t71 = qJ(2) + qJ(3);
	t11 = atan2(t17, t18) + t71;
	t10 = cos(t11);
	t65 = pkin(5) * cos(t71);
	t107 = -pkin(1) * t76 + t65;
	t140 = pkin(12) + r_i_i_C(3);
	t9 = sin(t11);
	t114 = t140 * t9;
	t150 = pkin(10) * t10 + pkin(16) + t107 + t114;
	t147 = pkin(4) * pkin(5);
	t145 = -t79 / 0.2e1;
	t137 = pkin(5) * sin(t71);
	t16 = 0.1e1 / t18 ^ 2;
	t133 = 0.1e1 / (t16 * t17 ^ 2 + 0.1e1) * t83;
	t132 = t16 * t17;
	t131 = t32 * t73;
	t116 = pkin(1) * t135;
	t129 = 0.2e1 / t90 * (t53 + t54) * t116;
	t74 = sin(qJ(4));
	t77 = sin(qJ(1));
	t126 = t74 * t77;
	t81 = cos(qJ(1));
	t125 = t74 * t81;
	t78 = cos(qJ(4));
	t122 = t77 * t78;
	t121 = t78 * t81;
	t115 = 0.1e1 / t34 ^ 2 * t147;
	t26 = 0.1e1 / t89;
	t113 = t26 * t35 / 0.2e1;
	t112 = -t26 * t101 / 0.2e1;
	t111 = t32 * t146;
	t110 = -t131 / 0.2e1;
	t109 = t131 / 0.2e1;
	t108 = t86 * t101 * t148;
	t106 = t35 * t148 - t31;
	t104 = 0.1e1 / t58 ^ 2 * t116;
	t103 = t72 * t115;
	t102 = t73 * t115;
	t41 = (t61 * t90 + (-t129 / 0.2e1 - t55 + t59 * t149) * t62) * pkin(7);
	t42 = t59 * t129 / 0.2e1 + t84 * t62 ^ 2 * t149 + (-t61 * t55 - t127) * pkin(7);
	t27 = ((-t75 * t41 / 0.2e1 + t42 * t145) * t56 + (-t119 - t124) * t104) * t87;
	t28 = ((t41 * t145 + t75 * t42 / 0.2e1) * t56 + (-t120 + t123) * t104) * t87;
	t23 = t27 * t68 + t28 * t67;
	t100 = t23 * t103;
	t99 = t101 * t103;
	t98 = t23 * t102;
	t97 = t101 * t102;
	t96 = r_i_i_C(1) * t78 - r_i_i_C(2) * t74 + pkin(10);
	t95 = 0.2e1 * (t29 + t30) * t147;
	t94 = t140 * t10 - t96 * t9;
	t93 = t96 * t10 + t114;
	t15 = 0.1e1 / t18;
	t19 = t23 * t95;
	t24 = -t27 * t67 + t28 * t68;
	t3 = (t106 * t23 + t19 * t112 - t24 * t89) * pkin(4);
	t4 = t19 * t113 + t23 * t108 + (-t23 * t89 + t24 * t31) * pkin(4);
	t1 = 0.1e1 + ((t20 * t100 + t4 * t109 + t3 * t111 + t21 * t98) * t15 - (t21 * t100 + t3 * t110 + t4 * t111 - t20 * t98) * t132) * t133;
	t92 = -pkin(1) * t80 + t94 * t1 - t137;
	t22 = t101 * t95;
	t13 = (t101 * t106 + t22 * t112 - t37 * t89) * pkin(4);
	t14 = t22 * t113 + t101 * t108 + (-t101 * t89 + t37 * t31) * pkin(4);
	t2 = 0.1e1 + ((t14 * t109 + t13 * t111 + t20 * t99 + t21 * t97) * t15 - (t13 * t110 + t14 * t111 - t20 * t97 + t21 * t99) * t132) * t133;
	t91 = t94 * t2 - t137;
	t8 = t10 * t121 + t126;
	t7 = -t10 * t125 + t122;
	t6 = -t10 * t122 + t125;
	t5 = t10 * t126 + t121;
	t12 = [r_i_i_C(1) * t6 + r_i_i_C(2) * t5 - t150 * t77, t92 * t81, t91 * t81, r_i_i_C(1) * t7 - r_i_i_C(2) * t8; r_i_i_C(1) * t8 + r_i_i_C(2) * t7 + t150 * t81, t92 * t77, t91 * t77, -r_i_i_C(1) * t5 + r_i_i_C(2) * t6; 0, t93 * t1 + t107, t93 * t2 + t65, (-r_i_i_C(1) * t74 - r_i_i_C(2) * t78) * t9;];
	Ja_transl = t12;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (5669->40), mult. (8752->77), div. (380->7), fcn. (5636->13), ass. (0->50)
	t57 = -2 * pkin(7);
	t32 = cos(pkin(18));
	t56 = t32 / 0.2e1;
	t55 = (-pkin(8) - pkin(3));
	t54 = (-pkin(8) + pkin(3));
	t27 = sin(qJ(2));
	t31 = cos(pkin(19));
	t50 = sin(pkin(19));
	t51 = cos(qJ(2));
	t24 = t27 * t31 - t51 * t50;
	t53 = pkin(1) * t24;
	t25 = t27 * t50 + t51 * t31;
	t52 = pkin(1) * t25;
	t35 = pkin(1) ^ 2;
	t45 = t53 * t57 + t35;
	t16 = ((pkin(7) - t55) * (pkin(7) + t55)) + t45;
	t17 = ((pkin(7) - t54) * (pkin(7) + t54)) + t45;
	t36 = sqrt(-t17 * t16);
	t44 = pkin(7) * t52;
	t49 = 0.1e1 / t36 * (t16 + t17) * t44;
	t21 = (pkin(7) ^ 2) + t45;
	t19 = 0.1e1 / t21;
	t29 = sin(pkin(18));
	t48 = t19 * t29;
	t33 = 0.1e1 / pkin(8);
	t47 = t19 * t33;
	t46 = t25 * t36;
	t43 = t19 * t56;
	t42 = 0.1e1 / t21 ^ 2 * t44;
	t18 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t21;
	t22 = -pkin(7) + t53;
	t11 = -pkin(1) * t46 - t18 * t22;
	t41 = t11 * t42;
	t12 = t18 * t52 - t22 * t36;
	t40 = t12 * t42;
	t10 = (t12 * t56 + t11 * t29 / 0.2e1) * t47;
	t9 = (t11 * t56 - t29 * t12 / 0.2e1) * t47;
	t5 = atan2(t10, t9);
	t2 = sin(t5);
	t3 = cos(t5);
	t39 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t38 = -pkin(15) + t39;
	t6 = (t24 * t36 + (0.2e1 * t22 * pkin(7) - t18 - t49) * t25) * pkin(1);
	t7 = -t22 * t49 + t35 * t25 ^ 2 * t57 + (-t18 * t24 - t46) * pkin(1);
	t8 = 0.1e1 / t9 ^ 2;
	t1 = ((t7 * t43 + t32 * t40 + t6 * t48 / 0.2e1 + t29 * t41) / t9 - (t6 * t43 + t32 * t41 - t7 * t48 / 0.2e1 - t29 * t40) * t10 * t8) / (t10 ^ 2 * t8 + 0.1e1) * t33;
	t37 = t1 * (-r_i_i_C(1) * t2 - r_i_i_C(2) * t3);
	t30 = cos(qJ(1));
	t28 = sin(qJ(1));
	t4 = [t30 * r_i_i_C(3) - t38 * t28, t30 * t37, 0, 0; t28 * r_i_i_C(3) + t38 * t30, t28 * t37, 0, 0; 0, t39 * t1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (5690->44), mult. (8759->78), div. (380->7), fcn. (5643->13), ass. (0->51)
	t59 = -2 * pkin(1);
	t28 = sin(pkin(23));
	t58 = t28 / 0.2e1;
	t57 = -pkin(8) - pkin(3);
	t56 = -pkin(8) + pkin(3);
	t30 = sin(qJ(2));
	t32 = cos(qJ(2));
	t34 = cos(pkin(19));
	t52 = sin(pkin(19));
	t24 = t30 * t34 - t32 * t52;
	t55 = pkin(7) * t24;
	t25 = t30 * t52 + t32 * t34;
	t54 = pkin(7) * t25;
	t53 = t30 * pkin(1);
	t47 = (pkin(1) ^ 2) + t55 * t59;
	t16 = (pkin(7) - t57) * (pkin(7) + t57) + t47;
	t17 = (pkin(7) - t56) * (pkin(7) + t56) + t47;
	t38 = sqrt(-t17 * t16);
	t46 = pkin(1) * t54;
	t51 = 0.1e1 / t38 * (t16 + t17) * t46;
	t35 = pkin(7) ^ 2;
	t21 = t35 + t47;
	t19 = 0.1e1 / t21;
	t36 = 0.1e1 / pkin(3);
	t50 = t19 * t36;
	t49 = t25 * t38;
	t29 = cos(pkin(23));
	t48 = t29 * t19;
	t45 = t19 * t58;
	t44 = 0.1e1 / t21 ^ 2 * t46;
	t43 = t28 * t44;
	t42 = t29 * t44;
	t18 = pkin(3) ^ 2 - pkin(8) ^ 2 + t21;
	t22 = pkin(1) - t55;
	t11 = -pkin(7) * t49 + t22 * t18;
	t12 = t18 * t54 + t22 * t38;
	t10 = (t29 * t12 / 0.2e1 + t11 * t58) * t50;
	t9 = (-t29 * t11 / 0.2e1 + t12 * t58) * t50;
	t4 = qJ(2) + atan2(t10, t9);
	t2 = sin(t4);
	t3 = cos(t4);
	t41 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t40 = -pkin(16) + t53 - t41;
	t6 = (t24 * t38 + (t22 * t59 - t18 - t51) * t25) * pkin(7);
	t7 = t22 * t51 + t35 * t25 ^ 2 * t59 + (-t24 * t18 - t49) * pkin(7);
	t8 = 0.1e1 / t9 ^ 2;
	t1 = 0.1e1 + ((t7 * t48 / 0.2e1 + t12 * t42 + t6 * t45 + t11 * t43) / t9 - (-t6 * t48 / 0.2e1 - t11 * t42 + t7 * t45 + t12 * t43) * t10 * t8) / (t10 ^ 2 * t8 + 0.1e1) * t36;
	t39 = -pkin(1) * t32 + (-r_i_i_C(1) * t3 + r_i_i_C(2) * t2) * t1;
	t33 = cos(qJ(1));
	t31 = sin(qJ(1));
	t5 = [t33 * r_i_i_C(3) + t40 * t31, t39 * t33, 0, 0; t31 * r_i_i_C(3) - t40 * t33, t39 * t31, 0, 0; 0, t41 * t1 - t53, 0, 0;];
	Ja_transl = t5;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (3064->35), mult. (4608->66), div. (148->7), fcn. (2956->11), ass. (0->42)
	t28 = pkin(6) ^ 2;
	t22 = sin(pkin(20));
	t23 = cos(pkin(20));
	t24 = sin(qJ(3));
	t26 = cos(qJ(3));
	t20 = t26 * t22 + t24 * t23;
	t43 = pkin(6) * t20;
	t39 = 0.2e1 * pkin(1) * t43 + t28;
	t17 = pkin(1) ^ 2 + t39;
	t14 = pkin(2) ^ 2 - pkin(13) ^ 2 + t17;
	t18 = -pkin(1) - t43;
	t46 = -pkin(2) - pkin(13);
	t12 = (pkin(1) - t46) * (pkin(1) + t46) + t39;
	t45 = -pkin(2) + pkin(13);
	t13 = (pkin(1) - t45) * (pkin(1) + t45) + t39;
	t31 = sqrt(-t13 * t12);
	t21 = t24 * t22 - t26 * t23;
	t41 = t21 * pkin(6);
	t8 = t14 * t41 - t18 * t31;
	t48 = pkin(6) * t8;
	t15 = 0.1e1 / t17;
	t47 = t15 / 0.2e1;
	t44 = pkin(1) * t21;
	t42 = 0.1e1 / t31 * (t12 + t13) * pkin(1) * t41;
	t40 = t21 * t31;
	t29 = 0.1e1 / pkin(2);
	t38 = t29 * t47;
	t7 = -pkin(6) * t40 - t18 * t14;
	t4 = qJ(2) + atan2(t8 * t38, t7 * t38);
	t2 = sin(t4);
	t3 = cos(t4);
	t36 = -r_i_i_C(1) * t3 + r_i_i_C(2) * t2;
	t35 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t34 = -pkin(16) - t35;
	t25 = sin(qJ(1));
	t33 = t36 * t25;
	t27 = cos(qJ(1));
	t32 = t36 * t27;
	t16 = 0.1e1 / t17 ^ 2;
	t6 = 0.1e1 / t7 ^ 2;
	t1 = 0.2e1 * (((-t18 * t42 + (t20 * t14 - t40) * pkin(6)) * t47 + (-t15 * t28 * t21 + t16 * t48) * t44) / t7 - ((-t20 * t31 + (-t14 - t42) * t21) * t47 + (t15 * t18 + t16 * t7) * t44) * t6 * t48) * pkin(2) * t17 * t29 / (t8 ^ 2 * t6 + 0.1e1);
	t5 = [t27 * r_i_i_C(3) + t34 * t25, t32, t1 * t32, 0; t25 * r_i_i_C(3) - t34 * t27, t33, t1 * t33, 0; 0, t35, t35 * t1, 0;];
	Ja_transl = t5;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (5534->46), mult. (8150->80), div. (322->11), fcn. (5102->14), ass. (0->58)
	t33 = sin(pkin(20));
	t34 = cos(pkin(20));
	t35 = sin(qJ(3));
	t37 = cos(qJ(3));
	t31 = t37 * t33 + t35 * t34;
	t63 = pkin(6) * t31;
	t56 = pkin(1) * t63;
	t30 = 0.2e1 * t56;
	t42 = pkin(2) ^ 2;
	t41 = pkin(6) ^ 2;
	t57 = pkin(1) ^ 2 + t41;
	t54 = -pkin(13) ^ 2 + t57;
	t25 = t30 + t42 + t54;
	t29 = -pkin(1) - t63;
	t32 = t35 * t33 - t37 * t34;
	t58 = t30 + t41;
	t68 = -pkin(2) - pkin(13);
	t21 = (pkin(1) - t68) * (pkin(1) + t68) + t58;
	t67 = -pkin(2) + pkin(13);
	t22 = (pkin(1) - t67) * (pkin(1) + t67) + t58;
	t60 = t22 * t21;
	t45 = sqrt(-t60);
	t59 = t32 * t45;
	t55 = pkin(6) * t59;
	t15 = -t29 * t25 - t55;
	t62 = pkin(6) * t32;
	t16 = t25 * t62 - t29 * t45;
	t28 = t30 + t57;
	t26 = 0.1e1 / t28;
	t43 = 0.1e1 / pkin(2);
	t69 = t43 / 0.2e1;
	t53 = t26 * t69;
	t11 = qJ(2) + atan2(t16 * t53, t15 * t53);
	t24 = t42 - t54 - 0.2e1 * t56;
	t52 = 0.1e1 / pkin(13) * t69;
	t7 = atan2(t45 * t52, t24 * t52) + t11;
	t5 = sin(t7);
	t6 = cos(t7);
	t49 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t71 = pkin(2) * sin(t11);
	t72 = t49 - t71;
	t70 = t26 / 0.2e1;
	t14 = 0.1e1 / t15 ^ 2;
	t27 = 0.1e1 / t28 ^ 2;
	t61 = 0.1e1 / t45 * (t21 + t22) * pkin(1) * t62;
	t64 = pkin(6) * t16;
	t66 = pkin(1) * t32;
	t2 = 0.2e1 * (((-t29 * t61 + (t31 * t25 - t59) * pkin(6)) * t70 + (-t26 * t41 * t32 + t27 * t64) * t66) / t15 - ((-t31 * t45 + (-t25 - t61) * t32) * t70 + (t15 * t27 + t26 * t29) * t66) * t14 * t64) * pkin(2) / (t16 ^ 2 * t14 + 0.1e1) * t28 * t43;
	t65 = pkin(2) * cos(t11);
	t50 = r_i_i_C(1) * t6 - r_i_i_C(2) * t5;
	t48 = pkin(16) + t72;
	t47 = t50 - t65;
	t23 = 0.1e1 / t24 ^ 2;
	t1 = (0.1e1 / t24 * t61 - 0.2e1 * pkin(1) * t23 * t55) / (-t23 * t60 + 0.1e1) + t2;
	t46 = t50 * t1 - t2 * t65;
	t38 = cos(qJ(1));
	t36 = sin(qJ(1));
	t3 = [t38 * r_i_i_C(3) - t48 * t36, t47 * t38, t46 * t38, 0; t36 * r_i_i_C(3) + t48 * t38, t47 * t36, t46 * t36, 0; 0, t72, t49 * t1 - t2 * t71, 0;];
	Ja_transl = t3;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:56
	% DurationCPUTime: 5.92s
	% Computational Cost: add. (163697->96), mult. (247221->187), div. (10738->14), fcn. (156847->22), ass. (0->110)
	t137 = 2 * pkin(4);
	t135 = -2 * pkin(1);
	t117 = sin(pkin(19));
	t73 = sin(qJ(2));
	t76 = cos(qJ(2));
	t78 = cos(pkin(19));
	t65 = t73 * t117 + t76 * t78;
	t118 = pkin(7) * t65;
	t64 = -t76 * t117 + t73 * t78;
	t119 = pkin(7) * t64;
	t104 = (pkin(1) ^ 2) + t119 * t135;
	t80 = pkin(7) ^ 2;
	t61 = t80 + t104;
	t58 = pkin(3) ^ 2 - pkin(8) ^ 2 + t61;
	t62 = pkin(1) - t119;
	t128 = -pkin(8) - pkin(3);
	t56 = (pkin(7) - t128) * (pkin(7) + t128) + t104;
	t127 = -pkin(8) + pkin(3);
	t57 = (pkin(7) - t127) * (pkin(7) + t127) + t104;
	t86 = sqrt(-t57 * t56);
	t52 = t58 * t118 + t62 * t86;
	t75 = cos(qJ(3));
	t106 = t75 * t52;
	t113 = t65 * t86;
	t51 = -pkin(7) * t113 + t58 * t62;
	t72 = sin(qJ(3));
	t109 = t72 * t51;
	t59 = 0.1e1 / t61;
	t83 = 0.1e1 / pkin(3);
	t114 = t59 * t83;
	t47 = (t109 / 0.2e1 + t106 / 0.2e1) * t114;
	t107 = t75 * t51;
	t108 = t72 * t52;
	t48 = (-t107 / 0.2e1 + t108 / 0.2e1) * t114;
	t69 = pkin(23) + pkin(22);
	t67 = sin(t69);
	t68 = cos(t69);
	t32 = t47 * t68 - t48 * t67;
	t121 = t32 * pkin(5);
	t81 = pkin(5) ^ 2;
	t105 = t121 * t137 + t81;
	t126 = (-pkin(9) - pkin(11));
	t23 = ((pkin(4) - t126) * (pkin(4) + t126)) + t105;
	t125 = (-pkin(9) + pkin(11));
	t24 = ((pkin(4) - t125) * (pkin(4) + t125)) + t105;
	t85 = sqrt(-t24 * t23);
	t134 = -0.1e1 / t85 / 0.2e1;
	t28 = (pkin(4) ^ 2) + t105;
	t26 = 0.1e1 / t28;
	t133 = t26 / 0.2e1;
	t103 = pkin(1) * t118;
	t115 = 0.2e1 / t86 * (t56 + t57) * t103;
	t40 = (t64 * t86 + (-t115 / 0.2e1 - t58 + t62 * t135) * t65) * pkin(7);
	t132 = -t40 / 0.2e1;
	t41 = t62 * t115 / 0.2e1 + t80 * t65 ^ 2 * t135 + (-t64 * t58 - t113) * pkin(7);
	t131 = t41 / 0.2e1;
	t70 = sin(pkin(23));
	t130 = t70 / 0.2e1;
	t129 = -t75 / 0.2e1;
	t124 = pkin(1) * t73;
	t71 = cos(pkin(23));
	t111 = t71 * t51;
	t112 = t70 * t52;
	t45 = (-t111 / 0.2e1 + t112 / 0.2e1) * t114;
	t110 = t71 * t52;
	t116 = t51 * t70;
	t46 = (t110 / 0.2e1 + t116 / 0.2e1) * t114;
	t37 = pkin(22) - atan2(t46, t45) - qJ(2);
	t123 = pkin(4) * sin(t37);
	t97 = t47 * t67 + t68 * t48;
	t120 = pkin(5) * t97;
	t25 = pkin(9) ^ 2 - pkin(11) ^ 2 + t28;
	t29 = -pkin(4) - t121;
	t13 = t25 * t120 - t29 * t85;
	t122 = pkin(5) * t13;
	t102 = t97 * t134;
	t79 = 0.1e1 / pkin(9);
	t101 = t79 * t133;
	t100 = t29 * t134;
	t99 = 0.1e1 / t61 ^ 2 * t103;
	t12 = -t85 * t120 - t25 * t29;
	t11 = 0.1e1 / t12 ^ 2;
	t98 = pkin(9) * t28 / (t11 * t13 ^ 2 + 0.1e1) * t79;
	t96 = 0.1e1 / t12 * t98;
	t5 = -atan2(t13 * t101, t12 * t101) + t37;
	t3 = sin(t5);
	t4 = cos(t5);
	t95 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
	t94 = r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t44 = 0.1e1 / t45 ^ 2;
	t9 = -0.1e1 + (-((t40 * t130 + t71 * t131) * t59 + (t110 + t116) * t99) / t45 + ((t41 * t130 + t71 * t132) * t59 + (-t111 + t112) * t99) * t46 * t44) / (t44 * t46 ^ 2 + 0.1e1) * t83;
	t93 = -pkin(16) - t123 + t124 + t94;
	t91 = pkin(5) * (t23 + t24) * t137;
	t14 = t97 * t91;
	t27 = 0.1e1 / t28 ^ 2;
	t87 = pkin(4) * (-t26 * t81 * t97 + t27 * t122);
	t89 = pkin(4) * (t12 * t27 + t26 * t29);
	t90 = t11 * t98 * t122;
	t2 = -0.2e1 * ((t14 * t100 + (t32 * t25 - t85 * t97) * pkin(5)) * t133 + t97 * t87) * t96 + 0.2e1 * ((t14 * t102 - t25 * t97 - t32 * t85) * t133 + t97 * t89) * t90;
	t92 = t2 * t95;
	t21 = ((t41 * t129 + t72 * t132) * t59 + (-t106 - t109) * t99) * t83;
	t22 = ((t40 * t129 + t72 * t131) * t59 + (-t107 + t108) * t99) * t83;
	t16 = t21 * t68 + t22 * t67;
	t17 = -t21 * t67 + t22 * t68;
	t7 = t16 * t91;
	t1 = -0.2e1 * ((t7 * t100 + (-t16 * t85 + t17 * t25) * pkin(5)) * t133 + t16 * t87) * t96 + 0.2e1 * ((t7 * t102 - t16 * t25 - t17 * t85) * t133 + t16 * t89) * t90 + t9;
	t88 = t95 * t1 + pkin(4) * t9 * cos(t37) - t76 * pkin(1);
	t77 = cos(qJ(1));
	t74 = sin(qJ(1));
	t6 = [r_i_i_C(3) * t77 + t93 * t74, t88 * t77, t77 * t92, 0; r_i_i_C(3) * t74 - t93 * t77, t88 * t74, t74 * t92, 0; 0, t94 * t1 - t9 * t123 - t124, t94 * t2, 0;];
	Ja_transl = t6;
end