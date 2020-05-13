% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m1DE2
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
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m1DE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m1DE2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_jacobia_transl_sym_varpar: pkin has to be [19x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:42
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
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(13) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) - t5 * t2, t6 * t4, 0, 0; t2 * r_i_i_C(3) + t5 * t4, t6 * t2, 0, 0; 0, t7, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->9), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->13)
	t10 = qJ(2) + qJ(3);
	t7 = sin(t10);
	t8 = cos(t10);
	t23 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t15 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t22 = t15 + cos(qJ(2)) * pkin(1);
	t12 = sin(qJ(1));
	t18 = t23 * t12;
	t13 = cos(qJ(1));
	t17 = t23 * t13;
	t16 = sin(qJ(2)) * pkin(1);
	t14 = -pkin(13) - t22;
	t1 = [t13 * r_i_i_C(3) + t14 * t12, -t13 * t16 + t17, t17, 0; t12 * r_i_i_C(3) - t14 * t13, -t12 * t16 + t18, t18, 0; 0, t22, t15, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:21:19
	% DurationCPUTime: 11.24s
	% Computational Cost: add. (322885->98), mult. (487477->179), div. (21464->11), fcn. (309017->21), ass. (0->116)
	t137 = -2 * pkin(1);
	t136 = -2 * pkin(4);
	t135 = pkin(3) * pkin(4);
	t68 = sin(pkin(17));
	t134 = t68 / 0.2e1;
	t70 = sin(qJ(3));
	t133 = t70 / 0.2e1;
	t132 = -pkin(6) - pkin(2);
	t131 = -pkin(6) + pkin(2);
	t130 = (-pkin(8) - pkin(10));
	t129 = (-pkin(8) + pkin(10));
	t121 = cos(pkin(16));
	t71 = sin(qJ(2));
	t73 = sin(pkin(16));
	t75 = cos(qJ(2));
	t56 = -t75 * t121 + t71 * t73;
	t124 = pkin(5) * t56;
	t109 = (pkin(1) ^ 2) + t124 * t137;
	t78 = pkin(5) ^ 2;
	t53 = t78 + t109;
	t51 = 0.1e1 / t53;
	t81 = 0.1e1 / pkin(2);
	t112 = t51 * t81;
	t57 = t71 * t121 + t75 * t73;
	t123 = pkin(5) * t57;
	t50 = pkin(2) ^ 2 - pkin(6) ^ 2 + t53;
	t54 = pkin(1) - t124;
	t48 = (pkin(5) - t132) * (pkin(5) + t132) + t109;
	t49 = (pkin(5) - t131) * (pkin(5) + t131) + t109;
	t84 = sqrt(-t49 * t48);
	t44 = t50 * t123 + t54 * t84;
	t115 = t44 * t70;
	t111 = t57 * t84;
	t43 = -pkin(5) * t111 + t54 * t50;
	t74 = cos(qJ(3));
	t116 = t43 * t74;
	t39 = (-t116 / 0.2e1 + t115 / 0.2e1) * t112;
	t114 = t44 * t74;
	t117 = t43 * t70;
	t40 = (t114 / 0.2e1 + t117 / 0.2e1) * t112;
	t66 = pkin(18) + pkin(19);
	t61 = sin(t66);
	t62 = cos(t66);
	t89 = t61 * t39 - t62 * t40;
	t128 = pkin(3) * t89;
	t34 = t62 * t39 + t61 * t40;
	t127 = pkin(3) * t34;
	t67 = qJ(2) + qJ(3);
	t126 = pkin(4) * sin(t67);
	t125 = pkin(4) * cos(t67);
	t110 = (pkin(4) ^ 2) - t127 * t136;
	t80 = pkin(3) ^ 2;
	t30 = t80 + t110;
	t28 = 0.1e1 / t30;
	t77 = 0.1e1 / pkin(10);
	t118 = t28 * t77;
	t27 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t30;
	t31 = pkin(4) + t127;
	t25 = ((pkin(3) - t130) * (pkin(3) + t130)) + t110;
	t26 = ((pkin(3) - t129) * (pkin(3) + t129)) + t110;
	t83 = sqrt(-t26 * t25);
	t16 = -t83 * t128 + t31 * t27;
	t17 = t27 * t128 + t31 * t83;
	t69 = cos(pkin(17));
	t13 = (-t16 * t69 / 0.2e1 + t17 * t134) * t118;
	t12 = 0.1e1 / t13 ^ 2;
	t14 = (t17 * t69 / 0.2e1 + t16 * t134) * t118;
	t122 = t77 / (t14 ^ 2 * t12 + 0.1e1);
	t120 = t12 * t14;
	t119 = t28 * t69;
	t108 = pkin(1) * t123;
	t113 = 0.2e1 / t84 * (t48 + t49) * t108;
	t107 = 0.1e1 / t30 ^ 2 * t135;
	t22 = 0.1e1 / t83;
	t106 = t22 * t31 / 0.2e1;
	t105 = -t22 * t89 / 0.2e1;
	t104 = t28 * t134;
	t103 = -t119 / 0.2e1;
	t102 = t119 / 0.2e1;
	t101 = t80 * t89 * t136;
	t100 = t31 * t136 - t27;
	t99 = t75 * pkin(1) - t125;
	t98 = 0.1e1 / t53 ^ 2 * t108;
	t97 = t68 * t107;
	t96 = t69 * t107;
	t37 = (t56 * t84 + (-t113 / 0.2e1 - t50 + t54 * t137) * t57) * pkin(5);
	t38 = t54 * t113 / 0.2e1 + t78 * t57 ^ 2 * t137 + (-t56 * t50 - t111) * pkin(5);
	t23 = ((t38 * t74 / 0.2e1 + t37 * t133) * t51 + (t114 + t117) * t98) * t81;
	t24 = ((-t37 * t74 / 0.2e1 + t38 * t133) * t51 + (t115 - t116) * t98) * t81;
	t20 = -t61 * t23 - t62 * t24;
	t95 = t20 * t97;
	t94 = t89 * t97;
	t93 = t20 * t96;
	t92 = t89 * t96;
	t7 = atan2(t14, t13) + t67;
	t5 = sin(t7);
	t6 = cos(t7);
	t91 = -r_i_i_C(1) * t6 + r_i_i_C(2) * t5;
	t90 = r_i_i_C(1) * t5 + r_i_i_C(2) * t6;
	t88 = -pkin(13) - t99 - t91;
	t87 = 0.2e1 * (t25 + t26) * t135;
	t11 = 0.1e1 / t13;
	t15 = t20 * t87;
	t19 = -t62 * t23 + t61 * t24;
	t3 = (t100 * t20 + t15 * t105 - t19 * t83) * pkin(3);
	t4 = t15 * t106 + t20 * t101 + (t19 * t27 - t20 * t83) * pkin(3);
	t1 = 0.1e1 + ((t4 * t102 + t3 * t104 + t16 * t95 + t17 * t93) * t11 - (t3 * t103 + t4 * t104 - t16 * t93 + t17 * t95) * t120) * t122;
	t86 = -t71 * pkin(1) + t90 * t1 + t126;
	t18 = t89 * t87;
	t10 = t18 * t106 + t89 * t101 + (t34 * t27 - t83 * t89) * pkin(3);
	t9 = (t100 * t89 + t18 * t105 - t34 * t83) * pkin(3);
	t2 = 0.1e1 + ((t10 * t102 + t9 * t104 + t16 * t94 + t17 * t92) * t11 - (t10 * t104 + t9 * t103 - t16 * t92 + t17 * t94) * t120) * t122;
	t85 = t90 * t2 + t126;
	t76 = cos(qJ(1));
	t72 = sin(qJ(1));
	t8 = [t76 * r_i_i_C(3) + t88 * t72, t86 * t76, t85 * t76, 0; t72 * r_i_i_C(3) - t88 * t76, t86 * t72, t85 * t72, 0; 0, t91 * t1 + t99, t91 * t2 - t125, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:21:54
	% DurationCPUTime: 27.09s
	% Computational Cost: add. (823996->112), mult. (1243899->199), div. (54860->11), fcn. (788475->23), ass. (0->129)
	t134 = cos(pkin(16));
	t76 = sin(qJ(2));
	t78 = sin(pkin(16));
	t81 = cos(qJ(2));
	t60 = -t134 * t81 + t76 * t78;
	t136 = pkin(5) * t60;
	t150 = -2 * pkin(1);
	t117 = (pkin(1) ^ 2) + t136 * t150;
	t84 = pkin(5) ^ 2;
	t57 = t84 + t117;
	t55 = 0.1e1 / t57;
	t87 = 0.1e1 / pkin(2);
	t124 = t55 * t87;
	t61 = t134 * t76 + t81 * t78;
	t135 = pkin(5) * t61;
	t54 = pkin(2) ^ 2 - pkin(6) ^ 2 + t57;
	t58 = pkin(1) - t136;
	t145 = -pkin(6) - pkin(2);
	t52 = (pkin(5) - t145) * (pkin(5) + t145) + t117;
	t144 = -pkin(6) + pkin(2);
	t53 = (pkin(5) - t144) * (pkin(5) + t144) + t117;
	t90 = sqrt(-t53 * t52);
	t48 = t135 * t54 + t58 * t90;
	t75 = sin(qJ(3));
	t127 = t48 * t75;
	t123 = t61 * t90;
	t47 = -pkin(5) * t123 + t58 * t54;
	t80 = cos(qJ(3));
	t128 = t47 * t80;
	t43 = (-t128 / 0.2e1 + t127 / 0.2e1) * t124;
	t126 = t48 * t80;
	t129 = t47 * t75;
	t44 = (t126 / 0.2e1 + t129 / 0.2e1) * t124;
	t70 = pkin(18) + pkin(19);
	t65 = sin(t70);
	t66 = cos(t70);
	t38 = t66 * t43 + t65 * t44;
	t139 = pkin(3) * t38;
	t149 = -2 * pkin(4);
	t118 = (pkin(4) ^ 2) - t139 * t149;
	t86 = pkin(3) ^ 2;
	t34 = t86 + t118;
	t32 = 0.1e1 / t34;
	t83 = 0.1e1 / pkin(10);
	t130 = t32 * t83;
	t72 = sin(pkin(17));
	t147 = t72 / 0.2e1;
	t97 = t65 * t43 - t66 * t44;
	t140 = pkin(3) * t97;
	t31 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t34;
	t35 = pkin(4) + t139;
	t143 = -pkin(8) - pkin(10);
	t29 = (pkin(3) - t143) * (pkin(3) + t143) + t118;
	t142 = -pkin(8) + pkin(10);
	t30 = (pkin(3) - t142) * (pkin(3) + t142) + t118;
	t89 = sqrt(-t30 * t29);
	t20 = -t140 * t89 + t35 * t31;
	t21 = t140 * t31 + t35 * t89;
	t73 = cos(pkin(17));
	t17 = (-t20 * t73 / 0.2e1 + t21 * t147) * t130;
	t18 = (t21 * t73 / 0.2e1 + t20 * t147) * t130;
	t71 = qJ(2) + qJ(3);
	t11 = atan2(t18, t17) + t71;
	t10 = cos(t11);
	t137 = pkin(4) * cos(t71);
	t106 = t81 * pkin(1) - t137;
	t141 = pkin(11) + r_i_i_C(3);
	t9 = sin(t11);
	t114 = t141 * t9;
	t151 = pkin(9) * t10 - pkin(13) - t106 + t114;
	t148 = pkin(3) * pkin(4);
	t146 = t75 / 0.2e1;
	t138 = pkin(4) * sin(t71);
	t16 = 0.1e1 / t17 ^ 2;
	t133 = 0.1e1 / (t18 ^ 2 * t16 + 0.1e1) * t83;
	t132 = t16 * t18;
	t131 = t32 * t73;
	t116 = pkin(1) * t135;
	t125 = 0.2e1 / t90 * (t52 + t53) * t116;
	t74 = sin(qJ(4));
	t77 = sin(qJ(1));
	t122 = t77 * t74;
	t79 = cos(qJ(4));
	t121 = t77 * t79;
	t82 = cos(qJ(1));
	t120 = t82 * t74;
	t119 = t82 * t79;
	t115 = 0.1e1 / t34 ^ 2 * t148;
	t26 = 0.1e1 / t89;
	t113 = t26 * t35 / 0.2e1;
	t112 = -t26 * t97 / 0.2e1;
	t111 = t32 * t147;
	t110 = -t131 / 0.2e1;
	t109 = t131 / 0.2e1;
	t108 = t86 * t97 * t149;
	t107 = t149 * t35 - t31;
	t104 = 0.1e1 / t57 ^ 2 * t116;
	t103 = t72 * t115;
	t102 = t73 * t115;
	t41 = (t60 * t90 + (-t125 / 0.2e1 - t54 + t58 * t150) * t61) * pkin(5);
	t42 = t58 * t125 / 0.2e1 + t84 * t61 ^ 2 * t150 + (-t60 * t54 - t123) * pkin(5);
	t27 = ((t42 * t80 / 0.2e1 + t41 * t146) * t55 + (t126 + t129) * t104) * t87;
	t28 = ((-t41 * t80 / 0.2e1 + t42 * t146) * t55 + (t127 - t128) * t104) * t87;
	t24 = -t65 * t27 - t66 * t28;
	t101 = t24 * t103;
	t100 = t97 * t103;
	t99 = t24 * t102;
	t98 = t97 * t102;
	t96 = r_i_i_C(1) * t79 - r_i_i_C(2) * t74 + pkin(9);
	t95 = 0.2e1 * (t29 + t30) * t148;
	t94 = -t10 * t141 + t9 * t96;
	t93 = -t10 * t96 - t114;
	t15 = 0.1e1 / t17;
	t19 = t24 * t95;
	t23 = -t66 * t27 + t65 * t28;
	t3 = (t107 * t24 + t112 * t19 - t23 * t89) * pkin(3);
	t4 = t19 * t113 + t24 * t108 + (t23 * t31 - t24 * t89) * pkin(3);
	t1 = 0.1e1 + ((t101 * t20 + t109 * t4 + t111 * t3 + t21 * t99) * t15 - (t101 * t21 + t110 * t3 + t111 * t4 - t20 * t99) * t132) * t133;
	t92 = -t76 * pkin(1) + t1 * t94 + t138;
	t22 = t97 * t95;
	t13 = (t107 * t97 + t112 * t22 - t38 * t89) * pkin(3);
	t14 = t22 * t113 + t97 * t108 + (t38 * t31 - t89 * t97) * pkin(3);
	t2 = 0.1e1 + ((t100 * t20 + t109 * t14 + t111 * t13 + t21 * t98) * t15 - (t100 * t21 + t110 * t13 + t111 * t14 - t20 * t98) * t132) * t133;
	t91 = t2 * t94 + t138;
	t8 = t10 * t119 - t122;
	t7 = t10 * t120 + t121;
	t6 = t10 * t121 + t120;
	t5 = t10 * t122 - t119;
	t12 = [t6 * r_i_i_C(1) - t5 * r_i_i_C(2) + t151 * t77, t92 * t82, t91 * t82, r_i_i_C(1) * t7 + r_i_i_C(2) * t8; -t8 * r_i_i_C(1) + t7 * r_i_i_C(2) - t151 * t82, t92 * t77, t91 * t77, r_i_i_C(1) * t5 + r_i_i_C(2) * t6; 0, t1 * t93 + t106, t2 * t93 - t137, (r_i_i_C(1) * t74 + r_i_i_C(2) * t79) * t9;];
	Ja_transl = t12;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:43
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (5669->40), mult. (8752->77), div. (380->7), fcn. (5636->13), ass. (0->50)
	t57 = -2 * pkin(5);
	t32 = cos(pkin(15));
	t56 = t32 / 0.2e1;
	t55 = (-pkin(6) - pkin(2));
	t54 = (-pkin(6) + pkin(2));
	t27 = sin(qJ(2));
	t29 = sin(pkin(16));
	t50 = cos(qJ(2));
	t51 = cos(pkin(16));
	t24 = t27 * t29 - t50 * t51;
	t53 = pkin(1) * t24;
	t25 = t27 * t51 + t50 * t29;
	t52 = pkin(1) * t25;
	t35 = pkin(1) ^ 2;
	t45 = t53 * t57 + t35;
	t16 = ((pkin(5) - t55) * (pkin(5) + t55)) + t45;
	t17 = ((pkin(5) - t54) * (pkin(5) + t54)) + t45;
	t36 = sqrt(-t17 * t16);
	t44 = pkin(5) * t52;
	t49 = 0.1e1 / t36 * (t16 + t17) * t44;
	t21 = (pkin(5) ^ 2) + t45;
	t19 = 0.1e1 / t21;
	t30 = sin(pkin(15));
	t48 = t19 * t30;
	t33 = 0.1e1 / pkin(6);
	t47 = t19 * t33;
	t46 = t25 * t36;
	t43 = t19 * t56;
	t42 = 0.1e1 / t21 ^ 2 * t44;
	t18 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t21;
	t22 = -pkin(5) + t53;
	t11 = -pkin(1) * t46 - t18 * t22;
	t41 = t11 * t42;
	t12 = t18 * t52 - t22 * t36;
	t40 = t12 * t42;
	t10 = (t12 * t56 - t11 * t30 / 0.2e1) * t47;
	t9 = (t11 * t56 + t12 * t30 / 0.2e1) * t47;
	t5 = atan2(t10, t9);
	t2 = sin(t5);
	t3 = cos(t5);
	t39 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t38 = -pkin(7) + t39;
	t6 = (t24 * t36 + (0.2e1 * t22 * pkin(5) - t18 - t49) * t25) * pkin(1);
	t7 = -t22 * t49 + t35 * t25 ^ 2 * t57 + (-t18 * t24 - t46) * pkin(1);
	t8 = 0.1e1 / t9 ^ 2;
	t1 = ((t7 * t43 + t32 * t40 - t6 * t48 / 0.2e1 - t30 * t41) / t9 - (t6 * t43 + t32 * t41 + t7 * t48 / 0.2e1 + t30 * t40) * t10 * t8) / (t10 ^ 2 * t8 + 0.1e1) * t33;
	t37 = t1 * (-r_i_i_C(1) * t2 - r_i_i_C(2) * t3);
	t31 = cos(qJ(1));
	t28 = sin(qJ(1));
	t4 = [t31 * r_i_i_C(3) - t38 * t28, t31 * t37, 0, 0; t28 * r_i_i_C(3) + t38 * t31, t28 * t37, 0, 0; 0, t39 * t1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:43
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (5690->44), mult. (8759->78), div. (380->7), fcn. (5643->13), ass. (0->51)
	t59 = -2 * pkin(1);
	t29 = sin(pkin(19));
	t58 = t29 / 0.2e1;
	t57 = -pkin(6) - pkin(2);
	t56 = -pkin(6) + pkin(2);
	t31 = sin(qJ(2));
	t33 = sin(pkin(16));
	t34 = cos(qJ(2));
	t53 = cos(pkin(16));
	t24 = t31 * t33 - t34 * t53;
	t55 = pkin(5) * t24;
	t25 = t31 * t53 + t34 * t33;
	t54 = pkin(5) * t25;
	t48 = (pkin(1) ^ 2) + t55 * t59;
	t16 = (pkin(5) - t57) * (pkin(5) + t57) + t48;
	t17 = (pkin(5) - t56) * (pkin(5) + t56) + t48;
	t39 = sqrt(-t17 * t16);
	t47 = pkin(1) * t54;
	t52 = 0.1e1 / t39 * (t16 + t17) * t47;
	t36 = pkin(5) ^ 2;
	t21 = t36 + t48;
	t19 = 0.1e1 / t21;
	t30 = cos(pkin(19));
	t51 = t19 * t30;
	t37 = 0.1e1 / pkin(2);
	t50 = t19 * t37;
	t49 = t25 * t39;
	t46 = t19 * t58;
	t45 = 0.1e1 / t21 ^ 2 * t47;
	t44 = t29 * t45;
	t43 = t30 * t45;
	t18 = pkin(2) ^ 2 - pkin(6) ^ 2 + t21;
	t22 = pkin(1) - t55;
	t11 = -pkin(5) * t49 + t22 * t18;
	t12 = t18 * t54 + t22 * t39;
	t10 = (t12 * t30 / 0.2e1 + t11 * t58) * t50;
	t9 = (-t11 * t30 / 0.2e1 + t12 * t58) * t50;
	t4 = qJ(2) + atan2(t10, t9);
	t2 = sin(t4);
	t3 = cos(t4);
	t42 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t28 = t34 * pkin(1);
	t41 = t28 + pkin(13) + t42;
	t6 = (t24 * t39 + (t22 * t59 - t18 - t52) * t25) * pkin(5);
	t7 = t22 * t52 + t36 * t25 ^ 2 * t59 + (-t24 * t18 - t49) * pkin(5);
	t8 = 0.1e1 / t9 ^ 2;
	t1 = 0.1e1 + ((t7 * t51 / 0.2e1 + t12 * t43 + t6 * t46 + t11 * t44) / t9 - (-t6 * t51 / 0.2e1 - t11 * t43 + t7 * t46 + t12 * t44) * t10 * t8) / (t10 ^ 2 * t8 + 0.1e1) * t37;
	t40 = -pkin(1) * t31 + (-r_i_i_C(1) * t2 - r_i_i_C(2) * t3) * t1;
	t35 = cos(qJ(1));
	t32 = sin(qJ(1));
	t5 = [t35 * r_i_i_C(3) - t41 * t32, t40 * t35, 0, 0; t32 * r_i_i_C(3) + t41 * t35, t40 * t32, 0, 0; 0, t42 * t1 + t28, 0, 0;];
	Ja_transl = t5;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:21:03
	% DurationCPUTime: 5.98s
	% Computational Cost: add. (163697->96), mult. (247221->187), div. (10738->14), fcn. (156847->22), ass. (0->110)
	t136 = 2 * pkin(3);
	t134 = -2 * pkin(1);
	t117 = cos(pkin(16));
	t73 = sin(qJ(2));
	t75 = sin(pkin(16));
	t77 = cos(qJ(2));
	t63 = -t117 * t77 + t73 * t75;
	t119 = pkin(5) * t63;
	t104 = (pkin(1) ^ 2) + t119 * t134;
	t80 = pkin(5) ^ 2;
	t60 = t80 + t104;
	t58 = 0.1e1 / t60;
	t83 = 0.1e1 / pkin(2);
	t107 = t58 * t83;
	t64 = t117 * t73 + t75 * t77;
	t118 = pkin(5) * t64;
	t57 = pkin(2) ^ 2 - pkin(6) ^ 2 + t60;
	t61 = pkin(1) - t119;
	t127 = -pkin(6) - pkin(2);
	t55 = (pkin(5) - t127) * (pkin(5) + t127) + t104;
	t126 = -pkin(6) + pkin(2);
	t56 = (pkin(5) - t126) * (pkin(5) + t126) + t104;
	t86 = sqrt(-t56 * t55);
	t51 = t118 * t57 + t61 * t86;
	t72 = sin(qJ(3));
	t110 = t51 * t72;
	t106 = t64 * t86;
	t50 = -pkin(5) * t106 + t57 * t61;
	t76 = cos(qJ(3));
	t113 = t50 * t76;
	t46 = (-t113 / 0.2e1 + t110 / 0.2e1) * t107;
	t109 = t51 * t76;
	t114 = t50 * t72;
	t47 = (t109 / 0.2e1 + t114 / 0.2e1) * t107;
	t69 = pkin(18) + pkin(19);
	t66 = sin(t69);
	t67 = cos(t69);
	t33 = t46 * t67 + t47 * t66;
	t120 = t33 * pkin(4);
	t81 = pkin(4) ^ 2;
	t105 = t120 * t136 + t81;
	t125 = (-pkin(8) - pkin(10));
	t23 = ((pkin(3) - t125) * (pkin(3) + t125)) + t105;
	t124 = (-pkin(8) + pkin(10));
	t24 = ((pkin(3) - t124) * (pkin(3) + t124)) + t105;
	t85 = sqrt(-t24 * t23);
	t133 = -0.1e1 / t85 / 0.2e1;
	t28 = (pkin(3) ^ 2) + t105;
	t26 = 0.1e1 / t28;
	t132 = t26 / 0.2e1;
	t103 = pkin(1) * t118;
	t108 = 0.2e1 / t86 * (t55 + t56) * t103;
	t40 = (t63 * t86 + (-t108 / 0.2e1 - t57 + t61 * t134) * t64) * pkin(5);
	t131 = -t40 / 0.2e1;
	t41 = t61 * t108 / 0.2e1 + t80 * t64 ^ 2 * t134 + (-t57 * t63 - t106) * pkin(5);
	t130 = t41 / 0.2e1;
	t70 = sin(pkin(19));
	t129 = t70 / 0.2e1;
	t128 = t72 / 0.2e1;
	t112 = t51 * t70;
	t71 = cos(pkin(19));
	t115 = t50 * t71;
	t44 = (-t115 / 0.2e1 + t112 / 0.2e1) * t107;
	t111 = t51 * t71;
	t116 = t50 * t70;
	t45 = (t111 / 0.2e1 + t116 / 0.2e1) * t107;
	t37 = pkin(18) - atan2(t45, t44) - qJ(2);
	t123 = pkin(3) * cos(t37);
	t94 = t46 * t66 - t47 * t67;
	t121 = pkin(4) * t94;
	t25 = pkin(8) ^ 2 - pkin(10) ^ 2 + t28;
	t29 = -pkin(3) - t120;
	t13 = t121 * t25 - t29 * t85;
	t122 = pkin(4) * t13;
	t102 = t94 * t133;
	t79 = 0.1e1 / pkin(8);
	t101 = t79 * t132;
	t100 = t29 * t133;
	t99 = 0.1e1 / t60 ^ 2 * t103;
	t12 = -t121 * t85 - t25 * t29;
	t11 = 0.1e1 / t12 ^ 2;
	t98 = pkin(8) * t28 / (t11 * t13 ^ 2 + 0.1e1) * t79;
	t97 = 0.1e1 / t12 * t98;
	t5 = -atan2(t13 * t101, t12 * t101) + t37;
	t3 = sin(t5);
	t4 = cos(t5);
	t96 = r_i_i_C(1) * t4 + r_i_i_C(2) * t3;
	t95 = r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t43 = 0.1e1 / t44 ^ 2;
	t9 = -0.1e1 + (-((t129 * t40 + t130 * t71) * t58 + (t111 + t116) * t99) / t44 + ((t129 * t41 + t131 * t71) * t58 + (t112 - t115) * t99) * t45 * t43) / (t43 * t45 ^ 2 + 0.1e1) * t83;
	t68 = t77 * pkin(1);
	t93 = -pkin(13) - t68 - t123 + t96;
	t91 = pkin(4) * (t23 + t24) * t136;
	t14 = t94 * t91;
	t27 = 0.1e1 / t28 ^ 2;
	t87 = pkin(3) * (-t26 * t81 * t94 + t122 * t27);
	t89 = pkin(3) * (t12 * t27 + t26 * t29);
	t90 = t11 * t98 * t122;
	t2 = -0.2e1 * ((t14 * t100 + (t25 * t33 - t85 * t94) * pkin(4)) * t132 + t94 * t87) * t97 + 0.2e1 * ((t102 * t14 - t25 * t94 - t33 * t85) * t132 + t94 * t89) * t90;
	t92 = t2 * t95;
	t21 = ((t128 * t40 + t130 * t76) * t58 + (t109 + t114) * t99) * t83;
	t22 = ((t128 * t41 + t131 * t76) * t58 + (t110 - t113) * t99) * t83;
	t16 = -t21 * t67 + t22 * t66;
	t17 = -t21 * t66 - t22 * t67;
	t7 = t17 * t91;
	t1 = -0.2e1 * ((t7 * t100 + (t16 * t25 - t17 * t85) * pkin(4)) * t132 + t17 * t87) * t97 + 0.2e1 * ((t102 * t7 - t16 * t85 - t17 * t25) * t132 + t17 * t89) * t90 + t9;
	t88 = t95 * t1 - pkin(3) * t9 * sin(t37) - t73 * pkin(1);
	t78 = cos(qJ(1));
	t74 = sin(qJ(1));
	t6 = [r_i_i_C(3) * t78 + t74 * t93, t88 * t78, t78 * t92, 0; r_i_i_C(3) * t74 - t78 * t93, t88 * t74, t74 * t92, 0; 0, t1 * t96 - t123 * t9 + t68, t96 * t2, 0;];
	Ja_transl = t6;
end