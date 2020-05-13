% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m1TE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m1TE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1TE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_jacobia_transl_sym_varpar: pkin has to be [23x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
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
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.05s
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
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->15), mult. (85->26), div. (0->0), fcn. (101->6), ass. (0->19)
	t13 = sin(qJ(1));
	t11 = sin(qJ(3));
	t12 = sin(qJ(2));
	t14 = cos(qJ(3));
	t15 = cos(qJ(2));
	t17 = t12 * t11 - t15 * t14;
	t5 = t17 * t13;
	t18 = t15 * t11 + t12 * t14;
	t6 = t18 * t13;
	t24 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t16 = cos(qJ(1));
	t7 = t17 * t16;
	t8 = t18 * t16;
	t23 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t22 = pkin(1) * t15;
	t21 = t12 * pkin(1);
	t20 = -t17 * r_i_i_C(1) - t18 * r_i_i_C(2);
	t19 = -pkin(16) + t21;
	t1 = [t5 * r_i_i_C(1) + t6 * r_i_i_C(2) + t16 * r_i_i_C(3) + t19 * t13, -t16 * t22 + t23, t23, 0; -t7 * r_i_i_C(1) - t8 * r_i_i_C(2) + t13 * r_i_i_C(3) - t19 * t16, -t13 * t22 + t24, t24, 0; 0, t20 - t21, t20, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:18:11
	% DurationCPUTime: 7.94s
	% Computational Cost: add. (205385->139), mult. (311257->229), div. (12992->8), fcn. (197929->16), ass. (0->120)
	t145 = sin(pkin(19));
	t79 = sin(qJ(2));
	t82 = cos(qJ(2));
	t84 = cos(pkin(19));
	t65 = -t145 * t82 + t79 * t84;
	t149 = pkin(7) * t65;
	t160 = -2 * pkin(1);
	t127 = (pkin(1) ^ 2) + t149 * t160;
	t86 = pkin(7) ^ 2;
	t56 = t86 + t127;
	t52 = 0.1e1 / t56;
	t89 = 0.1e1 / pkin(3);
	t139 = t52 * t89;
	t67 = t145 * t79 + t82 * t84;
	t148 = pkin(7) * t67;
	t51 = pkin(3) ^ 2 - pkin(8) ^ 2 + t56;
	t61 = pkin(1) - t149;
	t155 = -pkin(8) - pkin(3);
	t49 = (pkin(7) - t155) * (pkin(7) + t155) + t127;
	t156 = pkin(3) - pkin(8);
	t50 = (pkin(7) - t156) * (pkin(7) + t156) + t127;
	t91 = sqrt(-t50 * t49);
	t45 = t148 * t51 + t61 * t91;
	t81 = cos(qJ(3));
	t141 = t45 * t81;
	t138 = t67 * t91;
	t44 = -pkin(7) * t138 + t51 * t61;
	t78 = sin(qJ(3));
	t144 = t44 * t78;
	t40 = (t144 / 0.2e1 + t141 / 0.2e1) * t139;
	t142 = t45 * t78;
	t143 = t44 * t81;
	t41 = (-t143 / 0.2e1 + t142 / 0.2e1) * t139;
	t75 = pkin(23) + pkin(22);
	t73 = sin(t75);
	t74 = cos(t75);
	t111 = t40 * t74 - t41 * t73;
	t101 = pkin(4) * t111;
	t32 = pkin(5) + t101;
	t171 = t32 / 0.2e1;
	t76 = sin(pkin(21));
	t170 = t76 / 0.2e1;
	t77 = cos(pkin(21));
	t169 = -t77 / 0.2e1;
	t114 = t40 * t73 + t41 * t74;
	t151 = pkin(4) * t114;
	t168 = -t151 / 0.2e1;
	t158 = pkin(4) * pkin(5);
	t125 = -0.2e1 * t158;
	t128 = pkin(5) ^ 2 - t111 * t125;
	t88 = pkin(4) ^ 2;
	t31 = t88 + t128;
	t159 = 0.1e1 / t31;
	t126 = t159 / 0.2e1;
	t121 = t76 * t126;
	t153 = -pkin(9) + pkin(11);
	t154 = -pkin(9) - pkin(11);
	t92 = sqrt(-((pkin(4) - t153) * (pkin(4) + t153) + t128) * ((pkin(4) - t154) * (pkin(4) + t154) + t128));
	t150 = pkin(4) * t92;
	t161 = pkin(9) ^ 2;
	t29 = pkin(11) ^ 2 - t161 + t31;
	t25 = -t114 * t150 + t29 * t32;
	t26 = t151 * t29 + t32 * t92;
	t85 = 0.1e1 / pkin(11);
	t21 = (t126 * t26 * t77 + t121 * t25) * t85;
	t22 = (t159 * t169 * t25 + t121 * t26) * t85;
	t137 = t78 * t79;
	t66 = t81 * t82 - t137;
	t80 = sin(qJ(1));
	t57 = t66 * t80;
	t108 = t78 * t82 + t79 * t81;
	t58 = t108 * t80;
	t133 = t21 * t58 - t22 * t57;
	t83 = cos(qJ(1));
	t135 = t82 * t83;
	t59 = t135 * t81 - t137 * t83;
	t60 = t108 * t83;
	t167 = -t60 * t21 + t59 * t22;
	t131 = -t59 * t21 - t60 * t22;
	t115 = 0.1e1 / t31 ^ 2 * t85 * t158;
	t112 = t26 * t115;
	t113 = t25 * t115;
	t166 = t76 * t112 - t77 * t113;
	t165 = t77 * t112 + t76 * t113;
	t164 = -0.2e1 * t88 * t114 * pkin(5) - t150;
	t163 = -pkin(4) * t29 + t32 * t125;
	t162 = 0.4e1 / t92 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t128 - t161) * t158;
	t157 = -t81 / 0.2e1;
	t147 = t79 * pkin(1);
	t146 = -t21 * t57 - t22 * t58;
	t123 = pkin(1) * t148;
	t140 = 0.2e1 / t91 * (t49 + t50) * t123;
	t134 = t85 * t159;
	t130 = -t108 * t21 + t22 * t66;
	t129 = -t108 * t22 - t21 * t66;
	t120 = t85 * t126;
	t117 = -pkin(16) + t147;
	t116 = 0.1e1 / t56 ^ 2 * t123;
	t110 = t142 - t143;
	t109 = -t141 - t144;
	t37 = (t65 * t91 + (-t140 / 0.2e1 - t51 + t61 * t160) * t67) * pkin(7);
	t38 = t61 * t140 / 0.2e1 + t86 * t67 ^ 2 * t160 + (-t51 * t65 - t138) * pkin(7);
	t103 = -t78 * t37 / 0.2e1 + t38 * t157;
	t102 = t37 * t157 + t78 * t38 / 0.2e1;
	t27 = ((t102 * t73 + t103 * t74) * t52 + (t109 * t74 + t110 * t73) * t116) * t89;
	t99 = t27 * t162;
	t98 = t114 * t162;
	t97 = (-t92 * t101 + t114 * t163 + t168 * t98) * t134;
	t96 = (t29 * t101 + t114 * t164 + t171 * t98) * t120;
	t95 = pkin(4) * t89 * ((t102 * t74 - t103 * t73) * t52 + (-t109 * t73 + t110 * t74) * t116);
	t94 = (t163 * t27 + t168 * t99 - t92 * t95) * t134;
	t93 = (t164 * t27 + t171 * t99 + t29 * t95) * t120;
	t63 = t66 * pkin(5);
	t55 = t60 * pkin(5);
	t54 = t58 * pkin(5);
	t4 = t114 * t166 + t169 * t97 + t76 * t96;
	t3 = t114 * t165 + t170 * t97 + t77 * t96;
	t2 = t166 * t27 + t169 * t94 + t76 * t93;
	t1 = t165 * t27 + t170 * t94 + t77 * t93;
	t5 = [-t57 * pkin(5) + r_i_i_C(1) * t133 - r_i_i_C(2) * t146 + t83 * r_i_i_C(3) + t117 * t80, (-t1 * t60 + t2 * t59 + t131) * r_i_i_C(1) + (-t1 * t59 - t2 * t60 - t167) * r_i_i_C(2) - t55 - pkin(1) * t135, (-t3 * t60 + t4 * t59 + t131) * r_i_i_C(1) + (-t3 * t59 - t4 * t60 - t167) * r_i_i_C(2) - t55, 0; t59 * pkin(5) + r_i_i_C(1) * t167 + r_i_i_C(2) * t131 + t80 * r_i_i_C(3) - t117 * t83, (-t1 * t58 + t2 * t57 + t146) * r_i_i_C(1) + (-t1 * t57 - t2 * t58 + t133) * r_i_i_C(2) - t54 - t80 * t82 * pkin(1), (-t3 * t58 + t4 * t57 + t146) * r_i_i_C(1) + (-t3 * t57 - t4 * t58 + t133) * r_i_i_C(2) - t54, 0; 0, (t1 * t66 + t108 * t2 + t130) * r_i_i_C(1) + (-t1 * t108 + t2 * t66 + t129) * r_i_i_C(2) + t63 - t147, (t108 * t4 + t3 * t66 + t130) * r_i_i_C(1) + (-t108 * t3 + t4 * t66 + t129) * r_i_i_C(2) + t63, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:18:29
	% DurationCPUTime: 25.80s
	% Computational Cost: add. (530252->152), mult. (803355->245), div. (33680->8), fcn. (510749->18), ass. (0->129)
	t164 = -2 * pkin(5);
	t71 = pkin(3) ^ 2;
	t172 = -pkin(8) ^ 2 + t71;
	t147 = sin(qJ(2));
	t148 = sin(pkin(19));
	t150 = cos(qJ(2));
	t151 = cos(pkin(19));
	t119 = t147 * t151 - t150 * t148;
	t118 = pkin(7) * t119;
	t116 = (-0.2e1 * t118 + pkin(1)) * pkin(1);
	t168 = pkin(7) ^ 2;
	t50 = t116 + t168;
	t115 = t50 + t172;
	t114 = pkin(7) * t115;
	t117 = -t118 + pkin(1);
	t169 = 0.1e1 / pkin(3);
	t58 = t147 * t148 + t150 * t151;
	t159 = -pkin(8) + pkin(3);
	t160 = -pkin(8) - pkin(3);
	t70 = sqrt(-((pkin(7) - t159) * (pkin(7) + t159) + t116) * ((pkin(7) - t160) * (pkin(7) + t160) + t116));
	t113 = t169 * (t58 * t114 + t117 * t70);
	t162 = 0.1e1 / t50;
	t110 = t162 * t113 / 0.2e1;
	t153 = pkin(7) * t58;
	t112 = t169 * (t117 * t115 - t70 * t153);
	t111 = t162 * t112;
	t146 = sin(qJ(3));
	t149 = cos(qJ(3));
	t106 = t146 * t111 / 0.2e1 + t149 * t110;
	t178 = -t149 / 0.2e1;
	t107 = t146 * t110 + t111 * t178;
	t131 = pkin(23) + pkin(22);
	t125 = sin(t131);
	t126 = cos(t131);
	t99 = t126 * t106 - t125 * t107;
	t41 = (pkin(5) ^ 2) + (-t164 * t99 + pkin(4)) * pkin(4);
	t163 = 0.1e1 / t41;
	t135 = sin(pkin(21));
	t165 = 0.1e1 / pkin(11);
	t157 = -pkin(9) + pkin(11);
	t158 = -pkin(9) - pkin(11);
	t97 = pkin(4) * t99;
	t94 = (0.2e1 * t97 + pkin(5)) * pkin(5);
	t69 = sqrt(-((pkin(4) - t157) * (pkin(4) + t157) + t94) * ((pkin(4) - t158) * (pkin(4) + t158) + t94));
	t166 = pkin(9) ^ 2;
	t92 = pkin(11) ^ 2 - t166 + t41;
	t95 = t97 + pkin(5);
	t170 = t125 * t106 + t126 * t107;
	t98 = pkin(4) * t170;
	t85 = t165 * (-t69 * t98 + t95 * t92);
	t81 = t135 * t85;
	t136 = cos(pkin(21));
	t91 = pkin(4) * t92;
	t86 = t165 * (t170 * t91 + t95 * t69);
	t84 = t136 * t86;
	t35 = (t81 / 0.2e1 + t84 / 0.2e1) * t163;
	t82 = t136 * t85;
	t83 = t135 * t86;
	t36 = (-t82 / 0.2e1 + t83 / 0.2e1) * t163;
	t57 = -t147 * t146 + t150 * t149;
	t66 = sin(qJ(1));
	t51 = t57 * t66;
	t56 = -t150 * t146 - t147 * t149;
	t52 = t56 * t66;
	t22 = -t35 * t52 - t36 * t51;
	t65 = sin(qJ(4));
	t67 = cos(qJ(4));
	t68 = cos(qJ(1));
	t184 = t22 * t67 + t68 * t65;
	t183 = t22 * t65 - t67 * t68;
	t182 = t95 / 0.2e1;
	t181 = -t98 / 0.2e1;
	t180 = t135 / 0.2e1;
	t179 = -t136 / 0.2e1;
	t139 = -t51 * t35 + t52 * t36;
	t53 = t57 * t68;
	t54 = t56 * t68;
	t138 = -t53 * t35 + t54 * t36;
	t161 = pkin(4) * pkin(5);
	t174 = t161 / t41 ^ 2;
	t177 = (t81 + t84) * t174;
	t176 = (-t82 + t83) * t174;
	t175 = -0.2e1 * t95 * t161 - t91;
	t154 = pkin(4) * t69;
	t173 = pkin(4) ^ 2 * t164 * t170 - t154;
	t171 = 0.4e1 / t69 * (-t166 + (pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t94) * t161;
	t156 = -pkin(12) - r_i_i_C(3);
	t155 = t169 * t162;
	t152 = t165 * t163;
	t133 = pkin(1) * t153;
	t145 = 0.4e1 / t70 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t116 - t71) * t133;
	t108 = t112 * t133;
	t109 = t113 * t133;
	t121 = t155 * t178;
	t127 = t146 * t155;
	t43 = (t119 * t70 + (-t145 / 0.2e1 - t168 - t172) * t58) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t118) * t133;
	t44 = t117 * t145 / 0.2e1 - t119 * t114 + (-t70 - 0.2e1 * t133) * t153;
	t47 = 0.1e1 / t50 ^ 2;
	t102 = -t43 * t127 / 0.2e1 + t44 * t121 + (-t146 * t108 - t149 * t109) * t47;
	t103 = t43 * t121 + t44 * t127 / 0.2e1 + (-t149 * t108 + t146 * t109) * t47;
	t39 = t126 * t102 + t125 * t103;
	t130 = t152 / 0.2e1;
	t87 = pkin(4) * (-t125 * t102 + t126 * t103);
	t89 = t39 * t171;
	t73 = (t173 * t39 + t89 * t182 + t92 * t87) * t130;
	t74 = (t175 * t39 + t89 * t181 - t69 * t87) * t152;
	t13 = t136 * t73 + t177 * t39 + t74 * t180;
	t143 = t13 + t36;
	t14 = t135 * t73 + t176 * t39 + t74 * t179;
	t142 = t14 - t35;
	t88 = t170 * t171;
	t75 = (t170 * t173 + t88 * t182 + t99 * t91) * t130;
	t76 = (-t99 * t154 + t170 * t175 + t88 * t181) * t152;
	t15 = t136 * t75 + t170 * t177 + t76 * t180;
	t141 = t15 + t36;
	t16 = t135 * t75 + t170 * t176 + t76 * t179;
	t140 = t16 - t35;
	t137 = t56 * t35 + t57 * t36;
	t129 = t150 * pkin(1);
	t128 = t147 * pkin(1);
	t122 = t128 - pkin(16);
	t120 = r_i_i_C(1) * t67 - r_i_i_C(2) * t65 + pkin(10);
	t55 = t57 * pkin(5);
	t49 = t54 * pkin(5);
	t48 = t52 * pkin(5);
	t25 = t35 * t54 + t36 * t53;
	t20 = t25 * t67 + t65 * t66;
	t19 = -t25 * t65 + t66 * t67;
	t1 = [-t51 * pkin(5) + t22 * pkin(10) + t184 * r_i_i_C(1) - t183 * r_i_i_C(2) + t122 * t66 - t156 * t139, -t68 * t129 + t49 + t156 * (t142 * t54 - t143 * t53) + t120 * (t13 * t54 + t14 * t53 + t138), t49 + t156 * (t140 * t54 - t141 * t53) + t120 * (t15 * t54 + t16 * t53 + t138), r_i_i_C(1) * t19 - r_i_i_C(2) * t20; t53 * pkin(5) + t25 * pkin(10) + t20 * r_i_i_C(1) + t19 * r_i_i_C(2) - t122 * t68 + t156 * t138, -t66 * t129 + t48 + t156 * (t142 * t52 - t143 * t51) + t120 * (t13 * t52 + t14 * t51 + t139), t48 + t156 * (t140 * t52 - t141 * t51) + t120 * (t15 * t52 + t16 * t51 + t139), t183 * r_i_i_C(1) + t184 * r_i_i_C(2); 0, -t128 + t55 + t156 * (t142 * t57 + t143 * t56) + t120 * (t13 * t57 - t14 * t56 + t137), t55 + t156 * (t140 * t57 + t141 * t56) + t120 * (t15 * t57 - t16 * t56 + t137), (-r_i_i_C(1) * t65 - r_i_i_C(2) * t67) * (t35 * t57 - t36 * t56);];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (1345->39), mult. (2098->69), div. (64->4), fcn. (1356->10), ass. (0->42)
	t53 = -2 * pkin(7);
	t25 = sin(pkin(18));
	t52 = -t25 / 0.2e1;
	t28 = cos(pkin(18));
	t51 = t28 / 0.2e1;
	t50 = (-pkin(8) - pkin(3));
	t49 = (-pkin(8) + pkin(3));
	t23 = sin(qJ(2));
	t27 = cos(pkin(19));
	t44 = sin(pkin(19));
	t45 = cos(qJ(2));
	t20 = t23 * t27 - t45 * t44;
	t48 = pkin(1) * t20;
	t31 = pkin(1) ^ 2;
	t41 = t48 * t53 + t31;
	t12 = ((pkin(7) - t50) * (pkin(7) + t50)) + t41;
	t13 = ((pkin(7) - t49) * (pkin(7) + t49)) + t41;
	t32 = sqrt(-t13 * t12);
	t21 = t23 * t44 + t45 * t27;
	t46 = t21 * pkin(1);
	t40 = pkin(7) * t46;
	t47 = 0.1e1 / t32 * (t12 + t13) * t40;
	t17 = (pkin(7) ^ 2) + t41;
	t29 = 0.1e1 / pkin(8);
	t43 = 0.1e1 / t17 * t29;
	t42 = t21 * t32;
	t39 = 0.1e1 / t17 ^ 2 * t29 * t40;
	t14 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t17;
	t18 = -pkin(7) + t48;
	t33 = (-t18 * t47 + t31 * t21 ^ 2 * t53 + (-t20 * t14 - t42) * pkin(1)) * t43;
	t34 = pkin(1) * (t20 * t32 + (0.2e1 * t18 * pkin(7) - t14 - t47) * t21) * t43 / 0.2e1;
	t36 = t28 * t39;
	t37 = t25 * t39;
	t7 = -pkin(1) * t42 - t18 * t14;
	t8 = t14 * t46 - t18 * t32;
	t1 = t28 * t34 + t33 * t52 + t7 * t36 - t8 * t37;
	t2 = t25 * t34 + t33 * t51 + t8 * t36 + t7 * t37;
	t38 = r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t35 = -pkin(15) + (r_i_i_C(1) * (t7 * t51 + t8 * t52) + r_i_i_C(2) * (-t28 * t8 / 0.2e1 + t7 * t52)) * t43;
	t26 = cos(qJ(1));
	t24 = sin(qJ(1));
	t3 = [t26 * r_i_i_C(3) - t35 * t24, t38 * t26, 0, 0; t24 * r_i_i_C(3) + t35 * t26, t38 * t24, 0, 0; 0, t2 * r_i_i_C(1) + t1 * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (3368->53), mult. (5239->83), div. (176->4), fcn. (3405->10), ass. (0->52)
	t13 = sin(qJ(2));
	t15 = sin(pkin(19));
	t16 = cos(qJ(2));
	t18 = cos(pkin(19));
	t38 = t13 * t18 - t16 * t15;
	t37 = pkin(7) * t38;
	t34 = (-0.2e1 * t37 + pkin(1)) * pkin(1);
	t56 = pkin(7) ^ 2;
	t11 = t34 + t56;
	t55 = 0.1e1 / t11;
	t65 = t55 / 0.2e1;
	t57 = 0.1e1 / pkin(3);
	t64 = t57 * t65;
	t44 = sin(pkin(23));
	t45 = cos(pkin(23));
	t20 = pkin(3) ^ 2;
	t32 = -pkin(8) ^ 2 + t11 + t20;
	t36 = -t37 + pkin(1);
	t12 = t13 * t15 + t16 * t18;
	t53 = -pkin(8) + pkin(3);
	t54 = -pkin(8) - pkin(3);
	t19 = sqrt(-((pkin(7) - t53) * (pkin(7) + t53) + t34) * ((pkin(7) - t54) * (pkin(7) + t54) + t34));
	t48 = pkin(7) * t19;
	t42 = t12 * t48;
	t27 = t57 * (t36 * t32 - t42);
	t49 = pkin(7) * t12;
	t43 = pkin(1) * t49;
	t29 = 0.2e1 / t19 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t34 - t20) * t43;
	t31 = pkin(7) * t32;
	t30 = t12 * t31;
	t61 = 0.1e1 / t11 ^ 2 * t43;
	t59 = (-t29 * t49 - 0.2e1 * t36 * t43 + t38 * t48 - t30) * t64 + t27 * t61;
	t28 = t57 * (t36 * t19 + t30);
	t60 = (-0.2e1 * t56 * t12 ^ 2 * pkin(1) + t36 * t29 - t38 * t31 - t42) * t64 + t28 * t61;
	t25 = t28 * t65;
	t26 = t55 * t27;
	t7 = t45 * t25 + t44 * t26 / 0.2e1;
	t52 = t60 * t44 - t59 * t45 - t7;
	t41 = t52 * r_i_i_C(1);
	t2 = t59 * t44 + t60 * t45;
	t6 = -t45 * t26 / 0.2e1 + t44 * t25;
	t51 = t2 + t6;
	t63 = -t51 * r_i_i_C(2) + t41;
	t58 = (t13 * t7 - t16 * t6) * r_i_i_C(2) - t13 * pkin(1) + pkin(16);
	t47 = t13 * t6;
	t46 = t16 * t7;
	t35 = -t51 * r_i_i_C(1) - t52 * r_i_i_C(2) - pkin(1);
	t33 = t35 * t16;
	t17 = cos(qJ(1));
	t14 = sin(qJ(1));
	t5 = t14 * t47;
	t1 = [t5 * r_i_i_C(1) + t17 * r_i_i_C(3) + (r_i_i_C(1) * t46 - t58) * t14, (-t63 * t13 + t33) * t17, 0, 0; t14 * r_i_i_C(3) + ((-t46 - t47) * r_i_i_C(1) + t58) * t17, t5 * r_i_i_C(2) + ((t2 * r_i_i_C(2) - t41) * t13 + t33) * t14, 0, 0; 0, t35 * t13 + t63 * t16, 0, 0;];
	Ja_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (1678->46), mult. (2632->88), div. (88->4), fcn. (1704->10), ass. (0->53)
	t28 = cos(qJ(2));
	t30 = pkin(6) ^ 2;
	t22 = sin(pkin(20));
	t23 = cos(pkin(20));
	t24 = sin(qJ(3));
	t27 = cos(qJ(3));
	t20 = t27 * t22 + t24 * t23;
	t53 = pkin(6) * t20;
	t57 = 2 * pkin(1);
	t44 = t53 * t57 + t30;
	t17 = (pkin(1) ^ 2) + t44;
	t14 = pkin(2) ^ 2 - pkin(13) ^ 2 + t17;
	t18 = -pkin(1) - t53;
	t55 = -pkin(2) - pkin(13);
	t12 = (pkin(1) - t55) * (pkin(1) + t55) + t44;
	t54 = -pkin(2) + pkin(13);
	t13 = (pkin(1) - t54) * (pkin(1) + t54) + t44;
	t33 = sqrt(-t13 * t12);
	t21 = t24 * t22 - t27 * t23;
	t52 = pkin(6) * t21;
	t8 = t14 * t52 - t18 * t33;
	t47 = t28 * t8;
	t25 = sin(qJ(2));
	t45 = t21 * t33;
	t7 = -pkin(6) * t45 - t18 * t14;
	t50 = t25 * t7;
	t36 = -t47 / 0.2e1 - t50 / 0.2e1;
	t15 = 0.1e1 / t17;
	t31 = 0.1e1 / pkin(2);
	t46 = t15 * t31;
	t58 = t36 * t46;
	t56 = -t25 / 0.2e1;
	t43 = pkin(1) * t52;
	t51 = 0.1e1 / t33 * (t12 + t13) * t43;
	t49 = t25 * t8;
	t48 = t28 * t7;
	t42 = 0.1e1 / t17 ^ 2 * t43;
	t41 = -t48 + t49;
	t40 = -t47 - t50;
	t39 = t49 / 0.2e1 - t48 / 0.2e1;
	t1 = (-t20 * t33 + (t18 * t57 - t14 - t51) * t21) * pkin(6);
	t2 = -t18 * t51 - 0.2e1 * t30 * t21 ^ 2 * pkin(1) + (t20 * t14 - t45) * pkin(6);
	t38 = t28 * t1 / 0.2e1 + t2 * t56;
	t37 = -t28 * t2 / 0.2e1 + t1 * t56;
	t35 = t39 * t46;
	t34 = ((t37 * r_i_i_C(1) - t38 * r_i_i_C(2)) * t15 + (t40 * r_i_i_C(1) + t41 * r_i_i_C(2)) * t42) * t31;
	t29 = cos(qJ(1));
	t26 = sin(qJ(1));
	t6 = t29 * t35;
	t5 = t29 * t58;
	t4 = t26 * t35;
	t3 = t26 * t58;
	t9 = [-t26 * pkin(16) - t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + t29 * r_i_i_C(3), t6 * r_i_i_C(1) - t5 * r_i_i_C(2), t29 * t34, 0; t29 * pkin(16) + t5 * r_i_i_C(1) + t6 * r_i_i_C(2) + t26 * r_i_i_C(3), t4 * r_i_i_C(1) - t3 * r_i_i_C(2), t26 * t34, 0; 0, (t36 * r_i_i_C(1) + t39 * r_i_i_C(2)) * t46, ((t38 * r_i_i_C(1) + t37 * r_i_i_C(2)) * t15 + (-t41 * r_i_i_C(1) + t40 * r_i_i_C(2)) * t42) * t31, 0;];
	Ja_transl = t9;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (5391->80), mult. (8346->145), div. (356->5), fcn. (5322->10), ass. (0->71)
	t41 = sin(pkin(20));
	t42 = cos(pkin(20));
	t43 = sin(qJ(3));
	t46 = cos(qJ(3));
	t39 = t46 * t41 + t43 * t42;
	t84 = pkin(6) * t39;
	t72 = pkin(1) * t84;
	t38 = 0.2e1 * t72;
	t52 = pkin(2) ^ 2;
	t51 = pkin(6) ^ 2;
	t73 = pkin(1) ^ 2 + t51;
	t68 = -pkin(13) ^ 2 + t73;
	t33 = t38 + t52 + t68;
	t37 = -pkin(1) - t84;
	t40 = t43 * t41 - t46 * t42;
	t74 = t38 + t51;
	t86 = -pkin(2) - pkin(13);
	t30 = (pkin(1) - t86) * (pkin(1) + t86) + t74;
	t85 = -pkin(2) + pkin(13);
	t31 = (pkin(1) - t85) * (pkin(1) + t85) + t74;
	t55 = sqrt(-t31 * t30);
	t80 = t40 * t55;
	t83 = pkin(6) * t40;
	t71 = pkin(1) * t83;
	t82 = 0.1e1 / t55 * (t30 + t31) * t71;
	t10 = -t37 * t82 - 0.2e1 * t51 * t40 ^ 2 * pkin(1) + (t39 * t33 - t80) * pkin(6);
	t36 = t38 + t73;
	t34 = 0.1e1 / t36;
	t47 = cos(qJ(2));
	t53 = 0.1e1 / pkin(2);
	t66 = 0.1e1 / t36 ^ 2 * t71;
	t25 = -pkin(6) * t80 - t37 * t33;
	t77 = t47 * t25;
	t26 = t33 * t83 - t37 * t55;
	t44 = sin(qJ(2));
	t78 = t44 * t26;
	t89 = -t44 / 0.2e1;
	t9 = (-t39 * t55 + (0.2e1 * t37 * pkin(1) - t33 - t82) * t40) * pkin(6);
	t93 = ((t47 * t9 / 0.2e1 + t10 * t89) * t34 + (t77 - t78) * t66) * t53;
	t76 = t47 * t26;
	t79 = t44 * t25;
	t81 = t34 * t53;
	t16 = (-t76 / 0.2e1 - t79 / 0.2e1) * t81;
	t48 = cos(qJ(1));
	t13 = t48 * t16;
	t91 = t13 / 0.2e1;
	t32 = t52 - t68 - 0.2e1 * t72;
	t90 = -t32 / 0.2e1;
	t88 = -t55 / 0.2e1;
	t87 = t55 / 0.2e1;
	t75 = 0.1e1 / pkin(13) * t53;
	t70 = r_i_i_C(1) * t75;
	t69 = r_i_i_C(2) * t75;
	t67 = -t82 / 0.2e1;
	t65 = t75 * t90;
	t64 = t75 * t87;
	t45 = sin(qJ(1));
	t11 = t45 * t16;
	t15 = (t78 / 0.2e1 - t77 / 0.2e1) * t81;
	t12 = t45 * t15;
	t62 = -t11 * t65 + t12 * t64;
	t14 = t48 * t15;
	t61 = t13 * t64 + t14 * t65;
	t59 = -t11 * t88 + t12 * t90;
	t58 = t14 * t87 + t32 * t91;
	t6 = ((-t47 * t10 / 0.2e1 + t9 * t89) * t34 + (-t76 - t79) * t66) * t53;
	t4 = t48 * t93;
	t3 = t48 * t6;
	t2 = t45 * t93;
	t1 = t45 * t6;
	t5 = [-t11 * pkin(2) - t45 * pkin(16) + t62 * r_i_i_C(1) + t48 * r_i_i_C(3) - t59 * t69, t14 * pkin(2) + t61 * r_i_i_C(1) + t58 * t69, t3 * pkin(2) + ((-t13 * t71 + t14 * t67 + t3 * t90 - t4 * t88) * r_i_i_C(1) + (-t14 * t71 + t3 * t87 - t4 * t90 + t82 * t91) * r_i_i_C(2)) * t75, 0; t13 * pkin(2) + t48 * pkin(16) + t61 * r_i_i_C(2) + t45 * r_i_i_C(3) - t58 * t70, t12 * pkin(2) + t62 * r_i_i_C(2) + t59 * t70, t1 * pkin(2) + ((t1 * t90 - t11 * t71 + t12 * t67 - t2 * t88) * r_i_i_C(1) + (t1 * t87 - t11 * t67 - t12 * t71 - t2 * t90) * r_i_i_C(2)) * t75, 0; 0, (t15 * t88 + t16 * t90) * t70 + (t15 * t90 + t16 * t87) * t69 + t16 * pkin(2), t93 * pkin(2) + ((t15 * t71 + t16 * t67 + t6 * t88 + t90 * t93) * r_i_i_C(1) + (t15 * t67 - t16 * t71 + t6 * t90 + t87 * t93) * r_i_i_C(2)) * t75, 0;];
	Ja_transl = t5;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:18:16
	% DurationCPUTime: 12.02s
	% Computational Cost: add. (201249->154), mult. (305597->271), div. (12512->8), fcn. (194527->18), ass. (0->131)
	t124 = sin(qJ(2));
	t125 = sin(pkin(19));
	t127 = cos(qJ(2));
	t128 = cos(pkin(19));
	t102 = t124 * t128 - t127 * t125;
	t101 = pkin(7) * t102;
	t46 = t124 * t125 + t127 * t128;
	t131 = pkin(7) * t46;
	t118 = pkin(1) * t131;
	t136 = -pkin(8) + pkin(3);
	t137 = -pkin(8) - pkin(3);
	t98 = (-0.2e1 * t101 + pkin(1)) * pkin(1);
	t53 = sqrt(-((pkin(7) - t136) * (pkin(7) + t136) + t98) * ((pkin(7) - t137) * (pkin(7) + t137) + t98));
	t54 = pkin(3) ^ 2;
	t122 = 0.4e1 / t53 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t98 - t54) * t118;
	t144 = pkin(7) ^ 2;
	t149 = -pkin(8) ^ 2 + t54;
	t36 = (t102 * t53 + (-t122 / 0.2e1 - t144 - t149) * t46) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t101) * t118;
	t162 = -t36 / 0.2e1;
	t45 = t98 + t144;
	t97 = t45 + t149;
	t96 = pkin(7) * t97;
	t99 = -t101 + pkin(1);
	t37 = t99 * t122 / 0.2e1 - t102 * t96 + (-t53 - 0.2e1 * t118) * t131;
	t161 = t37 / 0.2e1;
	t48 = sin(pkin(22));
	t160 = t48 / 0.2e1;
	t49 = cos(pkin(22));
	t159 = -t49 / 0.2e1;
	t116 = pkin(23) + pkin(22);
	t110 = sin(t116);
	t111 = cos(t116);
	t123 = sin(qJ(3));
	t126 = cos(qJ(3));
	t139 = 0.1e1 / t45;
	t146 = 0.1e1 / pkin(3);
	t95 = t146 * (t46 * t96 + t99 * t53);
	t92 = t139 * t95 / 0.2e1;
	t94 = t146 * (-t53 * t131 + t99 * t97);
	t93 = t139 * t94;
	t84 = t123 * t93 / 0.2e1 + t126 * t92;
	t156 = -t126 / 0.2e1;
	t85 = t123 * t92 + t156 * t93;
	t77 = -t110 * t85 + t111 * t84;
	t75 = pkin(5) * t77;
	t73 = -t75 - pkin(4);
	t158 = -t73 / 0.2e1;
	t148 = t110 * t84 + t111 * t85;
	t76 = pkin(5) * t148;
	t157 = -t76 / 0.2e1;
	t138 = pkin(4) * pkin(5);
	t145 = pkin(5) ^ 2;
	t26 = t145 + (0.2e1 * t75 + pkin(4)) * pkin(4);
	t55 = pkin(9) ^ 2;
	t70 = -pkin(11) ^ 2 + t26 + t55;
	t69 = pkin(5) * t70;
	t155 = 0.2e1 * t73 * t138 - t69;
	t44 = 0.1e1 / t45 ^ 2;
	t154 = t118 * t44;
	t153 = t138 / t26 ^ 2;
	t142 = 0.1e1 / pkin(9);
	t134 = -pkin(9) + pkin(11);
	t135 = -pkin(9) - pkin(11);
	t141 = -0.2e1 * pkin(4);
	t72 = (-t141 * t77 + pkin(5)) * pkin(5);
	t52 = sqrt(-((pkin(4) - t134) * (pkin(4) + t134) + t72) * ((pkin(4) - t135) * (pkin(4) + t135) + t72));
	t15 = -t52 * t76 - t73 * t70;
	t129 = t15 * t142;
	t109 = t129 * t153;
	t64 = t142 * (t148 * t69 - t73 * t52);
	t62 = t64 * t153;
	t152 = -t49 * t109 - t48 * t62;
	t151 = t48 * t109 - t49 * t62;
	t132 = pkin(5) * t52;
	t150 = t141 * t145 * t148 - t132;
	t147 = 0.4e1 / t52 * (-t55 + (pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t72) * t138;
	t140 = 0.1e1 / t26;
	t133 = t146 * t139;
	t130 = t142 * t140;
	t121 = cos(pkin(23));
	t120 = sin(pkin(23));
	t115 = -t130 / 0.2e1;
	t114 = t127 * pkin(1);
	t113 = t124 * pkin(1);
	t112 = t123 * t133;
	t108 = t121 * t133;
	t107 = t113 - pkin(16);
	t104 = t133 * t156;
	t103 = t120 * t133 / 0.2e1;
	t87 = t121 * t94;
	t88 = t120 * t95;
	t38 = (-t87 / 0.2e1 + t88 / 0.2e1) * t139;
	t86 = t120 * t94;
	t89 = t121 * t95;
	t39 = (t89 / 0.2e1 + t86 / 0.2e1) * t139;
	t32 = -t124 * t39 + t127 * t38;
	t100 = t124 * t38 + t127 * t39;
	t23 = t108 * t162 + t37 * t103 + (-t87 + t88) * t154;
	t24 = t108 * t161 + t36 * t103 + (t86 + t89) * t154;
	t20 = -t124 * t24 + t127 * t23 - t100;
	t21 = -t124 * t23 - t127 * t24 - t32;
	t91 = t95 * t118;
	t90 = t94 * t118;
	t81 = t36 * t104 + t112 * t161 + (t123 * t91 - t126 * t90) * t44;
	t80 = t112 * t162 + t37 * t104 + (-t123 * t90 - t126 * t91) * t44;
	t22 = t110 * t81 + t111 * t80;
	t67 = t22 * t147;
	t66 = t148 * t147;
	t65 = pkin(5) * (-t110 * t80 + t111 * t81);
	t63 = -t140 * t64 / 0.2e1;
	t59 = (-t77 * t132 + t148 * t155 + t157 * t66) * t130;
	t58 = (t148 * t150 + t158 * t66 + t77 * t69) * t115;
	t57 = (t155 * t22 + t157 * t67 - t52 * t65) * t130;
	t56 = (t150 * t22 + t158 * t67 + t70 * t65) * t115;
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t30 = t32 * t51;
	t29 = t100 * t51;
	t28 = t32 * t50;
	t27 = t100 * t50;
	t19 = t20 * t51;
	t18 = t21 * t51;
	t17 = t20 * t50;
	t16 = t21 * t50;
	t10 = t49 * t15 * t115 + t48 * t63;
	t9 = t129 * t140 * t160 + t49 * t63;
	t4 = t148 * t152 + t159 * t59 + t48 * t58;
	t3 = t148 * t151 + t160 * t59 + t49 * t58;
	t2 = t152 * t22 + t159 * t57 + t48 * t56;
	t1 = t151 * t22 + t160 * t57 + t49 * t56;
	t5 = [(t10 * t27 + t28 * t9) * r_i_i_C(1) + (t10 * t28 - t27 * t9) * r_i_i_C(2) + t51 * r_i_i_C(3) + t107 * t50 + (t27 * t49 - t28 * t48) * pkin(4), (-t1 * t30 + t10 * t18 - t19 * t9 - t2 * t29) * r_i_i_C(1) + (t1 * t29 - t10 * t19 - t18 * t9 - t2 * t30) * r_i_i_C(2) - t51 * t114 + (t18 * t49 + t19 * t48) * pkin(4), (-t29 * t4 - t3 * t30) * r_i_i_C(1) + (t29 * t3 - t30 * t4) * r_i_i_C(2), 0; (-t10 * t29 - t30 * t9) * r_i_i_C(1) + (-t10 * t30 + t29 * t9) * r_i_i_C(2) + t50 * r_i_i_C(3) - t107 * t51 + (-t29 * t49 + t30 * t48) * pkin(4), (-t1 * t28 + t10 * t16 - t17 * t9 - t2 * t27) * r_i_i_C(1) + (t1 * t27 - t10 * t17 - t16 * t9 - t2 * t28) * r_i_i_C(2) - t50 * t114 + (t16 * t49 + t17 * t48) * pkin(4), (-t27 * t4 - t28 * t3) * r_i_i_C(1) + (t27 * t3 - t28 * t4) * r_i_i_C(2), 0; 0, (-t1 * t100 + t10 * t20 + t2 * t32 + t21 * t9) * r_i_i_C(1) + (-t1 * t32 + t10 * t21 - t100 * t2 - t20 * t9) * r_i_i_C(2) - t113 + (t20 * t49 - t21 * t48) * pkin(4), (-t100 * t3 + t32 * t4) * r_i_i_C(1) + (-t100 * t4 - t3 * t32) * r_i_i_C(2), 0;];
	Ja_transl = t5;
end