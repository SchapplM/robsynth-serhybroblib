% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m1TE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m1TE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobia_transl_sym_varpar: pkin has to be [19x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
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
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
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
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->15), mult. (85->26), div. (0->0), fcn. (101->6), ass. (0->19)
	t13 = sin(qJ(1));
	t11 = sin(qJ(3));
	t12 = sin(qJ(2));
	t14 = cos(qJ(3));
	t15 = cos(qJ(2));
	t18 = t15 * t11 + t12 * t14;
	t5 = t18 * t13;
	t17 = t12 * t11 - t15 * t14;
	t6 = t17 * t13;
	t24 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t16 = cos(qJ(1));
	t7 = t17 * t16;
	t8 = t18 * t16;
	t23 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t22 = pkin(1) * t12;
	t21 = t15 * pkin(1);
	t20 = t17 * r_i_i_C(1) + t18 * r_i_i_C(2);
	t19 = pkin(13) + t21;
	t1 = [-t6 * r_i_i_C(1) - t5 * r_i_i_C(2) + t16 * r_i_i_C(3) - t19 * t13, -t16 * t22 + t23, t23, 0; t7 * r_i_i_C(1) + t8 * r_i_i_C(2) + t13 * r_i_i_C(3) + t19 * t16, -t13 * t22 + t24, t24, 0; 0, t20 + t21, t20, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:20
	% DurationCPUTime: 7.89s
	% Computational Cost: add. (205385->138), mult. (311257->227), div. (12992->8), fcn. (197929->16), ass. (0->119)
	t141 = sin(qJ(2));
	t142 = cos(pkin(16));
	t77 = sin(pkin(16));
	t79 = cos(qJ(2));
	t62 = t141 * t77 - t142 * t79;
	t146 = pkin(5) * t62;
	t157 = -2 * pkin(1);
	t125 = (pkin(1) ^ 2) + t146 * t157;
	t82 = pkin(5) ^ 2;
	t53 = t82 + t125;
	t49 = 0.1e1 / t53;
	t85 = 0.1e1 / pkin(2);
	t135 = t49 * t85;
	t64 = t141 * t142 + t79 * t77;
	t145 = pkin(5) * t64;
	t48 = pkin(2) ^ 2 - pkin(6) ^ 2 + t53;
	t58 = pkin(1) - t146;
	t153 = -pkin(2) - pkin(6);
	t46 = (pkin(5) - t153) * (pkin(5) + t153) + t125;
	t152 = -pkin(6) + pkin(2);
	t47 = (pkin(5) - t152) * (pkin(5) + t152) + t125;
	t87 = sqrt(-t47 * t46);
	t42 = t145 * t48 + t58 * t87;
	t75 = sin(qJ(3));
	t138 = t42 * t75;
	t134 = t64 * t87;
	t41 = -pkin(5) * t134 + t48 * t58;
	t78 = cos(qJ(3));
	t139 = t41 * t78;
	t39 = (-t139 / 0.2e1 + t138 / 0.2e1) * t135;
	t137 = t42 * t78;
	t140 = t41 * t75;
	t40 = (t137 / 0.2e1 + t140 / 0.2e1) * t135;
	t72 = pkin(18) + pkin(19);
	t70 = sin(t72);
	t71 = cos(t72);
	t108 = t39 * t71 + t40 * t70;
	t98 = pkin(3) * t108;
	t32 = pkin(4) + t98;
	t167 = t32 / 0.2e1;
	t73 = sin(pkin(17));
	t166 = t73 / 0.2e1;
	t74 = cos(pkin(17));
	t165 = -t74 / 0.2e1;
	t107 = t39 * t70 - t40 * t71;
	t148 = pkin(3) * t107;
	t164 = -t148 / 0.2e1;
	t155 = pkin(4) * pkin(3);
	t123 = -0.2e1 * t155;
	t126 = pkin(4) ^ 2 - t108 * t123;
	t84 = pkin(3) ^ 2;
	t31 = t84 + t126;
	t156 = 0.1e1 / t31;
	t124 = t156 / 0.2e1;
	t118 = t73 * t124;
	t150 = -pkin(8) + pkin(10);
	t151 = -pkin(8) - pkin(10);
	t88 = sqrt(-((pkin(3) - t150) * (pkin(3) + t150) + t126) * ((pkin(3) - t151) * (pkin(3) + t151) + t126));
	t147 = pkin(3) * t88;
	t158 = pkin(8) ^ 2;
	t29 = pkin(10) ^ 2 - t158 + t31;
	t25 = -t107 * t147 + t29 * t32;
	t26 = t148 * t29 + t32 * t88;
	t81 = 0.1e1 / pkin(10);
	t21 = (t25 * t156 * t165 + t26 * t118) * t81;
	t22 = (t124 * t26 * t74 + t118 * t25) * t81;
	t76 = sin(qJ(1));
	t97 = t141 * t78 + t79 * t75;
	t54 = t97 * t76;
	t61 = t141 * t75 - t79 * t78;
	t55 = t61 * t76;
	t131 = -t55 * t21 - t54 * t22;
	t80 = cos(qJ(1));
	t56 = t61 * t80;
	t57 = t97 * t80;
	t130 = -t56 * t21 - t57 * t22;
	t129 = t57 * t21 - t56 * t22;
	t111 = 0.1e1 / t31 ^ 2 * t81 * t155;
	t109 = t26 * t111;
	t110 = t25 * t111;
	t163 = t74 * t109 + t73 * t110;
	t162 = t73 * t109 - t74 * t110;
	t161 = -0.2e1 * t84 * t107 * pkin(4) - t147;
	t160 = -pkin(3) * t29 + t32 * t123;
	t159 = 0.4e1 / t88 * (-t158 + (pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t126) * t155;
	t154 = t75 / 0.2e1;
	t144 = t79 * pkin(1);
	t143 = t54 * t21 - t55 * t22;
	t121 = pkin(1) * t145;
	t136 = 0.2e1 / t87 * (t46 + t47) * t121;
	t132 = t81 * t156;
	t128 = t61 * t21 + t22 * t97;
	t127 = t21 * t97 - t61 * t22;
	t119 = t141 * pkin(1);
	t117 = t81 * t124;
	t113 = pkin(13) + t144;
	t112 = 0.1e1 / t53 ^ 2 * t121;
	t106 = t138 - t139;
	t105 = t137 + t140;
	t37 = (t62 * t87 + (-t136 / 0.2e1 - t48 + t58 * t157) * t64) * pkin(5);
	t38 = t58 * t136 / 0.2e1 + t82 * t64 ^ 2 * t157 + (-t48 * t62 - t134) * pkin(5);
	t100 = -t37 * t78 / 0.2e1 + t38 * t154;
	t99 = t38 * t78 / 0.2e1 + t37 * t154;
	t27 = ((-t100 * t71 - t70 * t99) * t49 + (-t105 * t70 - t106 * t71) * t112) * t85;
	t95 = t27 * t159;
	t94 = t107 * t159;
	t93 = (t160 * t107 + t94 * t164 - t88 * t98) * t132;
	t92 = (t161 * t107 + t94 * t167 + t29 * t98) * t117;
	t91 = pkin(3) * t85 * ((t100 * t70 - t71 * t99) * t49 + (-t105 * t71 + t106 * t70) * t112);
	t90 = (t160 * t27 + t95 * t164 - t88 * t91) * t132;
	t89 = (t161 * t27 + t95 * t167 + t29 * t91) * t117;
	t60 = t61 * pkin(4);
	t52 = t57 * pkin(4);
	t51 = t54 * pkin(4);
	t4 = t163 * t107 + t93 * t166 + t74 * t92;
	t3 = t162 * t107 + t93 * t165 + t73 * t92;
	t2 = t163 * t27 + t90 * t166 + t74 * t89;
	t1 = t162 * t27 + t90 * t165 + t73 * t89;
	t5 = [-t55 * pkin(4) + t131 * r_i_i_C(1) - t143 * r_i_i_C(2) + t80 * r_i_i_C(3) - t113 * t76, (t1 * t56 + t2 * t57 + t129) * r_i_i_C(1) + (t1 * t57 - t2 * t56 + t130) * r_i_i_C(2) + t52 - t80 * t119, (t3 * t56 + t4 * t57 + t129) * r_i_i_C(1) + (t3 * t57 - t4 * t56 + t130) * r_i_i_C(2) + t52, 0; t56 * pkin(4) - t130 * r_i_i_C(1) + t129 * r_i_i_C(2) + t76 * r_i_i_C(3) + t113 * t80, (t1 * t55 + t2 * t54 + t143) * r_i_i_C(1) + (t1 * t54 - t2 * t55 + t131) * r_i_i_C(2) + t51 - t76 * t119, (t3 * t55 + t4 * t54 + t143) * r_i_i_C(1) + (t3 * t54 - t4 * t55 + t131) * r_i_i_C(2) + t51, 0; 0, (-t1 * t97 + t2 * t61 + t128) * r_i_i_C(1) + (t1 * t61 + t2 * t97 + t127) * r_i_i_C(2) + t60 + t144, (-t3 * t97 + t4 * t61 + t128) * r_i_i_C(1) + (t3 * t61 + t4 * t97 + t127) * r_i_i_C(2) + t60, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:49
	% DurationCPUTime: 25.58s
	% Computational Cost: add. (530252->152), mult. (803355->245), div. (33680->8), fcn. (510749->18), ass. (0->129)
	t71 = pkin(2) ^ 2;
	t172 = -pkin(6) ^ 2 + t71;
	t147 = sin(qJ(2));
	t148 = sin(pkin(16));
	t150 = cos(qJ(2));
	t151 = cos(pkin(16));
	t119 = t147 * t148 - t150 * t151;
	t118 = pkin(5) * t119;
	t116 = (-0.2e1 * t118 + pkin(1)) * pkin(1);
	t168 = pkin(5) ^ 2;
	t50 = t116 + t168;
	t115 = t50 + t172;
	t114 = pkin(5) * t115;
	t117 = -t118 + pkin(1);
	t169 = 0.1e1 / pkin(2);
	t58 = t147 * t151 + t150 * t148;
	t159 = -pkin(6) + pkin(2);
	t160 = -pkin(6) - pkin(2);
	t70 = sqrt(-((pkin(5) - t159) * (pkin(5) + t159) + t116) * ((pkin(5) - t160) * (pkin(5) + t160) + t116));
	t113 = t169 * (t58 * t114 + t117 * t70);
	t162 = 0.1e1 / t50;
	t110 = t162 * t113 / 0.2e1;
	t153 = pkin(5) * t58;
	t112 = t169 * (t117 * t115 - t70 * t153);
	t111 = t162 * t112;
	t146 = sin(qJ(3));
	t149 = cos(qJ(3));
	t106 = -t149 * t111 / 0.2e1 + t146 * t110;
	t178 = t146 / 0.2e1;
	t107 = t149 * t110 + t111 * t178;
	t132 = pkin(18) + pkin(19);
	t126 = sin(t132);
	t127 = cos(t132);
	t100 = t127 * t106 + t126 * t107;
	t164 = -2 * pkin(4);
	t41 = (pkin(4) ^ 2) + (-t100 * t164 + pkin(3)) * pkin(3);
	t163 = 0.1e1 / t41;
	t136 = cos(pkin(17));
	t165 = 0.1e1 / pkin(10);
	t157 = -pkin(8) + pkin(10);
	t158 = -pkin(8) - pkin(10);
	t98 = pkin(3) * t100;
	t94 = (0.2e1 * t98 + pkin(4)) * pkin(4);
	t69 = sqrt(-((pkin(3) - t157) * (pkin(3) + t157) + t94) * ((pkin(3) - t158) * (pkin(3) + t158) + t94));
	t166 = pkin(8) ^ 2;
	t92 = pkin(10) ^ 2 - t166 + t41;
	t95 = t98 + pkin(4);
	t170 = t126 * t106 - t127 * t107;
	t97 = pkin(3) * t170;
	t85 = t165 * (-t69 * t97 + t95 * t92);
	t82 = t136 * t85;
	t135 = sin(pkin(17));
	t91 = pkin(3) * t92;
	t86 = t165 * (t170 * t91 + t95 * t69);
	t83 = t135 * t86;
	t35 = (-t82 / 0.2e1 + t83 / 0.2e1) * t163;
	t81 = t135 * t85;
	t84 = t136 * t86;
	t36 = (t84 / 0.2e1 + t81 / 0.2e1) * t163;
	t57 = -t150 * t146 - t147 * t149;
	t66 = sin(qJ(1));
	t51 = t57 * t66;
	t56 = t147 * t146 - t150 * t149;
	t52 = t56 * t66;
	t21 = -t35 * t52 + t36 * t51;
	t65 = sin(qJ(4));
	t67 = cos(qJ(4));
	t68 = cos(qJ(1));
	t184 = t21 * t67 + t68 * t65;
	t183 = t21 * t65 - t67 * t68;
	t182 = t95 / 0.2e1;
	t181 = -t97 / 0.2e1;
	t180 = t135 / 0.2e1;
	t179 = -t136 / 0.2e1;
	t139 = -t51 * t35 - t52 * t36;
	t53 = t56 * t68;
	t54 = t57 * t68;
	t138 = -t54 * t35 - t53 * t36;
	t161 = pkin(3) * pkin(4);
	t174 = t161 / t41 ^ 2;
	t177 = (-t82 + t83) * t174;
	t176 = (t81 + t84) * t174;
	t175 = -0.2e1 * t95 * t161 - t91;
	t154 = pkin(3) * t69;
	t173 = pkin(3) ^ 2 * t164 * t170 - t154;
	t171 = 0.4e1 / t69 * (-t166 + (pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t94) * t161;
	t156 = -pkin(11) - r_i_i_C(3);
	t155 = t169 * t162;
	t152 = t165 * t163;
	t133 = pkin(1) * t153;
	t145 = 0.4e1 / t70 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t116 - t71) * t133;
	t108 = t112 * t133;
	t109 = t113 * t133;
	t121 = t155 * t178;
	t125 = t149 * t155;
	t43 = (t119 * t70 + (-t145 / 0.2e1 - t168 - t172) * t58) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t118) * t133;
	t44 = t117 * t145 / 0.2e1 - t119 * t114 + (-t70 - 0.2e1 * t133) * t153;
	t47 = 0.1e1 / t50 ^ 2;
	t102 = t44 * t125 / 0.2e1 + t43 * t121 + (t146 * t108 + t149 * t109) * t47;
	t103 = -t43 * t125 / 0.2e1 + t44 * t121 + (-t149 * t108 + t146 * t109) * t47;
	t39 = -t126 * t102 - t127 * t103;
	t130 = t152 / 0.2e1;
	t87 = pkin(3) * (-t127 * t102 + t126 * t103);
	t89 = t39 * t171;
	t73 = (t173 * t39 + t89 * t182 + t92 * t87) * t130;
	t74 = (t175 * t39 + t89 * t181 - t69 * t87) * t152;
	t13 = t135 * t73 + t177 * t39 + t74 * t179;
	t143 = t13 - t36;
	t14 = t136 * t73 + t176 * t39 + t74 * t180;
	t142 = t14 + t35;
	t88 = t170 * t171;
	t75 = (t100 * t91 + t170 * t173 + t88 * t182) * t130;
	t76 = (-t100 * t154 + t170 * t175 + t88 * t181) * t152;
	t15 = t135 * t75 + t170 * t177 + t76 * t179;
	t141 = t15 - t36;
	t16 = t136 * t75 + t170 * t176 + t76 * t180;
	t140 = t16 + t35;
	t137 = t56 * t35 - t57 * t36;
	t129 = t150 * pkin(1);
	t128 = t147 * pkin(1);
	t122 = t129 + pkin(13);
	t120 = r_i_i_C(1) * t67 - r_i_i_C(2) * t65 + pkin(9);
	t55 = t56 * pkin(4);
	t49 = t54 * pkin(4);
	t48 = t51 * pkin(4);
	t24 = t35 * t53 - t36 * t54;
	t20 = t24 * t67 + t65 * t66;
	t19 = -t24 * t65 + t66 * t67;
	t1 = [-t52 * pkin(4) + t21 * pkin(9) + t184 * r_i_i_C(1) - t183 * r_i_i_C(2) - t122 * t66 - t156 * t139, -t68 * t128 - t49 + t156 * (-t142 * t53 - t143 * t54) + t120 * (t13 * t53 - t14 * t54 + t138), -t49 + t156 * (-t140 * t53 - t141 * t54) + t120 * (t15 * t53 - t16 * t54 + t138), r_i_i_C(1) * t19 - r_i_i_C(2) * t20; t53 * pkin(4) + t24 * pkin(9) + t20 * r_i_i_C(1) + t19 * r_i_i_C(2) + t122 * t68 + t156 * t138, -t66 * t128 - t48 + t156 * (-t142 * t52 - t143 * t51) + t120 * (t13 * t52 - t14 * t51 + t139), -t48 + t156 * (-t140 * t52 - t141 * t51) + t120 * (t15 * t52 - t16 * t51 + t139), t183 * r_i_i_C(1) + t184 * r_i_i_C(2); 0, t129 + t55 + t156 * (-t142 * t57 + t143 * t56) + t120 * (t13 * t57 + t14 * t56 + t137), t55 + t156 * (-t140 * t57 + t141 * t56) + t120 * (t15 * t57 + t16 * t56 + t137), (-r_i_i_C(1) * t65 - r_i_i_C(2) * t67) * (t35 * t57 + t36 * t56);];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:03
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (1345->39), mult. (2098->70), div. (64->4), fcn. (1356->10), ass. (0->42)
	t54 = -2 * pkin(5);
	t27 = sin(pkin(15));
	t53 = t27 / 0.2e1;
	t29 = cos(pkin(15));
	t52 = t29 / 0.2e1;
	t51 = (-pkin(6) - pkin(2));
	t50 = (-pkin(6) + pkin(2));
	t24 = sin(qJ(2));
	t26 = sin(pkin(16));
	t46 = cos(qJ(2));
	t47 = cos(pkin(16));
	t21 = t24 * t26 - t46 * t47;
	t49 = pkin(1) * t21;
	t22 = t24 * t47 + t46 * t26;
	t48 = t22 * pkin(1);
	t32 = pkin(1) ^ 2;
	t42 = t49 * t54 + t32;
	t13 = ((pkin(5) - t51) * (pkin(5) + t51)) + t42;
	t14 = ((pkin(5) - t50) * (pkin(5) + t50)) + t42;
	t33 = sqrt(-t14 * t13);
	t41 = pkin(5) * t48;
	t45 = 0.1e1 / t33 * (t13 + t14) * t41;
	t18 = (pkin(5) ^ 2) + t42;
	t30 = 0.1e1 / pkin(6);
	t44 = 0.1e1 / t18 * t30;
	t43 = t22 * t33;
	t40 = 0.1e1 / t18 ^ 2 * t30 * t41;
	t15 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t18;
	t19 = -pkin(5) + t49;
	t34 = (-t19 * t45 + t32 * t22 ^ 2 * t54 + (-t21 * t15 - t43) * pkin(1)) * t44 / 0.2e1;
	t35 = pkin(1) * (t21 * t33 + (0.2e1 * t19 * pkin(5) - t15 - t45) * t22) * t44;
	t37 = t29 * t40;
	t38 = t27 * t40;
	t8 = -pkin(1) * t43 - t19 * t15;
	t9 = t15 * t48 - t19 * t33;
	t1 = t29 * t34 + t9 * t37 - t27 * t35 / 0.2e1 - t8 * t38;
	t2 = t27 * t34 + t35 * t52 + t8 * t37 + t9 * t38;
	t39 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t36 = -pkin(7) + (r_i_i_C(1) * (t8 * t52 + t9 * t53) + r_i_i_C(2) * (-t9 * t29 / 0.2e1 + t8 * t53)) * t44;
	t28 = cos(qJ(1));
	t25 = sin(qJ(1));
	t3 = [t28 * r_i_i_C(3) - t36 * t25, t39 * t28, 0, 0; t25 * r_i_i_C(3) + t36 * t28, t39 * t25, 0, 0; 0, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:04
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (3368->48), mult. (5239->77), div. (176->4), fcn. (3405->10), ass. (0->47)
	t12 = sin(qJ(2));
	t14 = sin(pkin(16));
	t15 = cos(qJ(2));
	t17 = cos(pkin(16));
	t37 = t12 * t14 - t15 * t17;
	t36 = pkin(5) * t37;
	t33 = (-0.2e1 * t36 + pkin(1)) * pkin(1);
	t53 = pkin(5) ^ 2;
	t10 = t33 + t53;
	t52 = 0.1e1 / t10;
	t61 = t52 / 0.2e1;
	t54 = 0.1e1 / pkin(2);
	t60 = t54 * t61;
	t11 = t12 * t17 + t15 * t14;
	t46 = pkin(5) * t11;
	t42 = pkin(1) * t46;
	t59 = t42 / t10 ^ 2;
	t50 = -pkin(6) + pkin(2);
	t51 = -pkin(6) - pkin(2);
	t18 = sqrt(-((pkin(5) - t50) * (pkin(5) + t50) + t33) * ((pkin(5) - t51) * (pkin(5) + t51) + t33));
	t19 = pkin(2) ^ 2;
	t32 = -pkin(6) ^ 2 + t10 + t19;
	t31 = pkin(5) * t32;
	t30 = t11 * t31;
	t34 = -t36 + pkin(1);
	t27 = t54 * (t34 * t18 + t30);
	t28 = 0.2e1 / t18 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t33 - t19) * t42;
	t45 = pkin(5) * t18;
	t41 = t11 * t45;
	t58 = (-0.2e1 * t53 * t11 ^ 2 * pkin(1) + t34 * t28 - t37 * t31 - t41) * t60 + t27 * t59;
	t26 = t54 * (t34 * t32 - t41);
	t57 = (-t28 * t46 - 0.2e1 * t34 * t42 + t37 * t45 - t30) * t60 + t26 * t59;
	t43 = sin(pkin(19));
	t44 = cos(pkin(19));
	t24 = t27 * t61;
	t25 = t52 * t26;
	t6 = t44 * t24 + t43 * t25 / 0.2e1;
	t48 = t58 * t43 - t57 * t44 - t6;
	t5 = -t44 * t25 / 0.2e1 + t43 * t24;
	t49 = t57 * t43 + t58 * t44 + t5;
	t56 = t49 * r_i_i_C(1) + t48 * r_i_i_C(2) + pkin(1);
	t55 = (t12 * t6 - t15 * t5) * r_i_i_C(1) + (t12 * t5 + t15 * t6) * r_i_i_C(2) - t15 * pkin(1) - pkin(13);
	t35 = t48 * r_i_i_C(1) - t49 * r_i_i_C(2);
	t29 = -t56 * t12 + t35 * t15;
	t16 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t16 * r_i_i_C(3) + t55 * t13, t29 * t16, 0, 0; t13 * r_i_i_C(3) - t55 * t16, t29 * t13, 0, 0; 0, t35 * t12 + t56 * t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:26
	% DurationCPUTime: 12.01s
	% Computational Cost: add. (201249->157), mult. (305597->270), div. (12512->8), fcn. (194527->18), ass. (0->132)
	t131 = sin(qJ(2));
	t134 = cos(qJ(2));
	t127 = sin(pkin(19));
	t128 = cos(pkin(19));
	t132 = sin(pkin(16));
	t135 = cos(pkin(16));
	t105 = t131 * t132 - t134 * t135;
	t104 = pkin(5) * t105;
	t102 = (-0.2e1 * t104 + pkin(1)) * pkin(1);
	t151 = pkin(5) ^ 2;
	t48 = t102 + t151;
	t146 = 0.1e1 / t48;
	t103 = -t104 + pkin(1);
	t49 = t131 * t135 + t132 * t134;
	t138 = pkin(5) * t49;
	t153 = 0.1e1 / pkin(2);
	t143 = -pkin(6) + pkin(2);
	t144 = -pkin(6) - pkin(2);
	t56 = sqrt(-((pkin(5) - t143) * (pkin(5) + t143) + t102) * ((pkin(5) - t144) * (pkin(5) + t144) + t102));
	t57 = pkin(2) ^ 2;
	t156 = -pkin(6) ^ 2 + t57;
	t99 = t48 + t156;
	t95 = t153 * (t103 * t99 - t138 * t56);
	t94 = t146 * t95;
	t92 = -t94 / 0.2e1;
	t97 = pkin(5) * t99;
	t96 = t153 * (t103 * t56 + t49 * t97);
	t93 = t146 * t96 / 0.2e1;
	t41 = t127 * t93 + t128 * t92;
	t91 = t94 / 0.2e1;
	t42 = t127 * t91 + t128 * t93;
	t33 = -t131 * t42 + t134 * t41;
	t125 = pkin(1) * t138;
	t129 = 0.4e1 / t56 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t102 - t57) * t125;
	t36 = (t105 * t56 + (-t129 / 0.2e1 - t151 - t156) * t49) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t104) * t125;
	t167 = -t36 / 0.2e1;
	t37 = t103 * t129 / 0.2e1 - t105 * t97 + (-t56 - 0.2e1 * t125) * t138;
	t166 = t37 / 0.2e1;
	t51 = sin(pkin(18));
	t165 = t51 / 0.2e1;
	t52 = cos(pkin(18));
	t164 = -t52 / 0.2e1;
	t124 = pkin(18) + pkin(19);
	t114 = sin(t124);
	t115 = cos(t124);
	t130 = sin(qJ(3));
	t133 = cos(qJ(3));
	t87 = t130 * t93 + t133 * t92;
	t88 = t130 * t91 + t133 * t93;
	t81 = t114 * t88 + t115 * t87;
	t79 = pkin(4) * t81;
	t76 = -t79 - pkin(3);
	t163 = -t76 / 0.2e1;
	t155 = t114 * t87 - t115 * t88;
	t78 = pkin(4) * t155;
	t162 = -t78 / 0.2e1;
	t145 = pkin(3) * pkin(4);
	t152 = pkin(4) ^ 2;
	t26 = t152 + (0.2e1 * t79 + pkin(3)) * pkin(3);
	t58 = pkin(8) ^ 2;
	t73 = -pkin(10) ^ 2 + t26 + t58;
	t72 = pkin(4) * t73;
	t161 = 0.2e1 * t145 * t76 - t72;
	t160 = t145 / t26 ^ 2;
	t149 = 0.1e1 / pkin(8);
	t141 = -pkin(8) + pkin(10);
	t142 = -pkin(8) - pkin(10);
	t148 = -0.2e1 * pkin(3);
	t75 = (-t148 * t81 + pkin(4)) * pkin(4);
	t55 = sqrt(-((pkin(3) - t141) * (pkin(3) + t141) + t75) * ((pkin(3) - t142) * (pkin(3) + t142) + t75));
	t15 = -t55 * t78 - t73 * t76;
	t136 = t15 * t149;
	t111 = t136 * t160;
	t67 = t149 * (t155 * t72 - t55 * t76);
	t65 = t67 * t160;
	t159 = -t111 * t52 - t51 * t65;
	t158 = t111 * t51 - t52 * t65;
	t139 = pkin(4) * t55;
	t157 = t148 * t152 * t155 - t139;
	t140 = t153 * t146;
	t122 = t140 / 0.2e1;
	t106 = t127 * t122;
	t112 = t128 * t140;
	t47 = 0.1e1 / t48 ^ 2;
	t89 = t95 * t125;
	t90 = t96 * t125;
	t23 = t112 * t166 + t36 * t106 + (t127 * t89 + t128 * t90) * t47;
	t24 = t112 * t167 + t37 * t106 + (t127 * t90 - t128 * t89) * t47;
	t20 = t131 * t24 + t134 * t23 + t33;
	t154 = 0.4e1 / t55 * (-t58 + (pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t75) * t145;
	t147 = 0.1e1 / t26;
	t137 = t149 * t147;
	t121 = -t137 / 0.2e1;
	t120 = t134 * pkin(1);
	t119 = t131 * pkin(1);
	t118 = t134 * t42;
	t117 = t131 * t41;
	t113 = t133 * t140;
	t110 = t120 + pkin(13);
	t107 = t130 * t122;
	t32 = t117 + t118;
	t100 = -t131 * t23 + t134 * t24 - t118;
	t21 = -t117 + t100;
	t84 = t113 * t167 + t37 * t107 + (t130 * t90 - t133 * t89) * t47;
	t83 = t113 * t166 + t36 * t107 + (t130 * t89 + t133 * t90) * t47;
	t22 = -t114 * t83 - t115 * t84;
	t70 = t22 * t154;
	t69 = t155 * t154;
	t68 = pkin(4) * (t114 * t84 - t115 * t83);
	t66 = -t147 * t67 / 0.2e1;
	t62 = (-t81 * t139 + t155 * t161 + t162 * t69) * t137;
	t61 = (t155 * t157 + t163 * t69 + t81 * t72) * t121;
	t60 = (t161 * t22 + t162 * t70 - t55 * t68) * t137;
	t59 = (t157 * t22 + t163 * t70 + t73 * t68) * t121;
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t38 = t53 * t117;
	t30 = t32 * t54;
	t29 = t33 * t54;
	t28 = -t118 * t53 - t38;
	t27 = t33 * t53;
	t19 = t20 * t54;
	t18 = t21 * t54;
	t17 = t20 * t53;
	t16 = t100 * t53 - t38;
	t10 = t121 * t15 * t52 + t51 * t66;
	t9 = t136 * t147 * t165 + t52 * t66;
	t4 = t155 * t159 + t164 * t62 + t51 * t61;
	t3 = t155 * t158 + t165 * t62 + t52 * t61;
	t2 = t159 * t22 + t164 * t60 + t51 * t59;
	t1 = t158 * t22 + t165 * t60 + t52 * t59;
	t5 = [(-t10 * t27 - t28 * t9) * r_i_i_C(1) + (-t10 * t28 + t27 * t9) * r_i_i_C(2) + t54 * r_i_i_C(3) - t110 * t53 + (-t27 * t52 + t28 * t51) * pkin(3), (-t1 * t30 + t10 * t18 - t19 * t9 + t2 * t29) * r_i_i_C(1) + (-t1 * t29 - t10 * t19 - t18 * t9 - t2 * t30) * r_i_i_C(2) - t54 * t119 + (t18 * t52 + t19 * t51) * pkin(3), (t29 * t4 - t3 * t30) * r_i_i_C(1) + (-t29 * t3 - t30 * t4) * r_i_i_C(2), 0; (t10 * t29 - t30 * t9) * r_i_i_C(1) + (-t10 * t30 - t29 * t9) * r_i_i_C(2) + t53 * r_i_i_C(3) + t110 * t54 + (t29 * t52 + t30 * t51) * pkin(3), (t1 * t28 + t10 * t16 - t17 * t9 + t2 * t27) * r_i_i_C(1) + (-t1 * t27 - t10 * t17 - t16 * t9 + t2 * t28) * r_i_i_C(2) - t53 * t119 + (t16 * t52 + t17 * t51) * pkin(3), (t27 * t4 + t28 * t3) * r_i_i_C(1) + (-t27 * t3 + t28 * t4) * r_i_i_C(2), 0; 0, (t1 * t33 + t10 * t20 + t2 * t32 + t21 * t9) * r_i_i_C(1) + (-t1 * t32 + t10 * t21 + t2 * t33 - t20 * t9) * r_i_i_C(2) + t120 + (t20 * t52 - t21 * t51) * pkin(3), (t3 * t33 + t32 * t4) * r_i_i_C(1) + (-t3 * t32 + t33 * t4) * r_i_i_C(2), 0;];
	Ja_transl = t5;
end