% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh3m1DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m1DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_jacobia_transl_sym_varpar: pkin has to be [19x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
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
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
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
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
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
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:19:02
	% DurationCPUTime: 23.44s
	% Computational Cost: add. (712897->115), mult. (1075969->208), div. (47728->11), fcn. (682025->19), ass. (0->126)
	t146 = sin(qJ(3));
	t147 = sin(qJ(2));
	t148 = cos(pkin(16));
	t87 = sin(pkin(16));
	t89 = cos(qJ(2));
	t75 = t147 * t148 + t89 * t87;
	t151 = pkin(5) * t75;
	t73 = t147 * t87 - t89 * t148;
	t152 = pkin(5) * t73;
	t163 = -2 * pkin(1);
	t131 = (pkin(1) ^ 2) + t152 * t163;
	t92 = pkin(5) ^ 2;
	t64 = t92 + t131;
	t59 = pkin(2) ^ 2 - pkin(6) ^ 2 + t64;
	t69 = pkin(1) - t152;
	t159 = -pkin(6) - pkin(2);
	t57 = (pkin(5) - t159) * (pkin(5) + t159) + t131;
	t158 = -pkin(6) + pkin(2);
	t58 = (pkin(5) - t158) * (pkin(5) + t158) + t131;
	t98 = sqrt(-t58 * t57);
	t53 = t59 * t151 + t69 * t98;
	t120 = t53 * t146;
	t60 = 0.1e1 / t64;
	t95 = 0.1e1 / pkin(2);
	t135 = t60 * t95;
	t134 = t75 * t98;
	t52 = -pkin(5) * t134 + t59 * t69;
	t88 = cos(qJ(3));
	t138 = t52 * t88;
	t48 = (-t138 / 0.2e1 + t120 / 0.2e1) * t135;
	t121 = t52 * t146;
	t137 = t53 * t88;
	t49 = (t137 / 0.2e1 + t121 / 0.2e1) * t135;
	t83 = pkin(18) + pkin(19);
	t81 = sin(t83);
	t82 = cos(t83);
	t43 = t48 * t82 + t49 * t81;
	t153 = pkin(3) * t43;
	t162 = -2 * pkin(4);
	t132 = (pkin(4) ^ 2) - t153 * t162;
	t94 = pkin(3) ^ 2;
	t39 = t94 + t132;
	t37 = 0.1e1 / t39;
	t91 = 0.1e1 / pkin(10);
	t139 = t37 * t91;
	t84 = sin(pkin(17));
	t160 = t84 / 0.2e1;
	t106 = t48 * t81 - t49 * t82;
	t154 = pkin(3) * t106;
	t36 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t39;
	t40 = pkin(4) + t153;
	t157 = -pkin(8) - pkin(10);
	t34 = (pkin(3) - t157) * (pkin(3) + t157) + t132;
	t156 = -pkin(8) + pkin(10);
	t35 = (pkin(3) - t156) * (pkin(3) + t156) + t132;
	t97 = sqrt(-t35 * t34);
	t25 = -t97 * t154 + t36 * t40;
	t26 = t36 * t154 + t40 * t97;
	t85 = cos(pkin(17));
	t22 = (-t25 * t85 / 0.2e1 + t26 * t160) * t139;
	t23 = (t26 * t85 / 0.2e1 + t25 * t160) * t139;
	t17 = atan2(t23, t22);
	t14 = sin(t17);
	t15 = cos(t17);
	t72 = t147 * t146 - t89 * t88;
	t90 = cos(qJ(1));
	t67 = t72 * t90;
	t74 = -t89 * t146 - t147 * t88;
	t68 = t74 * t90;
	t149 = -t67 * t14 - t68 * t15;
	t164 = t14 * t68 - t67 * t15;
	t100 = t149 * r_i_i_C(1) + t164 * r_i_i_C(2);
	t86 = sin(qJ(1));
	t65 = t74 * t86;
	t66 = t72 * t86;
	t116 = t14 * t65 - t66 * t15;
	t155 = -t66 * t14 - t65 * t15;
	t101 = t155 * r_i_i_C(1) + t116 * r_i_i_C(2);
	t99 = (-t74 * t14 + t72 * t15) * r_i_i_C(1) + (-t14 * t72 - t74 * t15) * r_i_i_C(2);
	t161 = pkin(3) * pkin(4);
	t150 = t89 * pkin(1);
	t21 = 0.1e1 / t22 ^ 2;
	t142 = 0.1e1 / (t21 * t23 ^ 2 + 0.1e1) * t91;
	t141 = t21 * t23;
	t140 = t37 * t85;
	t130 = pkin(1) * t151;
	t136 = 0.2e1 / t98 * (t57 + t58) * t130;
	t129 = 0.1e1 / t39 ^ 2 * t161;
	t128 = t146 / 0.2e1;
	t127 = t147 * pkin(1);
	t31 = 0.1e1 / t97;
	t126 = t31 * t40 / 0.2e1;
	t125 = -t31 * t106 / 0.2e1;
	t124 = t37 * t160;
	t123 = -t140 / 0.2e1;
	t122 = t140 / 0.2e1;
	t119 = t94 * t106 * t162;
	t118 = pkin(13) + t150;
	t117 = t40 * t162 - t36;
	t115 = 0.1e1 / t64 ^ 2 * t130;
	t46 = (t73 * t98 + (-t136 / 0.2e1 - t59 + t69 * t163) * t75) * pkin(5);
	t47 = t69 * t136 / 0.2e1 + t92 * t75 ^ 2 * t163 + (-t73 * t59 - t134) * pkin(5);
	t32 = ((t47 * t88 / 0.2e1 + t46 * t128) * t60 + (t121 + t137) * t115) * t95;
	t33 = ((-t46 * t88 / 0.2e1 + t47 * t128) * t60 + (t120 - t138) * t115) * t95;
	t29 = -t32 * t81 - t33 * t82;
	t114 = t29 * t129;
	t113 = t106 * t129;
	t112 = t25 * t113;
	t111 = t26 * t113;
	t110 = t84 * t114;
	t109 = t85 * t114;
	t105 = 0.2e1 * (t34 + t35) * t161;
	t104 = -t65 * pkin(4) + t101;
	t103 = -t68 * pkin(4) + t100;
	t102 = t72 * pkin(4) + t99;
	t28 = -t32 * t82 + t33 * t81;
	t27 = t106 * t105;
	t24 = t29 * t105;
	t20 = 0.1e1 / t22;
	t19 = t27 * t126 + t106 * t119 + (-t106 * t97 + t43 * t36) * pkin(3);
	t18 = (t106 * t117 + t27 * t125 - t43 * t97) * pkin(3);
	t4 = t24 * t126 + t29 * t119 + (t28 * t36 - t29 * t97) * pkin(3);
	t3 = (t117 * t29 + t24 * t125 - t28 * t97) * pkin(3);
	t2 = ((t85 * t111 + t84 * t112 + t19 * t122 + t18 * t124) * t20 - (t84 * t111 - t85 * t112 + t18 * t123 + t19 * t124) * t141) * t142;
	t1 = ((t26 * t109 + t25 * t110 + t4 * t122 + t3 * t124) * t20 - (-t25 * t109 + t26 * t110 + t3 * t123 + t4 * t124) * t141) * t142;
	t5 = [-t66 * pkin(4) + t116 * r_i_i_C(1) - t155 * r_i_i_C(2) + t90 * r_i_i_C(3) - t118 * t86, t100 * t1 - t90 * t127 + t103, t100 * t2 + t103, 0; t67 * pkin(4) - r_i_i_C(1) * t164 + t149 * r_i_i_C(2) + t86 * r_i_i_C(3) + t118 * t90, t101 * t1 - t86 * t127 + t104, t101 * t2 + t104, 0; 0, t99 * t1 + t102 + t150, t99 * t2 + t102, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:20:10
	% DurationCPUTime: 66.65s
	% Computational Cost: add. (1815820->129), mult. (2740335->232), div. (121720->11), fcn. (1736925->21), ass. (0->138)
	t106 = cos(pkin(17));
	t187 = pkin(3) * pkin(4);
	t117 = pkin(3) ^ 2;
	t104 = pkin(18) + pkin(19);
	t102 = sin(t104);
	t103 = cos(t104);
	t174 = sin(qJ(3));
	t109 = sin(pkin(16));
	t112 = cos(qJ(2));
	t175 = sin(qJ(2));
	t176 = cos(pkin(16));
	t93 = t175 * t109 - t112 * t176;
	t179 = pkin(5) * t93;
	t189 = -2 * pkin(1);
	t156 = (pkin(1) ^ 2) + t179 * t189;
	t186 = -pkin(6) - pkin(2);
	t77 = (pkin(5) - t186) * (pkin(5) + t186) + t156;
	t185 = -pkin(6) + pkin(2);
	t78 = (pkin(5) - t185) * (pkin(5) + t185) + t156;
	t121 = sqrt(-t78 * t77);
	t95 = t112 * t109 + t175 * t176;
	t178 = pkin(5) * t95;
	t115 = pkin(5) ^ 2;
	t84 = t115 + t156;
	t79 = pkin(2) ^ 2 - pkin(6) ^ 2 + t84;
	t89 = pkin(1) - t179;
	t73 = t121 * t89 + t79 * t178;
	t146 = t73 * t174;
	t111 = cos(qJ(3));
	t161 = t121 * t95;
	t72 = -pkin(5) * t161 + t89 * t79;
	t159 = t72 * t111;
	t118 = 0.1e1 / pkin(2);
	t80 = 0.1e1 / t84;
	t162 = t118 * t80;
	t68 = (-t159 / 0.2e1 + t146 / 0.2e1) * t162;
	t147 = t72 * t174;
	t158 = t73 * t111;
	t69 = (t158 / 0.2e1 + t147 / 0.2e1) * t162;
	t63 = t102 * t69 + t103 * t68;
	t180 = pkin(3) * t63;
	t188 = -2 * pkin(4);
	t157 = (pkin(4) ^ 2) - t180 * t188;
	t59 = t117 + t157;
	t152 = 0.1e1 / t59 ^ 2 * t187;
	t134 = t106 * t152;
	t153 = pkin(1) * t178;
	t141 = 0.1e1 / t84 ^ 2 * t153;
	t151 = t174 / 0.2e1;
	t171 = 0.2e1 / t121 * (t77 + t78) * t153;
	t66 = (t93 * t121 + (-t171 / 0.2e1 - t79 + t89 * t189) * t95) * pkin(5);
	t67 = t89 * t171 / 0.2e1 + t115 * t95 ^ 2 * t189 + (-t93 * t79 - t161) * pkin(5);
	t52 = ((t67 * t111 / 0.2e1 + t66 * t151) * t80 + (t147 + t158) * t141) * t118;
	t53 = ((-t66 * t111 / 0.2e1 + t67 * t151) * t80 + (t146 - t159) * t141) * t118;
	t49 = -t102 * t52 - t103 * t53;
	t131 = t49 * t134;
	t105 = sin(pkin(17));
	t135 = t105 * t152;
	t133 = t49 * t135;
	t57 = 0.1e1 / t59;
	t160 = t57 * t106;
	t143 = t160 / 0.2e1;
	t144 = -t160 / 0.2e1;
	t177 = t105 / 0.2e1;
	t145 = t57 * t177;
	t114 = 0.1e1 / pkin(10);
	t163 = t114 * t57;
	t184 = -pkin(8) - pkin(10);
	t54 = (pkin(3) - t184) * (pkin(3) + t184) + t157;
	t183 = -pkin(8) + pkin(10);
	t55 = (pkin(3) - t183) * (pkin(3) + t183) + t157;
	t120 = sqrt(-t55 * t54);
	t124 = t102 * t68 - t103 * t69;
	t181 = pkin(3) * t124;
	t56 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t59;
	t60 = pkin(4) + t180;
	t45 = -t120 * t181 + t60 * t56;
	t46 = t120 * t60 + t56 * t181;
	t42 = (-t45 * t106 / 0.2e1 + t46 * t177) * t163;
	t41 = 0.1e1 / t42 ^ 2;
	t43 = (t46 * t106 / 0.2e1 + t45 * t177) * t163;
	t164 = t114 / (t41 * t43 ^ 2 + 0.1e1);
	t172 = t41 * t43;
	t142 = t60 * t188 - t56;
	t51 = 0.1e1 / t120;
	t148 = -t51 * t124 / 0.2e1;
	t123 = 0.2e1 * (t54 + t55) * t187;
	t44 = t49 * t123;
	t48 = t102 * t53 - t103 * t52;
	t23 = (-t48 * t120 + t142 * t49 + t44 * t148) * pkin(3);
	t140 = t117 * t124 * t188;
	t149 = t51 * t60 / 0.2e1;
	t24 = t44 * t149 + t49 * t140 + (-t120 * t49 + t48 * t56) * pkin(3);
	t40 = 0.1e1 / t42;
	t203 = ((t46 * t131 + t45 * t133 + t24 * t143 + t23 * t145) * t40 - (-t45 * t131 + t46 * t133 + t23 * t144 + t24 * t145) * t172) * t164 + 0.1e1;
	t130 = t124 * t134;
	t132 = t124 * t135;
	t47 = t124 * t123;
	t38 = (-t63 * t120 + t124 * t142 + t47 * t148) * pkin(3);
	t39 = t47 * t149 + t124 * t140 + (-t120 * t124 + t63 * t56) * pkin(3);
	t202 = ((t46 * t130 + t45 * t132 + t39 * t143 + t38 * t145) * t40 - (-t45 * t130 + t46 * t132 + t38 * t144 + t39 * t145) * t172) * t164 + 0.1e1;
	t107 = sin(qJ(4));
	t110 = cos(qJ(4));
	t122 = r_i_i_C(1) * t110 - r_i_i_C(2) * t107 + pkin(9);
	t37 = atan2(t43, t42);
	t34 = sin(t37);
	t35 = cos(t37);
	t113 = cos(qJ(1));
	t92 = -t112 * t111 + t175 * t174;
	t87 = t92 * t113;
	t94 = -t175 * t111 - t112 * t174;
	t88 = t94 * t113;
	t194 = t87 * t34 + t88 * t35;
	t204 = t122 * t194;
	t201 = t122 * t203;
	t200 = t122 * t202;
	t182 = pkin(11) + r_i_i_C(3);
	t199 = t182 * t203;
	t198 = t182 * t202;
	t108 = sin(qJ(1));
	t85 = t94 * t108;
	t86 = t92 * t108;
	t192 = -t85 * t34 + t86 * t35;
	t197 = t113 * t107 - t110 * t192;
	t196 = t107 * t192 + t113 * t110;
	t195 = -t86 * t34 - t85 * t35;
	t193 = -t94 * t34 + t92 * t35;
	t191 = -t88 * t34 + t87 * t35;
	t190 = t92 * t34 + t94 * t35;
	t173 = t112 * pkin(1);
	t150 = t175 * pkin(1);
	t139 = pkin(13) + t173;
	t91 = t92 * pkin(4);
	t83 = t88 * pkin(4);
	t82 = t85 * pkin(4);
	t16 = t108 * t107 + t110 * t191;
	t15 = -t107 * t191 + t108 * t110;
	t1 = [-t86 * pkin(4) - t192 * pkin(9) + t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t139 * t108 + t182 * t195, -t113 * t150 + t191 * t199 - t203 * t204 - t83, t191 * t198 - t202 * t204 - t83, r_i_i_C(1) * t15 - r_i_i_C(2) * t16; t87 * pkin(4) + pkin(9) * t191 + t16 * r_i_i_C(1) + t15 * r_i_i_C(2) + t139 * t113 + t182 * t194, -t108 * t150 + t192 * t199 + t195 * t201 - t82, t192 * t198 + t195 * t200 - t82, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2); 0, t190 * t199 + t193 * t201 + t173 + t91, t190 * t198 + t193 * t200 + t91, (-t107 * r_i_i_C(1) - t110 * r_i_i_C(2)) * t190;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:15
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
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:16
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (12676->58), mult. (19555->99), div. (856->7), fcn. (12637->13), ass. (0->59)
	t35 = pkin(5) ^ 2;
	t30 = sin(qJ(2));
	t32 = sin(pkin(16));
	t33 = cos(qJ(2));
	t55 = cos(pkin(16));
	t25 = t30 * t32 - t33 * t55;
	t60 = pkin(5) * t25;
	t64 = -2 * pkin(1);
	t50 = (pkin(1) ^ 2) + t60 * t64;
	t22 = t35 + t50;
	t19 = pkin(2) ^ 2 - pkin(6) ^ 2 + t22;
	t23 = pkin(1) - t60;
	t26 = t30 * t55 + t33 * t32;
	t62 = -pkin(6) - pkin(2);
	t17 = (pkin(5) - t62) * (pkin(5) + t62) + t50;
	t61 = -pkin(6) + pkin(2);
	t18 = (pkin(5) - t61) * (pkin(5) + t61) + t50;
	t38 = sqrt(-t18 * t17);
	t51 = t26 * t38;
	t12 = -pkin(5) * t51 + t23 * t19;
	t59 = pkin(5) * t26;
	t13 = t19 * t59 + t23 * t38;
	t29 = cos(pkin(19));
	t20 = 0.1e1 / t22;
	t36 = 0.1e1 / pkin(2);
	t52 = t20 * t36;
	t28 = sin(pkin(19));
	t63 = t28 / 0.2e1;
	t10 = (-t12 * t29 / 0.2e1 + t13 * t63) * t52;
	t11 = (t13 * t29 / 0.2e1 + t12 * t63) * t52;
	t49 = pkin(1) * t59;
	t45 = 0.1e1 / t22 ^ 2 * t49;
	t43 = t29 * t45;
	t44 = t28 * t45;
	t46 = t20 * t63;
	t53 = t20 * t29;
	t54 = 0.1e1 / t38 * (t17 + t18) * t49;
	t7 = (t25 * t38 + (t23 * t64 - t19 - t54) * t26) * pkin(5);
	t8 = t23 * t54 + t35 * t26 ^ 2 * t64 + (-t25 * t19 - t51) * pkin(5);
	t9 = 0.1e1 / t10 ^ 2;
	t1 = ((t8 * t53 / 0.2e1 + t13 * t43 + t7 * t46 + t12 * t44) / t10 - (-t7 * t53 / 0.2e1 - t12 * t43 + t8 * t46 + t13 * t44) * t11 * t9) / (t11 ^ 2 * t9 + 0.1e1) * t36;
	t66 = t1 + 0.1e1;
	t6 = atan2(t11, t10);
	t4 = cos(t6);
	t47 = t66 * t4;
	t3 = sin(t6);
	t48 = t66 * t3;
	t65 = t47 * r_i_i_C(1) - t48 * r_i_i_C(2) + pkin(1);
	t58 = t30 * t3;
	t57 = t30 * t4;
	t56 = t33 * t4;
	t42 = -t33 * t3 - t57;
	t41 = -t56 + t58;
	t40 = t42 * r_i_i_C(1);
	t39 = -t48 * r_i_i_C(1) - t47 * r_i_i_C(2);
	t34 = cos(qJ(1));
	t31 = sin(qJ(1));
	t2 = t31 * t58;
	t5 = [t2 * r_i_i_C(1) + t34 * r_i_i_C(3) + (r_i_i_C(2) * t57 - pkin(13) + (-t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - pkin(1)) * t33) * t31, (-t65 * t30 + t39 * t33) * t34, 0, 0; t31 * r_i_i_C(3) + (t33 * pkin(1) - t41 * r_i_i_C(1) + t42 * r_i_i_C(2) + pkin(13)) * t34, t2 * r_i_i_C(2) + (t40 - r_i_i_C(2) * t56 - t30 * pkin(1) + (t41 * r_i_i_C(2) + t40) * t1) * t31, 0, 0; 0, t39 * t30 + t65 * t33, 0, 0;];
	Ja_transl = t5;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:19:11
	% DurationCPUTime: 28.55s
	% Computational Cost: add. (694749->162), mult. (1049593->286), div. (45832->15), fcn. (666143->24), ass. (0->146)
	t51 = pkin(2) ^ 2;
	t153 = -pkin(6) ^ 2 + t51;
	t133 = sin(qJ(2));
	t134 = sin(pkin(16));
	t136 = cos(qJ(2));
	t137 = cos(pkin(16));
	t118 = t133 * t134 - t136 * t137;
	t117 = pkin(5) * t118;
	t114 = (-0.2e1 * t117 + pkin(1)) * pkin(1);
	t150 = pkin(5) ^ 2;
	t40 = t114 + t150;
	t112 = t40 + t153;
	t115 = -t117 + pkin(1);
	t41 = t133 * t137 + t136 * t134;
	t138 = pkin(5) * t41;
	t152 = 0.1e1 / pkin(2);
	t142 = -pkin(6) + pkin(2);
	t143 = -pkin(6) - pkin(2);
	t50 = sqrt(-((pkin(5) - t142) * (pkin(5) + t142) + t114) * ((pkin(5) - t143) * (pkin(5) + t143) + t114));
	t109 = t152 * (t115 * t112 - t50 * t138);
	t146 = 0.1e1 / t40;
	t108 = t146 * t109;
	t106 = -t108 / 0.2e1;
	t111 = pkin(5) * t112;
	t110 = t152 * (t41 * t111 + t115 * t50);
	t107 = t146 * t110 / 0.2e1;
	t132 = sin(qJ(3));
	t135 = cos(qJ(3));
	t101 = t135 * t106 + t132 * t107;
	t105 = t108 / 0.2e1;
	t102 = t132 * t105 + t135 * t107;
	t129 = pkin(18) + pkin(19);
	t123 = sin(t129);
	t124 = cos(t129);
	t93 = t123 * t101 - t124 * t102;
	t163 = pkin(4) * t93;
	t141 = (-pkin(8) - pkin(10));
	t147 = -2 * pkin(3);
	t94 = t124 * t101 + t123 * t102;
	t86 = (-t147 * t94 + pkin(4)) * pkin(4);
	t77 = ((pkin(3) - t141) * (pkin(3) + t141)) + t86;
	t140 = (-pkin(8) + pkin(10));
	t78 = ((pkin(3) - t140) * (pkin(3) + t140)) + t86;
	t49 = sqrt(-t78 * t77);
	t164 = t163 * t49;
	t91 = pkin(4) * t94;
	t52 = pkin(8) ^ 2;
	t151 = pkin(4) ^ 2;
	t85 = t151 + (0.2e1 * t91 + pkin(3)) * pkin(3);
	t81 = -pkin(10) ^ 2 + t52 + t85;
	t79 = pkin(4) * t81;
	t162 = t93 * t79;
	t87 = -t91 - pkin(3);
	t161 = 0.2e1 * t87;
	t130 = pkin(1) * t138;
	t131 = 0.4e1 / t50 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t114 - t51) * t130;
	t32 = (t118 * t50 + (-t131 / 0.2e1 - t150 - t153) * t41) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t117) * t130;
	t160 = -t32 / 0.2e1;
	t83 = 0.1e1 / t85;
	t159 = t83 / 0.2e1;
	t158 = -t87 / 0.2e1;
	t157 = -t163 / 0.2e1;
	t43 = sin(pkin(19));
	t45 = cos(pkin(19));
	t35 = t45 * t106 + t43 * t107;
	t36 = t43 * t105 + t45 * t107;
	t31 = atan2(t36, t35);
	t28 = sin(t31);
	t125 = t133 * t28;
	t29 = cos(t31);
	t27 = t136 * t29;
	t154 = -t125 + t27;
	t103 = t109 * t130;
	t104 = t110 * t130;
	t139 = t152 * t146;
	t128 = t139 / 0.2e1;
	t122 = t43 * t128;
	t33 = t115 * t131 / 0.2e1 - t118 * t111 + (-t50 - 0.2e1 * t130) * t138;
	t34 = 0.1e1 / t35 ^ 2;
	t39 = 0.1e1 / t40 ^ 2;
	t19 = ((t33 * t45 * t128 + t32 * t122 + (t43 * t103 + t45 * t104) * t39) / t35 - (t45 * t139 * t160 + t33 * t122 + (-t45 * t103 + t43 * t104) * t39) * t36 * t34) / (t34 * t36 ^ 2 + 0.1e1);
	t156 = t154 * t19;
	t144 = pkin(3) * pkin(4);
	t155 = t144 / t85 ^ 2;
	t148 = 0.1e1 / pkin(8);
	t145 = 0.1e1 / t49;
	t127 = t136 * pkin(1);
	t126 = t133 * pkin(1);
	t121 = t135 * t139;
	t120 = t127 + pkin(13);
	t119 = t132 * t128;
	t25 = t133 * t29 + t136 * t28;
	t113 = -t27 - t156;
	t15 = (-t19 - 0.1e1) * t25;
	t98 = t121 * t160 + t33 * t119 + (-t135 * t103 + t132 * t104) * t39;
	t97 = t33 * t121 / 0.2e1 + t32 * t119 + (t132 * t103 + t135 * t104) * t39;
	t89 = pkin(3) * t163;
	t88 = t151 * t93 * t147;
	t82 = t148 * t83;
	t80 = t82 / 0.2e1;
	t76 = -t123 * t97 - t124 * t98;
	t75 = pkin(4) * t76;
	t74 = pkin(4) * (t123 * t98 - t124 * t97);
	t73 = pkin(3) * t75;
	t72 = -t87 * t49 + t162;
	t71 = -t87 * t81 - t164;
	t70 = 0.1e1 / t71 ^ 2;
	t69 = t148 * t72;
	t68 = t148 * t71;
	t67 = 0.4e1 * t145 * (((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t86) * t89 - t93 * t52 * t144);
	t66 = t69 * t155;
	t65 = t68 * t155;
	t64 = 0.2e1 * t145 * (t77 + t78) * t73;
	t63 = 0.1e1 / (t72 ^ 2 * t70 + 0.1e1) / t82;
	t62 = atan2(t69 * t159, t68 * t159);
	t61 = sin(t62);
	t60 = 0.2e1 / t71 * t63;
	t59 = -0.2e1 * t70 * t72 * t63;
	t58 = ((t67 * t158 + t94 * t79 + t88 * t93 - t164) * t80 + t93 * t66) * t60 + ((t67 * t157 + t89 * t161 - t49 * t91 - t162) * t80 + t93 * t65) * t59;
	t18 = cos(t62);
	t57 = t18 * t58;
	t56 = ((t64 * t158 - t49 * t75 + t81 * t74 + t76 * t88) * t80 + t76 * t66) * t60 + ((t64 * t157 + t73 * t161 - t49 * t74 - t81 * t75) * t80 + t76 * t65) * t59;
	t55 = t18 * t56;
	t54 = t58 * t61;
	t53 = t56 * t61;
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = cos(pkin(18));
	t44 = sin(pkin(18));
	t26 = t47 * t125;
	t23 = t154 * t48;
	t22 = t25 * t48;
	t21 = t47 * t27 - t26;
	t20 = t25 * t47;
	t16 = t154 + t156;
	t14 = t15 * t48;
	t13 = (t125 + t113) * t48;
	t12 = t15 * t47;
	t11 = t113 * t47 + t26;
	t10 = -t46 * t18 - t44 * t61;
	t9 = t18 * t44 - t46 * t61;
	t4 = -t44 * t57 + t46 * t54;
	t3 = -t44 * t54 - t46 * t57;
	t2 = -t44 * t55 + t46 * t53;
	t1 = -t44 * t53 - t46 * t55;
	t5 = [(-t10 * t21 + t20 * t9) * r_i_i_C(1) + (t10 * t20 + t21 * t9) * r_i_i_C(2) + t48 * r_i_i_C(3) - t120 * t47 + (-t20 * t44 - t21 * t46) * pkin(3), (-t1 * t22 + t10 * t14 + t13 * t9 + t2 * t23) * r_i_i_C(1) + (-t1 * t23 + t10 * t13 - t14 * t9 - t2 * t22) * r_i_i_C(2) - t48 * t126 + (-t13 * t44 + t14 * t46) * pkin(3), (-t22 * t3 + t23 * t4) * r_i_i_C(1) + (-t22 * t4 - t23 * t3) * r_i_i_C(2), 0; (t10 * t23 - t22 * t9) * r_i_i_C(1) + (-t10 * t22 - t23 * t9) * r_i_i_C(2) + t47 * r_i_i_C(3) + t120 * t48 + (t22 * t44 + t23 * t46) * pkin(3), (-t1 * t20 + t10 * t12 + t11 * t9 + t2 * t21) * r_i_i_C(1) + (-t1 * t21 + t10 * t11 - t12 * t9 - t2 * t20) * r_i_i_C(2) - t47 * t126 + (-t11 * t44 + t12 * t46) * pkin(3), (-t20 * t3 + t21 * t4) * r_i_i_C(1) + (-t20 * t4 - t21 * t3) * r_i_i_C(2), 0; 0, (t1 * t154 + t10 * t16 + t15 * t9 + t2 * t25) * r_i_i_C(1) + (-t1 * t25 + t10 * t15 + t154 * t2 - t16 * t9) * r_i_i_C(2) + t127 + (-t15 * t44 + t16 * t46) * pkin(3), (t154 * t3 + t25 * t4) * r_i_i_C(1) + (t154 * t4 - t25 * t3) * r_i_i_C(2), 0;];
	Ja_transl = t5;
end