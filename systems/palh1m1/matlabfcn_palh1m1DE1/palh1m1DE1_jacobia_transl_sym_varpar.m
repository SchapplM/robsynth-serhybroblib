% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh1m1DE1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m1DE1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_jacobia_transl_sym_varpar: pkin has to be [23x1] (double)');
Ja_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
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
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:58
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
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:58
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
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:44:06
	% DurationCPUTime: 23.77s
	% Computational Cost: add. (712897->116), mult. (1075969->210), div. (47728->11), fcn. (682025->19), ass. (0->127)
	t149 = sin(pkin(19));
	t88 = sin(qJ(2));
	t91 = cos(qJ(2));
	t93 = cos(pkin(19));
	t74 = -t91 * t149 + t88 * t93;
	t152 = pkin(7) * t74;
	t165 = -2 * pkin(1);
	t131 = (pkin(1) ^ 2) + t152 * t165;
	t160 = -pkin(8) - pkin(3);
	t58 = (pkin(7) - t160) * (pkin(7) + t160) + t131;
	t159 = -pkin(8) + pkin(3);
	t59 = (pkin(7) - t159) * (pkin(7) + t159) + t131;
	t101 = sqrt(-t59 * t58);
	t76 = t88 * t149 + t91 * t93;
	t151 = pkin(7) * t76;
	t95 = pkin(7) ^ 2;
	t65 = t95 + t131;
	t60 = pkin(3) ^ 2 - pkin(8) ^ 2 + t65;
	t70 = pkin(1) - t152;
	t54 = t101 * t70 + t60 * t151;
	t90 = cos(qJ(3));
	t135 = t90 * t54;
	t130 = t101 * t76;
	t53 = -pkin(7) * t130 + t70 * t60;
	t87 = sin(qJ(3));
	t139 = t87 * t53;
	t61 = 0.1e1 / t65;
	t98 = 0.1e1 / pkin(3);
	t143 = t61 * t98;
	t49 = (t139 / 0.2e1 + t135 / 0.2e1) * t143;
	t136 = t90 * t53;
	t138 = t87 * t54;
	t50 = (-t136 / 0.2e1 + t138 / 0.2e1) * t143;
	t84 = pkin(23) + pkin(22);
	t82 = sin(t84);
	t83 = cos(t84);
	t42 = t83 * t49 - t82 * t50;
	t154 = pkin(4) * t42;
	t164 = -2 * pkin(5);
	t132 = (pkin(5) ^ 2) - t154 * t164;
	t97 = pkin(4) ^ 2;
	t39 = t97 + t132;
	t37 = 0.1e1 / t39;
	t94 = 0.1e1 / pkin(11);
	t145 = t37 * t94;
	t85 = sin(pkin(21));
	t162 = t85 / 0.2e1;
	t158 = -pkin(9) - pkin(11);
	t34 = (pkin(4) - t158) * (pkin(4) + t158) + t132;
	t157 = -pkin(9) + pkin(11);
	t35 = (pkin(4) - t157) * (pkin(4) + t157) + t132;
	t100 = sqrt(-t35 * t34);
	t115 = t82 * t49 + t83 * t50;
	t153 = pkin(4) * t115;
	t36 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t39;
	t40 = pkin(5) + t154;
	t25 = -t100 * t153 + t40 * t36;
	t26 = t100 * t40 + t36 * t153;
	t86 = cos(pkin(21));
	t22 = (t25 * t162 + t26 * t86 / 0.2e1) * t145;
	t23 = (-t25 * t86 / 0.2e1 + t26 * t162) * t145;
	t17 = atan2(t22, t23);
	t14 = sin(t17);
	t15 = cos(t17);
	t137 = t88 * t87;
	t75 = t91 * t90 - t137;
	t89 = sin(qJ(1));
	t66 = t75 * t89;
	t109 = t91 * t87 + t88 * t90;
	t67 = t109 * t89;
	t119 = t67 * t14 - t66 * t15;
	t156 = -t66 * t14 - t67 * t15;
	t104 = t156 * r_i_i_C(1) + t119 * r_i_i_C(2);
	t92 = cos(qJ(1));
	t133 = t91 * t92;
	t68 = t90 * t133 - t92 * t137;
	t69 = t109 * t92;
	t110 = -t69 * t14 + t68 * t15;
	t155 = -t68 * t14 - t69 * t15;
	t103 = t155 * r_i_i_C(1) - t110 * r_i_i_C(2);
	t102 = (-t109 * t14 + t75 * t15) * r_i_i_C(1) + (-t109 * t15 - t75 * t14) * r_i_i_C(2);
	t163 = pkin(4) * pkin(5);
	t161 = -t90 / 0.2e1;
	t150 = t88 * pkin(1);
	t21 = 0.1e1 / t23 ^ 2;
	t148 = 0.1e1 / (t22 ^ 2 * t21 + 0.1e1) * t94;
	t147 = t21 * t22;
	t146 = t37 * t86;
	t129 = pkin(1) * t151;
	t144 = 0.2e1 / t101 * (t58 + t59) * t129;
	t128 = 0.1e1 / t39 ^ 2 * t163;
	t31 = 0.1e1 / t100;
	t127 = t31 * t40 / 0.2e1;
	t126 = -t31 * t115 / 0.2e1;
	t125 = t37 * t162;
	t124 = -t146 / 0.2e1;
	t123 = t146 / 0.2e1;
	t122 = t115 * t97 * t164;
	t121 = -pkin(16) + t150;
	t120 = t40 * t164 - t36;
	t118 = 0.1e1 / t65 ^ 2 * t129;
	t46 = (t74 * t101 + (-t144 / 0.2e1 - t60 + t70 * t165) * t76) * pkin(7);
	t47 = t70 * t144 / 0.2e1 + t95 * t76 ^ 2 * t165 + (-t74 * t60 - t130) * pkin(7);
	t32 = ((-t87 * t46 / 0.2e1 + t47 * t161) * t61 + (-t135 - t139) * t118) * t98;
	t33 = ((t46 * t161 + t87 * t47 / 0.2e1) * t61 + (-t136 + t138) * t118) * t98;
	t28 = t83 * t32 + t82 * t33;
	t117 = t28 * t128;
	t116 = t115 * t128;
	t114 = t25 * t117;
	t113 = t25 * t116;
	t112 = t26 * t117;
	t111 = t26 * t116;
	t108 = 0.2e1 * (t34 + t35) * t163;
	t107 = -t67 * pkin(5) + t104;
	t106 = -t69 * pkin(5) + t103;
	t105 = t75 * pkin(5) + t102;
	t29 = -t82 * t32 + t83 * t33;
	t27 = t115 * t108;
	t24 = t28 * t108;
	t20 = 0.1e1 / t23;
	t19 = t27 * t127 + t115 * t122 + (-t100 * t115 + t42 * t36) * pkin(4);
	t18 = (-t42 * t100 + t115 * t120 + t27 * t126) * pkin(4);
	t4 = t24 * t127 + t28 * t122 + (-t100 * t28 + t29 * t36) * pkin(4);
	t3 = (-t29 * t100 + t120 * t28 + t24 * t126) * pkin(4);
	t2 = ((t86 * t111 + t85 * t113 + t19 * t123 + t18 * t125) * t20 - (t85 * t111 - t86 * t113 + t18 * t124 + t19 * t125) * t147) * t148;
	t1 = ((t86 * t112 + t85 * t114 + t4 * t123 + t3 * t125) * t20 - (t85 * t112 - t86 * t114 + t3 * t124 + t4 * t125) * t147) * t148;
	t5 = [-t66 * pkin(5) + t119 * r_i_i_C(1) - t156 * r_i_i_C(2) + t92 * r_i_i_C(3) + t121 * t89, -pkin(1) * t133 + t103 * t1 + t106, t103 * t2 + t106, 0; t68 * pkin(5) + t110 * r_i_i_C(1) + t155 * r_i_i_C(2) + t89 * r_i_i_C(3) - t121 * t92, -t89 * t91 * pkin(1) + t104 * t1 + t107, t104 * t2 + t107, 0; 0, t102 * t1 + t105 - t150, t102 * t2 + t105, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:45:13
	% DurationCPUTime: 65.62s
	% Computational Cost: add. (1815820->129), mult. (2740335->232), div. (121720->11), fcn. (1736925->21), ass. (0->138)
	t107 = cos(pkin(21));
	t188 = pkin(4) * pkin(5);
	t119 = pkin(4) ^ 2;
	t105 = pkin(23) + pkin(22);
	t103 = sin(t105);
	t104 = cos(t105);
	t120 = 0.1e1 / pkin(3);
	t117 = pkin(7) ^ 2;
	t110 = sin(qJ(2));
	t115 = cos(pkin(19));
	t175 = sin(pkin(19));
	t176 = cos(qJ(2));
	t94 = t110 * t115 - t176 * t175;
	t180 = pkin(7) * t94;
	t190 = -2 * pkin(1);
	t155 = (pkin(1) ^ 2) + t180 * t190;
	t85 = t117 + t155;
	t81 = 0.1e1 / t85;
	t159 = t120 * t81;
	t113 = cos(qJ(3));
	t187 = -pkin(8) - pkin(3);
	t78 = (pkin(7) - t187) * (pkin(7) + t187) + t155;
	t186 = -pkin(8) + pkin(3);
	t79 = (pkin(7) - t186) * (pkin(7) + t186) + t155;
	t123 = sqrt(-t79 * t78);
	t96 = t110 * t175 + t176 * t115;
	t179 = pkin(7) * t96;
	t80 = pkin(3) ^ 2 - pkin(8) ^ 2 + t85;
	t90 = pkin(1) - t180;
	t74 = t123 * t90 + t80 * t179;
	t162 = t113 * t74;
	t109 = sin(qJ(3));
	t158 = t123 * t96;
	t73 = -pkin(7) * t158 + t90 * t80;
	t165 = t109 * t73;
	t69 = (t165 / 0.2e1 + t162 / 0.2e1) * t159;
	t163 = t113 * t73;
	t164 = t109 * t74;
	t70 = (-t163 / 0.2e1 + t164 / 0.2e1) * t159;
	t63 = t103 * t70 - t104 * t69;
	t182 = pkin(4) * t63;
	t189 = -2 * pkin(5);
	t156 = (pkin(5) ^ 2) + t182 * t189;
	t59 = t119 + t156;
	t151 = 0.1e1 / t59 ^ 2 * t188;
	t136 = t107 * t151;
	t152 = pkin(1) * t179;
	t142 = 0.1e1 / t85 ^ 2 * t152;
	t177 = -t113 / 0.2e1;
	t172 = 0.2e1 / t123 * (t78 + t79) * t152;
	t66 = (t94 * t123 + (-t172 / 0.2e1 - t80 + t90 * t190) * t96) * pkin(7);
	t67 = t90 * t172 / 0.2e1 + t117 * t96 ^ 2 * t190 + (-t94 * t80 - t158) * pkin(7);
	t52 = ((-t109 * t66 / 0.2e1 + t67 * t177) * t81 + (-t162 - t165) * t142) * t120;
	t53 = ((t66 * t177 + t109 * t67 / 0.2e1) * t81 + (-t163 + t164) * t142) * t120;
	t48 = t103 * t53 + t104 * t52;
	t132 = t48 * t136;
	t106 = sin(pkin(21));
	t137 = t106 * t151;
	t134 = t48 * t137;
	t57 = 0.1e1 / t59;
	t157 = t57 * t107;
	t145 = t157 / 0.2e1;
	t146 = -t157 / 0.2e1;
	t178 = t106 / 0.2e1;
	t147 = t57 * t178;
	t116 = 0.1e1 / pkin(11);
	t160 = t116 * t57;
	t185 = -pkin(9) - pkin(11);
	t54 = (pkin(4) - t185) * (pkin(4) + t185) + t156;
	t184 = -pkin(9) + pkin(11);
	t55 = (pkin(4) - t184) * (pkin(4) + t184) + t156;
	t122 = sqrt(-t55 * t54);
	t135 = t103 * t69 + t104 * t70;
	t181 = pkin(4) * t135;
	t56 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t59;
	t60 = pkin(5) - t182;
	t45 = -t122 * t181 + t60 * t56;
	t46 = t122 * t60 + t56 * t181;
	t43 = (-t45 * t107 / 0.2e1 + t46 * t178) * t160;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = (t45 * t178 + t46 * t107 / 0.2e1) * t160;
	t161 = t116 / (t41 * t42 ^ 2 + 0.1e1);
	t173 = t41 * t42;
	t144 = t60 * t189 - t56;
	t51 = 0.1e1 / t122;
	t148 = -t51 * t135 / 0.2e1;
	t126 = 0.2e1 * (t54 + t55) * t188;
	t44 = t48 * t126;
	t49 = -t103 * t52 + t104 * t53;
	t23 = (-t49 * t122 + t144 * t48 + t148 * t44) * pkin(4);
	t141 = t119 * t135 * t189;
	t149 = t51 * t60 / 0.2e1;
	t24 = t44 * t149 + t48 * t141 + (-t122 * t48 + t49 * t56) * pkin(4);
	t40 = 0.1e1 / t43;
	t203 = ((t132 * t46 + t134 * t45 + t145 * t24 + t147 * t23) * t40 - (-t132 * t45 + t134 * t46 + t146 * t23 + t147 * t24) * t173) * t161 + 0.1e1;
	t131 = t135 * t136;
	t133 = t135 * t137;
	t47 = t135 * t126;
	t38 = (t122 * t63 + t135 * t144 + t148 * t47) * pkin(4);
	t39 = t47 * t149 + t135 * t141 + (-t122 * t135 - t56 * t63) * pkin(4);
	t202 = ((t131 * t46 + t133 * t45 + t145 * t39 + t147 * t38) * t40 - (-t131 * t45 + t133 * t46 + t146 * t38 + t147 * t39) * t173) * t161 + 0.1e1;
	t37 = atan2(t42, t43);
	t34 = sin(t37);
	t35 = cos(t37);
	t111 = sin(qJ(1));
	t95 = -t110 * t109 + t176 * t113;
	t86 = t95 * t111;
	t124 = t176 * t109 + t110 * t113;
	t87 = t124 * t111;
	t18 = t87 * t34 - t86 * t35;
	t183 = pkin(12) + r_i_i_C(3);
	t204 = t183 * t18;
	t108 = sin(qJ(4));
	t112 = cos(qJ(4));
	t125 = r_i_i_C(1) * t112 - r_i_i_C(2) * t108 + pkin(10);
	t201 = t125 * t203;
	t200 = t125 * t202;
	t199 = t183 * t203;
	t198 = t183 * t202;
	t114 = cos(qJ(1));
	t197 = t114 * t108 + t18 * t112;
	t196 = t18 * t108 - t114 * t112;
	t195 = -t86 * t34 - t87 * t35;
	t88 = t95 * t114;
	t89 = t124 * t114;
	t194 = -t88 * t34 - t89 * t35;
	t193 = -t124 * t34 + t95 * t35;
	t192 = -t89 * t34 + t88 * t35;
	t191 = t124 * t35 + t95 * t34;
	t174 = t110 * pkin(1);
	t150 = t176 * pkin(1);
	t140 = -pkin(16) + t174;
	t92 = t95 * pkin(5);
	t84 = t89 * pkin(5);
	t83 = t87 * pkin(5);
	t16 = t111 * t108 + t112 * t192;
	t15 = -t108 * t192 + t111 * t112;
	t1 = [-t86 * pkin(5) + t18 * pkin(10) + t197 * r_i_i_C(1) - t196 * r_i_i_C(2) + t140 * t111 + t183 * t195, -t114 * t150 + t192 * t199 + t194 * t201 - t84, t192 * t198 + t194 * t200 - t84, r_i_i_C(1) * t15 - r_i_i_C(2) * t16; t88 * pkin(5) + pkin(10) * t192 + t16 * r_i_i_C(1) + t15 * r_i_i_C(2) - t140 * t114 - t183 * t194, -t111 * t150 + t195 * t201 - t203 * t204 - t83, t195 * t200 - t202 * t204 - t83, t196 * r_i_i_C(1) + t197 * r_i_i_C(2); 0, t191 * t199 + t193 * t201 - t174 + t92, t191 * t198 + t193 * t200 + t92, (-t108 * r_i_i_C(1) - t112 * r_i_i_C(2)) * t191;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:59
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
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:43:01
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (12676->58), mult. (19555->99), div. (856->7), fcn. (12637->13), ass. (0->60)
	t32 = sin(qJ(2));
	t34 = cos(qJ(2));
	t37 = pkin(7) ^ 2;
	t36 = cos(pkin(19));
	t58 = sin(pkin(19));
	t27 = t32 * t36 - t34 * t58;
	t64 = pkin(7) * t27;
	t70 = -2 * pkin(1);
	t53 = (pkin(1) ^ 2) + t64 * t70;
	t24 = t37 + t53;
	t21 = pkin(3) ^ 2 - pkin(8) ^ 2 + t24;
	t25 = pkin(1) - t64;
	t28 = t32 * t58 + t34 * t36;
	t67 = -pkin(8) - pkin(3);
	t19 = (pkin(7) - t67) * (pkin(7) + t67) + t53;
	t66 = -pkin(8) + pkin(3);
	t20 = (pkin(7) - t66) * (pkin(7) + t66) + t53;
	t40 = sqrt(-t20 * t19);
	t55 = t28 * t40;
	t14 = -pkin(7) * t55 + t25 * t21;
	t63 = pkin(7) * t28;
	t15 = t21 * t63 + t25 * t40;
	t31 = cos(pkin(23));
	t22 = 0.1e1 / t24;
	t38 = 0.1e1 / pkin(3);
	t56 = t22 * t38;
	t30 = sin(pkin(23));
	t68 = t30 / 0.2e1;
	t12 = (-t31 * t14 / 0.2e1 + t15 * t68) * t56;
	t13 = (t31 * t15 / 0.2e1 + t14 * t68) * t56;
	t8 = atan2(t13, t12);
	t5 = sin(t8);
	t6 = cos(t8);
	t44 = t32 * t6 + t34 * t5;
	t52 = pkin(1) * t63;
	t57 = 0.1e1 / t40 * (t19 + t20) * t52;
	t10 = t25 * t57 + t37 * t28 ^ 2 * t70 + (-t27 * t21 - t55) * pkin(7);
	t11 = 0.1e1 / t12 ^ 2;
	t47 = 0.1e1 / t24 ^ 2 * t52;
	t45 = t31 * t47;
	t46 = t30 * t47;
	t49 = t22 * t68;
	t54 = t31 * t22;
	t9 = (t27 * t40 + (t25 * t70 - t21 - t57) * t28) * pkin(7);
	t1 = ((t10 * t54 / 0.2e1 + t15 * t45 + t9 * t49 + t14 * t46) / t12 - (-t9 * t54 / 0.2e1 - t14 * t45 + t10 * t49 + t15 * t46) * t13 * t11) / (t13 ^ 2 * t11 + 0.1e1) * t38;
	t50 = (t1 + 0.1e1) * t6;
	t71 = t50 * r_i_i_C(2);
	t69 = t1 * t5;
	t33 = sin(qJ(1));
	t65 = t44 * t33;
	t62 = t32 * t5;
	t59 = t34 * t6;
	t51 = t5 + t69;
	t48 = t32 * pkin(1) - pkin(16);
	t43 = -t59 + t62;
	t42 = t43 * r_i_i_C(1);
	t41 = -t50 * r_i_i_C(1) + t51 * r_i_i_C(2) - pkin(1);
	t35 = cos(qJ(1));
	t4 = t35 * t62;
	t2 = [t65 * r_i_i_C(1) + t35 * r_i_i_C(3) + (-t43 * r_i_i_C(2) + t48) * t33, t4 * r_i_i_C(1) + ((r_i_i_C(1) * t69 + t71) * t32 + t41 * t34) * t35, 0, 0; t4 * r_i_i_C(2) + t33 * r_i_i_C(3) + (-t44 * r_i_i_C(1) - r_i_i_C(2) * t59 - t48) * t35, t65 * r_i_i_C(2) + (t42 - t34 * pkin(1) + (t44 * r_i_i_C(2) + t42) * t1) * t33, 0, 0; 0, (-t51 * r_i_i_C(1) - t71) * t34 + t41 * t32, 0, 0;];
	Ja_transl = t2;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:59
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (6100->41), mult. (9228->81), div. (296->7), fcn. (5940->13), ass. (0->47)
	t34 = pkin(6) ^ 2;
	t26 = sin(pkin(20));
	t27 = cos(pkin(20));
	t28 = sin(qJ(3));
	t31 = cos(qJ(3));
	t24 = t31 * t26 + t28 * t27;
	t48 = pkin(6) * t24;
	t44 = 0.2e1 * pkin(1) * t48 + t34;
	t21 = pkin(1) ^ 2 + t44;
	t19 = 0.1e1 / t21;
	t53 = t19 / 0.2e1;
	t52 = -pkin(2) - pkin(13);
	t51 = -pkin(2) + pkin(13);
	t25 = t28 * t26 - t31 * t27;
	t50 = pkin(1) * t25;
	t18 = pkin(2) ^ 2 - pkin(13) ^ 2 + t21;
	t22 = -pkin(1) - t48;
	t16 = (pkin(1) - t52) * (pkin(1) + t52) + t44;
	t17 = (pkin(1) - t51) * (pkin(1) + t51) + t44;
	t37 = sqrt(-t17 * t16);
	t47 = t25 * pkin(6);
	t12 = t18 * t47 - t22 * t37;
	t49 = pkin(6) * t12;
	t46 = 0.1e1 / t37 * (t16 + t17) * pkin(1) * t47;
	t45 = t25 * t37;
	t35 = 0.1e1 / pkin(2);
	t43 = t35 * t53;
	t29 = sin(qJ(2));
	t32 = cos(qJ(2));
	t11 = -pkin(6) * t45 - t22 * t18;
	t8 = atan2(t12 * t43, t11 * t43);
	t6 = sin(t8);
	t7 = cos(t8);
	t41 = t29 * t7 + t32 * t6;
	t40 = t29 * t6 - t32 * t7;
	t39 = -t41 * r_i_i_C(1) + t40 * r_i_i_C(2);
	t10 = 0.1e1 / t11 ^ 2;
	t20 = 0.1e1 / t21 ^ 2;
	t1 = 0.2e1 * (((-t22 * t46 + (t24 * t18 - t45) * pkin(6)) * t53 + (-t19 * t34 * t25 + t20 * t49) * t50) / t11 - ((-t24 * t37 + (-t18 - t46) * t25) * t53 + (t11 * t20 + t19 * t22) * t50) * t10 * t49) * pkin(2) * t21 * t35 / (t12 ^ 2 * t10 + 0.1e1);
	t38 = t1 * (t40 * r_i_i_C(1) + t41 * r_i_i_C(2));
	t33 = cos(qJ(1));
	t30 = sin(qJ(1));
	t5 = t41 * t33;
	t4 = t40 * t33;
	t3 = t41 * t30;
	t2 = t40 * t30;
	t9 = [-t30 * pkin(16) + t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + t33 * r_i_i_C(3), t4 * r_i_i_C(1) + t5 * r_i_i_C(2), t33 * t38, 0; t33 * pkin(16) - t5 * r_i_i_C(1) + t4 * r_i_i_C(2) + t30 * r_i_i_C(3), t2 * r_i_i_C(1) + t3 * r_i_i_C(2), t30 * t38, 0; 0, t39, t39 * t1, 0;];
	Ja_transl = t9;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:43:02
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (18644->71), mult. (27920->131), div. (1084->11), fcn. (17748->16), ass. (0->73)
	t45 = sin(pkin(20));
	t46 = cos(pkin(20));
	t47 = sin(qJ(3));
	t50 = cos(qJ(3));
	t43 = t50 * t45 + t47 * t46;
	t84 = pkin(6) * t43;
	t73 = pkin(1) * t84;
	t42 = 0.2e1 * t73;
	t56 = pkin(2) ^ 2;
	t55 = pkin(6) ^ 2;
	t74 = pkin(1) ^ 2 + t55;
	t71 = -pkin(13) ^ 2 + t74;
	t37 = t42 + t56 + t71;
	t41 = -pkin(1) - t84;
	t44 = t47 * t45 - t50 * t46;
	t75 = t42 + t55;
	t88 = -pkin(2) - pkin(13);
	t33 = (pkin(1) - t88) * (pkin(1) + t88) + t75;
	t87 = -pkin(2) + pkin(13);
	t34 = (pkin(1) - t87) * (pkin(1) + t87) + t75;
	t81 = t34 * t33;
	t59 = sqrt(-t81);
	t80 = t44 * t59;
	t72 = pkin(6) * t80;
	t24 = -t41 * t37 - t72;
	t83 = t44 * pkin(6);
	t25 = t37 * t83 - t41 * t59;
	t40 = t42 + t74;
	t38 = 0.1e1 / t40;
	t57 = 0.1e1 / pkin(2);
	t89 = t57 / 0.2e1;
	t70 = t38 * t89;
	t21 = atan2(t25 * t70, t24 * t70);
	t19 = sin(t21);
	t20 = cos(t21);
	t48 = sin(qJ(2));
	t51 = cos(qJ(2));
	t63 = t51 * t19 + t48 * t20;
	t23 = 0.1e1 / t24 ^ 2;
	t39 = 0.1e1 / t40 ^ 2;
	t82 = 0.1e1 / t59 * (t33 + t34) * pkin(1) * t83;
	t85 = pkin(6) * t25;
	t86 = pkin(1) * t44;
	t90 = t38 / 0.2e1;
	t7 = 0.2e1 * (((-t41 * t82 + (t43 * t37 - t80) * pkin(6)) * t90 + (-t38 * t55 * t44 + t39 * t85) * t86) / t24 - ((-t43 * t59 + (-t37 - t82) * t44) * t90 + (t24 * t39 + t38 * t41) * t86) * t23 * t85) * pkin(2) / (t25 ^ 2 * t23 + 0.1e1) * t40 * t57;
	t5 = t63 * t7;
	t12 = t48 * t19 - t51 * t20;
	t69 = 0.1e1 / pkin(13) * t89;
	t36 = t56 - t71 - 0.2e1 * t73;
	t28 = atan2(t59 * t69, t36 * t69);
	t26 = sin(t28);
	t27 = cos(t28);
	t49 = sin(qJ(1));
	t8 = t12 * t49;
	t9 = t63 * t49;
	t67 = t9 * t26 + t8 * t27;
	t66 = t8 * t26 - t9 * t27;
	t52 = cos(qJ(1));
	t10 = t12 * t52;
	t11 = t63 * t52;
	t65 = -t10 * t27 - t11 * t26;
	t64 = t10 * t26 - t11 * t27;
	t6 = t12 * t7;
	t62 = -t67 * r_i_i_C(1) + t66 * r_i_i_C(2);
	t61 = t65 * r_i_i_C(1) + t64 * r_i_i_C(2);
	t60 = (-t12 * t26 + t27 * t63) * r_i_i_C(1) + (-t12 * t27 - t26 * t63) * r_i_i_C(2);
	t35 = 0.1e1 / t36 ^ 2;
	t14 = (0.1e1 / t36 * t82 - 0.2e1 * pkin(1) * t35 * t72) / (-t35 * t81 + 0.1e1);
	t4 = t52 * t5;
	t3 = t7 * t10;
	t2 = t7 * t9;
	t1 = t49 * t6;
	t13 = [t9 * pkin(2) - t49 * pkin(16) + t66 * r_i_i_C(1) + t67 * r_i_i_C(2) + t52 * r_i_i_C(3), t10 * pkin(2) + t61, (-t4 * t26 - t3 * t27) * r_i_i_C(1) + (t3 * t26 - t4 * t27) * r_i_i_C(2) + t3 * pkin(2) + t61 * t14, 0; -t11 * pkin(2) + t52 * pkin(16) - t64 * r_i_i_C(1) + t65 * r_i_i_C(2) + t49 * r_i_i_C(3), t8 * pkin(2) + t62, (-t1 * t27 - t2 * t26) * r_i_i_C(1) + (t1 * t26 - t2 * t27) * r_i_i_C(2) + t1 * pkin(2) + t62 * t14, 0; 0, -pkin(2) * t63 + t60, (-t6 * t26 + t5 * t27) * r_i_i_C(1) + (-t5 * t26 - t6 * t27) * r_i_i_C(2) - t5 * pkin(2) + t60 * t14, 0;];
	Ja_transl = t13;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:44:14
	% DurationCPUTime: 28.58s
	% Computational Cost: add. (694749->164), mult. (1049593->290), div. (45832->15), fcn. (666143->24), ass. (0->147)
	t51 = pkin(3) ^ 2;
	t158 = -pkin(8) ^ 2 + t51;
	t137 = sin(qJ(2));
	t138 = sin(pkin(19));
	t140 = cos(qJ(2));
	t141 = cos(pkin(19));
	t120 = t137 * t141 - t140 * t138;
	t119 = pkin(7) * t120;
	t115 = (-0.2e1 * t119 + pkin(1)) * pkin(1);
	t154 = pkin(7) ^ 2;
	t40 = t115 + t154;
	t113 = t40 + t158;
	t117 = -t119 + pkin(1);
	t41 = t137 * t138 + t140 * t141;
	t142 = pkin(7) * t41;
	t156 = 0.1e1 / pkin(3);
	t146 = -pkin(8) + pkin(3);
	t147 = -pkin(8) - pkin(3);
	t50 = sqrt(-((pkin(7) - t146) * (pkin(7) + t146) + t115) * ((pkin(7) - t147) * (pkin(7) + t147) + t115));
	t109 = t156 * (t117 * t113 - t50 * t142);
	t150 = 0.1e1 / t40;
	t107 = t150 * t109;
	t105 = t107 / 0.2e1;
	t112 = pkin(7) * t113;
	t110 = t156 * (t41 * t112 + t117 * t50);
	t108 = t150 * t110;
	t106 = t108 / 0.2e1;
	t136 = sin(qJ(3));
	t139 = cos(qJ(3));
	t101 = t136 * t105 + t139 * t106;
	t102 = t139 * t105 - t136 * t108 / 0.2e1;
	t133 = pkin(23) + pkin(22);
	t123 = sin(t133);
	t124 = cos(t133);
	t93 = t124 * t101 + t123 * t102;
	t90 = pkin(5) * t93;
	t157 = t123 * t101 - t124 * t102;
	t167 = pkin(5) * t157;
	t145 = (-pkin(9) - pkin(11));
	t151 = -2 * pkin(4);
	t86 = (-t151 * t93 + pkin(5)) * pkin(5);
	t77 = ((pkin(4) - t145) * (pkin(4) + t145)) + t86;
	t144 = (-pkin(9) + pkin(11));
	t78 = ((pkin(4) - t144) * (pkin(4) + t144)) + t86;
	t49 = sqrt(-t78 * t77);
	t168 = t167 * t49;
	t52 = pkin(9) ^ 2;
	t155 = pkin(5) ^ 2;
	t85 = t155 + (0.2e1 * t90 + pkin(4)) * pkin(4);
	t81 = -pkin(11) ^ 2 + t52 + t85;
	t79 = pkin(5) * t81;
	t166 = t157 * t79;
	t43 = sin(pkin(23));
	t45 = cos(pkin(23));
	t35 = -t45 * t107 / 0.2e1 + t43 * t106;
	t36 = t43 * t105 + t45 * t106;
	t31 = atan2(t36, t35);
	t28 = sin(t31);
	t29 = cos(t31);
	t118 = t137 * t29 + t140 * t28;
	t134 = pkin(1) * t142;
	t103 = t109 * t134;
	t104 = t110 * t134;
	t143 = t156 * t150;
	t131 = t143 / 0.2e1;
	t126 = t43 * t131;
	t132 = -t143 / 0.2e1;
	t135 = 0.4e1 / t50 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t115 - t51) * t134;
	t32 = (t120 * t50 + (-t135 / 0.2e1 - t154 - t158) * t41) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t119) * t134;
	t33 = t117 * t135 / 0.2e1 - t120 * t112 + (-t50 - 0.2e1 * t134) * t142;
	t34 = 0.1e1 / t35 ^ 2;
	t39 = 0.1e1 / t40 ^ 2;
	t19 = ((t45 * t33 * t131 + t32 * t126 + (t43 * t103 + t45 * t104) * t39) / t35 - (t45 * t32 * t132 + t33 * t126 + (-t45 * t103 + t43 * t104) * t39) * t36 * t34) / (t34 * t36 ^ 2 + 0.1e1);
	t165 = (t19 + 0.1e1) * t118;
	t87 = -t90 - pkin(4);
	t164 = 0.2e1 * t87;
	t83 = 0.1e1 / t85;
	t163 = t83 / 0.2e1;
	t162 = -t87 / 0.2e1;
	t161 = -t167 / 0.2e1;
	t148 = pkin(4) * pkin(5);
	t159 = t148 / t85 ^ 2;
	t152 = 0.1e1 / pkin(9);
	t149 = 0.1e1 / t49;
	t130 = t140 * pkin(1);
	t129 = t137 * pkin(1);
	t128 = t140 * t29;
	t127 = t137 * t28;
	t125 = t136 * t143;
	t122 = t129 - pkin(16);
	t121 = t139 * t132;
	t25 = t128 - t127;
	t116 = t25 * t19;
	t114 = t127 - t116;
	t98 = t32 * t121 + t33 * t125 / 0.2e1 + (-t139 * t103 + t136 * t104) * t39;
	t97 = -t32 * t125 / 0.2e1 + t33 * t121 + (-t136 * t103 - t139 * t104) * t39;
	t89 = pkin(4) * t167;
	t88 = t155 * t157 * t151;
	t82 = t152 * t83;
	t80 = t82 / 0.2e1;
	t76 = t123 * t98 + t124 * t97;
	t75 = pkin(5) * (-t123 * t97 + t124 * t98);
	t74 = pkin(5) * t76;
	t73 = pkin(4) * t74;
	t72 = -t87 * t49 + t166;
	t71 = -t87 * t81 - t168;
	t70 = 0.1e1 / t71 ^ 2;
	t69 = t152 * t72;
	t68 = t152 * t71;
	t67 = 0.4e1 * t149 * (((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t86) * t89 - t157 * t52 * t148);
	t66 = t69 * t159;
	t65 = t68 * t159;
	t64 = 0.2e1 * t149 * (t77 + t78) * t73;
	t63 = 0.1e1 / (t72 ^ 2 * t70 + 0.1e1) / t82;
	t62 = atan2(t69 * t163, t68 * t163);
	t61 = sin(t62);
	t60 = 0.2e1 / t71 * t63;
	t59 = -0.2e1 * t70 * t72 * t63;
	t58 = ((t157 * t88 + t67 * t162 + t93 * t79 - t168) * t80 + t157 * t66) * t60 + ((t67 * t161 + t89 * t164 - t49 * t90 - t166) * t80 + t157 * t65) * t59;
	t18 = cos(t62);
	t57 = t18 * t58;
	t56 = ((t64 * t162 - t49 * t74 + t81 * t75 + t76 * t88) * t80 + t76 * t66) * t60 + ((t64 * t161 + t73 * t164 - t49 * t75 - t81 * t74) * t80 + t76 * t65) * t59;
	t55 = t18 * t56;
	t54 = t58 * t61;
	t53 = t56 * t61;
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = cos(pkin(22));
	t44 = sin(pkin(22));
	t27 = t48 * t127;
	t26 = t47 * t128;
	t23 = t118 * t48;
	t22 = -t48 * t128 + t27;
	t21 = t118 * t47;
	t20 = t47 * t127 - t26;
	t15 = -t128 + t114;
	t14 = t27 + (-t128 - t116) * t48;
	t13 = t165 * t48;
	t12 = t114 * t47 - t26;
	t11 = t165 * t47;
	t10 = -t46 * t18 - t44 * t61;
	t9 = t18 * t44 - t46 * t61;
	t4 = -t44 * t57 + t46 * t54;
	t3 = -t44 * t54 - t46 * t57;
	t2 = -t44 * t55 + t46 * t53;
	t1 = -t44 * t53 - t46 * t55;
	t5 = [(t10 * t21 - t20 * t9) * r_i_i_C(1) + (-t10 * t20 - t21 * t9) * r_i_i_C(2) + t48 * r_i_i_C(3) + t122 * t47 + (t20 * t44 + t21 * t46) * pkin(4), (t1 * t22 + t10 * t14 + t13 * t9 - t2 * t23) * r_i_i_C(1) + (t1 * t23 + t10 * t13 - t14 * t9 + t2 * t22) * r_i_i_C(2) - t48 * t130 + (-t13 * t44 + t14 * t46) * pkin(4), (t22 * t3 - t23 * t4) * r_i_i_C(1) + (t22 * t4 + t23 * t3) * r_i_i_C(2), 0; (-t10 * t23 + t22 * t9) * r_i_i_C(1) + (t10 * t22 + t23 * t9) * r_i_i_C(2) + t47 * r_i_i_C(3) - t122 * t48 + (-t22 * t44 - t23 * t46) * pkin(4), (t1 * t20 + t10 * t12 + t11 * t9 - t2 * t21) * r_i_i_C(1) + (t1 * t21 + t10 * t11 - t12 * t9 + t2 * t20) * r_i_i_C(2) - t47 * t130 + (-t11 * t44 + t12 * t46) * pkin(4), (t20 * t3 - t21 * t4) * r_i_i_C(1) + (t20 * t4 + t21 * t3) * r_i_i_C(2), 0; 0, (-t1 * t118 - t10 * t165 + t15 * t9 + t2 * t25) * r_i_i_C(1) + (-t1 * t25 + t10 * t15 - t118 * t2 + t165 * t9) * r_i_i_C(2) - t129 + (-t15 * t44 - t165 * t46) * pkin(4), (-t118 * t3 + t25 * t4) * r_i_i_C(1) + (-t118 * t4 - t25 * t3) * r_i_i_C(2), 0;];
	Ja_transl = t5;
end