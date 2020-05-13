% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh3m1TE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energykin_fixb_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 00:29:44
% EndTime: 2020-04-18 00:30:20
% DurationCPUTime: 35.77s
% Computational Cost: add. (982060->178), mult. (1522615->445), div. (62952->22), fcn. (949994->22), ass. (0->195)
t224 = -2 * pkin(1);
t140 = pkin(4) ^ 2;
t115 = pkin(18) + pkin(19);
t112 = sin(t115);
t113 = cos(t115);
t123 = sin(qJ(3));
t124 = sin(qJ(2));
t125 = sin(pkin(16));
t128 = cos(qJ(2));
t129 = cos(pkin(16));
t108 = t124 * t129 + t125 * t128;
t144 = pkin(1) ^ 2;
t107 = t124 * t125 - t128 * t129;
t203 = pkin(5) * t107;
t176 = t203 * t224 + t144;
t219 = -pkin(6) - pkin(2);
t92 = (pkin(5) - t219) * (pkin(5) + t219) + t176;
t218 = -pkin(6) + pkin(2);
t93 = (pkin(5) - t218) * (pkin(5) + t218) + t176;
t146 = sqrt(-t93 * t92);
t174 = pkin(2) ^ 2 - pkin(6) ^ 2;
t139 = pkin(5) ^ 2;
t98 = t139 + t176;
t95 = t98 + t174;
t99 = pkin(1) - t203;
t86 = pkin(5) * t108 * t95 + t146 * t99;
t184 = t86 * t123;
t127 = cos(qJ(3));
t179 = t108 * t146;
t85 = -pkin(5) * t179 + t95 * t99;
t187 = t85 * t127;
t154 = -t184 / 0.2e1 + t187 / 0.2e1;
t143 = 0.1e1 / pkin(2);
t96 = 0.1e1 / t98;
t193 = t143 * t96;
t81 = t154 * t193;
t183 = t86 * t127;
t188 = t85 * t123;
t153 = t183 / 0.2e1 + t188 / 0.2e1;
t82 = t153 * t193;
t64 = -t112 * t82 + t113 * t81;
t211 = pkin(4) * t64;
t182 = -0.2e1 * pkin(3) * t211 + t140;
t217 = -pkin(8) - pkin(10);
t53 = (pkin(3) - t217) * (pkin(3) + t217) + t182;
t216 = -pkin(8) + pkin(10);
t54 = (pkin(3) - t216) * (pkin(3) + t216) + t182;
t145 = sqrt(-t54 * t53);
t63 = -t112 * t81 - t113 * t82;
t191 = t145 * t63;
t175 = pkin(8) ^ 2 - pkin(10) ^ 2;
t141 = pkin(3) ^ 2;
t59 = t141 + t182;
t55 = t59 + t175;
t60 = -pkin(3) + t211;
t42 = -pkin(4) * t191 - t55 * t60;
t223 = t42 / 0.2e1;
t57 = 0.1e1 / t59;
t222 = t57 / 0.2e1;
t105 = t108 * qJD(2);
t106 = t107 * qJD(2);
t171 = pkin(1) * pkin(5) * t105;
t201 = 0.2e1 * (t92 + t93) * t171 / t146;
t164 = -t201 / 0.2e1;
t151 = t106 * t146 + t108 * t164;
t70 = ((t99 * t224 - t95) * t105 + t151) * pkin(5);
t221 = -t70 / 0.2e1;
t170 = -0.2e1 * t105 * t108;
t180 = t105 * t146;
t71 = t99 * t201 / 0.2e1 + t139 * pkin(1) * t170 + (-t106 * t95 - t180) * pkin(5);
t220 = t71 / 0.2e1;
t120 = cos(pkin(19));
t160 = 0.1e1 / t98 ^ 2 * t171;
t185 = t86 * t120;
t118 = sin(pkin(19));
t186 = t86 * t118;
t189 = t85 * t120;
t190 = t85 * t118;
t209 = t118 / 0.2e1;
t78 = (-t189 / 0.2e1 + t186 / 0.2e1) * t193;
t147 = t78 ^ 2;
t74 = 0.1e1 / t147;
t79 = (t185 / 0.2e1 + t190 / 0.2e1) * t193;
t75 = t79 ^ 2;
t35 = qJD(2) + (((t120 * t220 + t70 * t209) * t96 + (t185 + t190) * t160) / t78 - ((t120 * t221 + t71 * t209) * t96 + (t186 - t189) * t160) * t79 * t74) * t143 / (t74 * t75 + 0.1e1);
t215 = pkin(3) * t35;
t207 = t123 / 0.2e1;
t47 = ((t183 + t188) * t160 + (t154 * qJD(3) + t127 * t220 + t70 * t207) * t96) * t143;
t48 = ((t184 - t187) * t160 + (t153 * qJD(3) + t127 * t221 + t71 * t207) * t96) * t143;
t45 = -t112 * t47 - t113 * t48;
t214 = pkin(3) * t45;
t212 = pkin(4) * t63;
t43 = -t145 * t60 + t55 * t212;
t213 = pkin(4) * t43;
t116 = sin(pkin(17));
t210 = t116 / 0.2e1;
t121 = cos(pkin(18));
t208 = -t121 / 0.2e1;
t130 = cos(pkin(15));
t206 = t130 / 0.2e1;
t132 = qJD(1) ^ 2;
t205 = t132 / 0.2e1;
t204 = cos(qJ(4));
t172 = pkin(4) * t214;
t202 = 0.2e1 * (t53 + t54) * t172 / t145;
t200 = pkin(1) * qJD(2);
t117 = cos(pkin(17));
t199 = t117 * t57;
t126 = sin(pkin(15));
t198 = t126 * t96;
t138 = 0.1e1 / pkin(6);
t194 = t138 * t96;
t100 = pkin(1) * t107 - pkin(5);
t94 = t98 - t174;
t84 = -pkin(1) * t179 - t100 * t94;
t87 = pkin(1) * t108 * t94 - t100 * t146;
t83 = (t87 * t206 - t84 * t126 / 0.2e1) * t194;
t197 = t132 * t83;
t134 = 0.1e1 / pkin(10);
t196 = t134 * t57;
t136 = 0.1e1 / pkin(8);
t195 = t136 * t57;
t192 = t145 * t45;
t157 = t130 * t160;
t158 = t126 * t160;
t162 = t96 * t206;
t69 = ((0.2e1 * pkin(5) * t100 - t94) * t105 + t151) * pkin(1);
t72 = t100 * t164 + t144 * pkin(5) * t170 + (-t106 * t94 - t180) * pkin(1);
t80 = (t84 * t206 + t87 * t126 / 0.2e1) * t194;
t148 = t80 ^ 2;
t76 = 0.1e1 / t148;
t77 = t83 ^ 2;
t36 = ((t72 * t162 + t87 * t157 - t69 * t198 / 0.2e1 - t84 * t158) / t80 - (t69 * t162 + t84 * t157 + t72 * t198 / 0.2e1 + t87 * t158) * t83 * t76) / (t76 * t77 + 0.1e1) * t138;
t181 = qJD(1) * t36;
t178 = t128 * t132;
t131 = qJD(2) ^ 2;
t177 = t131 * t144;
t173 = qJD(1) * qJD(2);
t114 = qJD(2) + qJD(3);
t169 = t78 * t200;
t168 = t79 * t200;
t167 = t123 * t200;
t166 = t127 * t200;
t165 = -t202 / 0.2e1;
t163 = t57 * t210;
t58 = 0.1e1 / t59 ^ 2;
t161 = t58 * t172;
t111 = (-pkin(1) * t128 - pkin(13)) * qJD(1);
t102 = (t123 * t124 - t127 * t128) * qJD(1);
t103 = (-t123 * t128 - t124 * t127) * qJD(1);
t56 = t59 - t175;
t61 = -pkin(3) * t64 + pkin(4);
t155 = -pkin(3) * t191 + t56 * t61;
t156 = pkin(3) * t56 * t63 + t145 * t61;
t27 = (-t155 * t117 / 0.2e1 + t156 * t210) * t196;
t28 = (t156 * t117 / 0.2e1 + t155 * t210) * t196;
t21 = t102 * t27 - t103 * t28;
t110 = pkin(4) * t114 - t166;
t22 = t28 * t110 - t27 * t167;
t44 = t112 * t48 - t113 * t47;
t152 = -t44 * t145 + t63 * t165;
t23 = t110 * t27 + t28 * t167;
t91 = -pkin(4) * t102 + t111;
t150 = t156 * t161;
t149 = t155 * t161;
t122 = sin(qJ(4));
t119 = sin(pkin(18));
t109 = t111 ^ 2 / 0.2e1;
t66 = (-t124 * t79 + t128 * t78) * qJD(1);
t65 = (t124 * t78 + t128 * t79) * qJD(1);
t51 = t111 + (-t119 * t65 - t121 * t66) * pkin(3);
t41 = 0.1e1 / t42 ^ 2;
t34 = t119 * t215 + t168;
t33 = t121 * t215 + t169;
t30 = (t42 * t208 - t119 * t43 / 0.2e1) * t195;
t29 = (t119 * t223 + t43 * t208) * t195;
t26 = 0.1e1 / t27 ^ 2;
t20 = t102 * t28 + t103 * t27;
t18 = -qJD(4) + t21;
t17 = -t29 * t65 + t30 * t66;
t16 = t29 * t66 + t30 * t65;
t15 = t61 * t202 / 0.2e1 - 0.2e1 * t141 * t45 * t212 + (t44 * t56 - t192) * pkin(3);
t14 = ((-0.2e1 * pkin(4) * t61 - t56) * t45 + t152) * pkin(3);
t13 = -t29 * t34 + t30 * t33;
t12 = t29 * t33 + t30 * t34;
t11 = -pkin(9) * t21 - pkin(11) * t20 + t91;
t10 = t35 + (((t60 * t165 + (t44 * t55 - t192) * pkin(4)) * t222 + (-t140 * t57 * t63 + t58 * t213) * t214) / t223 - 0.2e1 * ((-t45 * t55 + t152) * t222 + (t42 * t58 + t57 * t60) * t214) * t41 * t213) * pkin(8) * t136 / (t41 * t43 ^ 2 + 0.1e1) * t59;
t9 = ((t15 * t199 / 0.2e1 + t117 * t150 + t14 * t163 + t116 * t149) / t27 - (-t14 * t199 / 0.2e1 - t117 * t149 + t15 * t163 + t116 * t150) * t28 * t26) / (t26 * t28 ^ 2 + 0.1e1) * t134 + t114;
t7 = -pkin(9) * t9 - t23;
t6 = pkin(11) * t9 + t22;
t5 = t122 * t9 + t204 * t20;
t3 = t122 * t20 - t204 * t9;
t2 = t122 * t11 + t204 * t6;
t1 = t204 * t11 - t122 * t6;
t4 = [0, 0, 0, 0, 0, t205, 0, 0, 0, 0, t124 ^ 2 * t205, t124 * t178, t124 * t173, t128 ^ 2 * t205, t128 * t173, t131 / 0.2e1, pkin(13) * t178, -t132 * pkin(13) * t124, 0, pkin(13) ^ 2 * t205, t103 ^ 2 / 0.2e1, t103 * t102, t114 * t103, t102 ^ 2 / 0.2e1, t114 * t102, t114 ^ 2 / 0.2e1, -t102 * t111 - t114 * t166, t103 * t111 + t114 * t167, (-t102 * t123 + t103 * t127) * t200, t109 + (t123 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1) * t177, t20 ^ 2 / 0.2e1, t21 * t20, t9 * t20, t21 ^ 2 / 0.2e1, t21 * t9, t9 ^ 2 / 0.2e1, -t21 * t91 + t23 * t9, t20 * t91 - t22 * t9, -t20 * t23 + t21 * t22, t22 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t5 * t3, -t5 * t18, t3 ^ 2 / 0.2e1, t18 * t3, t18 ^ 2 / 0.2e1, -t1 * t18 + t3 * t7, t18 * t2 + t5 * t7, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t77 * t205, t80 * t197, t83 * t181, t148 * t205, t80 * t181, t36 ^ 2 / 0.2e1, -t132 * pkin(7) * t80, pkin(7) * t197, 0, pkin(7) ^ 2 * t205, t65 ^ 2 / 0.2e1, t66 * t65, t35 * t65, t66 ^ 2 / 0.2e1, t35 * t66, t35 ^ 2 / 0.2e1, -t111 * t66 + t169 * t35, t111 * t65 - t168 * t35, (-t65 * t78 + t66 * t79) * t200, t109 + (t75 / 0.2e1 + t147 / 0.2e1) * t177, t16 ^ 2 / 0.2e1, t17 * t16, t16 * t10, t17 ^ 2 / 0.2e1, t17 * t10, t10 ^ 2 / 0.2e1, t10 * t13 - t17 * t51, -t10 * t12 + t16 * t51, t12 * t17 - t13 * t16, t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1;];
T_reg = t4;
