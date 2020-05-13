% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh3m1DE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energykin_fixb_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 05:39:24
% EndTime: 2020-04-20 05:40:08
% DurationCPUTime: 42.18s
% Computational Cost: add. (1177682->182), mult. (1816383->446), div. (76944->22), fcn. (1136090->34), ass. (0->201)
t230 = -2 * pkin(1);
t152 = pkin(3) ^ 2;
t151 = pkin(4) ^ 2;
t126 = pkin(18) + pkin(19);
t123 = sin(t126);
t124 = cos(t126);
t134 = sin(qJ(3));
t150 = pkin(5) ^ 2;
t155 = pkin(1) ^ 2;
t135 = sin(qJ(2));
t136 = sin(pkin(16));
t139 = cos(qJ(2));
t140 = cos(pkin(16));
t118 = t135 * t136 - t139 * t140;
t210 = pkin(5) * t118;
t184 = t210 * t230 + t155;
t109 = t150 + t184;
t182 = pkin(2) ^ 2 - pkin(6) ^ 2;
t106 = t109 + t182;
t110 = pkin(1) - t210;
t119 = t135 * t140 + t139 * t136;
t226 = -pkin(6) - pkin(2);
t103 = (pkin(5) - t226) * (pkin(5) + t226) + t184;
t225 = -pkin(6) + pkin(2);
t104 = (pkin(5) - t225) * (pkin(5) + t225) + t184;
t156 = sqrt(-t104 * t103);
t97 = pkin(5) * t119 * t106 + t110 * t156;
t195 = t97 * t134;
t138 = cos(qJ(3));
t187 = t119 * t156;
t96 = -pkin(5) * t187 + t110 * t106;
t198 = t96 * t138;
t161 = -t195 / 0.2e1 + t198 / 0.2e1;
t107 = 0.1e1 / t109;
t154 = 0.1e1 / pkin(2);
t189 = t107 * t154;
t92 = t161 * t189;
t194 = t97 * t138;
t199 = t96 * t134;
t160 = t194 / 0.2e1 + t199 / 0.2e1;
t93 = t160 * t189;
t72 = -t123 * t93 + t124 * t92;
t217 = pkin(4) * t72;
t193 = -0.2e1 * pkin(3) * t217 + t151;
t67 = t152 + t193;
t65 = 0.1e1 / t67;
t229 = t65 / 0.2e1;
t116 = t119 * qJD(2);
t117 = t118 * qJD(2);
t179 = pkin(1) * pkin(5) * t116;
t207 = 0.2e1 / t156 * (t103 + t104) * t179;
t172 = -t207 / 0.2e1;
t158 = t117 * t156 + t119 * t172;
t83 = ((t110 * t230 - t106) * t116 + t158) * pkin(5);
t228 = -t83 / 0.2e1;
t178 = -0.2e1 * t116 * t119;
t188 = t116 * t156;
t84 = t110 * t207 / 0.2e1 + t150 * pkin(1) * t178 + (-t117 * t106 - t188) * pkin(5);
t227 = t84 / 0.2e1;
t224 = -pkin(8) - pkin(10);
t223 = -pkin(8) + pkin(10);
t147 = 0.1e1 / pkin(8);
t170 = t147 * t229;
t61 = (pkin(3) - t224) * (pkin(3) + t224) + t193;
t62 = (pkin(3) - t223) * (pkin(3) + t223) + t193;
t157 = sqrt(-t62 * t61);
t71 = -t123 * t92 - t124 * t93;
t202 = t157 * t71;
t183 = pkin(8) ^ 2 - pkin(10) ^ 2;
t63 = t67 + t183;
t68 = -pkin(3) + t217;
t45 = -pkin(4) * t202 - t68 * t63;
t218 = pkin(4) * t71;
t47 = -t157 * t68 + t63 * t218;
t37 = atan2(t47 * t170, t45 * t170);
t222 = sin(t37);
t131 = cos(pkin(19));
t167 = 0.1e1 / t109 ^ 2 * t179;
t196 = t97 * t131;
t129 = sin(pkin(19));
t197 = t97 * t129;
t200 = t96 * t131;
t201 = t96 * t129;
t215 = t129 / 0.2e1;
t89 = (-t200 / 0.2e1 + t197 / 0.2e1) * t189;
t87 = 0.1e1 / t89 ^ 2;
t90 = (t196 / 0.2e1 + t201 / 0.2e1) * t189;
t42 = qJD(2) + (((t131 * t227 + t83 * t215) * t107 + (t196 + t201) * t167) / t89 - ((t131 * t228 + t84 * t215) * t107 + (t197 - t200) * t167) * t90 * t87) * t154 / (t90 ^ 2 * t87 + 0.1e1);
t221 = pkin(3) * t42;
t214 = t134 / 0.2e1;
t53 = ((t194 + t199) * t167 + (t161 * qJD(3) + t138 * t227 + t83 * t214) * t107) * t154;
t54 = ((t195 - t198) * t167 + (t160 * qJD(3) + t138 * t228 + t84 * t214) * t107) * t154;
t50 = -t123 * t53 - t124 * t54;
t220 = pkin(3) * t50;
t219 = pkin(4) * t47;
t127 = sin(pkin(17));
t216 = t127 / 0.2e1;
t141 = cos(pkin(15));
t213 = t141 / 0.2e1;
t143 = qJD(1) ^ 2;
t212 = t143 / 0.2e1;
t211 = cos(qJ(4));
t180 = pkin(4) * t220;
t209 = 0.2e1 * (t61 + t62) * t180 / t157;
t208 = pkin(1) * qJD(2);
t128 = cos(pkin(17));
t206 = t128 * t65;
t137 = sin(pkin(15));
t149 = 0.1e1 / pkin(6);
t190 = t107 * t149;
t105 = t109 - t182;
t111 = pkin(1) * t118 - pkin(5);
t95 = -pkin(1) * t187 - t111 * t105;
t98 = pkin(1) * t119 * t105 - t111 * t156;
t91 = (t95 * t213 + t98 * t137 / 0.2e1) * t190;
t94 = (t98 * t213 - t95 * t137 / 0.2e1) * t190;
t81 = atan2(t94, t91);
t78 = cos(t81);
t205 = t143 * t78;
t145 = 0.1e1 / pkin(10);
t204 = t145 * t65;
t203 = t157 * t50;
t162 = t141 * t167;
t163 = t137 * t167;
t168 = t107 * t213;
t191 = t107 * t137;
t82 = ((0.2e1 * t111 * pkin(5) - t105) * t116 + t158) * pkin(1);
t85 = t111 * t172 + t155 * pkin(5) * t178 + (-t117 * t105 - t188) * pkin(1);
t88 = 0.1e1 / t91 ^ 2;
t43 = ((t85 * t168 + t98 * t162 - t82 * t191 / 0.2e1 - t95 * t163) / t91 - (t82 * t168 + t95 * t162 + t85 * t191 / 0.2e1 + t98 * t163) * t94 * t88) / (t94 ^ 2 * t88 + 0.1e1) * t149;
t192 = qJD(1) * t43;
t186 = t139 * t143;
t142 = qJD(2) ^ 2;
t185 = t142 * t155;
t181 = qJD(1) * qJD(2);
t125 = qJD(2) + qJD(3);
t80 = atan2(t90, t89);
t74 = sin(t80);
t177 = t74 * t208;
t75 = cos(t80);
t176 = t75 * t208;
t175 = t134 * t208;
t174 = t138 * t208;
t173 = -t209 / 0.2e1;
t171 = t65 * t216;
t66 = 0.1e1 / t67 ^ 2;
t169 = t66 * t180;
t113 = (t134 * t135 - t138 * t139) * qJD(1);
t114 = (-t134 * t139 - t135 * t138) * qJD(1);
t64 = t67 - t183;
t69 = -pkin(3) * t72 + pkin(4);
t46 = -pkin(3) * t202 + t69 * t64;
t48 = pkin(3) * t71 * t64 + t157 * t69;
t34 = (-t46 * t128 / 0.2e1 + t48 * t216) * t204;
t35 = (t48 * t128 / 0.2e1 + t46 * t216) * t204;
t31 = atan2(t35, t34);
t28 = sin(t31);
t29 = cos(t31);
t17 = -t29 * t113 + t28 * t114;
t165 = t127 * t169;
t164 = t128 * t169;
t122 = (-pkin(1) * t139 - pkin(13)) * qJD(1);
t121 = t125 * pkin(4) - t174;
t21 = t28 * t121 - t29 * t175;
t49 = t123 * t54 - t124 * t53;
t159 = -t49 * t157 + t71 * t173;
t20 = t29 * t121 + t28 * t175;
t102 = -t113 * pkin(4) + t122;
t133 = sin(qJ(4));
t132 = cos(pkin(18));
t130 = sin(pkin(18));
t120 = t122 ^ 2 / 0.2e1;
t77 = sin(t81);
t59 = (t135 * t75 + t139 * t74) * qJD(1);
t57 = (t135 * t74 - t139 * t75) * qJD(1);
t51 = t122 + (-t130 * t59 + t132 * t57) * pkin(3);
t44 = 0.1e1 / t45 ^ 2;
t41 = t132 * t221 + t176;
t40 = t130 * t221 + t177;
t36 = cos(t37);
t33 = 0.1e1 / t34 ^ 2;
t27 = -t130 * t222 - t132 * t36;
t26 = t130 * t36 - t132 * t222;
t23 = t69 * t209 / 0.2e1 - 0.2e1 * t152 * t50 * t218 + (t49 * t64 - t203) * pkin(3);
t22 = ((-0.2e1 * t69 * pkin(4) - t64) * t50 + t159) * pkin(3);
t19 = t28 * t113 + t29 * t114;
t16 = qJD(4) + t17;
t15 = -t26 * t57 + t27 * t59;
t14 = -t26 * t59 - t27 * t57;
t13 = t26 * t41 + t27 * t40;
t12 = -t26 * t40 + t27 * t41;
t11 = t17 * pkin(9) - t19 * pkin(11) + t102;
t10 = t42 + 0.2e1 * (((t68 * t173 + (t49 * t63 - t203) * pkin(4)) * t229 + (-t151 * t65 * t71 + t66 * t219) * t220) / t45 - ((-t50 * t63 + t159) * t229 + (t45 * t66 + t65 * t68) * t220) * t44 * t219) * pkin(8) * t147 / (t44 * t47 ^ 2 + 0.1e1) * t67;
t9 = ((t23 * t206 / 0.2e1 + t48 * t164 + t22 * t171 + t46 * t165) / t34 - (-t22 * t206 / 0.2e1 - t46 * t164 + t23 * t171 + t48 * t165) * t35 * t33) / (t33 * t35 ^ 2 + 0.1e1) * t145 + t125;
t7 = t9 * pkin(11) + t21;
t6 = -t9 * pkin(9) - t20;
t5 = t133 * t9 + t211 * t19;
t3 = t133 * t19 - t211 * t9;
t2 = t133 * t11 + t211 * t7;
t1 = t211 * t11 - t133 * t7;
t4 = [0, 0, 0, 0, 0, t212, 0, 0, 0, 0, t135 ^ 2 * t212, t135 * t186, t135 * t181, t139 ^ 2 * t212, t139 * t181, t142 / 0.2e1, pkin(13) * t186, -t143 * pkin(13) * t135, 0, pkin(13) ^ 2 * t212, t114 ^ 2 / 0.2e1, t113 * t114, t125 * t114, t113 ^ 2 / 0.2e1, t113 * t125, t125 ^ 2 / 0.2e1, -t122 * t113 - t125 * t174, t122 * t114 + t125 * t175, (-t113 * t134 + t114 * t138) * t208, t120 + (t134 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1) * t185, t19 ^ 2 / 0.2e1, -t19 * t17, t9 * t19, t17 ^ 2 / 0.2e1, -t9 * t17, t9 ^ 2 / 0.2e1, t102 * t17 + t20 * t9, t102 * t19 - t21 * t9, -t17 * t21 - t19 * t20, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t3 * t5, t16 * t5, t3 ^ 2 / 0.2e1, -t16 * t3, t16 ^ 2 / 0.2e1, t1 * t16 + t3 * t6, -t16 * t2 + t5 * t6, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t77 ^ 2 * t212, t77 * t205, t77 * t192, t78 ^ 2 * t212, t78 * t192, t43 ^ 2 / 0.2e1, -pkin(7) * t205, t143 * pkin(7) * t77, 0, pkin(7) ^ 2 * t212, t59 ^ 2 / 0.2e1, -t57 * t59, t42 * t59, t57 ^ 2 / 0.2e1, -t42 * t57, t42 ^ 2 / 0.2e1, t122 * t57 + t176 * t42, t122 * t59 - t177 * t42, (-t57 * t74 - t59 * t75) * t208, t120 + (t74 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1) * t185, t15 ^ 2 / 0.2e1, t14 * t15, t10 * t15, t14 ^ 2 / 0.2e1, t10 * t14, t10 ^ 2 / 0.2e1, t10 * t12 - t14 * t51, -t10 * t13 + t15 * t51, -t12 * t15 + t13 * t14, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1;];
T_reg = t4;
