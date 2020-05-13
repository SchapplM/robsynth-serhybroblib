% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m2DE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energykin_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:06:20
% EndTime: 2020-05-02 21:06:22
% DurationCPUTime: 1.90s
% Computational Cost: add. (632->223), mult. (1611->429), div. (0->0), fcn. (1141->22), ass. (0->190)
t118 = cos(pkin(18));
t95 = t118 ^ 2;
t164 = t95 / 0.2e1 - 0.1e1 / 0.4e1;
t100 = sin(pkin(20));
t104 = cos(pkin(20));
t113 = sin(pkin(18));
t36 = -t118 * t100 + t113 * t104;
t87 = pkin(22) + pkin(21);
t74 = cos(t87);
t210 = t74 * t36;
t41 = t113 * t100 + t118 * t104;
t73 = sin(t87);
t211 = t73 * t41;
t12 = t210 - t211;
t231 = t12 ^ 2;
t122 = qJD(1) ^ 2;
t112 = sin(qJ(2));
t111 = sin(qJ(3));
t213 = pkin(5) * t111;
t67 = pkin(1) + t213;
t195 = t67 * t112;
t116 = cos(qJ(3));
t117 = cos(qJ(2));
t179 = t117 * t116;
t59 = pkin(5) * t179;
t208 = t59 + pkin(15);
t31 = -t195 + t208;
t201 = t122 * t31;
t230 = t201 * t73;
t191 = t100 * t104;
t72 = t95 - 0.1e1 / 0.2e1;
t144 = t72 * pkin(11);
t177 = t118 * t113;
t34 = -pkin(9) * t177 + t144;
t160 = pkin(11) * t177;
t35 = t72 * pkin(9) + t160;
t89 = t104 ^ 2;
t138 = -t34 * t191 + t35 * t89;
t66 = t74 ^ 2;
t199 = t122 * t66;
t50 = -t113 * pkin(9) + t118 * pkin(11);
t51 = t118 * pkin(9) + t113 * pkin(11);
t26 = -t100 * t50 + t51 * t104;
t157 = t35 * t191;
t56 = t177 / 0.2e1;
t158 = pkin(9) * t56 - pkin(11) * t164;
t4 = t34 * t89 + t157 + t158;
t140 = t164 * pkin(9);
t57 = -t177 / 0.2e1;
t5 = pkin(11) * t57 + t138 - t140;
t229 = (-qJD(1) * qJD(4) * t26 + (t31 * t41 + 0.4e1 * t4 * t73) * t122) * t74 - 0.4e1 * t5 * t199 + t36 * t230 + (t31 * qJD(4) + (-t118 * t51 + 0.2e1 * t138) * qJD(1)) * qJD(1);
t197 = t36 * t116;
t13 = -pkin(5) * t197 + t67 * t41;
t135 = t111 * t41 - t197;
t196 = t41 * t116;
t14 = pkin(5) * t196 + t67 * t36;
t20 = t111 * t36 + t196;
t207 = pkin(5) * qJD(3);
t228 = qJD(4) * ((t14 * t112 + t13 * t117) * qJD(2) + (t20 * t112 + t135 * t117) * t207) * t74;
t83 = t112 * pkin(15);
t227 = t83 * t122;
t226 = qJD(4) * t73;
t110 = sin(qJ(4));
t171 = qJD(1) * t110;
t115 = cos(qJ(4));
t170 = qJD(1) * t115;
t225 = -t89 * t177 + t72 * t191;
t221 = 0.2e1 * pkin(15);
t220 = -0.4e1 * t122;
t219 = t31 / 0.4e1;
t91 = qJD(2) + qJD(3);
t71 = t91 ^ 2 / 0.2e1;
t99 = sin(pkin(21));
t218 = pkin(4) * t99;
t92 = qJD(2) ^ 2 / 0.2e1;
t217 = t122 / 0.2e1;
t216 = t83 - pkin(1);
t215 = pkin(2) * t112;
t103 = cos(pkin(21));
t214 = pkin(4) * t103;
t84 = t112 * pkin(1);
t212 = t122 * t231;
t124 = pkin(11) ^ 2;
t125 = pkin(9) ^ 2;
t86 = -t124 + t125;
t58 = t86 * t95;
t142 = pkin(9) * t160;
t209 = 0.2e1 * t142 + t58;
t64 = t84 - pkin(15);
t206 = pkin(14) * t122;
t205 = pkin(15) * t122;
t102 = cos(pkin(22));
t98 = sin(pkin(22));
t203 = t102 * t98;
t202 = t118 * t98;
t200 = t122 * t64;
t97 = qJ(2) + qJ(3);
t79 = cos(t97);
t198 = t122 * t79;
t126 = pkin(5) ^ 2;
t93 = t116 ^ 2;
t194 = t93 * t126;
t193 = qJD(1) * t91;
t192 = qJD(2) * t91;
t101 = sin(pkin(19));
t105 = cos(pkin(19));
t190 = t105 * t101;
t189 = t111 * t105;
t188 = t112 * t111;
t187 = t112 * t116;
t186 = t112 * t117;
t185 = t113 * t102;
t114 = sin(pkin(17));
t184 = t114 * t113;
t119 = cos(pkin(17));
t183 = t114 * t119;
t182 = t116 * t101;
t181 = t116 * t111;
t180 = t117 * t111;
t178 = t117 * t122;
t176 = t119 * t118;
t175 = t125 / 0.4e1 - t124 / 0.4e1;
t75 = -t125 / 0.2e1 + t124 / 0.2e1;
t129 = pkin(1) ^ 2;
t174 = t126 + t129;
t172 = t93 - 0.1e1 / 0.2e1;
t169 = qJD(1) * qJD(2);
t168 = pkin(1) * t213;
t167 = t118 * t218;
t165 = t73 * (t177 * t191 + t72 * t89 - t164) * t74;
t69 = 0.2e1 * t168;
t163 = (t69 + t174) * t92 + (pkin(1) * t111 + pkin(5)) * qJD(2) * t207 + t126 * qJD(3) ^ 2 / 0.2e1;
t162 = pkin(1) * t192;
t161 = pkin(4) * t185;
t156 = pkin(15) * t178;
t29 = pkin(9) * t144 + t75 * t177;
t155 = t29 * t191;
t88 = t102 ^ 2;
t154 = t88 * t177;
t43 = t180 + t187;
t44 = t179 - t188;
t15 = t44 * t101 + t43 * t105;
t152 = t15 * t193;
t151 = t93 * t190;
t150 = t212 / 0.2e1;
t149 = -t211 / 0.4e1;
t148 = t105 * t188;
t147 = t112 * t178;
t145 = -t177 / 0.4e1;
t143 = t112 * t169;
t63 = t117 * t169;
t141 = t98 * t167;
t139 = pkin(1) * t63;
t123 = pkin(15) ^ 2;
t94 = t117 ^ 2;
t136 = -t94 * t129 + t123 + t129;
t10 = t36 * t73 + t41 * t74;
t42 = t113 * t119 - t118 * t114;
t45 = t176 + t184;
t21 = t112 * t45 - t42 * t117;
t134 = t112 * t42 + t45 * t117;
t132 = t123 + t174 + (-0.2e1 * t168 - t174 + 0.2e1 * t194) * t94 - t194;
t38 = t102 * t118 + t113 * t98;
t131 = t99 * t161 + t38 * t214 - pkin(15) - t141;
t90 = t105 ^ 2;
t68 = t90 - 0.1e1 / 0.2e1;
t130 = -t68 * t93 + t181 * t190 - 0.1e1 / 0.4e1;
t96 = t119 ^ 2;
t82 = t129 * t92;
t80 = t123 * t217;
t78 = sin(t97);
t76 = t90 / 0.2e1;
t65 = -0.2e1 * pkin(1) * t83;
t62 = t94 * t217;
t61 = t112 ^ 2 * t217;
t53 = t111 * t162;
t52 = t116 * t162;
t47 = -pkin(1) * t94 - t216;
t39 = -t185 + t202;
t37 = -t111 * t101 + t116 * t105;
t28 = t75 + t209;
t27 = (t136 + t65) * t217 + t82;
t25 = t100 * t51 + t50 * t104;
t19 = (pkin(5) * t187 + t117 * t67) * qJD(2) + t43 * t207;
t16 = -t43 * t101 + t44 * t105;
t8 = t10 * qJD(1) + qJD(4);
t7 = t131 + t84;
t3 = (t210 / 0.4e1 + t149) * qJD(4) + (-t165 + (t66 - 0.1e1 / 0.2e1) * (-t225 + t57)) * qJD(1);
t1 = (-t13 * t112 + t14 * t117) * qJD(2) + (-t112 * t135 + t20 * t117) * t207;
t2 = [0, 0, 0, 0, 0, t217, 0, 0, 0, 0, t62, -t147, t63, t61, -t143, t92, -t227, -t156, 0, t80, t43 ^ 2 * t217, 0.2e1 * t122 * (t172 * t186 + (t94 - 0.1e1 / 0.2e1) * t181), t43 * t193, t44 ^ 2 * t217, t44 * t193, t71, (t47 * t111 - t64 * t179) * t122 + t53, (t47 * t116 + t64 * t180) * t122 + t52, -t139, t27, t231 * t217, ((t225 + t56) * t66 + t165 + t89 * t56 - t164 * t191 + t145) * t220, 0, t10 ^ 2 * t217, 0, 0, -t10 * t201, t12 * t201, -t19 * qJD(1), (-0.2e1 * (pkin(5) * t188 + t64) * t59 - 0.2e1 * t216 * t213 + t65 + t132) * t217 + t163, t115 ^ 2 * t150, -t110 * t115 * t212, -0.4e1 * t3 * t170, t110 ^ 2 * t150, 0.4e1 * t3 * t171, t8 ^ 2 / 0.2e1, -t110 * t228 + (-t1 * t110 - t25 * t170) * t226 - t19 * t171 - t229 * t115, -t115 * t228 + (-t1 * t115 + t25 * t171) * t226 - t19 * t170 + t229 * t110, (t4 * t66 + (t36 * t219 + t5 * t73) * t74 + t31 * t149 + t158 * t89 - t157 / 0.2e1 + pkin(9) * t145 + (-0.1e1 / 0.4e1 + t95 / 0.4e1) * pkin(11)) * t220, 0.2e1 * (t28 * t89 + t75 * t95 - t142 - 0.2e1 * t155 + t175) * t199 + ((t29 * t89 + t28 * t191 / 0.2e1 + t175 * t177 - pkin(11) * t140) * t73 + t26 * t219) * t74 * t220 + t25 * t230 + (-0.2e1 * t208 * t195 + t59 * t221 + t69 + (-0.4e1 * t142 + t86 - 0.2e1 * t58) * t89 + 0.4e1 * t155 + t124 + t132 + t209) * t217 + t163, t134 ^ 2 * t217, 0.4e1 * t122 * ((t96 * t177 - t72 * t183 + t57) * t94 - (t176 * t184 + t72 * t96 - t164) * t186 + t96 * t57 + t164 * t183 + t177 / 0.4e1), t134 * t169, t21 ^ 2 * t217, -t21 * t169, t92, t21 * t206, t134 * t206, 0, pkin(14) ^ 2 * t217, t39 ^ 2 * t217, -0.2e1 * t122 * (t72 * t203 - t154 + t56), 0, t38 ^ 2 * t217, 0, 0, t38 * t200, t39 * t200, -t139, t27, t15 ^ 2 * t217, ((-t172 * t190 - t68 * t181) * t94 + (t76 + t130) * t186 + t151 / 0.2e1 + (t76 - 0.1e1 / 0.4e1) * t181 - t190 / 0.4e1) * t220, t152, t16 ^ 2 * t217, t16 * t193, t71, t16 * t205, -t15 * t205, 0, t80, t62, -t147, t63, t61, -t143, t92, -t227 + ((-(t94 * t101 + t105 * t186 - t101) * t116 - (-t101 * t186 + t94 * t105 - t105) * t111) * t122 + (t182 + t189) * t192) * pkin(2), -t156 + (-(-t112 * t182 + t117 * t37 - t148) * t178 - t37 * t192) * pkin(2), -pkin(2) * t152, (-0.4e1 * pkin(2) * (t151 * t215 + (pkin(2) * t68 * t188 - pkin(15) * t105 / 0.2e1) * t116 - t101 * (-t111 * pkin(15) + t105 * t215) / 0.2e1) * t117 + t123 - 0.2e1 * ((-pkin(2) * t189 + t83) * t182 + pkin(15) * t148) * pkin(2)) * t217 + ((0.4e1 * (-t90 / 0.2e1 - t130) * t94 + (-0.2e1 * t90 + 0.1e1) * t93 + t90) * t217 + t71) * pkin(2) ^ 2, t78 ^ 2 * t217, t78 * t198, t78 * t193, t79 ^ 2 * t217, t79 * t193, t71, -t198 * t7 + t53, t7 * t122 * t78 + t52, -t139, t82 + (0.2e1 * t131 * t84 - 0.4e1 * (-t154 * t218 + (pkin(15) * t118 / 0.2e1 + t72 * t98 * t218) * t102 + t113 * (t98 * pkin(15) + t167) / 0.2e1) * t214 - 0.2e1 * (pkin(4) * t202 + t99 * pkin(15)) * t161 + t141 * t221 + t136 + (0.4e1 * (t177 * t203 + t72 * t88 - t164) * t103 ^ 2 + (-0.2e1 * t95 + 0.1e1) * t88 + t95) * pkin(4) ^ 2) * t217;];
T_reg = t2;
