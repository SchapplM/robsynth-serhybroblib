% Calculate minimal parameter regressor of coriolis matrix for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% cmat_reg [(2*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1turnDE1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:26
% EndTime: 2020-06-27 16:35:44
% DurationCPUTime: 3.37s
% Computational Cost: add. (49115->151), mult. (69258->371), div. (2823->21), fcn. (19133->8), ass. (0->178)
t88 = cos(qJ(2));
t201 = pkin(2) * t88;
t97 = pkin(1) ^ 2;
t168 = -0.2e1 * pkin(1) * t201 + t97;
t209 = -pkin(3) - pkin(4);
t74 = (pkin(2) - t209) * (pkin(2) + t209) + t168;
t208 = -pkin(3) + pkin(4);
t75 = (pkin(2) - t208) * (pkin(2) + t208) + t168;
t179 = t74 * t75;
t100 = sqrt(-t179);
t87 = sin(qJ(2));
t165 = t100 * t87;
t66 = pkin(2) * t165;
t218 = pkin(3) ^ 2;
t219 = pkin(4) ^ 2;
t167 = t218 - t219;
t96 = pkin(2) ^ 2;
t81 = t96 + t168;
t77 = t81 - t167;
t83 = pkin(1) - t201;
t59 = t77 * t83 - t66;
t54 = 0.1e1 / t59 ^ 2;
t202 = pkin(2) * t87;
t70 = t77 * t202;
t60 = t100 * t83 + t70;
t187 = t54 * t60;
t206 = pkin(1) * t87;
t79 = 0.1e1 / t81 ^ 2;
t150 = t79 * t206;
t133 = pkin(2) * t150;
t116 = t60 * t133;
t86 = t87 ^ 2;
t146 = t96 * t86 * pkin(1);
t169 = t77 * t201 + t66;
t68 = 0.1e1 / t100;
t181 = t68 * t83;
t211 = pkin(1) * pkin(2);
t125 = (-t74 - t75) * t211;
t64 = t87 * t125;
t34 = t64 * t181 + 0.2e1 * t146 + t169;
t78 = 0.1e1 / t81;
t189 = t34 * t78;
t91 = 0.1e1 / pkin(4);
t27 = (t189 / 0.2e1 - t116) * t91;
t53 = 0.1e1 / t59;
t193 = t27 * t53;
t117 = t59 * t133;
t161 = t88 * t100;
t182 = t68 * t64;
t32 = t70 + (-t161 + (0.2e1 * t83 * pkin(1) - t182) * t87) * pkin(2);
t190 = t32 * t78;
t25 = (-t190 / 0.2e1 + t117) * t91;
t221 = -0.2e1 * t25 * t187 - 0.2e1 * t193;
t52 = t59 ^ 2;
t191 = t32 * t53 / t52;
t56 = t60 ^ 2;
t46 = t54 * t56 + 0.1e1;
t217 = (t34 * t187 - t56 * t191) / t46 ^ 2;
t65 = pkin(1) * t165;
t76 = t81 + t167;
t205 = pkin(1) * t88;
t84 = -pkin(2) + t205;
t58 = -t76 * t84 - t65;
t50 = 0.1e1 / t58 ^ 2;
t72 = t76 * t206;
t61 = -t100 * t84 + t72;
t188 = t50 * t61;
t214 = -0.2e1 * t84;
t31 = t72 + (-t161 + (pkin(2) * t214 - t182) * t87) * pkin(1);
t48 = t58 ^ 2;
t49 = 0.1e1 / t58;
t192 = t31 * t49 / t48;
t175 = t97 * t86;
t145 = pkin(2) * t175;
t170 = t76 * t205 + t65;
t180 = t68 * t84;
t33 = -t64 * t180 + 0.2e1 * t145 + t170;
t57 = t61 ^ 2;
t47 = t50 * t57 + 0.1e1;
t216 = (t33 * t188 - t57 * t192) / t47 ^ 2;
t212 = 0.2e1 * t87;
t215 = t182 * t212 + t161;
t171 = t52 + t56;
t92 = 0.1e1 / t219;
t40 = t171 * t92 * t79;
t39 = 0.1e1 / t40;
t35 = t40 ^ (-0.1e1 / 0.2e1);
t172 = t48 + t57;
t95 = 0.1e1 / t218;
t41 = t172 * t95 * t79;
t37 = t41 ^ (-0.1e1 / 0.2e1);
t213 = -0.4e1 * t87;
t210 = t78 / 0.2e1;
t207 = pkin(1) * t78;
t204 = pkin(2) * t78;
t203 = pkin(2) * t79;
t200 = pkin(3) * t81;
t199 = pkin(4) * t81;
t94 = 0.1e1 / pkin(3);
t141 = t37 * t78 * t94;
t177 = t87 * t58;
t183 = t61 * t88;
t22 = (-t177 - t183) * t141;
t151 = 0.2e1 * t79;
t121 = t151 * t211;
t152 = pkin(1) * t202;
t80 = t78 * t79;
t132 = t80 * t152;
t112 = 0.4e1 * t132;
t194 = ((t31 * t58 + t33 * t61) * t151 - t172 * t112) * t95 * t37 / t41;
t143 = t78 * t194;
t173 = t33 + t58;
t174 = -t31 + t61;
t176 = t87 * t88;
t6 = ((t177 / 0.2e1 + t183 / 0.2e1) * t143 + ((t61 * t176 + t58 * t86) * t121 + (-t173 * t88 + t174 * t87) * t78) * t37) * t94;
t198 = t22 * t6;
t184 = t61 * t87;
t186 = t58 * t88;
t23 = (t184 - t186) * t141;
t7 = ((-t184 / 0.2e1 + t186 / 0.2e1) * t143 + ((t58 * t176 - t61 * t86) * t121 + (t173 * t87 + t174 * t88) * t78) * t37) * t94;
t197 = t23 * t7;
t140 = 0.4e1 * t68 / t179 * t64 ^ 2;
t123 = -t140 / 0.4e1;
t139 = t96 * t175;
t62 = t88 * t125 - 0.4e1 * t139;
t109 = (-0.2e1 * t88 * t182 + (-t62 * t68 + t123) * t87) * t78;
t44 = 0.1e1 / t47;
t129 = t44 * t50 * t200;
t114 = t61 * t129;
t126 = t80 * t139;
t115 = 0.8e1 * t126;
t120 = pkin(3) * t44 * t152;
t149 = t49 * t200;
t130 = t44 * t149;
t144 = t97 * t202;
t178 = t78 * t88;
t122 = -0.2e1 * t133;
t26 = (t58 * t122 + t31 * t78) * t94;
t28 = (t61 * t122 + t33 * t78) * t94;
t3 = (((t84 * t123 + 0.6e1 * t88 * t144 - t62 * t180 - t72) * t78 + t61 * t115 + (t215 * t78 + (t33 * t213 - 0.2e1 * t183) * t203) * pkin(1)) * t130 - ((0.4e1 * t145 + t170) * t78 + t58 * t115 + (t109 + (t178 * t214 + (t31 * t213 - 0.2e1 * t186) * t79) * pkin(2)) * pkin(1)) * t114) * t94 + (0.2e1 * t49 * t120 - t31 * t129 - 0.2e1 * t149 * t216) * t28 + (-0.2e1 * t120 * t188 - t33 * t129 + 0.2e1 * (t44 * t192 + t50 * t216) * t61 * t200) * t26;
t196 = t37 * t3;
t17 = ((t32 * t59 + t34 * t60) * t151 - t171 * t112) * t92;
t195 = t17 * t35 * t39;
t185 = t60 * t79;
t1 = t22 * t7 + t23 * t6;
t166 = t1 * qJD(1);
t4 = (t22 * t87 - t6 * t88) * pkin(2);
t164 = t4 * qJD(1);
t5 = (-t23 * t87 + t7 * t88) * pkin(2);
t163 = t5 * qJD(1);
t138 = t17 * t79 / t40 ^ 2;
t8 = (-t60 * t59 * t138 + (t32 * t185 + (-0.4e1 * t60 * t132 + t34 * t79) * t59) * t39) * t92;
t162 = t8 * qJD(1);
t9 = (-t56 * t138 / 0.2e1 + (-0.2e1 * t56 * t132 + t34 * t185) * t39) * t92;
t160 = t9 * qJD(1);
t159 = qJD(2) * t91;
t158 = qJD(2) * t94;
t111 = -0.2e1 * t79 * t35 * t144;
t134 = -t195 / 0.2e1;
t10 = (t59 * t111 + (t59 * t134 + t32 * t35) * t207) * t91;
t157 = t10 * qJD(1);
t124 = t60 * t134;
t11 = (t60 * t111 + (t34 * t35 + t124) * t207) * t91;
t156 = t11 * qJD(1);
t82 = t88 ^ 2 - t86;
t155 = t82 * qJD(1);
t154 = t88 * qJD(2);
t42 = 0.1e1 / t46;
t148 = t42 * t199;
t2 = 0.2e1 * pkin(4) * t42 * t152 * t221 + 0.4e1 * (t193 * t217 + (t42 * t191 + t54 * t217) * t25 * t60) * t199 + 0.2e1 * ((-t25 * t34 + t27 * t32) * t54 + (-((t83 * t140 / 0.4e1 + t62 * t181 - t70 + t215 * pkin(2)) * t210 + 0.4e1 * t60 * t126 + (0.3e1 * t96 * t78 * t176 + (-0.2e1 * t34 * t87 - t60 * t88) * t203) * pkin(1)) * t53 - (-(0.4e1 * t146 + t169) * t78 / 0.2e1 - 0.4e1 * t59 * t126 + (-t109 / 0.2e1 + (-t83 * t178 + (t32 * t212 + t59 * t88) * t79) * pkin(1)) * pkin(2)) * t187) * t91) * t148;
t147 = t78 * t35 * t2;
t137 = qJD(1) * t198;
t136 = qJD(1) * t197;
t135 = qJD(1) * t176;
t13 = -t26 * t114 + t28 * t130 + 0.1e1;
t110 = t37 * t96 * t13 * t150;
t12 = t148 * t221;
t14 = [0, 0, 0, t87 * t154, t82 * qJD(2), 0, 0, 0, 0, 0, qJD(2) * t198, t1 * qJD(2), 0, qJD(2) * t197, 0, 0, t5 * qJD(2), t4 * qJD(2), t9 * qJD(2), -t8 * qJD(2), 0, 0, 0, -t10 * qJD(2), -t11 * qJD(2); 0, 0, 0, t135, t155, t154, -t87 * qJD(2), 0, 0, 0, t137, t166, (t13 * t6 + t22 * t3) * qJD(2), t136, (t13 * t7 + t23 * t3) * qJD(2), 0, t163, t164, t160, -t162, (t60 * t147 + (t78 * t124 + (-0.2e1 * t116 + t189) * t35) * t12) * t159, (-t59 * t147 + (t59 * t195 * t210 + (0.2e1 * t117 - t190) * t35) * t12) * t159, 0, -t157, -t156; 0, 0, 0, -t135, -t155, 0, 0, 0, 0, 0, -t137, -t166, 0, -t136, 0, 0, -t163, -t164, -t160, t162, 0, 0, 0, t157, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t3 * qJD(2), (0.2e1 * t58 * t110 + (-t58 * t196 + (-t31 * t37 + t58 * t194 / 0.2e1) * t13) * t204) * t158, (-0.2e1 * t61 * t110 + (t61 * t196 + (t33 * t37 - t61 * t194 / 0.2e1) * t13) * t204) * t158, 0, 0, 0, 0, t12 * t2 * qJD(2), 0, 0;];
cmat_reg = t14;
