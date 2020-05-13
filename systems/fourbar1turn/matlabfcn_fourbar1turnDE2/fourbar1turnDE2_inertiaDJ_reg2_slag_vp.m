% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% fourbar1turnDE2
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
% MMD_reg [((2+1)*2/2)x(2*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbar1turnDE2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:34:57
% EndTime: 2020-04-12 19:35:11
% DurationCPUTime: 3.63s
% Computational Cost: add. (39004->167), mult. (55784->431), div. (2060->22), fcn. (15023->8), ass. (0->178)
t86 = cos(qJ(2));
t202 = pkin(2) * t86;
t96 = pkin(1) ^ 2;
t167 = -0.2e1 * pkin(1) * t202 + t96;
t95 = pkin(2) ^ 2;
t80 = t95 + t167;
t78 = 0.1e1 / t80 ^ 2;
t85 = sin(qJ(2));
t176 = t78 * t85;
t158 = pkin(1) * t176;
t142 = pkin(2) * t158;
t77 = 0.1e1 / t80;
t210 = -t77 / 0.2e1;
t81 = pkin(1) - t202;
t161 = 0.2e1 * t81 * pkin(1);
t208 = -pkin(3) - pkin(4);
t73 = (pkin(2) - t208) * (pkin(2) + t208) + t167;
t207 = pkin(4) - pkin(3);
t74 = (pkin(2) - t207) * (pkin(2) + t207) + t167;
t178 = t73 * t74;
t99 = sqrt(-t178);
t172 = t86 * t99;
t130 = pkin(1) * pkin(2) * (-t73 - t74);
t63 = t85 * t130;
t67 = 0.1e1 / t99;
t181 = t67 * t63;
t221 = pkin(3) ^ 2;
t222 = pkin(4) ^ 2;
t166 = t221 - t222;
t76 = t80 - t166;
t204 = pkin(2) * t76;
t71 = t85 * t204;
t28 = t71 + (-t172 + (t161 - t181) * t85) * pkin(2);
t173 = t85 * t99;
t58 = -pkin(2) * t173 + t76 * t81;
t89 = 0.1e1 / pkin(4);
t16 = (t58 * t142 + t28 * t210) * t89;
t53 = 0.1e1 / t58 ^ 2;
t59 = t81 * t99 + t71;
t188 = t53 * t59;
t209 = t77 / 0.2e1;
t84 = t85 ^ 2;
t153 = t95 * t84 * pkin(1);
t180 = t67 * t81;
t30 = t63 * t180 + 0.2e1 * t153 + (t76 * t86 + t173) * pkin(2);
t18 = (-t59 * t142 + t30 * t209) * t89;
t52 = 0.1e1 / t58;
t194 = t18 * t52;
t223 = -0.2e1 * t16 * t188 - 0.2e1 * t194;
t220 = 0.2e1 * t58;
t159 = 0.2e1 * t78;
t75 = t80 + t166;
t82 = pkin(1) * t86 - pkin(2);
t57 = -pkin(1) * t173 - t75 * t82;
t49 = 0.1e1 / t57 ^ 2;
t206 = pkin(1) * t75;
t72 = t85 * t206;
t60 = -t82 * t99 + t72;
t189 = t49 * t60;
t62 = qJD(2) * t63;
t182 = t67 * t62;
t148 = t85 * t182;
t211 = -0.2e1 * t82;
t163 = pkin(2) * t211;
t23 = (-t148 + (-t172 + (t75 + t163) * t85) * qJD(2)) * pkin(1);
t47 = t57 ^ 2;
t48 = 0.1e1 / t57;
t193 = t23 * t48 / t47;
t175 = t84 * t96;
t146 = qJD(2) * t175;
t131 = pkin(2) * t146;
t165 = qJD(2) * t85;
t145 = t99 * t165;
t164 = qJD(2) * t86;
t168 = pkin(1) * t145 + t164 * t206;
t179 = t67 * t82;
t26 = -t62 * t179 + 0.2e1 * t131 + t168;
t56 = t60 ^ 2;
t46 = t49 * t56 + 0.1e1;
t219 = (t26 * t189 - t56 * t193) / t46 ^ 2;
t171 = t47 + t56;
t217 = t171 * t78;
t93 = 0.1e1 / t221;
t40 = t93 * t217;
t33 = t40 ^ (-0.1e1 / 0.2e1);
t92 = 0.1e1 / pkin(3);
t218 = t33 * t92;
t147 = pkin(2) * t165;
t136 = pkin(1) * t147;
t79 = t77 * t78;
t122 = t79 * t136;
t114 = 0.4e1 * t122;
t216 = -t171 * t114 + (t23 * t57 + t26 * t60) * t159;
t214 = 0.2e1 * pkin(2);
t51 = t58 ^ 2;
t55 = t59 ^ 2;
t170 = t51 + t55;
t90 = 0.1e1 / t222;
t39 = t170 * t90 * t78;
t35 = 0.1e1 / t39;
t213 = 0.1e1 / t40;
t31 = t39 ^ (-0.1e1 / 0.2e1);
t212 = -0.2e1 * t31;
t205 = pkin(1) * t77;
t203 = pkin(2) * t77;
t201 = pkin(3) * t80;
t200 = pkin(4) * t80;
t135 = 0.4e1 / t178 * t62 * t181;
t120 = -t135 / 0.4e1;
t61 = (t86 * t130 - 0.4e1 * t95 * t175) * qJD(2);
t108 = (-0.2e1 * t86 * t182 + (-t67 * t61 + t120) * t85) * t77;
t121 = t79 * t95 * t146;
t162 = 0.2e1 * t200;
t125 = t16 * t59 * t162;
t144 = t85 * t164;
t129 = t95 * t144;
t132 = qJD(2) * t153;
t149 = t77 * t182;
t45 = t53 * t55 + 0.1e1;
t41 = 0.1e1 / t45;
t156 = t41 * t200;
t169 = pkin(2) * t145 + t164 * t204;
t177 = t77 * t86;
t25 = t62 * t180 + 0.2e1 * t132 + t169;
t190 = t25 * t78;
t24 = (-t148 + (-t172 + (t76 + t161) * t85) * qJD(2)) * pkin(2);
t192 = t24 * t52 / t51;
t1 = 0.2e1 * (t125 * t53 + t162 * t194) / t45 ^ 2 * (t25 * t188 - t55 * t192) + 0.2e1 * (pkin(4) * t136 * t223 + t125 * t192) * t41 + 0.2e1 * ((-t16 * t25 + t18 * t24) * t53 + (-((t81 * t135 / 0.4e1 + t61 * t180 + 0.6e1 * pkin(1) * t129) * t209 + 0.4e1 * t59 * t121 + ((t149 / 0.2e1 - pkin(1) * t190) * t85 + ((t172 + (-t76 + t181) * t85) * t209 + (-t30 * t85 - t86 * t59) * t78 * pkin(1)) * qJD(2)) * pkin(2)) * t52 - ((0.4e1 * t132 + t169) * t210 - 0.4e1 * t58 * t121 + (-t108 / 0.2e1 + (t24 * t176 + (-t81 * t177 + (t28 * t85 + t58 * t86) * t78) * qJD(2)) * pkin(1)) * pkin(2)) * t188) * t89) * t156;
t199 = t1 * t77;
t191 = t24 * t58;
t9 = ((t25 * t59 + t191) * t159 - t170 * t114) * t90;
t198 = t31 * t35 * t9;
t10 = t216 * t93;
t197 = t10 * t33 * t213;
t183 = t60 * t86;
t187 = t57 * t85;
t112 = (t183 + t187) * t77;
t14 = t112 * t218;
t196 = t14 * t57;
t150 = t77 * t218;
t184 = t60 * t85;
t186 = t57 * t86;
t15 = (t184 - t186) * t150;
t195 = t15 * t60;
t185 = t59 * t78;
t174 = t85 * t86;
t160 = 0.2e1 * t33;
t157 = t48 * t201;
t43 = 0.1e1 / t46;
t139 = t43 * t49 * t201;
t124 = t60 * t139;
t140 = t43 * t157;
t128 = -0.2e1 * t142;
t27 = t72 + (-t172 + (t163 - t181) * t85) * pkin(1);
t17 = (t57 * t128 + t27 * t77) * t92;
t29 = -t63 * t179 + t175 * t214 + (t75 * t86 + t173) * pkin(1);
t19 = (t60 * t128 + t29 * t77) * t92;
t6 = -t17 * t124 + t19 * t140 + 0.1e1;
t155 = t6 * t197;
t154 = 0.1e1 / t39 ^ 2 * t78 * t9;
t152 = t77 * t197;
t143 = pkin(2) * t159;
t5 = t156 * t223;
t141 = t5 * t77 * t198;
t134 = qJD(2) * t158;
t133 = t96 * t147;
t127 = pkin(1) * t143;
t123 = pkin(2) * t134;
t119 = -0.4e1 * t122;
t118 = pkin(3) * t43 * t136;
t117 = 0.8e1 * t121;
t116 = t33 * t95 * t134;
t113 = 0.4e1 * t31 * t78 * t133;
t111 = t6 * t116;
t4 = ((-t184 / 0.2e1 + t186 / 0.2e1) * t152 + ((-t23 * t86 + t26 * t85) * t77 + (t112 + (t57 * t174 - t60 * t84) * t127) * qJD(2)) * t33) * t92;
t3 = -t60 * t150 * t165 + ((-t187 / 0.2e1 - t183 / 0.2e1) * t152 + ((t23 * t85 + t26 * t86) * t77 + (t57 * t177 + (-t60 * t174 - t57 * t84) * t127) * qJD(2)) * t33) * t92;
t2 = (((t82 * t120 + 0.6e1 * t86 * t133 - t61 * t179) * t77 + t60 * t117 + ((-0.2e1 * t26 * t78 * pkin(2) + t149) * t85 + ((t172 + (-t75 + t181) * t85) * t77 + (-t29 * t85 - t183) * t143) * qJD(2)) * pkin(1)) * t140 - ((0.4e1 * t131 + t168) * t77 + t57 * t117 + (t108 + (-0.2e1 * t23 * t176 + (t177 * t211 + (-t27 * t85 - t186) * t159) * qJD(2)) * pkin(2)) * pkin(1)) * t124) * t92 + (0.2e1 * t48 * t118 - t23 * t139 - 0.2e1 * t157 * t219) * t19 + (-0.2e1 * t118 * t189 - t26 * t139 + 0.2e1 * (t43 * t193 + t49 * t219) * t60 * t201) * t17;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t144, 0.2e1 * (t86 ^ 2 - t84) * qJD(2), 0, -0.2e1 * t144, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t3, -0.2e1 * t14 * t4 - 0.2e1 * t15 * t3, 0, 0.2e1 * t15 * t4, 0, 0, (-t15 * t165 + t4 * t86) * t214, (-t14 * t165 + t3 * t86) * t214, 0, -0.2e1 * t129, (-t55 * t154 + (t55 * t119 + 0.2e1 * t25 * t185) * t35) * t90, (t59 * t154 * t220 + (-0.2e1 * t24 * t185 + (0.8e1 * t59 * t122 - 0.2e1 * t190) * t58) * t35) * t90, 0, (-t51 * t154 + (t51 * t119 + t191 * t159) * t35) * t90, 0, 0, (t58 * t113 + (t58 * t198 + t24 * t212) * t205) * t89, (t59 * t113 + (t59 * t198 + t25 * t212) * t205) * t89, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, 0, -t165, 0, 0, 0, 0, 0, 0, 0, -t14 * t2 - t3 * t6, 0, t15 * t2 + t4 * t6, 0, 0, 0, (0.2e1 * (t195 + t196) * t116 + ((t195 / 0.2e1 + t196 / 0.2e1) * t197 + (-t14 * t23 - t15 * t26 - t3 * t57 - t4 * t60) * t33) * t203) * t92, 0, 0, 0, (-t59 * t141 / 0.2e1 + (t59 * t199 + (-0.2e1 * t59 * t123 + t25 * t77) * t5) * t31) * t89, 0, (t58 * t141 / 0.2e1 + (-t58 * t199 + (t123 * t220 - t24 * t77) * t5) * t31) * t89, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t6 * t2, (0.4e1 * t57 * t111 + (t57 * t155 + (-t2 * t57 - t23 * t6) * t160) * t203) * t92, (-0.4e1 * t60 * t111 + (-t60 * t155 + (t2 * t60 + t26 * t6) * t160) * t203) * t92, 0, (-0.1e1 / t40 ^ 2 * t10 * t217 + t216 * t213) * t95 * t93, 0, 0, 0, 0, 0, 0.2e1 * t5 * t1, 0, 0, 0, 0;];
MMD_reg = t7;
