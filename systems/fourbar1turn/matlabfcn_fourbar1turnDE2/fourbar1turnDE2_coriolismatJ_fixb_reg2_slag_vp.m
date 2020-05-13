% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(2*2)x(2*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1turnDE2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:34:59
% EndTime: 2020-04-12 19:35:14
% DurationCPUTime: 4.48s
% Computational Cost: add. (54815->174), mult. (77530->423), div. (3232->22), fcn. (21389->8), ass. (0->198)
t106 = pkin(1) ^ 2;
t96 = cos(qJ(2));
t227 = pkin(2) * t96;
t184 = -0.2e1 * pkin(1) * t227 + t106;
t235 = -pkin(3) - pkin(4);
t81 = (pkin(2) - t235) * (pkin(2) + t235) + t184;
t234 = -pkin(3) + pkin(4);
t82 = (pkin(2) - t234) * (pkin(2) + t234) + t184;
t203 = t81 * t82;
t109 = sqrt(-t203);
t95 = sin(qJ(2));
t190 = t109 * t95;
t73 = pkin(2) * t190;
t247 = pkin(3) ^ 2;
t248 = pkin(4) ^ 2;
t183 = -t247 + t248;
t105 = pkin(2) ^ 2;
t88 = t105 + t184;
t84 = t88 + t183;
t91 = pkin(1) - t227;
t66 = t84 * t91 - t73;
t61 = 0.1e1 / t66 ^ 2;
t228 = pkin(2) * t95;
t77 = t84 * t228;
t67 = t109 * t91 + t77;
t210 = t61 * t67;
t232 = pkin(1) * t95;
t86 = 0.1e1 / t88 ^ 2;
t172 = t86 * t232;
t156 = pkin(2) * t172;
t134 = t67 * t156;
t94 = t95 ^ 2;
t167 = t105 * t94 * pkin(1);
t196 = t84 * t227 + t73;
t75 = 0.1e1 / t109;
t205 = t75 * t91;
t237 = pkin(1) * pkin(2);
t148 = (-t81 - t82) * t237;
t71 = t95 * t148;
t39 = t71 * t205 + 0.2e1 * t167 + t196;
t85 = 0.1e1 / t88;
t212 = t39 * t85;
t99 = 0.1e1 / pkin(4);
t29 = (t212 / 0.2e1 - t134) * t99;
t60 = 0.1e1 / t66;
t220 = t29 * t60;
t135 = t66 * t156;
t185 = t96 * t109;
t206 = t75 * t71;
t37 = t77 + (-t185 + (0.2e1 * t91 * pkin(1) - t206) * t95) * pkin(2);
t214 = t37 * t85;
t27 = (-t214 / 0.2e1 + t135) * t99;
t250 = -0.2e1 * t27 * t210 - 0.2e1 * t220;
t103 = 0.1e1 / t247;
t192 = t103 * t86;
t72 = pkin(1) * t190;
t83 = t88 - t183;
t231 = pkin(1) * t96;
t92 = -pkin(2) + t231;
t65 = -t83 * t92 - t72;
t55 = t65 ^ 2;
t79 = t83 * t232;
t68 = -t109 * t92 + t79;
t64 = t68 ^ 2;
t199 = t55 + t64;
t48 = t199 * t192;
t246 = 0.1e1 / t48;
t42 = t48 ^ (-0.1e1 / 0.2e1);
t174 = pkin(1) * t228;
t87 = t85 * t86;
t155 = t87 * t174;
t245 = -0.4e1 * t155;
t193 = t103 * t85;
t244 = pkin(2) * t193;
t59 = t66 ^ 2;
t216 = t37 * t60 / t59;
t63 = t67 ^ 2;
t53 = t61 * t63 + 0.1e1;
t243 = (t39 * t210 - t63 * t216) / t53 ^ 2;
t57 = 0.1e1 / t65 ^ 2;
t211 = t57 * t68;
t240 = -0.2e1 * t92;
t36 = t79 + (-t185 + (pkin(2) * t240 - t206) * t95) * pkin(1);
t56 = 0.1e1 / t65;
t219 = t36 * t56 / t55;
t191 = t106 * t94;
t166 = pkin(2) * t191;
t197 = t83 * t231 + t72;
t204 = t75 * t92;
t38 = -t71 * t204 + 0.2e1 * t166 + t197;
t54 = t57 * t64 + 0.1e1;
t242 = (t38 * t211 - t64 * t219) / t54 ^ 2;
t126 = (t36 * t65 + t38 * t68) * t86;
t239 = 0.2e1 * t95;
t241 = t206 * t239 + t185;
t100 = 0.1e1 / t248;
t198 = t59 + t63;
t47 = t198 * t86 * t100;
t43 = 0.1e1 / t47;
t40 = t47 ^ (-0.1e1 / 0.2e1);
t238 = t42 * t246;
t236 = t85 / 0.2e1;
t233 = pkin(1) * t85;
t230 = pkin(2) * t85;
t229 = pkin(2) * t86;
t226 = pkin(3) * t88;
t225 = pkin(4) * t88;
t102 = 0.1e1 / pkin(3);
t194 = t102 * t42;
t162 = t85 * t194;
t207 = t68 * t96;
t24 = (t65 * t95 + t207) * t162;
t136 = t64 * t155;
t142 = -0.2e1 * t155;
t123 = 0.2e1 * t238 * (t55 * t142 + t126 - 0.2e1 * t136);
t122 = -t123 / 0.2e1;
t119 = t95 * t122;
t173 = 0.2e1 * t86;
t141 = t173 * t237;
t200 = t38 + t65;
t201 = t95 * t96;
t217 = t36 * t95;
t35 = t68 * t95 * t162;
t6 = -t35 + ((-t68 * t201 - t65 * t94) * t141 + (t200 * t96 + t217) * t85) * t194 + (t65 * t119 + t122 * t207) * t102 * t193;
t224 = t24 * t6;
t209 = t65 * t96;
t25 = -t162 * t209 + t35;
t121 = t123 / 0.2e1;
t120 = t65 * t121;
t7 = (((t65 * t201 - t68 * t94) * t141 + ((-t36 + t68) * t96 + t200 * t95) * t85) * t42 + (t68 * t119 + t96 * t120) * t193) * t102;
t223 = t25 * t7;
t215 = t37 * t66;
t20 = ((t39 * t67 + t215) * t173 + t198 * t245) * t100;
t222 = t20 * t40 * t43;
t221 = t25 * t68;
t208 = t67 * t86;
t202 = t85 * t96;
t1 = -t24 * t7 - t25 * t6;
t195 = t1 * qJD(1);
t4 = (-t24 * t95 + t6 * t96) * pkin(2);
t189 = t4 * qJD(1);
t5 = (-t25 * t95 + t7 * t96) * pkin(2);
t188 = t5 * qJD(1);
t165 = t20 / t47 ^ 2 * t86;
t8 = (-t67 * t66 * t165 + (t37 * t208 + (t67 * t245 + t39 * t86) * t66) * t43) * t100;
t187 = t8 * qJD(1);
t146 = -t165 / 0.2e1;
t9 = (t59 * t146 + (t59 * t142 + t86 * t215) * t43) * t100;
t186 = t9 * qJD(1);
t182 = qJD(2) * t96;
t181 = qJD(2) * t99;
t10 = (t63 * t146 + (t63 * t142 + t39 * t208) * t43) * t100;
t180 = t10 * qJD(1);
t168 = t106 * t228;
t128 = -0.2e1 * t40 * t86 * t168;
t157 = -t222 / 0.2e1;
t11 = (t66 * t128 + (t66 * t157 + t37 * t40) * t233) * t99;
t179 = t11 * qJD(1);
t145 = t67 * t157;
t12 = (t67 * t128 + (t39 * t40 + t145) * t233) * t99;
t178 = t12 * qJD(1);
t89 = t96 ^ 2 - t94;
t177 = t89 * qJD(1);
t176 = qJD(2) * t102;
t171 = t56 * t226;
t49 = 0.1e1 / t53;
t170 = t49 * t225;
t163 = 0.4e1 * t75 / t203 * t71 ^ 2;
t144 = -t163 / 0.4e1;
t159 = t105 * t191;
t69 = t96 * t148 - 0.4e1 * t159;
t124 = (-0.2e1 * t96 * t206 + (-t69 * t75 + t144) * t95) * t85;
t147 = t87 * t159;
t2 = 0.2e1 * pkin(4) * t49 * t174 * t250 + 0.4e1 * (t220 * t243 + (t49 * t216 + t61 * t243) * t27 * t67) * t225 + 0.2e1 * ((-t27 * t39 + t29 * t37) * t61 + (-((t91 * t163 / 0.4e1 + t69 * t205 - t77 + t241 * pkin(2)) * t236 + 0.4e1 * t67 * t147 + (0.3e1 * t105 * t85 * t201 + (-0.2e1 * t39 * t95 - t67 * t96) * t229) * pkin(1)) * t60 - (-(0.4e1 * t167 + t196) * t85 / 0.2e1 - 0.4e1 * t66 * t147 + (-t124 / 0.2e1 + (-t91 * t202 + (t37 * t239 + t66 * t96) * t86) * pkin(1)) * pkin(2)) * t210) * t99) * t170;
t169 = t2 * t40 * t85;
t161 = qJD(1) * t224;
t160 = qJD(1) * t223;
t158 = t95 * t182;
t90 = qJD(1) * t201;
t51 = 0.1e1 / t54;
t153 = t51 * t171;
t152 = t51 * t57 * t226;
t149 = t105 * t172;
t143 = -0.2e1 * t156;
t140 = t105 * t90;
t139 = pkin(3) * t51 * t174;
t133 = t68 * t152;
t28 = (t65 * t143 + t36 * t85) * t102;
t30 = (t68 * t143 + t38 * t85) * t102;
t14 = -t28 * t133 + t30 * t153 + 0.1e1;
t132 = t14 * t149;
t131 = 0.8e1 * t147;
t127 = t199 * t155;
t125 = 0.2e1 * t126;
t118 = t14 * t238 * (t125 - 0.4e1 * t127) * t244;
t13 = t170 * t250;
t3 = (((t92 * t144 + 0.6e1 * t96 * t168 - t69 * t204 - t79) * t85 + t68 * t131 + (t241 * t85 + (-0.4e1 * t38 * t95 - 0.2e1 * t207) * t229) * pkin(1)) * t153 - ((0.4e1 * t166 + t197) * t85 + t65 * t131 + (t124 + (t202 * t240 + (-0.2e1 * t209 - 0.4e1 * t217) * t86) * pkin(2)) * pkin(1)) * t133) * t102 + (0.2e1 * t56 * t139 - t36 * t152 - 0.2e1 * t171 * t242) * t30 + (-0.2e1 * t139 * t211 - t38 * t152 + 0.2e1 * (t51 * t219 + t57 * t242) * t68 * t226) * t28;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t89 * qJD(2), 0, -t158, 0, 0, 0, 0, 0, 0, qJD(2) * t224, t1 * qJD(2), 0, qJD(2) * t223, 0, 0, t5 * qJD(2), t4 * qJD(2), 0, -t105 * t158, t10 * qJD(2), -t8 * qJD(2), 0, t9 * qJD(2), 0, 0, -t11 * qJD(2), -t12 * qJD(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t177, t182, -t90, -qJD(2) * t95, 0, 0, 0, 0, 0, t161, t195, (-t14 * t6 - t24 * t3) * qJD(2), t160, (t14 * t7 + t25 * t3) * qJD(2), 0, t188, t189, ((0.2e1 * (t24 * t65 + t221) * t149 + (-t24 * t36 - t25 * t38 - t6 * t65 - t68 * t7) * t230) * t42 + (t120 * t24 + t121 * t221) * t244) * t176, -t140, t180, -t187, (t67 * t169 + (t85 * t145 + (-0.2e1 * t134 + t212) * t40) * t13) * t181, t186, (-t66 * t169 + (t66 * t222 * t236 + (0.2e1 * t135 - t214) * t40) * t13) * t181, 0, -t179, -t178, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t177, 0, t90, 0, 0, 0, 0, 0, 0, -t161, -t195, 0, -t160, 0, 0, -t188, -t189, 0, t140, -t180, t187, 0, -t186, 0, 0, t179, t178, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t3 * qJD(2), ((0.2e1 * t65 * t132 + (-t14 * t36 - t3 * t65) * t230) * t42 + t65 * t118 / 0.2e1) * t176, ((-0.2e1 * t68 * t132 + (t14 * t38 + t3 * t68) * t230) * t42 - t68 * t118 / 0.2e1) * t176, 0, (((-t125 / 0.2e1 + 0.2e1 * t136) * t64 + (-t126 + (0.2e1 * t55 + 0.4e1 * t64) * t155) * t55) / t48 ^ 2 * t192 + (t126 - 0.2e1 * t127) * t246) * t105 * t103 * qJD(2), 0, 0, 0, 0, 0, t13 * t2 * qJD(2), 0, 0, 0, 0;];
cmat_reg = t15;
