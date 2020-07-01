% Calculate minimal parameter regressor of coriolis matrix for
% fourbar1turnTE
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
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1turnTE_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:36
% EndTime: 2020-06-27 16:22:53
% DurationCPUTime: 2.16s
% Computational Cost: add. (25595->131), mult. (35160->320), div. (1166->15), fcn. (9788->4), ass. (0->156)
t83 = 0.1e1 / pkin(4);
t79 = sin(qJ(2));
t180 = pkin(1) * t79;
t80 = cos(qJ(2));
t175 = pkin(2) * t80;
t88 = pkin(1) ^ 2;
t150 = -0.2e1 * pkin(1) * t175 + t88;
t87 = pkin(2) ^ 2;
t73 = t87 + t150;
t71 = 0.1e1 / t73 ^ 2;
t134 = t71 * t180;
t119 = pkin(2) * t134;
t70 = 0.1e1 / t73;
t184 = -t70 / 0.2e1;
t182 = -pkin(3) - pkin(4);
t66 = (pkin(2) - t182) * (pkin(2) + t182) + t150;
t181 = -pkin(3) + pkin(4);
t67 = (pkin(2) - t181) * (pkin(2) + t181) + t150;
t160 = t66 * t67;
t89 = sqrt(-t160);
t154 = t80 * t89;
t110 = pkin(1) * pkin(2) * (-t66 - t67);
t56 = t79 * t110;
t60 = 0.1e1 / t89;
t163 = t60 * t56;
t176 = pkin(2) * t79;
t195 = pkin(4) ^ 2;
t149 = pkin(3) ^ 2 - t195;
t69 = t73 - t149;
t62 = t69 * t176;
t75 = pkin(1) - t175;
t26 = t62 + (-t154 + (0.2e1 * t75 * pkin(1) - t163) * t79) * pkin(2);
t157 = t79 * t89;
t58 = pkin(2) * t157;
t51 = t69 * t75 - t58;
t99 = t51 * t119 + t26 * t184;
t16 = t99 * t83;
t46 = 0.1e1 / t51 ^ 2;
t52 = t75 * t89 + t62;
t165 = t46 * t52;
t183 = t70 / 0.2e1;
t78 = t79 ^ 2;
t131 = t87 * t78 * pkin(1);
t151 = t69 * t175 + t58;
t162 = t60 * t75;
t28 = t56 * t162 + 0.2e1 * t131 + t151;
t98 = -t52 * t119 + t28 * t183;
t19 = t98 * t83;
t45 = 0.1e1 / t51;
t169 = t19 * t45;
t197 = -0.2e1 * t16 * t165 - 0.2e1 * t169;
t167 = t26 * t45 * t46;
t48 = t52 ^ 2;
t36 = t46 * t48 + 0.1e1;
t194 = (t28 * t165 - t48 * t167) / t36 ^ 2;
t57 = pkin(1) * t157;
t68 = t73 + t149;
t179 = pkin(1) * t80;
t76 = -pkin(2) + t179;
t50 = -t68 * t76 - t57;
t43 = 0.1e1 / t50 ^ 2;
t64 = t68 * t180;
t53 = -t76 * t89 + t64;
t166 = t43 * t53;
t191 = -0.2e1 * t76;
t25 = t64 + (-t154 + (pkin(2) * t191 - t163) * t79) * pkin(1);
t42 = 0.1e1 / t50;
t168 = t25 * t42 * t43;
t153 = t88 * t78;
t130 = pkin(2) * t153;
t152 = t68 * t179 + t57;
t161 = t60 * t76;
t27 = -t56 * t161 + 0.2e1 * t130 + t152;
t49 = t53 ^ 2;
t37 = t43 * t49 + 0.1e1;
t193 = (t27 * t166 - t49 * t168) / t37 ^ 2;
t189 = 0.2e1 * t79;
t192 = t163 * t189 + t154;
t190 = -0.4e1 * t79;
t188 = -t25 / 0.2e1;
t187 = -t50 / 0.2e1;
t186 = t52 / 0.4e1;
t185 = t53 / 0.2e1;
t178 = pkin(2) * t70;
t177 = pkin(2) * t71;
t174 = pkin(3) * t73;
t173 = pkin(4) * t73;
t86 = 0.1e1 / pkin(3);
t159 = t70 * t86;
t164 = t53 * t80;
t29 = -(t79 * t50 + t164) * t159 / 0.2e1;
t120 = t187 - t27 / 0.2e1;
t121 = t188 + t185;
t136 = pkin(1) * t177;
t158 = t79 * t80;
t8 = ((t53 * t158 + t50 * t78) * t136 + (t120 * t80 + t121 * t79) * t70) * t86;
t172 = t29 * t8;
t156 = t80 * t50;
t30 = (t79 * t185 - t156 / 0.2e1) * t159;
t9 = ((t79 * t156 - t53 * t78) * t136 + (-t120 * t79 + t121 * t80) * t70) * t86;
t171 = t30 * t9;
t126 = t87 * t153;
t72 = t70 * t71;
t112 = t72 * t126;
t127 = 0.4e1 * t60 / t160 * t56 ^ 2;
t32 = 0.1e1 / t36;
t132 = t32 * t173;
t135 = pkin(1) * t176;
t155 = t80 * t70;
t54 = t110 * t80 - 0.4e1 * t126;
t109 = -t127 / 0.4e1;
t97 = (-0.2e1 * t80 * t163 + (-t54 * t60 + t109) * t79) * t70;
t1 = 0.2e1 * pkin(4) * t32 * t135 * t197 + 0.4e1 * (t169 * t194 + (t32 * t167 + t46 * t194) * t16 * t52) * t173 + 0.2e1 * ((-t16 * t28 + t19 * t26) * t46 + (-((t75 * t127 / 0.4e1 + t54 * t162 - t62 + t192 * pkin(2)) * t183 + 0.4e1 * t52 * t112 + (0.3e1 * t87 * t79 * t155 + (-0.2e1 * t28 * t79 - t52 * t80) * t177) * pkin(1)) * t45 - ((0.4e1 * t131 + t151) * t184 - 0.4e1 * t51 * t112 + (-t97 / 0.2e1 + (-t75 * t155 + (t26 * t189 + t51 * t80) * t71) * pkin(1)) * pkin(2)) * t165) * t83) * t132;
t170 = t70 * t1;
t3 = t29 * t9 + t30 * t8;
t148 = t3 * qJD(1);
t6 = (t29 * t79 - t8 * t80) * pkin(2);
t147 = t6 * qJD(1);
t7 = (-t30 * t79 + t80 * t9) * pkin(2);
t146 = t7 * qJD(1);
t145 = qJD(2) * t83;
t144 = qJD(2) * t86;
t118 = t72 * t135;
t84 = 0.1e1 / t195;
t10 = (-t52 * t51 * t118 + (t28 * t51 / 0.4e1 + t26 * t186) * t71) * t84;
t143 = t10 * qJD(1);
t14 = (t71 * t28 * t186 - t48 * t118 / 0.2e1) * t84;
t142 = t14 * qJD(1);
t129 = t88 * t176;
t102 = t71 * t83 * t129;
t111 = pkin(1) * t83 * t183;
t15 = -t102 * t51 + t111 * t26;
t141 = t15 * qJD(1);
t18 = -t102 * t52 + t111 * t28;
t140 = t18 * qJD(1);
t74 = t80 ^ 2 - t78;
t139 = t74 * qJD(1);
t138 = t80 * qJD(2);
t133 = t42 * t174;
t125 = qJD(1) * t172;
t124 = qJD(1) * t171;
t123 = qJD(1) * t158;
t34 = 0.1e1 / t37;
t116 = t34 * t133;
t115 = t34 * t43 * t174;
t108 = -0.2e1 * t119;
t107 = pkin(3) * t34 * t135;
t101 = t53 * t115;
t17 = (t108 * t50 + t25 * t70) * t86;
t20 = (t108 * t53 + t27 * t70) * t86;
t5 = -t101 * t17 + t116 * t20 + 0.1e1;
t104 = t5 * t87 * t134;
t103 = 0.8e1 * t112;
t4 = t132 * t197;
t2 = (((t109 * t76 + 0.6e1 * t129 * t80 - t54 * t161 - t64) * t70 + t53 * t103 + (t192 * t70 + (t27 * t190 - 0.2e1 * t164) * t177) * pkin(1)) * t116 - ((0.4e1 * t130 + t152) * t70 + t50 * t103 + (t97 + (t155 * t191 + (t25 * t190 - 0.2e1 * t156) * t71) * pkin(2)) * pkin(1)) * t101) * t86 + (0.2e1 * t42 * t107 - t25 * t115 - 0.2e1 * t133 * t193) * t20 + (-0.2e1 * t107 * t166 - t27 * t115 + 0.2e1 * (t34 * t168 + t43 * t193) * t53 * t174) * t17;
t11 = [0, 0, 0, t79 * t138, t74 * qJD(2), 0, 0, 0, 0, 0, qJD(2) * t172, t3 * qJD(2), 0, qJD(2) * t171, 0, 0, t7 * qJD(2), t6 * qJD(2), t14 * qJD(2), -t10 * qJD(2), 0, 0, 0, -t15 * qJD(2), -t18 * qJD(2); 0, 0, 0, t123, t139, t138, -t79 * qJD(2), 0, 0, 0, t125, t148, (t2 * t29 + t5 * t8) * qJD(2), t124, (t2 * t30 + t5 * t9) * qJD(2), 0, t146, t147, t142, -t143, (t52 * t170 / 0.2e1 + t98 * t4) * t145, (-t51 * t170 / 0.2e1 + t99 * t4) * t145, 0, -t141, -t140; 0, 0, 0, -t123, -t139, 0, 0, 0, 0, 0, -t125, -t148, 0, -t124, 0, 0, -t146, -t147, -t142, t143, 0, 0, 0, t141, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t2 * qJD(2), (t50 * t104 + (t2 * t187 + t5 * t188) * t178) * t144, (-t53 * t104 + (t27 * t5 / 0.2e1 + t2 * t185) * t178) * t144, 0, 0, 0, 0, t4 * t1 * qJD(2), 0, 0;];
cmat_reg = t11;
