% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% palh2m2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 14:45:26
% EndTime: 2019-06-06 14:45:27
% DurationCPUTime: 1.34s
% Computational Cost: add. (1845->207), mult. (5936->323), div. (0->0), fcn. (3426->8), ass. (0->145)
t139 = qJD(2) * qJD(3);
t95 = qJD(2) ^ 2;
t121 = -t95 + t139;
t92 = cos(qJ(2));
t142 = qJDD(2) * t92;
t88 = sin(qJ(2));
t143 = qJDD(2) * t88;
t87 = sin(qJ(3));
t91 = cos(qJ(3));
t23 = ((t121 * t92 - t143) * t91 + (t121 * t88 + t142) * t87) * pkin(4);
t51 = (t87 * t88 + t91 * t92) * pkin(4);
t32 = (-qJD(2) + qJD(3)) * t51;
t47 = qJD(2) * t51;
t184 = t32 + t47;
t82 = t87 ^ 2;
t84 = t91 ^ 2;
t183 = t82 + t84;
t89 = sin(qJ(1));
t93 = cos(qJ(1));
t178 = g(1) * t89 - g(2) * t93;
t150 = qJD(2) * t88;
t131 = t91 * t150;
t123 = pkin(4) * t131;
t161 = t87 * t92;
t136 = pkin(4) * t161;
t68 = qJD(2) * t136;
t48 = -t68 + t123;
t43 = t91 * t48;
t44 = qJD(3) * pkin(5) + t47;
t27 = -t87 * t44 - t43;
t182 = qJD(3) * t27;
t65 = g(1) * t93 + g(2) * t89;
t83 = t88 ^ 2;
t85 = t92 ^ 2;
t75 = t83 + t85;
t180 = pkin(4) * t75;
t132 = pkin(4) * t150;
t62 = t75 * qJD(1);
t179 = t62 * t132 - t178;
t168 = t91 * pkin(5);
t126 = -pkin(2) - t168;
t167 = t92 * pkin(4);
t77 = pkin(1) + t167;
t63 = t77 * qJD(1);
t30 = t126 * t62 - t63;
t39 = t183 * t62;
t18 = -t39 * pkin(3) + t30;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t110 = t90 * t18 + t86 * t27;
t148 = qJD(3) * t87;
t130 = t62 * t148;
t140 = qJD(1) * qJD(2);
t128 = t88 * t140;
t46 = pkin(4) * t128 - t77 * qJDD(1);
t79 = t83 * qJDD(1);
t80 = t85 * qJDD(1);
t61 = t79 + t80;
t19 = pkin(5) * t130 + t126 * t61 + t46;
t53 = t82 * t61;
t54 = t84 * t61;
t38 = t53 + t54;
t13 = -t38 * pkin(3) + t19;
t162 = t87 * t48;
t26 = t91 * t44 - t162;
t149 = qJD(3) * t26;
t160 = t91 * t23;
t125 = qJD(3) * t123 + t95 * t136 + (t142 * t91 + t143 * t87) * pkin(4);
t127 = t92 * t139;
t158 = t95 * t88;
t24 = (-t87 * t127 - t91 * t158) * pkin(4) + t125;
t20 = qJDD(3) * pkin(5) + t24;
t5 = -t87 * t20 - t149 + t160;
t1 = -t110 * qJD(4) - t86 * t13 + t90 * t5;
t9 = -t86 * t18 + t90 * t27;
t114 = -t110 * t90 + t9 * t86;
t146 = qJD(4) * t86;
t2 = -t86 * t5 + t18 * t146 + (-qJD(4) * t27 - t13) * t90;
t175 = -qJD(4) * t114 + t1 * t90 - t2 * t86;
t174 = pkin(4) * t88;
t173 = pkin(5) * t87;
t169 = g(3) * t92;
t35 = qJDD(4) + t38;
t166 = t35 * t86;
t165 = t35 * t90;
t37 = qJD(4) + t39;
t45 = t75 * pkin(5) * t148 + t132;
t164 = t37 * t45;
t159 = t91 * t61;
t108 = -t88 * t91 + t161;
t52 = t108 * pkin(4);
t109 = t91 * t51 + t87 * t52;
t33 = t68 + (-t108 * qJD(3) - t131) * pkin(4);
t10 = -t109 * qJD(3) + t91 * t32 - t87 * t33;
t29 = -t91 * t47 + t162;
t157 = t10 - t29;
t155 = -t44 + t47;
t154 = t65 * t174;
t153 = t82 - t84;
t152 = t83 - t85;
t151 = qJD(1) * t88;
t147 = qJD(3) * t91;
t145 = qJD(4) * t90;
t141 = t92 * qJDD(1);
t138 = qJDD(1) * pkin(1);
t137 = t62 * t173;
t60 = t62 ^ 2;
t135 = t87 * t60 * t91;
t96 = qJD(1) ^ 2;
t134 = t88 * t96 * t92;
t133 = pkin(4) * t151;
t124 = t62 * t133;
t122 = t91 * t130;
t120 = t92 * t128;
t119 = -0.2e1 * pkin(1) * t140;
t57 = t89 * t86 - t93 * t90;
t58 = -t93 * t86 - t89 * t90;
t118 = -g(1) * t58 + g(2) * t57;
t117 = -g(1) * t57 - g(2) * t58;
t113 = -t110 * t86 - t9 * t90;
t111 = pkin(2) + t77;
t31 = -t87 * t51 + t91 * t52;
t6 = t91 * t20 + t87 * t23 + t182;
t106 = -t26 * t155 * t87 + t6 * t168 + t65 * t173;
t105 = t95 * t92 + t143;
t104 = t6 * t109 + t154 + (t31 * qJD(3) + t184 * t87 + t91 * t33 + t43) * t26;
t41 = -t62 * pkin(2) - t63;
t49 = -t75 * pkin(2) - t77;
t103 = qJD(3) * (t41 * t75 + t49 * t62);
t101 = -t41 * t62 + t65;
t100 = pkin(1) * t96 + t65;
t99 = t178 + 0.2e1 * t138;
t40 = t126 * t75 - t77;
t36 = -t61 * pkin(2) + t46;
t98 = t36 * t75 + t49 * t61 + t179;
t94 = qJD(3) ^ 2;
t59 = t111 + t168;
t42 = t183 * t75;
t25 = -t42 * pkin(3) + t40;
t22 = -t86 * t133 + t90 * t29;
t21 = -t90 * t133 - t86 * t29;
t14 = t155 * t91;
t12 = -t86 * t137 + t90 * t14;
t11 = -t90 * t137 - t86 * t14;
t3 = [0, 0, 0, 0, 0, qJDD(1), t178, t65, 0, 0, t79 + 0.2e1 * t120, -0.2e1 * t152 * t140 + 0.2e1 * t88 * t141, t105, t80 - 0.2e1 * t120, t142 - t158, 0, t88 * t119 + t99 * t92, t92 * t119 - t99 * t88, -t65, (t178 + t138) * pkin(1), 0, 0, 0, t61 * t75, 0, 0, -t46 * t75 + t77 * t61 - t179, 0, -t105 * t180 - t65, -t63 * t132 + (t178 - t46) * t77, (0.2e1 * t122 + t53) * t75, 0.2e1 * (-t153 * t62 * qJD(3) + t87 * t159) * t75, (qJDD(3) * t87 + t91 * t94) * t75, (-0.2e1 * t122 + t54) * t75, (qJDD(3) * t91 - t87 * t94) * t75, 0, t103 * t87 - t91 * t98, t103 * t91 + t87 * t98, (qJD(3) * t29 - t24 * t87 + t160) * t75 - t65, t111 * t178 + t132 * t41 + t36 * t49, 0, 0, 0, t38 * t42, 0, 0, -t19 * t42 - t40 * t38 - t45 * t39 + t178, 0, t5 * t42 - t65, t178 * t59 + t19 * t40 + t30 * t45, 0, 0, 0, 0, 0, t35 * t42, -t90 * t164 + t2 * t42 + (t37 * t146 - t165) * t25 + t118, t86 * t164 - t1 * t42 + (t37 * t145 + t166) * t25 + t117, 0, -t114 * t45 + (qJD(4) * t113 - t1 * t86 - t2 * t90) * t25 + t178 * (pkin(3) + t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t152 * t96, t88 * qJDD(1), t134, t141, qJDD(2), t100 * t88 - t169, g(3) * t88 + t100 * t92, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t61 * t174, t154 + (qJDD(2) * t180 + t63 * t151 - t169) * pkin(4), 0, 0, 0, 0, 0, 0, t91 * t124 + t51 * qJDD(3) + (t33 + t48) * qJD(3), -qJD(3) * t184 - t52 * qJDD(3) - t87 * t124, t31 * t61 + (-t162 + (-qJD(3) * t52 - t33) * t87 + (-qJD(3) * t51 + t184) * t91) * t62, t23 * t52 + t24 * t51 - t48 * t32 + t47 * t33 + (-t41 * t151 - t169) * pkin(4) + t154, 0, 0, 0, 0, 0, 0, t39 * t133, 0, t157 * t39 + t31 * t38, t5 * t31 + t157 * t27 + (-t30 * t151 - t169) * pkin(4) + t104, 0, 0, 0, 0, 0, 0, -t31 * t166 + (-t10 * t86 - t31 * t145 - t21) * t37, -t31 * t165 + (-t10 * t90 + t31 * t146 + t22) * t37, 0, -g(3) * t167 - t113 * t10 + t110 * t21 + t175 * t31 - t9 * t22 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t153 * t60, t87 * t61, t135, t159, qJDD(3), -t48 * qJD(3) + (-pkin(4) * t158 - g(3)) * t91 + (-pkin(4) * t127 + t101) * t87 + t125, g(3) * t87 + t47 * qJD(3) + t101 * t91 - t23, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t137, 0, -t14 * t39 + (-t39 * t147 - t38 * t87) * pkin(5), -t27 * t14 + ((-g(3) - t182) * t91 + (-t30 * t62 - t149 - t5) * t87) * pkin(5) + t106, 0, 0, 0, 0, 0, 0, -t11 * t37 + (t87 * t166 + (t87 * t145 + t86 * t147) * t37) * pkin(5), t12 * t37 + (t87 * t165 + (-t87 * t146 + t90 * t147) * t37) * pkin(5), 0, t110 * t11 - t9 * t12 + ((qJD(3) * t113 - g(3)) * t91 + (-t149 - t175) * t87) * pkin(5) + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t9 * t37 + t118 + t2, -t110 * t37 - t1 + t117, 0, 0;];
tau_reg  = t3;
