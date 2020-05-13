% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((2+1)*2/2)x(2*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbar1turnTE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:19:52
% EndTime: 2020-04-12 19:20:02
% DurationCPUTime: 2.68s
% Computational Cost: add. (24448->148), mult. (34317->379), div. (1034->16), fcn. (9255->4), ass. (0->155)
t193 = pkin(4) ^ 2;
t192 = pkin(3) ^ 2;
t75 = cos(qJ(2));
t173 = pkin(2) * t75;
t85 = pkin(1) ^ 2;
t145 = -0.2e1 * pkin(1) * t173 + t85;
t84 = pkin(2) ^ 2;
t69 = t84 + t145;
t67 = 0.1e1 / t69 ^ 2;
t137 = 0.2e1 * t67;
t74 = sin(qJ(2));
t73 = t74 ^ 2;
t152 = t73 * t85;
t128 = qJD(2) * t152;
t115 = pkin(2) * t128;
t143 = t74 * qJD(2);
t180 = -pkin(3) - pkin(4);
t62 = (pkin(2) - t180) * (pkin(2) + t180) + t145;
t179 = pkin(4) - pkin(3);
t63 = (pkin(2) - t179) * (pkin(2) + t179) + t145;
t156 = t62 * t63;
t86 = sqrt(-t156);
t127 = t86 * t143;
t142 = t75 * qJD(2);
t144 = t192 - t193;
t64 = t69 + t144;
t178 = pkin(1) * t64;
t146 = pkin(1) * t127 + t142 * t178;
t56 = 0.1e1 / t86;
t71 = pkin(1) * t75 - pkin(2);
t157 = t56 * t71;
t114 = pkin(1) * pkin(2) * (-t62 - t63);
t52 = t74 * t114;
t51 = qJD(2) * t52;
t19 = -t157 * t51 + 0.2e1 * t115 + t146;
t191 = t19 / 0.2e1;
t151 = t74 * t86;
t65 = t69 - t144;
t70 = pkin(1) - t173;
t47 = -pkin(2) * t151 + t65 * t70;
t90 = t47 ^ 2;
t42 = 0.1e1 / t90;
t176 = pkin(2) * t65;
t60 = t74 * t176;
t48 = t70 * t86 + t60;
t163 = t42 * t48;
t160 = t56 * t51;
t131 = t74 * t160;
t139 = 0.2e1 * t70 * pkin(1);
t148 = t75 * t86;
t17 = (-t131 + (-t148 + (t65 + t139) * t74) * qJD(2)) * pkin(2);
t41 = 0.1e1 / t47;
t167 = t17 * t41 * t42;
t134 = t84 * t73 * pkin(1);
t116 = qJD(2) * t134;
t147 = pkin(2) * t127 + t142 * t176;
t158 = t56 * t70;
t18 = t158 * t51 + 0.2e1 * t116 + t147;
t44 = t48 ^ 2;
t30 = t42 * t44 + 0.1e1;
t190 = 0.1e1 / t30 ^ 2 * (t163 * t18 - t167 * t44);
t61 = t74 * t178;
t49 = -t71 * t86 + t61;
t161 = t49 * t19;
t186 = -0.2e1 * t71;
t141 = pkin(2) * t186;
t16 = (-t131 + (-t148 + (t64 + t141) * t74) * qJD(2)) * pkin(1);
t46 = -pkin(1) * t151 - t64 * t71;
t38 = 0.1e1 / t46;
t88 = t46 ^ 2;
t39 = 0.1e1 / t88;
t168 = t16 * t38 * t39;
t45 = t49 ^ 2;
t31 = t39 * t45 + 0.1e1;
t189 = 0.1e1 / t31 ^ 2 * (t161 * t39 - t168 * t45);
t159 = t56 * t52;
t187 = 0.2e1 * pkin(2);
t66 = 0.1e1 / t69;
t185 = -t16 / 0.2e1;
t184 = -t66 / 0.2e1;
t183 = t66 / 0.2e1;
t182 = -t74 / 0.2e1;
t181 = t74 / 0.2e1;
t177 = pkin(1) * t67;
t175 = pkin(2) * t66;
t174 = pkin(2) * t67;
t172 = pkin(3) * t69;
t171 = pkin(4) * t69;
t138 = pkin(1) * t174;
t125 = t74 * t138;
t21 = t60 + (-t148 + (t139 - t159) * t74) * pkin(2);
t78 = 0.1e1 / pkin(4);
t10 = (t125 * t47 + t184 * t21) * t78;
t68 = t66 * t67;
t104 = t68 * t84 * t128;
t126 = t74 * t142;
t113 = t84 * t126;
t119 = 0.4e1 / t156 * t51 * t159;
t23 = t52 * t158 + 0.2e1 * t134 + (t65 * t75 + t151) * pkin(2);
t12 = (-t125 * t48 + t183 * t23) * t78;
t132 = t66 * t160;
t133 = t10 * t163;
t26 = 0.1e1 / t30;
t135 = t26 * t171;
t140 = 0.2e1 * t171;
t153 = t67 * t74;
t155 = t66 * t75;
t164 = t18 * t67;
t169 = t12 * t41;
t50 = (t114 * t75 - 0.4e1 * t152 * t84) * qJD(2);
t103 = -t119 / 0.4e1;
t95 = (t74 * t103 + 0.2e1 * (t50 * t182 - t75 * t51) * t56) * t66;
t129 = pkin(2) * t143;
t120 = pkin(1) * t129;
t98 = -0.2e1 * pkin(4) * t26 * t120;
t1 = 0.2e1 * t98 * t133 + 0.2e1 * (t140 * t190 + t98) * t169 + 0.2e1 * (t26 * t167 + t42 * t190) * t10 * t48 * t140 + 0.2e1 * ((-t10 * t18 + t12 * t17) * t42 + (-((t70 * t119 / 0.4e1 + t50 * t158 + 0.6e1 * pkin(1) * t113) * t183 + 0.4e1 * t48 * t104 + ((t132 / 0.2e1 - pkin(1) * t164) * t74 + ((t148 + (-t65 + t159) * t74) * t183 + (-t23 * t74 - t75 * t48) * t177) * qJD(2)) * pkin(2)) * t41 - ((0.4e1 * t116 + t147) * t184 - 0.4e1 * t47 * t104 + (-t95 / 0.2e1 + (t17 * t153 + (-t70 * t155 + (t21 * t74 + t47 * t75) * t67) * qJD(2)) * pkin(1)) * pkin(2)) * t163) * t78) * t135;
t170 = t1 * t66;
t166 = t17 * t66;
t165 = t18 * t66;
t162 = t48 * t67;
t81 = 0.1e1 / pkin(3);
t154 = t66 * t81;
t150 = t75 * t46;
t149 = t75 * t49;
t136 = t38 * t172;
t130 = qJD(2) * t177;
t28 = 0.1e1 / t31;
t124 = t28 * t136;
t123 = t28 * t39 * t172;
t118 = t74 * t130;
t117 = t85 * t129;
t112 = -0.2e1 * t125;
t109 = t49 * t123;
t108 = pkin(2) * t118;
t106 = t68 * t120;
t105 = t84 * t118;
t102 = pkin(3) * t28 * t120;
t101 = t117 * t137;
t20 = t61 + (-t148 + (t141 - t159) * t74) * pkin(1);
t11 = (t112 * t46 + t20 * t66) * t81;
t22 = -t52 * t157 + t152 * t187 + (t64 * t75 + t151) * pkin(1);
t13 = (t112 * t49 + t22 * t66) * t81;
t4 = -t109 * t11 + t124 * t13 + 0.1e1;
t100 = t4 * t105;
t99 = 0.8e1 * t104;
t97 = t150 / 0.2e1 + t49 * t182;
t96 = (t46 * t181 + t149 / 0.2e1) * t66;
t79 = 0.1e1 / t193;
t25 = t97 * t154;
t24 = t81 * t96;
t6 = ((t181 * t19 + t185 * t75) * t66 + (t96 + (t150 * t74 - t49 * t73) * t138) * qJD(2)) * t81;
t5 = (t97 * qJD(2) + t16 * t181 + t75 * t191) * t154 + (-t149 * t74 - t46 * t73) * pkin(2) * t81 * t130;
t3 = 0.2e1 * (-t133 - t169) * t135;
t2 = (((t71 * t103 + 0.6e1 * t75 * t117 - t50 * t157) * t66 + t49 * t99 + ((-0.2e1 * t19 * t174 + t132) * t74 + ((t148 + (-t64 + t159) * t74) * t66 + (-t22 * t74 - t149) * pkin(2) * t137) * qJD(2)) * pkin(1)) * t124 - ((0.4e1 * t115 + t146) * t66 + t46 * t99 + (t95 + (-0.2e1 * t16 * t153 + (t155 * t186 + (-t20 * t74 - t150) * t137) * qJD(2)) * pkin(2)) * pkin(1)) * t109) * t81 + (0.2e1 * t38 * t102 - t16 * t123 - 0.2e1 * t136 * t189) * t13 + (-t19 * t123 + (-0.2e1 * t39 * t102 + 0.2e1 * (t28 * t168 + t39 * t189) * t172) * t49) * t11;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t126, 0.2e1 * (t75 ^ 2 - t73) * qJD(2), 0, -0.2e1 * t126, 0, 0, 0, 0, 0, 0, 0.2e1 * t24 * t5, -0.2e1 * t24 * t6 + 0.2e1 * t25 * t5, 0, -0.2e1 * t25 * t6, 0, 0, (t143 * t25 + t6 * t75) * t187, (-t143 * t24 + t5 * t75) * t187, 0, -0.2e1 * t113, (t18 * t162 / 0.2e1 - t44 * t106) * t79, (-t17 * t162 / 0.2e1 + (-t164 / 0.2e1 + 0.2e1 * t48 * t106) * t47) * t79, 0, (t47 * t67 * t17 / 0.2e1 - t90 * t106) * t79, 0, 0, (-pkin(1) * t166 + t101 * t47) * t78, (-pkin(1) * t165 + t101 * t48) * t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, -t143, 0, 0, 0, 0, 0, 0, 0, -t2 * t24 - t4 * t5, 0, -t2 * t25 + t4 * t6, 0, 0, 0, ((t24 * t46 - t25 * t49) * t105 + (t25 * t191 - t49 * t6 / 0.2e1 + t24 * t185 - t46 * t5 / 0.2e1) * t175) * t81, 0, 0, 0, (t48 * t170 / 0.2e1 + (t165 / 0.2e1 - t48 * t108) * t3) * t78, 0, (-t47 * t170 / 0.2e1 + (-t166 / 0.2e1 + t47 * t108) * t3) * t78, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t4 * t2, (0.2e1 * t46 * t100 + (-t16 * t4 - t2 * t46) * t175) * t81, (-0.2e1 * t49 * t100 + (t19 * t4 + t2 * t49) * t175) * t81, 0, ((t161 / 0.2e1 + t46 * t16 / 0.2e1) * t67 + (-t45 - t88) * t106) * t84 / t192, 0, 0, 0, 0, 0, 0.2e1 * t3 * t1, 0, 0, 0, 0;];
MMD_reg = t7;
