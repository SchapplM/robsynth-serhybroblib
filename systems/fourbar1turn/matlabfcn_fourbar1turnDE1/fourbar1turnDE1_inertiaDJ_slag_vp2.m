% Calculate time derivative of joint inertia matrix for
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
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:35
% EndTime: 2020-04-12 19:25:50
% DurationCPUTime: 4.34s
% Computational Cost: add. (51068->200), mult. (73237->468), div. (2760->21), fcn. (19739->8), ass. (0->179)
t80 = cos(qJ(2));
t214 = pkin(2) * t80;
t167 = -0.2e1 * t214;
t90 = pkin(1) ^ 2;
t172 = pkin(1) * t167 + t90;
t89 = pkin(2) ^ 2;
t74 = t89 + t172;
t72 = 0.1e1 / t74 ^ 2;
t79 = sin(qJ(2));
t182 = t72 * t79;
t163 = pkin(1) * t182;
t145 = pkin(2) * t163;
t131 = -0.2e1 * t145;
t76 = pkin(1) * t80 - pkin(2);
t226 = -0.2e1 * t76;
t168 = pkin(2) * t226;
t221 = -pkin(3) - pkin(4);
t67 = (pkin(2) - t221) * (pkin(2) + t221) + t172;
t220 = pkin(4) - pkin(3);
t68 = (pkin(2) - t220) * (pkin(2) + t220) + t172;
t185 = t67 * t68;
t93 = sqrt(-t185);
t177 = t80 * t93;
t134 = pkin(1) * pkin(2) * (-t67 - t68);
t57 = t79 * t134;
t61 = 0.1e1 / t93;
t188 = t61 * t57;
t236 = pkin(3) ^ 2;
t237 = pkin(4) ^ 2;
t171 = t236 - t237;
t69 = t74 + t171;
t217 = pkin(1) * t69;
t66 = t79 * t217;
t23 = t66 + (-t177 + (t168 - t188) * t79) * pkin(1);
t178 = t79 * t93;
t51 = -pkin(1) * t178 - t69 * t76;
t71 = 0.1e1 / t74;
t86 = 0.1e1 / pkin(3);
t16 = (t131 * t51 + t23 * t71) * t86;
t78 = t79 ^ 2;
t181 = t78 * t90;
t186 = t61 * t76;
t25 = -t57 * t186 + 0.2e1 * pkin(2) * t181 + (t69 * t80 + t178) * pkin(1);
t54 = -t76 * t93 + t66;
t18 = (t131 * t54 + t25 * t71) * t86;
t43 = 0.1e1 / t51 ^ 2;
t195 = t43 * t54;
t42 = 0.1e1 / t51;
t238 = -t16 * t195 + t18 * t42;
t157 = t89 * t78 * pkin(1);
t136 = qJD(2) * t157;
t170 = qJD(2) * t79;
t147 = t93 * t170;
t169 = qJD(2) * t80;
t70 = t74 - t171;
t215 = pkin(2) * t70;
t174 = pkin(2) * t147 + t169 * t215;
t75 = pkin(1) - t214;
t187 = t61 * t75;
t56 = qJD(2) * t57;
t21 = t187 * t56 + 0.2e1 * t136 + t174;
t65 = t79 * t215;
t53 = t75 * t93 + t65;
t192 = t53 * t21;
t189 = t61 * t56;
t151 = t79 * t189;
t165 = 0.2e1 * t75 * pkin(1);
t20 = (-t151 + (-t177 + (t70 + t165) * t79) * qJD(2)) * pkin(2);
t52 = -pkin(2) * t178 + t70 * t75;
t45 = t52 ^ 2;
t46 = 0.1e1 / t52;
t197 = t20 * t46 / t45;
t47 = 0.1e1 / t52 ^ 2;
t49 = t53 ^ 2;
t39 = t47 * t49 + 0.1e1;
t235 = 0.1e1 / t39 ^ 2 * (t192 * t47 - t197 * t49);
t19 = (-t151 + (-t177 + (t69 + t168) * t79) * qJD(2)) * pkin(1);
t41 = t51 ^ 2;
t198 = t19 * t42 / t41;
t149 = qJD(2) * t181;
t135 = pkin(2) * t149;
t173 = pkin(1) * t147 + t169 * t217;
t22 = -t186 * t56 + 0.2e1 * t135 + t173;
t50 = t54 ^ 2;
t40 = t43 * t50 + 0.1e1;
t234 = (t195 * t22 - t198 * t50) / t40 ^ 2;
t180 = t79 * t51;
t190 = t54 * t80;
t118 = t180 + t190;
t183 = t71 * t86;
t176 = t41 + t50;
t87 = 0.1e1 / t236;
t34 = t176 * t87 * t72;
t29 = t34 ^ (-0.1e1 / 0.2e1);
t11 = t118 * t29 * t183;
t233 = t11 * t51;
t175 = t45 + t49;
t84 = 0.1e1 / t237;
t33 = t175 * t84 * t72;
t27 = t33 ^ (-0.1e1 / 0.2e1);
t232 = t27 * t71;
t231 = Ifges(3,1) - Ifges(3,2);
t223 = -t71 / 0.2e1;
t24 = t65 + (-t177 + (t165 - t188) * t79) * pkin(2);
t83 = 0.1e1 / pkin(4);
t15 = (t145 * t52 + t223 * t24) * t83;
t194 = t47 * t53;
t222 = t71 / 0.2e1;
t26 = t57 * t187 + 0.2e1 * t157 + (t70 * t80 + t178) * pkin(2);
t17 = (-t145 * t53 + t222 * t26) * t83;
t200 = t17 * t46;
t230 = t15 * t194 + t200;
t228 = 0.1e1 / t34;
t227 = 0.2e1 * t54;
t216 = pkin(1) * t72;
t213 = pkin(3) * t74;
t212 = pkin(4) * t74;
t150 = pkin(2) * t170;
t139 = pkin(1) * t150;
t73 = t71 * t72;
t123 = t73 * t139;
t108 = t176 * t123;
t111 = t72 * (t19 * t51 + t22 * t54);
t8 = (-0.4e1 * t108 + 0.2e1 * t111) * t87;
t211 = t29 * t228 * t8;
t210 = Ifges(5,1) * t53;
t207 = Ifges(5,4) * t52;
t206 = Ifges(5,4) * t53;
t204 = Ifges(5,5) * t53;
t203 = Ifges(5,2) * t52;
t201 = Ifges(5,6) * t52;
t196 = t27 * t83;
t193 = t51 * t80;
t191 = t54 * t79;
t184 = t71 * t80;
t179 = t79 * t80;
t164 = 0.2e1 * t72;
t37 = 0.1e1 / t40;
t162 = t37 * t213;
t35 = 0.1e1 / t39;
t160 = t35 * t212;
t159 = 0.1e1 / t33 * ((t52 * t20 + t192) * t164 - 0.4e1 * t175 * t123) * t84 * t232;
t158 = t71 * t211;
t154 = t71 * t196;
t152 = t71 * t189;
t148 = t79 * t169;
t146 = pkin(2) * t164;
t140 = t83 * t159;
t138 = 0.4e1 / t185 * t56 * t188;
t137 = qJD(2) * t163;
t133 = t89 * t148;
t132 = t154 / 0.2e1;
t130 = pkin(1) * t146;
t126 = t52 * t140;
t125 = t53 * t140;
t124 = pkin(2) * t137;
t122 = t73 * t89 * t149;
t121 = Ifges(5,5) * t132;
t120 = -Ifges(5,6) * t154 / 0.2e1;
t119 = -t138 / 0.4e1;
t115 = 0.8e1 * t122;
t114 = t29 * t89 * t137;
t6 = t238 * t162 + 0.1e1;
t112 = t6 * (mrSges(4,1) * t51 - mrSges(4,2) * t54);
t110 = (t191 - t193) * t71;
t109 = t124 * t196;
t107 = t52 * t109;
t106 = t53 * t109;
t104 = 0.2e1 * t124 * t52 - t20 * t71;
t103 = -0.2e1 * t124 * t53 + t21 * t71;
t55 = (t134 * t80 - 0.4e1 * t181 * t89) * qJD(2);
t102 = (-0.2e1 * t80 * t189 + (-t55 * t61 + t119) * t79) * t71;
t12 = t86 * t29 * t110;
t5 = -0.2e1 * t230 * t160;
t4 = ((-t191 / 0.2e1 + t193 / 0.2e1) * t158 + ((-t19 * t80 + t22 * t79) * t71 + (t118 * t71 + (t179 * t51 - t54 * t78) * t130) * qJD(2)) * t29) * t86;
t3 = ((t180 / 0.2e1 + t190 / 0.2e1) * t158 + ((-t19 * t79 - t22 * t80) * t71 + (t110 + (t179 * t54 + t51 * t78) * t130) * qJD(2)) * t29) * t86;
t2 = 0.2e1 * t238 * pkin(3) * t37 * t139 + (-t16 * t22 - t18 * t19 - ((0.4e1 * t135 + t173) * t71 + t51 * t115 + (t102 + (-0.2e1 * t19 * t182 + (t184 * t226 + (-t23 * t79 - t193) * t164) * qJD(2)) * pkin(2)) * pkin(1)) * t54 * t86) * t43 * t162 + ((((0.6e1 * t90 * pkin(2) * t148 + t76 * t119 - t55 * t186) * t71 + t54 * t115 + ((-0.2e1 * t22 * t72 * pkin(2) + t152) * t79 + ((t177 + (-t69 + t188) * t79) * t71 + (-t25 * t79 - t190) * t146) * qJD(2)) * pkin(1)) * t86 * t37 - 0.2e1 * t18 * t234) * t42 + (t37 * t198 + t43 * t234) * t16 * t227) * t213;
t1 = -0.4e1 * t230 * pkin(4) * t35 * t139 + 0.4e1 * (t200 * t235 + (t197 * t35 + t47 * t235) * t15 * t53) * t212 + 0.2e1 * ((-t15 * t21 + t17 * t20) * t47 + (-((t75 * t138 / 0.4e1 + t55 * t187 + 0.6e1 * pkin(1) * t133) * t222 + 0.4e1 * t53 * t122 + ((t152 / 0.2e1 - t21 * t216) * t79 + ((t177 + (-t70 + t188) * t79) * t222 + (-t26 * t79 - t53 * t80) * t216) * qJD(2)) * pkin(2)) * t46 - ((0.4e1 * t136 + t174) * t223 - 0.4e1 * t52 * t122 + (-t102 / 0.2e1 + (t20 * t182 + (-t75 * t184 + (t24 * t79 + t52 * t80) * t72) * qJD(2)) * pkin(1)) * pkin(2)) * t194) * t83) * t160;
t7 = [-0.2e1 * m(4) * t133 + 0.2e1 * (-mrSges(4,1) * t12 - mrSges(4,2) * t11) * t150 + (-mrSges(4,1) * t4 + mrSges(4,2) * t3) * t167 - 0.2e1 * t3 * t11 * Ifges(4,1) + 0.2e1 * t4 * Ifges(4,2) * t12 + (-0.2e1 * Ifges(3,4) * t79 + t231 * t80) * t170 + (0.2e1 * Ifges(3,4) * t80 + t231 * t79) * t169 + (t21 * t154 - 0.2e1 * t106 - t125 / 0.2e1) * (-t207 + t210) * t154 + (-t20 * t154 + 0.2e1 * t107 + t126 / 0.2e1) * (-t203 + t206) * t154 + 0.2e1 * (-t11 * t4 + t3 * t12) * Ifges(4,4) + (-0.2e1 * pkin(1) * ((-t53 * mrSges(5,2) / 0.2e1 - t52 * mrSges(5,1) / 0.2e1) * t159 + (-mrSges(5,1) * t104 + mrSges(5,2) * t103) * t27) + (t53 * ((-t210 / 0.2e1 + t207 / 0.2e1) * t159 + (Ifges(5,1) * t103 + Ifges(5,4) * t104) * t27) - t52 * ((-t206 / 0.2e1 + t203 / 0.2e1) * t159 + (Ifges(5,4) * t103 + Ifges(5,2) * t104) * t27)) * t232 * t83) * t83; Ifges(3,5) * t169 - Ifges(3,6) * t170 + ((-t201 + t204) * t132 + t52 * t120 + t53 * t121) * t1 + (-Ifges(5,5) * t106 + ((-t204 / 0.2e1 + t201 / 0.2e1) * t159 + (Ifges(5,5) * t103 + Ifges(5,6) * t104) * t27) * t83 / 0.2e1 + Ifges(5,6) * t107 + t20 * t120 + t21 * t121 - Ifges(5,5) * t125 / 0.4e1 + Ifges(5,6) * t126 / 0.4e1) * t5 + (-Ifges(4,5) * t11 + Ifges(4,6) * t12) * t2 + t6 * (Ifges(4,5) * t3 + Ifges(4,6) * t4) + ((t12 * t227 + 0.2e1 * t233) * t86 * t114 + ((t233 / 0.2e1 + t12 * t54 / 0.2e1) * t211 + (-t19 * t11 - t22 * t12 + t51 * t3 - t54 * t4) * t29) * pkin(2) * t183) * mrSges(4,3); 0.2e1 * t6 * Ifges(4,3) * t2 + 0.2e1 * t5 * Ifges(5,3) * t1 + m(4) * (-t176 * t8 / t34 ^ 2 * t164 + (-0.8e1 * t108 + 0.4e1 * t111) * t228) * t89 * t87 / 0.2e1 + (0.4e1 * t112 * t114 + (t112 * t211 + 0.2e1 * ((t2 * t54 + t22 * t6) * mrSges(4,2) + (-t19 * t6 - t2 * t51) * mrSges(4,1)) * t29) * t71 * pkin(2)) * t86;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
