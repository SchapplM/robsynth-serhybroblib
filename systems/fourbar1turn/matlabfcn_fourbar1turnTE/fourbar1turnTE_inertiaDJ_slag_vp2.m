% Calculate time derivative of joint inertia matrix for
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnTE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:35
% EndTime: 2020-04-12 19:18:45
% DurationCPUTime: 3.01s
% Computational Cost: add. (30796->177), mult. (43267->408), div. (1326->15), fcn. (11657->4), ass. (0->151)
t68 = cos(qJ(2));
t176 = pkin(2) * t68;
t141 = -0.2e1 * t176;
t77 = pkin(1) ^ 2;
t146 = pkin(1) * t141 + t77;
t76 = pkin(2) ^ 2;
t62 = t76 + t146;
t59 = 0.1e1 / t62;
t71 = 0.1e1 / pkin(4);
t157 = t59 * t71;
t60 = 0.1e1 / t62 ^ 2;
t177 = pkin(2) * t60;
t138 = pkin(1) * t177;
t67 = sin(qJ(2));
t122 = t67 * t138;
t108 = -0.2e1 * t122;
t64 = pkin(1) * t68 - pkin(2);
t193 = -0.2e1 * t64;
t142 = pkin(2) * t193;
t183 = -pkin(3) - pkin(4);
t55 = (pkin(2) - t183) * (pkin(2) + t183) + t146;
t182 = pkin(4) - pkin(3);
t56 = (pkin(2) - t182) * (pkin(2) + t182) + t146;
t159 = t55 * t56;
t78 = sqrt(-t159);
t149 = t68 * t78;
t112 = pkin(1) * pkin(2) * (-t55 - t56);
t45 = t67 * t112;
t49 = 0.1e1 / t78;
t162 = t49 * t45;
t198 = pkin(3) ^ 2;
t145 = -pkin(4) ^ 2 + t198;
t57 = t62 + t145;
t179 = pkin(1) * t57;
t54 = t67 * t179;
t17 = t54 + (-t149 + (t142 - t162) * t67) * pkin(1);
t152 = t67 * t78;
t39 = -pkin(1) * t152 - t57 * t64;
t73 = 0.1e1 / pkin(3);
t10 = (t108 * t39 + t17 * t59) * t73;
t66 = t67 ^ 2;
t153 = t66 * t77;
t160 = t49 * t64;
t19 = -t45 * t160 + 0.2e1 * pkin(2) * t153 + (t57 * t68 + t152) * pkin(1);
t42 = -t64 * t78 + t54;
t12 = (t108 * t42 + t19 * t59) * t73;
t80 = t39 ^ 2;
t32 = 0.1e1 / t80;
t165 = t32 * t42;
t31 = 0.1e1 / t39;
t199 = -t10 * t165 + t12 * t31;
t132 = t76 * t66 * pkin(1);
t114 = qJD(2) * t132;
t144 = qJD(2) * t67;
t124 = t78 * t144;
t143 = qJD(2) * t68;
t127 = pkin(2) * t143;
t58 = t62 - t145;
t148 = pkin(2) * t124 + t58 * t127;
t63 = pkin(1) - t176;
t161 = t49 * t63;
t44 = qJD(2) * t45;
t15 = t161 * t44 + 0.2e1 * t114 + t148;
t40 = -pkin(2) * t152 + t58 * t63;
t35 = 0.1e1 / t40 ^ 2;
t53 = pkin(2) * t67 * t58;
t41 = t63 * t78 + t53;
t164 = t35 * t41;
t163 = t44 * t49;
t129 = t67 * t163;
t139 = 0.2e1 * t63 * pkin(1);
t14 = (-t129 + (-t149 + (t58 + t139) * t67) * qJD(2)) * pkin(2);
t34 = 0.1e1 / t40;
t166 = t14 * t34 * t35;
t37 = t41 ^ 2;
t29 = t35 * t37 + 0.1e1;
t197 = 0.1e1 / t29 ^ 2 * (t15 * t164 - t166 * t37);
t126 = qJD(2) * t153;
t113 = pkin(2) * t126;
t147 = pkin(1) * t124 + t143 * t179;
t16 = -t160 * t44 + 0.2e1 * t113 + t147;
t13 = (-t129 + (-t149 + (t57 + t142) * t67) * qJD(2)) * pkin(1);
t167 = t13 * t31 * t32;
t38 = t42 ^ 2;
t30 = t32 * t38 + 0.1e1;
t196 = 0.1e1 / t30 ^ 2 * (t16 * t165 - t167 * t38);
t195 = Ifges(3,1) - Ifges(3,2);
t190 = -t40 / 0.2e1;
t189 = t41 / 0.2e1;
t188 = -t59 / 0.2e1;
t187 = t59 / 0.2e1;
t186 = -t67 / 0.2e1;
t185 = t67 / 0.2e1;
t184 = -t68 / 0.2e1;
t178 = pkin(1) * t60;
t175 = pkin(3) * t62;
t174 = pkin(4) * t62;
t20 = t45 * t161 + 0.2e1 * t132 + (t58 * t68 + t152) * pkin(2);
t11 = (-t122 * t41 + t187 * t20) * t71;
t169 = t11 * t34;
t158 = t59 * t68;
t155 = t59 * t73;
t154 = t60 * t67;
t151 = t68 * t39;
t150 = t68 * t42;
t140 = 0.2e1 * t174;
t137 = 0.2e1 * t60;
t27 = 0.1e1 / t30;
t136 = t27 * t175;
t25 = 0.1e1 / t29;
t134 = t25 * t174;
t18 = t53 + (-t149 + (t139 - t162) * t67) * pkin(2);
t9 = (t122 * t40 + t18 * t188) * t71;
t133 = t9 * t164;
t130 = t59 * t163;
t128 = pkin(2) * t144;
t125 = t76 * t144;
t123 = t157 / 0.2e1;
t121 = mrSges(4,3) * pkin(2) * t155;
t117 = pkin(1) * t128;
t116 = 0.4e1 / t159 * t44 * t162;
t115 = pkin(1) * t125;
t111 = t68 * t125;
t110 = Ifges(5,5) * t157 / 0.4e1;
t109 = -Ifges(5,6) * t157 / 0.4e1;
t105 = t60 * t117;
t61 = t59 * t60;
t102 = t61 * t76 * t126;
t101 = -t116 / 0.4e1;
t99 = 0.2e1 * t115;
t98 = t40 * t105;
t97 = t41 * t105;
t96 = 0.8e1 * t102;
t95 = -0.2e1 * pkin(4) * t25 * t117;
t94 = t150 / 0.2e1 + t39 * t185;
t93 = t71 * t98;
t92 = t71 * t97;
t90 = (t42 * t185 - t151 / 0.2e1) * t59;
t89 = t14 * t187 - t98;
t88 = t15 * t187 - t97;
t43 = (t112 * t68 - 0.4e1 * t153 * t76) * qJD(2);
t87 = (t67 * t101 + 0.2e1 * (t43 * t186 - t68 * t44) * t49) * t59;
t22 = t73 * t90;
t21 = t94 * t155;
t6 = ((t13 * t184 + t16 * t185) * t59 + (t94 * t59 + (t151 * t67 - t42 * t66) * t138) * qJD(2)) * t73;
t5 = ((t13 * t186 + t16 * t184) * t59 + (t90 + (t150 * t67 + t39 * t66) * t138) * qJD(2)) * t73;
t4 = t199 * t136 + 0.1e1;
t3 = 0.2e1 * (-t133 - t169) * t134;
t2 = 0.2e1 * t199 * pkin(3) * t27 * t117 + (-t10 * t16 - t12 * t13 - ((0.4e1 * t113 + t147) * t59 + t39 * t96 + (t87 + (-0.2e1 * t13 * t154 + (t158 * t193 + (-t17 * t67 - t151) * t137) * qJD(2)) * pkin(2)) * pkin(1)) * t42 * t73) * t32 * t136 + ((((0.6e1 * t127 * t67 * t77 + t101 * t64 - t160 * t43) * t59 + t42 * t96 + ((-0.2e1 * t16 * t177 + t130) * t67 + ((t149 + (-t57 + t162) * t67) * t59 + (-t19 * t67 - t150) * pkin(2) * t137) * qJD(2)) * pkin(1)) * t73 * t27 - 0.2e1 * t12 * t196) * t31 + 0.2e1 * (t27 * t167 + t32 * t196) * t10 * t42) * t175;
t1 = 0.2e1 * t95 * t133 + 0.2e1 * (t140 * t197 + t95) * t169 + 0.2e1 * (t25 * t166 + t35 * t197) * t41 * t9 * t140 + 0.2e1 * ((t11 * t14 - t9 * t15) * t35 + (-((t63 * t116 / 0.4e1 + t43 * t161 + 0.6e1 * pkin(1) * t111) * t187 + 0.4e1 * t41 * t102 + ((t130 / 0.2e1 - t15 * t178) * t67 + ((t149 + (-t58 + t162) * t67) * t187 + (-t20 * t67 - t41 * t68) * t178) * qJD(2)) * pkin(2)) * t34 - ((0.4e1 * t114 + t148) * t188 - 0.4e1 * t40 * t102 + (-t87 / 0.2e1 + (t14 * t154 + (-t63 * t158 + (t18 * t67 + t40 * t68) * t60) * qJD(2)) * pkin(1)) * pkin(2)) * t164) * t71) * t134;
t7 = [-0.2e1 * m(4) * t111 + 0.2e1 * (-mrSges(4,1) * t22 - mrSges(4,2) * t21) * t128 + (-mrSges(4,1) * t6 + mrSges(4,2) * t5) * t141 + 0.2e1 * t6 * Ifges(4,2) * t22 - 0.2e1 * t5 * t21 * Ifges(4,1) + (-0.2e1 * Ifges(3,4) * t67 + t195 * t68) * t144 + (0.2e1 * Ifges(3,4) * t68 + t195 * t67) * t143 + (t15 * t123 - t92) * (Ifges(5,1) * t189 + Ifges(5,4) * t190) * t157 + (-t14 * t157 / 0.2e1 + t93) * (Ifges(5,4) * t189 + Ifges(5,2) * t190) * t157 + 0.2e1 * (-t6 * t21 + t22 * t5) * Ifges(4,4) + (-0.2e1 * pkin(1) * (mrSges(5,1) * t89 + mrSges(5,2) * t88) + ((Ifges(5,1) * t88 - Ifges(5,4) * t89) * t189 + (Ifges(5,4) * t88 - Ifges(5,2) * t89) * t190) * t157) * t71; Ifges(3,5) * t143 - Ifges(3,6) * t144 + (t39 * t21 + t42 * t22) * mrSges(4,3) * t60 * t73 * t115 - (t16 * t22 + t42 * t6) * t121 / 0.2e1 + (-t13 * t21 + t39 * t5) * t121 / 0.2e1 + (t41 * t110 + t40 * t109 + (Ifges(5,5) * t189 + Ifges(5,6) * t190) * t123) * t1 + (t15 * t110 - Ifges(5,5) * t92 / 0.2e1 + t14 * t109 + Ifges(5,6) * t93 / 0.2e1 + (Ifges(5,5) * t88 - Ifges(5,6) * t89) * t71 / 0.2e1) * t3 + (-Ifges(4,5) * t21 + Ifges(4,6) * t22) * t2 + t4 * (Ifges(4,5) * t5 + Ifges(4,6) * t6); 0.2e1 * t4 * Ifges(4,3) * t2 + 0.2e1 * t3 * Ifges(5,3) * t1 + m(4) * ((t13 * t39 + t16 * t42) * t76 * t60 + (-t38 - t80) * pkin(2) * t61 * t99) / t198 / 0.2e1 + ((mrSges(4,1) * t39 - mrSges(4,2) * t42) * t60 * t4 * t99 + ((t16 * t4 + t2 * t42) * mrSges(4,2) + (-t13 * t4 - t2 * t39) * mrSges(4,1)) * t59 * pkin(2)) * t73;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
