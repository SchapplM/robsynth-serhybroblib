% Calculate time derivative of joint inertia matrix for
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
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:29
% EndTime: 2020-05-02 20:58:39
% DurationCPUTime: 2.28s
% Computational Cost: add. (4191->195), mult. (8273->350), div. (0->0), fcn. (10472->22), ass. (0->123)
t124 = sin(pkin(20));
t82 = cos(pkin(20));
t87 = sin(pkin(18));
t91 = cos(pkin(18));
t54 = -t91 * t124 + t82 * t87;
t58 = t87 * t124 + t91 * t82;
t85 = sin(qJ(3));
t89 = cos(qJ(3));
t38 = t54 * t89 - t58 * t85;
t41 = t54 * t85 + t58 * t89;
t86 = sin(qJ(2));
t90 = cos(qJ(2));
t105 = t38 * t86 + t41 * t90;
t165 = -t38 * t90 + t41 * t86;
t74 = pkin(22) + pkin(21);
t69 = sin(t74);
t70 = cos(t74);
t12 = t69 * t105 + t165 * t70;
t31 = t41 * qJD(3);
t167 = t165 * qJD(2) + t86 * t31;
t34 = t38 * qJD(3);
t166 = t105 * qJD(2) + t86 * t34;
t61 = -t85 * t86 + t89 * t90;
t164 = -0.2e1 * t61;
t78 = qJ(2) + qJ(3);
t71 = sin(t78);
t163 = 0.2e1 * t71;
t162 = (mrSges(11,2) + mrSges(4,2)) * t85;
t13 = -t105 * t70 + t165 * t69;
t123 = qJD(3) * t90;
t132 = t69 * (t34 * t90 - t167);
t7 = (t41 * t123 + t166) * t70 + t132;
t14 = -t31 * t90 - t166;
t8 = (-t38 * t123 + t167) * t70 - t14 * t69;
t2 = t7 * pkin(5) + (t7 * t85 - t8 * t89 + (t12 * t89 + t13 * t85) * qJD(3)) * pkin(1);
t116 = pkin(1) * t85 + pkin(5);
t146 = pkin(1) * t89;
t5 = t12 * t116 - t13 * t146;
t161 = t12 * t2 + t5 * t7;
t67 = t86 * pkin(1) - pkin(15);
t159 = 0.2e1 * t67;
t84 = sin(qJ(4));
t88 = cos(qJ(4));
t127 = t84 ^ 2 + t88 ^ 2;
t158 = mrSges(11,1) + mrSges(4,1);
t75 = qJD(2) + qJD(3);
t81 = sin(pkin(19));
t83 = cos(pkin(19));
t53 = t81 * t89 + t83 * t85;
t51 = t53 * qJD(3);
t55 = -t85 * t81 + t83 * t89;
t52 = t55 * qJD(3);
t156 = (mrSges(10,1) * t52 + mrSges(10,2) * t51) * pkin(2);
t154 = 2 * m(5);
t153 = 2 * m(6);
t152 = -0.2e1 * pkin(15);
t125 = cos(pkin(22));
t79 = sin(pkin(22));
t56 = t125 * t91 + t87 * t79;
t57 = -t87 * t125 + t91 * t79;
t151 = 0.2e1 * (cos(pkin(21)) * t56 - t57 * sin(pkin(21))) * pkin(4) + 0.2e1 * t67;
t60 = t85 * t90 + t86 * t89;
t150 = 0.2e1 * t60;
t149 = t2 * t5;
t126 = pkin(1) * qJD(3);
t117 = t89 * t126;
t96 = pkin(1) * t12;
t1 = (t14 * t70 - t132) * t146 + qJD(3) * t85 * t96 - t8 * t116 - t13 * t117;
t4 = -t13 * t116 - t89 * t96;
t148 = t4 * t1;
t145 = pkin(2) * mrSges(10,3);
t143 = sin(pkin(17));
t27 = -t54 * t70 + t58 * t69;
t140 = mrSges(6,3) * t27;
t139 = t27 * t88;
t134 = t51 * t86;
t133 = t53 * t86;
t128 = -Ifges(10,4) - Ifges(3,4);
t122 = qJD(4) * t84;
t121 = qJD(4) * t88;
t120 = 0.2e1 * t86;
t119 = 0.2e1 * t88;
t93 = pkin(5) ^ 2;
t118 = t12 * t7 * t93;
t109 = mrSges(6,1) * t84 + mrSges(6,2) * t88;
t115 = (mrSges(5,3) + t109) * t27;
t112 = -pkin(5) * t61 + t67;
t111 = -t1 * t13 - t4 * t8;
t110 = mrSges(6,1) * t88 - mrSges(6,2) * t84;
t28 = -t54 * t69 - t58 * t70;
t19 = t28 * mrSges(6,2) - t84 * t140;
t20 = -t28 * mrSges(6,1) - mrSges(6,3) * t139;
t108 = t88 * t19 - t84 * t20;
t107 = -t19 * t84 - t20 * t88;
t104 = -t52 * t90 + t134;
t37 = t53 * t90 + t55 * t86;
t40 = t55 * t90 - t133;
t103 = t56 * t90 - t57 * t86;
t102 = -t56 * t86 - t57 * t90;
t92 = cos(pkin(17));
t59 = -t91 * t143 + t87 * t92;
t62 = t143 * t87 + t92 * t91;
t43 = t59 * t90 - t62 * t86;
t44 = t59 * t86 + t62 * t90;
t99 = mrSges(5,3) * t28 + t108;
t23 = -t37 * qJD(2) - t51 * t90 - t52 * t86;
t24 = t40 * qJD(2) - t104;
t47 = t75 * t61;
t48 = t75 * t60;
t72 = cos(t78);
t98 = qJD(2) * t133 * t145 + Ifges(4,5) * t47 + Ifges(9,5) * t24 - Ifges(4,6) * t48 + Ifges(9,6) * t23 + (Ifges(11,5) * t72 - Ifges(11,6) * t71) * t75;
t97 = qJD(4) * (-t127 * t140 + t107);
t36 = t43 * qJD(2);
t35 = t44 * qJD(2);
t33 = t102 * qJD(2);
t32 = t103 * qJD(2);
t30 = -pkin(2) * t40 - pkin(15);
t29 = pkin(1) * qJD(2) * t90 + t48 * pkin(5);
t25 = 0.2e1 * m(10) * (-t51 * t55 + t52 * t53) * pkin(2) ^ 2;
t17 = t110 * t27 * qJD(4);
t16 = -t28 * pkin(9) - t27 * pkin(11) + t112;
t3 = t13 * t93 * t8;
t6 = [t47 * Ifges(4,1) * t150 + t48 * Ifges(4,2) * t164 + (mrSges(4,1) * t48 + mrSges(4,2) * t47) * t159 + 0.2e1 * t24 * t37 * Ifges(9,1) - 0.2e1 * t43 * Ifges(7,2) * t35 + 0.2e1 * t36 * t44 * Ifges(7,1) + 0.2e1 * pkin(14) * (mrSges(7,1) * t35 + mrSges(7,2) * t36) + (-mrSges(9,1) * t23 + mrSges(9,2) * t24) * t152 + 0.2e1 * (t23 * t37 + t24 * t40) * Ifges(9,4) + 0.2e1 * (-t35 * t44 + t43 * t36) * Ifges(7,4) + 0.2e1 * (t47 * t61 - t60 * t48) * Ifges(4,4) + 0.2e1 * (m(6) * t127 * t16 + m(5) * t112 - mrSges(5,1) * t28 + mrSges(5,2) * t27 - t107) * t29 + ((mrSges(11,1) * t71 + mrSges(11,2) * t72) * t151 + 0.2e1 * (-t71 ^ 2 + t72 ^ 2) * Ifges(11,4) + (-Ifges(11,2) + Ifges(11,1)) * t72 * t163) * t75 + (0.2e1 * t108 * t16 + ((-Ifges(6,4) * t139 + Ifges(6,6) * t28) * t119 + (0.2e1 * Ifges(6,5) * t28 + (0.2e1 * Ifges(6,4) * t84 + (-Ifges(6,1) + Ifges(6,2)) * t119) * t27) * t84) * t27) * qJD(4) + ((0.2e1 * pkin(15) * mrSges(3,2) - 0.2e1 * t30 * mrSges(10,2) - t128 * t120) * t86 + (mrSges(3,1) * t152 + 0.2e1 * t30 * mrSges(10,1) + (Ifges(3,2) + Ifges(10,2) - Ifges(10,1) - Ifges(3,1)) * t120 + 0.2e1 * t128 * t90 + (m(11) * t151 - 0.2e1 * mrSges(11,1) * t72 + mrSges(4,1) * t164 + 0.2e1 * mrSges(8,1) * t56 + mrSges(11,2) * t163 + mrSges(4,2) * t150 + 0.2e1 * mrSges(8,2) * t57 + (m(8) + m(4)) * t159) * pkin(1)) * t90) * qJD(2) + 0.2e1 * (Ifges(9,2) * t40 + (-m(10) * t30 - mrSges(10,1) * t86 - mrSges(10,2) * t90) * pkin(2)) * t23; Ifges(7,5) * t36 - Ifges(7,6) * t35 + t4 * t17 + t115 * t1 + t104 * t145 + t99 * t2 + t5 * t97 + ((-Ifges(3,5) - Ifges(10,5)) * t86 + (-t55 * t145 - Ifges(3,6) - Ifges(10,6)) * t90) * qJD(2) + ((-t32 * t57 - t33 * t56) * mrSges(8,3) + (-t47 * t85 + t48 * t89 + (-t60 * t89 + t61 * t85) * qJD(3)) * mrSges(4,3) + (-qJD(3) + t75) * mrSges(11,3) * (t71 * t89 - t72 * t85)) * pkin(1) + t98; (t127 * t149 + t148) * t153 + (t148 + t149) * t154 + t25 + 0.2e1 * m(8) * (-t102 * t32 + t103 * t33) * pkin(1) ^ 2 + 0.2e1 * t156 - 0.2e1 * t126 * t162 + 0.2e1 * t158 * t117; (t134 + (-qJD(2) * t55 - t52) * t90) * t145 + (-t115 * t8 + t12 * t97 - t13 * t17 + t99 * t7) * pkin(5) + t98; t25 + t156 + (t158 * t89 - t162) * t126 + 0.2e1 * (m(6) * (t161 * t127 + t111) / 0.2e1 + m(5) * (t111 + t161) / 0.2e1) * pkin(5); t25 + (t127 * t118 + t3) * t153 + (t3 + t118) * t154; t110 * t29 + ((-Ifges(6,5) * t84 - Ifges(6,6) * t88) * t27 - t109 * t16) * qJD(4); (t5 * t122 - t2 * t88) * mrSges(6,2) + (-t5 * t121 - t2 * t84) * mrSges(6,1); ((t12 * t122 - t7 * t88) * mrSges(6,2) + (-t12 * t121 - t7 * t84) * mrSges(6,1)) * pkin(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
