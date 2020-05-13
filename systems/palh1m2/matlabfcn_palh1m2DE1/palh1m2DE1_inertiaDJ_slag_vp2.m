% Calculate time derivative of joint inertia matrix for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:38
% EndTime: 2020-05-01 20:57:45
% DurationCPUTime: 1.32s
% Computational Cost: add. (684->227), mult. (872->272), div. (0->0), fcn. (410->90), ass. (0->121)
t79 = pkin(22) + pkin(21);
t122 = (pkin(18) - t79);
t117 = (-pkin(20) + t122);
t126 = 2 * pkin(15);
t80 = qJD(3) + qJD(2);
t83 = qJD(2) - qJD(4);
t142 = t83 / 0.2e1;
t107 = m(5) + m(6);
t66 = -qJD(4) + t80;
t141 = pkin(5) * t66;
t140 = pkin(9) * mrSges(6,2);
t99 = sin(pkin(18));
t139 = pkin(9) * t99;
t112 = 2 * qJ(2);
t138 = sin(t112);
t103 = cos(pkin(18));
t137 = pkin(9) * t103;
t94 = mrSges(11,2) + mrSges(4,2);
t97 = sin(qJ(3));
t136 = t94 * t97;
t101 = cos(qJ(3));
t53 = pkin(5) * t107 + mrSges(11,1) + mrSges(4,1);
t134 = t101 * t53;
t98 = sin(qJ(2));
t133 = t101 * t98;
t74 = pkin(11) * mrSges(6,1) - Ifges(6,5);
t73 = pkin(11) * mrSges(6,2) - Ifges(6,6);
t132 = qJD(2) * t98;
t102 = cos(qJ(2));
t131 = t101 * t102;
t130 = qJ(4) - qJ(2);
t85 = qJ(4) + qJ(2);
t129 = qJD(2) * t102;
t127 = pkin(18) - pkin(22);
t84 = qJ(3) + pkin(19);
t125 = (pkin(1) * t134 + (mrSges(10,1) * cos(t84) + mrSges(10,2) * sin(t84)) * pkin(2)) * qJD(3);
t65 = qJD(4) + t80;
t124 = -pkin(5) * t65 / 0.2e1;
t123 = qJD(3) * t136;
t121 = 2 * t117;
t59 = pkin(5) * t97 + pkin(1);
t15 = -pkin(5) * t131 + t59 * t98;
t120 = -qJ(2) + t122;
t119 = qJ(2) + t122;
t100 = cos(qJ(4));
t96 = sin(qJ(4));
t34 = mrSges(6,1) * t96 + mrSges(6,2) * t100;
t118 = mrSges(6,1) * t100 - mrSges(6,2) * t96;
t25 = t102 * t97 + t133;
t24 = -t97 * t98 + t131;
t116 = mrSges(6,1) * t139 - t103 * t74;
t47 = -qJ(2) + t117;
t46 = qJ(2) + t117;
t45 = -qJ(4) + t117;
t44 = qJ(4) + t117;
t115 = t24 * qJD(3);
t38 = -qJ(3) + t47;
t35 = qJ(3) + t46;
t68 = qJ(2) + t84;
t55 = sin(t68);
t56 = cos(t68);
t87 = qJ(3) + qJ(2);
t71 = sin(t87);
t72 = cos(t87);
t75 = qJ(3) + t85;
t76 = qJ(3) - t130;
t114 = (cos(t75) * t124 - cos(t76) * t141 / 0.2e1) * mrSges(6,2) + (sin(t75) * t124 + sin(t76) * t141 / 0.2e1) * mrSges(6,1) + ((-pkin(2) * mrSges(10,3) + Ifges(9,5)) * t56 + (-pkin(5) * mrSges(5,3) + Ifges(11,5) + Ifges(4,5)) * t72 - Ifges(9,6) * t55 - t71 * (Ifges(11,6) + Ifges(4,6))) * t80;
t111 = 2 * qJ(4);
t106 = pkin(9) * mrSges(6,1);
t92 = cos(pkin(19));
t91 = cos(pkin(20));
t90 = sin(pkin(19));
t89 = sin(pkin(20));
t88 = t112 + qJ(3);
t82 = qJD(2) + qJD(4);
t81 = 0.2e1 * qJD(2) + qJD(3);
t78 = pkin(9) * m(6) + mrSges(5,1);
t77 = qJ(2) + pkin(17) - pkin(18);
t70 = -qJ(2) + t127;
t69 = qJ(2) + t127;
t67 = t112 + t84;
t63 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t62 = cos(t79);
t61 = sin(t79);
t60 = 0.2e1 * t87;
t58 = cos(t77);
t57 = sin(t77);
t54 = 0.2e1 * t77;
t52 = 0.2e1 * t68;
t51 = qJ(3) + t119;
t50 = -qJ(3) + t120;
t49 = -qJ(4) + t121;
t48 = qJ(4) + t121;
t40 = -qJ(2) + t45;
t39 = -qJ(2) + t44;
t37 = qJ(2) + t45;
t36 = qJ(2) + t44;
t33 = -qJ(4) + t38;
t32 = qJ(4) + t38;
t31 = -qJ(4) + t35;
t30 = qJ(4) + t35;
t29 = 2 * t45;
t28 = 2 * t44;
t27 = t118 * qJD(4);
t23 = -mrSges(6,2) * t139 + t103 * t73;
t22 = mrSges(6,1) * t137 + t74 * t99;
t21 = mrSges(6,2) * t137 + t73 * t99;
t16 = pkin(5) * t133 + t102 * t59;
t12 = t80 * t25;
t11 = t24 * qJD(2) + t115;
t10 = t103 * t25 - t24 * t99;
t9 = t103 * t24 + t25 * t99;
t8 = t59 * t129 + (t25 * qJD(3) + t101 * t132) * pkin(5);
t7 = -t59 * t132 + (t101 * t129 + t115) * pkin(5);
t6 = t103 * t16 + t15 * t99;
t5 = t103 * t15 - t16 * t99;
t4 = -t103 * t12 + t11 * t99;
t3 = t103 * t11 + t12 * t99;
t2 = t103 * t8 - t7 * t99;
t1 = t103 * t7 + t8 * t99;
t13 = [(-mrSges(10,1) * cos(t67) + mrSges(10,2) * sin(t67)) * t81 * pkin(2) + (((-cos(t31) / 0.2e1 - cos(t32) / 0.2e1) * t66 + (cos(t30) / 0.2e1 + cos(t33) / 0.2e1) * t65) * mrSges(6,2) + ((sin(t31) / 0.2e1 - sin(t32) / 0.2e1) * t66 + (sin(t30) / 0.2e1 - sin(t33) / 0.2e1) * t65) * mrSges(6,1)) * pkin(5) + (-t123 + (-t53 * cos(t88) + t94 * sin(t88)) * t81 + ((sin(t37) / 0.2e1 - sin(t39) / 0.2e1) * t83 + (-sin(t36) / 0.2e1 + sin(t40) / 0.2e1) * t82) * mrSges(6,2) + ((cos(t37) / 0.2e1 + cos(t39) / 0.2e1) * t83 + (cos(t36) / 0.2e1 + cos(t40) / 0.2e1) * t82) * mrSges(6,1)) * pkin(1) + ((-t94 * t72 - (pkin(2) * m(10) + mrSges(9,1)) * t55 - mrSges(9,2) * t56 - t53 * t71) * t126 + ((-sin(t38) + sin(t35)) * t78 + (cos(t38) - cos(t35)) * t63) * pkin(5) + ((cos(t50) + cos(t51)) * mrSges(11,2) + (-sin(t50) + sin(t51)) * mrSges(11,1)) * pkin(4) - (pkin(2) ^ 2 * m(10) - Ifges(9,1) + Ifges(9,2)) * sin(t52) + 0.2e1 * (Ifges(11,4) + Ifges(4,4)) * cos(t60) + 0.2e1 * Ifges(9,4) * cos(t52) - (pkin(5) ^ 2 * t107 - Ifges(11,1) - Ifges(4,1) + Ifges(11,2) + Ifges(4,2)) * sin(t60)) * t80 + (-t34 * pkin(9) + (-sin(t28) / 0.4e1 + sin(t29) / 0.4e1 + sin(t111) / 0.2e1) * (Ifges(6,2) - Ifges(6,1)) + (-cos(t111) + cos(t28) / 0.2e1 + cos(t29) / 0.2e1) * Ifges(6,4) + ((cos(t44) + cos(t45)) * mrSges(6,2) + (sin(t44) - sin(t45)) * mrSges(6,1)) * pkin(15) + (t74 - t140) * cos(t48) / 0.2e1 - (t74 + t140) * cos(t49) / 0.2e1 - (t106 + t73) * sin(t48) / 0.2e1 + (t106 - t73) * sin(t49) / 0.2e1) * qJD(4) + (-(Ifges(3,1) + Ifges(10,1) - Ifges(3,2) - Ifges(10,2)) * t138 - 0.2e1 * Ifges(7,4) * cos(t54) - (-Ifges(7,2) + Ifges(7,1)) * sin(t54) + 0.2e1 * (-Ifges(3,4) - Ifges(10,4)) * cos(t112) + ((mrSges(3,2) + mrSges(10,2)) * t98 - (mrSges(3,1) + mrSges(10,1)) * t102) * t126 + 0.2e1 * (mrSges(7,1) * t58 - mrSges(7,2) * t57) * pkin(14) + ((cos(t46) + cos(t47)) * t78 + (-t138 * pkin(1) + t102 * t126) * (-m(4) - m(8) - m(11) - t107) + (sin(t46) + sin(t47)) * t63 + (-sin(t69) - sin(t70)) * mrSges(8,2) + (cos(t69) + cos(t70)) * mrSges(8,1) + (cos(t119) + cos(t120)) * m(11) * pkin(4)) * pkin(1)) * qJD(2) + t125; ((sin(-t130) * t142 + t82 * sin(t85) / 0.2e1) * mrSges(6,2) + (cos(-t130) * t142 - t82 * cos(t85) / 0.2e1) * mrSges(6,1)) * pkin(1) + (-Ifges(7,5) * t57 - Ifges(7,6) * t58 - ((-mrSges(11,3) - mrSges(4,3) - mrSges(5,3) - mrSges(8,3)) * pkin(1) + Ifges(3,5) + Ifges(10,5)) * t98 - t102 * (Ifges(3,6) + Ifges(10,6))) * qJD(2) + t114; -0.2e1 * pkin(1) * t123 + 0.2e1 * t125; t114; ((t134 - t136) * pkin(1) + ((mrSges(10,1) * t92 + mrSges(10,2) * t90) * t101 - (mrSges(10,1) * t90 - mrSges(10,2) * t92) * t97) * pkin(2)) * qJD(3); 0; t8 * t118 + ((-(t89 * t116 + t22 * t91) * t96 + (-t21 * t91 + t23 * t89) * t100) * t62 + (-(t116 * t91 - t89 * t22) * t96 + (t21 * t89 + t23 * t91) * t100) * t61 - (-pkin(15) + t15) * t34) * qJD(4); -t27 * ((-t5 * t89 + t6 * t91) * t62 - (t5 * t91 + t6 * t89) * t61) - t34 * ((t1 * t91 - t2 * t89) * t62 - (t1 * t89 + t2 * t91) * t61); (-((t3 * t91 + t4 * t89) * t62 + t61 * (-t3 * t89 + t4 * t91)) * t34 - ((t10 * t91 + t89 * t9) * t62 + t61 * (-t10 * t89 + t9 * t91)) * t27) * pkin(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq = res;
