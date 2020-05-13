% Calculate time derivative of joint inertia matrix for
% palh1m2TE
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
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2TE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:15:48
% EndTime: 2020-05-01 20:15:52
% DurationCPUTime: 1.16s
% Computational Cost: add. (666->220), mult. (709->263), div. (0->0), fcn. (374->90), ass. (0->118)
t75 = pkin(22) + pkin(21);
t116 = (pkin(18) - t75);
t111 = (-pkin(20) + t116);
t74 = qJD(2) + qJD(3) / 0.2e1;
t132 = 0.2e1 * t74;
t91 = sin(qJ(4));
t95 = cos(qJ(4));
t30 = mrSges(6,1) * t95 - mrSges(6,2) * t91;
t131 = qJD(4) * t30;
t76 = qJD(3) + qJD(2);
t78 = qJD(2) - qJD(4);
t130 = t78 / 0.2e1;
t102 = m(5) + m(6);
t62 = -qJD(4) + t76;
t129 = pkin(5) * t62;
t92 = sin(qJ(3));
t128 = pkin(5) * t92;
t96 = cos(qJ(3));
t127 = pkin(5) * t96;
t126 = pkin(9) * mrSges(6,2);
t29 = mrSges(6,1) * t91 + mrSges(6,2) * t95;
t125 = pkin(9) * t29;
t107 = 2 * qJ(2);
t124 = sin(t107);
t97 = cos(qJ(2));
t123 = t96 * t97;
t122 = pkin(11) * mrSges(6,2) - Ifges(6,6);
t121 = pkin(11) * mrSges(6,1) - Ifges(6,5);
t80 = qJ(3) + qJ(2);
t120 = (pkin(18) - pkin(22));
t79 = qJ(3) + pkin(19);
t119 = 2 * pkin(15);
t118 = t76 * t127;
t61 = qJD(4) + t76;
t117 = -pkin(5) * t61 / 0.2e1;
t115 = 2 * t111;
t114 = -qJ(2) + t116;
t113 = qJ(2) + t116;
t46 = pkin(5) * t102 + mrSges(11,1) + mrSges(4,1);
t89 = mrSges(11,2) + mrSges(4,2);
t112 = t46 * t96 - t89 * t92;
t93 = sin(qJ(2));
t21 = t92 * t97 + t93 * t96;
t20 = -t92 * t93 + t123;
t110 = t112 * pkin(1);
t42 = -qJ(2) + t111;
t41 = qJ(2) + t111;
t38 = -qJ(4) + t111;
t37 = qJ(4) + t111;
t32 = -qJ(3) + t42;
t31 = qJ(3) + t41;
t65 = qJ(2) + t79;
t48 = sin(t65);
t49 = cos(t65);
t67 = sin(t80);
t68 = cos(t80);
t70 = qJ(4) + t80;
t71 = -qJ(4) + t80;
t109 = (-Ifges(9,6) * t48 - t49 * (pkin(2) * mrSges(10,3) - Ifges(9,5)) - (pkin(5) * mrSges(5,3) - Ifges(11,5) - Ifges(4,5)) * t68 - t67 * (Ifges(11,6) + Ifges(4,6))) * t76 + (-cos(t71) * t129 / 0.2e1 + cos(t70) * t117) * mrSges(6,2) + (sin(t71) * t129 / 0.2e1 + sin(t70) * t117) * mrSges(6,1);
t106 = 0.2e1 * qJ(4);
t101 = pkin(9) * mrSges(6,1);
t98 = cos(pkin(18));
t94 = sin(pkin(18));
t87 = cos(pkin(19));
t86 = cos(pkin(20));
t85 = sin(pkin(19));
t84 = sin(pkin(20));
t83 = qJ(2) - qJ(4);
t82 = qJ(2) + qJ(4);
t81 = t107 + qJ(3);
t77 = qJD(2) + qJD(4);
t73 = pkin(9) * m(6) + mrSges(5,1);
t72 = pkin(17) + qJ(2) - pkin(18);
t66 = t107 + t79;
t64 = -qJ(2) + t120;
t63 = qJ(2) + t120;
t60 = cos(t79);
t59 = sin(t79);
t57 = -pkin(11) * m(6) + mrSges(5,2) - mrSges(6,3);
t56 = cos(t75);
t55 = sin(t75);
t53 = 0.2e1 * t80;
t52 = pkin(1) + t128;
t51 = cos(t72);
t50 = sin(t72);
t47 = 0.2e1 * t72;
t45 = 0.2e1 * t65;
t44 = qJ(3) + t113;
t43 = -qJ(3) + t114;
t40 = -qJ(4) + t115;
t39 = qJ(4) + t115;
t36 = -qJ(2) + t38;
t35 = -qJ(2) + t37;
t34 = qJ(2) + t38;
t33 = qJ(2) + t37;
t28 = -qJ(4) + t32;
t27 = qJ(4) + t32;
t26 = -qJ(4) + t31;
t25 = qJ(4) + t31;
t24 = 0.2e1 * t38;
t23 = 0.2e1 * t37;
t22 = pkin(1) * qJD(2) + t76 * t128;
t19 = t84 * t94 + t86 * t98;
t18 = t84 * t98 - t86 * t94;
t13 = t91 * t121 + t122 * t95;
t12 = t76 * t21;
t11 = t76 * t20;
t10 = -t20 * t94 + t21 * t98;
t9 = t20 * t98 + t21 * t94;
t8 = t84 * t125 + t13 * t86;
t7 = t86 * t125 - t13 * t84;
t6 = -t19 * t127 + t18 * t52;
t5 = t18 * t127 + t19 * t52;
t4 = -t19 * t118 + t18 * t22;
t3 = t18 * t118 + t19 * t22;
t2 = t11 * t94 - t12 * t98;
t1 = t11 * t98 + t12 * t94;
t14 = [((sin(t66) * t132 + qJD(3) * t59) * mrSges(10,2) + (-0.2e1 * cos(t66) * t74 + qJD(3) * t60) * mrSges(10,1)) * pkin(2) + (((-cos(t26) / 0.2e1 - cos(t27) / 0.2e1) * t62 + (cos(t25) / 0.2e1 + cos(t28) / 0.2e1) * t61) * mrSges(6,2) + ((sin(t26) / 0.2e1 - sin(t27) / 0.2e1) * t62 + (sin(t25) / 0.2e1 - sin(t28) / 0.2e1) * t61) * mrSges(6,1)) * pkin(5) + ((-t46 * cos(t81) + t89 * sin(t81)) * t132 + t112 * qJD(3) + ((sin(t34) / 0.2e1 - sin(t35) / 0.2e1) * t78 + (-sin(t33) / 0.2e1 + sin(t36) / 0.2e1) * t77) * mrSges(6,2) + ((cos(t34) / 0.2e1 + cos(t35) / 0.2e1) * t78 + (cos(t33) / 0.2e1 + cos(t36) / 0.2e1) * t77) * mrSges(6,1)) * pkin(1) + (-t125 + (t121 - t126) * cos(t39) / 0.2e1 - (t121 + t126) * cos(t40) / 0.2e1 - (t101 + t122) * sin(t39) / 0.2e1 + (t101 - t122) * sin(t40) / 0.2e1 + (sin(t23) / 0.4e1 - sin(t24) / 0.4e1 - sin(t106) / 0.2e1) * (-Ifges(6,2) + Ifges(6,1)) + (cos(t23) / 0.2e1 + cos(t24) / 0.2e1 - cos(t106)) * Ifges(6,4) + ((cos(t37) + cos(t38)) * mrSges(6,2) + (sin(t37) - sin(t38)) * mrSges(6,1)) * pkin(15)) * qJD(4) + (0.2e1 * Ifges(9,4) * cos(t45) - (pkin(2) ^ 2 * m(10) - Ifges(9,1) + Ifges(9,2)) * sin(t45) + 0.2e1 * (Ifges(4,4) + Ifges(11,4)) * cos(t53) + (-pkin(5) ^ 2 * t102 + Ifges(11,1) + Ifges(4,1) - Ifges(11,2) - Ifges(4,2)) * sin(t53) + (sin(t31) - sin(t32)) * pkin(5) * t73 + (cos(t31) - cos(t32)) * pkin(5) * t57 + (-(pkin(2) * m(10) + mrSges(9,1)) * t48 - mrSges(9,2) * t49 - t46 * t67 - t89 * t68) * t119 + ((cos(t43) + cos(t44)) * mrSges(11,2) + (sin(t44) - sin(t43)) * mrSges(11,1)) * pkin(4)) * t76 + (-(Ifges(7,1) - Ifges(7,2)) * sin(t47) - 0.2e1 * Ifges(7,4) * cos(t47) - 0.2e1 * (Ifges(3,4) + Ifges(10,4)) * cos(t107) + (-Ifges(3,1) - Ifges(10,1) + Ifges(3,2) + Ifges(10,2)) * t124 + ((mrSges(3,2) + mrSges(10,2)) * t93 - (mrSges(3,1) + mrSges(10,1)) * t97) * t119 + 0.2e1 * (mrSges(7,1) * t51 - mrSges(7,2) * t50) * pkin(14) + ((cos(t41) + cos(t42)) * t73 + (pkin(1) * t124 - 0.2e1 * pkin(15) * t97) * (m(11) + m(4) + m(8) + t102) + (-sin(t41) - sin(t42)) * t57 + (-sin(t63) - sin(t64)) * mrSges(8,2) + (cos(t63) + cos(t64)) * mrSges(8,1) + (cos(t114) + cos(t113)) * m(11) * pkin(4)) * pkin(1)) * qJD(2); (-Ifges(7,6) * t51 - Ifges(7,5) * t50 + (-Ifges(3,5) - Ifges(10,5)) * t93 - t97 * (Ifges(3,6) + Ifges(10,6))) * qJD(2) + ((mrSges(4,3) + mrSges(5,3) + mrSges(8,3) + mrSges(11,3)) * t93 * qJD(2) + (sin(t83) * t130 + t77 * sin(t82) / 0.2e1) * mrSges(6,2) + (cos(t83) * t130 - t77 * cos(t82) / 0.2e1) * mrSges(6,1)) * pkin(1) + t109; 0.2e1 * qJD(3) * (t110 + (mrSges(10,1) * t60 + mrSges(10,2) * t59) * pkin(2)); t109; qJD(3) * (t110 + ((mrSges(10,1) * t87 + mrSges(10,2) * t85) * t96 - t92 * (mrSges(10,1) * t85 - mrSges(10,2) * t87)) * pkin(2)); 0; t30 * (t93 * t118 + t22 * t97) + (-t56 * (t7 * t98 + t8 * t94) + t55 * (-t7 * t94 + t8 * t98) - (-pkin(5) * t123 + t52 * t93 - pkin(15)) * t29) * qJD(4); -((t5 * t97 - t6 * t93) * t56 - t55 * (t5 * t93 + t6 * t97)) * t131 + ((t3 * t93 + t4 * t97) * t56 + (t3 * t97 - t4 * t93) * t55) * t29; (-t29 * ((t1 * t86 + t2 * t84) * t56 + t55 * (-t1 * t84 + t2 * t86)) - ((t10 * t86 + t84 * t9) * t56 + t55 * (-t10 * t84 + t86 * t9)) * t131) * pkin(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t14(1), t14(2), t14(4), t14(7); t14(2), t14(3), t14(5), t14(8); t14(4), t14(5), t14(6), t14(9); t14(7), t14(8), t14(9), t14(10);];
Mq = res;
