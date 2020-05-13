% Calculate joint inertia matrix for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_inertiaJ_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:28
% EndTime: 2020-05-01 20:57:32
% DurationCPUTime: 1.98s
% Computational Cost: add. (643->267), mult. (678->273), div. (0->0), fcn. (258->100), ass. (0->110)
t73 = pkin(22) + pkin(21);
t57 = pkin(18) - t73;
t45 = (-pkin(20) + t57);
t128 = m(11) * pkin(4);
t98 = m(5) + m(6);
t132 = pkin(9) * mrSges(6,2);
t90 = sin(pkin(18));
t131 = pkin(9) * t90;
t94 = cos(pkin(18));
t130 = pkin(9) * t94;
t129 = (pkin(11) * mrSges(6,3));
t91 = cos(qJ(4));
t127 = mrSges(6,1) * t91;
t87 = sin(qJ(4));
t126 = mrSges(6,2) * t87;
t89 = sin(qJ(2));
t92 = cos(qJ(3));
t125 = t89 * t92;
t93 = cos(qJ(2));
t124 = t92 * t93;
t67 = pkin(11) * mrSges(6,1) - Ifges(6,5);
t66 = pkin(11) * mrSges(6,2) - Ifges(6,6);
t56 = -m(4) - m(8) - m(11) - t98;
t123 = pkin(1) ^ 2 * t56;
t72 = t98 * pkin(5) ^ 2;
t122 = qJ(4) - qJ(2);
t76 = qJ(4) + qJ(2);
t75 = pkin(18) - pkin(22);
t74 = qJ(3) + pkin(19);
t119 = pkin(11) * m(6) + mrSges(6,3);
t88 = sin(qJ(3));
t51 = pkin(5) * t88 + pkin(1);
t5 = -pkin(5) * t124 + t51 * t89;
t118 = 2 * t45;
t99 = pkin(2) ^ 2 * m(10);
t117 = Ifges(11,3) + Ifges(4,3) + Ifges(9,3) + t72 + t99;
t116 = -qJ(2) + t57;
t115 = qJ(2) + t57;
t114 = mrSges(6,1) * t131 - t67 * t94;
t43 = pkin(5) * t98 + mrSges(11,1) + mrSges(4,1);
t85 = mrSges(11,2) + mrSges(4,2);
t113 = (t43 * t88 + t85 * t92) * pkin(1);
t36 = -qJ(2) + t45;
t35 = qJ(2) + t45;
t34 = -qJ(4) + t45;
t33 = qJ(4) + t45;
t29 = -qJ(3) + t36;
t26 = qJ(3) + t35;
t61 = qJ(2) + t74;
t46 = sin(t61);
t47 = cos(t61);
t78 = qJ(3) + qJ(2);
t64 = sin(t78);
t65 = cos(t78);
t68 = qJ(3) + t76;
t69 = qJ(3) - t122;
t112 = t65 * (Ifges(11,6) + Ifges(4,6)) + Ifges(9,6) * t47 + (-pkin(2) * mrSges(10,3) + Ifges(9,5)) * t46 + (Ifges(11,5) + Ifges(4,5)) * t64 + (-mrSges(5,3) * t64 + (cos(t68) / 0.2e1 - cos(t69) / 0.2e1) * mrSges(6,1) - (sin(t68) + sin(t69)) * mrSges(6,2) / 0.2e1) * pkin(5);
t107 = pkin(9) ^ 2;
t106 = pkin(11) ^ 2;
t104 = 0.2e1 * qJ(2);
t103 = 0.2e1 * qJ(4);
t97 = pkin(9) * mrSges(6,1);
t83 = cos(pkin(19));
t82 = cos(pkin(20));
t81 = sin(pkin(19));
t80 = sin(pkin(20));
t79 = t104 + qJ(3);
t71 = pkin(9) * m(6) + mrSges(5,1);
t70 = qJ(2) + pkin(17) - pkin(18);
t63 = -qJ(2) + t75;
t62 = qJ(2) + t75;
t60 = t104 + t74;
t59 = cos(t74);
t58 = sin(t74);
t55 = -mrSges(5,2) + t119;
t54 = cos(t73);
t53 = sin(t73);
t52 = 0.2e1 * t78;
t50 = cos(t70);
t49 = sin(t70);
t48 = 0.2e1 * t75;
t44 = 0.2e1 * t70;
t42 = 0.2e1 * t61;
t41 = qJ(3) + t115;
t40 = -qJ(3) + t116;
t38 = -qJ(4) + t118;
t37 = qJ(4) + t118;
t31 = -qJ(2) + t34;
t30 = -qJ(2) + t33;
t28 = qJ(2) + t34;
t27 = qJ(2) + t33;
t25 = 2 * t45;
t24 = mrSges(6,1) * t87 + mrSges(6,2) * t91;
t19 = -qJ(4) + t29;
t18 = qJ(4) + t29;
t17 = -qJ(4) + t26;
t16 = qJ(4) + t26;
t15 = 0.2e1 * t34;
t14 = 0.2e1 * t33;
t13 = t88 * t93 + t125;
t12 = -t88 * t89 + t124;
t11 = -mrSges(6,2) * t131 + t66 * t94;
t10 = mrSges(6,1) * t130 + t67 * t90;
t9 = mrSges(6,2) * t130 + t66 * t90;
t6 = pkin(5) * t125 + t51 * t93;
t4 = -t12 * t90 + t13 * t94;
t3 = t12 * t94 + t13 * t90;
t2 = t5 * t90 + t6 * t94;
t1 = t5 * t94 - t6 * t90;
t7 = [(-sin(t103) / 0.2e1 + sin(t14) / 0.4e1 - sin(t15) / 0.4e1) * Ifges(6,4) + pkin(9) * t127 + ((-cos(t60) - t59) * mrSges(10,2) + (-sin(t60) + t58) * mrSges(10,1)) * pkin(2) + (pkin(9) * t119 + Ifges(5,4)) * sin(t25) + Ifges(8,4) * sin(t48) + (t97 - t66) * cos(t38) / 0.2e1 + (t97 + t66) * cos(t37) / 0.2e1 - pkin(9) * t126 + ((-sin(t40) + sin(t41)) * mrSges(11,2) + (-cos(t40) - cos(t41)) * mrSges(11,1) + (cos(0.2e1 * t57) / 0.2e1 + 0.1e1 / 0.2e1) * t128) * pkin(4) + (pkin(14) * m(7) + 0.2e1 * mrSges(7,1) * t49 + 0.2e1 * mrSges(7,2) * t50) * pkin(14) + (-cos(t103) / 0.4e1 + cos(t14) / 0.8e1 + cos(t15) / 0.8e1) * (Ifges(6,2) - Ifges(6,1)) + (t106 + t107) * m(6) / 0.2e1 + (-Ifges(9,1) + Ifges(9,2) + t99) * cos(t42) / 0.2e1 + (-Ifges(7,2) + Ifges(7,1)) * cos(t44) / 0.2e1 + (Ifges(8,2) - Ifges(8,1)) * cos(t48) / 0.2e1 - t123 / 0.2e1 + (0.2e1 * pkin(15) * t56 * t89 + (-cos(t79) + t92) * t85 + (sin(t35) - sin(t36)) * t71 + (-cos(t35) + cos(t36)) * t55 + (-sin(t79) + t88) * t43 + (cos(t62) - cos(t63)) * mrSges(8,2) + (sin(t62) - sin(t63)) * mrSges(8,1) + (sin(t115) - sin(t116)) * t128 + (cos(t27) / 0.2e1 - cos(t28) / 0.2e1 - cos(t30) / 0.2e1 + cos(t31) / 0.2e1) * mrSges(6,2) + (sin(t27) / 0.2e1 + sin(t28) / 0.2e1 - sin(t30) / 0.2e1 - sin(t31) / 0.2e1) * mrSges(6,1)) * pkin(1) + (0.2e1 * mrSges(8,2) * sin(t75) - 0.2e1 * t85 * t64 - 0.2e1 * (mrSges(3,2) + mrSges(10,2)) * t93 - 0.2e1 * t55 * sin(t45) - 0.2e1 * t71 * cos(t45) - 0.2e1 * mrSges(8,1) * cos(t75) - 0.2e1 * (mrSges(3,1) + mrSges(10,1)) * t89 - 0.2e1 * cos(t57) * t128 + 0.2e1 * t43 * t65 - 0.2e1 * mrSges(9,2) * t46 + 0.2e1 * (pkin(2) * m(10) + mrSges(9,1)) * t47 + (m(9) + m(3) + m(10) - t56) * pkin(15) + (sin(t33) - sin(t34)) * mrSges(6,2) + (-cos(t33) - cos(t34)) * mrSges(6,1)) * pkin(15) + (t67 - t132) * sin(t37) / 0.2e1 + (t67 + t132) * sin(t38) / 0.2e1 + Ifges(9,4) * sin(t42) - Ifges(7,4) * sin(t44) + (-Ifges(11,1) - Ifges(4,1) + Ifges(11,2) + Ifges(4,2) + t72) * cos(t52) / 0.2e1 + (Ifges(3,1) + Ifges(10,1) - Ifges(3,2) - Ifges(10,2) + t123) * cos(t104) / 0.2e1 + ((-cos(t26) - cos(t29)) * t71 + (-sin(t26) - sin(t29)) * t55 + (sin(t16) / 0.2e1 - sin(t17) / 0.2e1 + sin(t18) / 0.2e1 - sin(t19) / 0.2e1) * mrSges(6,2) + (-cos(t16) / 0.2e1 - cos(t17) / 0.2e1 - cos(t18) / 0.2e1 - cos(t19) / 0.2e1) * mrSges(6,1)) * pkin(5) + (0.2e1 * (-t106 + t107) * m(6) - (4 * t129) - (2 * Ifges(5,1)) - Ifges(6,1) + (2 * Ifges(5,2)) - Ifges(6,2) + (2 * Ifges(6,3))) * cos(t25) / 0.4e1 + t99 / 0.2e1 + t72 / 0.2e1 + (-Ifges(3,4) - Ifges(10,4)) * sin(t104) + (Ifges(11,4) + Ifges(4,4)) * sin(t52) + t129 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.2e1 + Ifges(8,1) / 0.2e1 + Ifges(9,1) / 0.2e1 + Ifges(10,1) / 0.2e1 + Ifges(11,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.2e1 + Ifges(8,2) / 0.2e1 + Ifges(9,2) / 0.2e1 + Ifges(10,2) / 0.2e1 + Ifges(11,2) / 0.2e1 + Ifges(2,3) + Ifges(6,3) / 0.2e1; Ifges(7,5) * t50 - Ifges(7,6) * t49 + (Ifges(3,5) + Ifges(10,5)) * t93 - t89 * (Ifges(3,6) + Ifges(10,6)) + ((-mrSges(11,3) - mrSges(4,3) - mrSges(5,3) - mrSges(8,3)) * t93 + (-cos(-t122) / 0.2e1 - cos(t76) / 0.2e1) * mrSges(6,2) + (sin(-t122) / 0.2e1 - sin(t76) / 0.2e1) * mrSges(6,1)) * pkin(1) + t112; -t123 + Ifges(3,3) + Ifges(7,3) + Ifges(10,3) + t117 + 0.2e1 * (mrSges(10,1) * t58 - mrSges(10,2) * t59) * pkin(2) + 0.2e1 * t113; t112; t113 + ((mrSges(10,1) * t83 + mrSges(10,2) * t81) * t88 + (mrSges(10,1) * t81 - mrSges(10,2) * t83) * t92) * pkin(2) + t117; t117; ((t10 * t82 + t114 * t80) * t91 + (t11 * t80 - t82 * t9) * t87 + Ifges(6,3) * (t80 * t90 + t82 * t94)) * t54 + ((-t80 * t10 + t114 * t82) * t91 + (t11 * t82 + t80 * t9) * t87 - Ifges(6,3) * (t80 * t94 - t82 * t90)) * t53 + (-pkin(15) + t5) * (-t126 + t127); -t24 * ((-t1 * t80 + t2 * t82) * t54 - (t1 * t82 + t2 * t80) * t53); -((t3 * t80 + t4 * t82) * t54 + t53 * (t3 * t82 - t4 * t80)) * t24 * pkin(5); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
