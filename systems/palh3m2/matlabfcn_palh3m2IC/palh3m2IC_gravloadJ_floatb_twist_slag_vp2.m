% Calculate Gravitation load on the joints for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 05:00:10
% EndTime: 2020-05-07 05:00:11
% DurationCPUTime: 1.56s
% Computational Cost: add. (767->160), mult. (624->153), div. (14->8), fcn. (492->39), ass. (0->91)
t59 = cos(qJ(5));
t133 = mrSges(6,1) * t59;
t53 = qJ(2) + qJ(3);
t50 = qJ(4) + t53;
t42 = sin(t50);
t43 = cos(t50);
t171 = -mrSges(5,2) * t43 + (-pkin(8) * m(6) - mrSges(5,1) - t133) * t42;
t52 = qJ(2) + qJ(7);
t45 = sin(t52);
t47 = cos(t52);
t111 = pkin(15) - qJ(7);
t49 = -qJ(2) + t111;
t101 = (-qJ(8) + t111);
t44 = -qJ(2) + t101;
t35 = sin(t44);
t36 = cos(t44);
t87 = t36 * mrSges(9,1) + t35 * mrSges(9,2);
t169 = -t47 * mrSges(8,1) + t45 * mrSges(8,2) - m(9) * pkin(3) * cos(t49) + t87;
t61 = cos(qJ(1));
t152 = t171 * t61;
t46 = sin(t53);
t48 = cos(t53);
t157 = mrSges(4,1) * t46 + mrSges(4,2) * t48;
t167 = -t157 * t61 + t152;
t57 = sin(qJ(1));
t151 = t171 * t57;
t166 = -t157 * t57 + t151;
t150 = g(1) * t61 + g(2) * t57;
t100 = -t48 * mrSges(4,1) + t46 * mrSges(4,2);
t56 = sin(qJ(2));
t60 = cos(qJ(2));
t99 = -t43 * mrSges(5,1) + t42 * mrSges(5,2);
t165 = -t60 * mrSges(3,1) + t56 * mrSges(3,2) - t100 + t169 - t99;
t142 = pkin(4) * t48;
t95 = pkin(8) * t43 + pkin(10) * t42;
t164 = m(6) * (-t95 - t142);
t161 = -m(4) - m(8) - m(9);
t132 = mrSges(9,1) * t35;
t141 = pkin(10) * t43;
t144 = pkin(3) * sin(t49);
t145 = pkin(1) * t56;
t143 = pkin(4) * t46;
t17 = t143 - t145;
t55 = sin(qJ(5));
t129 = mrSges(6,2) * t55;
t105 = t42 * t129;
t85 = -t43 * mrSges(6,3) - t105;
t160 = m(4) * t145 - m(5) * t17 - m(6) * (t17 - t141) - t85 - m(9) * (t144 - t145) + t132;
t112 = qJ(2) - qJ(6);
t107 = (qJ(4) + pkin(14));
t97 = (qJ(3) + t107);
t94 = pkin(15) + t97;
t82 = t94 - t112;
t69 = -(2 * qJ(7)) - pkin(16) + t82;
t83 = t94 + t112;
t71 = pkin(16) + t83;
t155 = -cos(qJ(8) - t69) + cos(qJ(8) - t71);
t127 = t42 * mrSges(6,3);
t84 = -t127 + (t129 - t133) * t43;
t154 = -t84 - t99;
t153 = -m(5) * t143 - m(6) * (-t141 + t143) - t85;
t149 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t51 = t60 * pkin(1);
t104 = t51 - t142;
t14 = pkin(12) + t104;
t54 = sin(qJ(6));
t58 = cos(qJ(6));
t90 = t58 * mrSges(7,1) - t54 * mrSges(7,2);
t148 = pkin(6) * m(7) - m(3) * pkin(12) - m(5) * t14 + m(6) * (-t14 + t95) - mrSges(2,1) + t127 - t90 + t161 * (t51 + pkin(12)) + t165;
t70 = -t105 + (-m(6) * pkin(10) - mrSges(6,3)) * t43;
t137 = ((m(6) * t95 + t154) * g(3) + (-t57 * t70 + t151) * g(2) + (-t61 * t70 + t152) * g(1)) / pkin(9);
t128 = mrSges(9,2) * t36;
t15 = t57 * t128;
t16 = t61 * t128;
t136 = (g(3) * t87 - g(1) * (-t132 * t61 + t16) - g(2) * (-t132 * t57 + t15)) / pkin(7);
t121 = t55 * t57;
t120 = t55 * t61;
t119 = t57 * t59;
t118 = t59 * t61;
t110 = qJ(7) + pkin(16);
t96 = t110 + t112;
t88 = -mrSges(8,1) * t45 - mrSges(8,2) * t47;
t74 = -qJ(7) + t82;
t73 = -qJ(7) + t83;
t65 = 0.1e1 / pkin(2);
t34 = sin(t96);
t7 = t118 * t43 - t121;
t6 = t120 * t43 + t119;
t5 = t119 * t43 + t120;
t4 = t121 * t43 - t118;
t1 = [(t7 * mrSges(6,1) - t6 * mrSges(6,2) + t148 * t61 + t149 * t57) * g(2) + (-t5 * mrSges(6,1) + t4 * mrSges(6,2) - t148 * t57 + t149 * t61) * g(1); (t160 * t57 - t15 + t166) * g(2) + (t160 * t61 - t16 + t167) * g(1) + (-m(5) * t104 - t164 - t84 + (-m(6) + t161) * t51 + t165) * g(3) + (pkin(1) * sin(t110) / pkin(5) * (-g(3) * t90 + t150 * (mrSges(7,1) * t54 + mrSges(7,2) * t58)) + (-pkin(2) * t34 - pkin(1) * sin(t112)) * t65 * (-g(1) * t16 - g(2) * t15 + t169 * g(3) + t150 * (-m(9) * t144 + t132 - t88))) / t34 + (pkin(3) * (pkin(1) * (cos(-qJ(8) + t112) - cos(qJ(8) + t112)) + (cos(-qJ(8) + t96) - cos(qJ(8) + t96)) * pkin(2)) * t137 + (((cos(t74) - cos(t73)) * pkin(1) + (cos(t69) - cos(t71)) * pkin(2)) * pkin(3) + ((-cos(-qJ(8) + t74) + cos(-qJ(8) + t73)) * pkin(1) + t155 * pkin(2)) * pkin(7)) * t136) / t155 * t65 + t150 * (m(8) * t145 + mrSges(3,1) * t56 + mrSges(3,2) * t60 - t88); -pkin(9) * t137 + (-(cos((2 * qJ(3) + t107)) - cos((2 * qJ(7)) - (2 * pkin(15)) + 0.2e1 * qJ(8) + t107)) / (cos((2 * t97)) - cos((2 * t101))) * t137 + sin(t107) / sin(-qJ(7) - qJ(8) + t94) * t136) * pkin(4) + (m(5) * t142 - t100 + t154 - t164) * g(3) + (t153 * t57 + t166) * g(2) + (t153 * t61 + t167) * g(1); -g(1) * (mrSges(6,1) * t6 + mrSges(6,2) * t7) - g(2) * (mrSges(6,1) * t4 + mrSges(6,2) * t5) - g(3) * (mrSges(6,1) * t55 + mrSges(6,2) * t59) * t42;];
taug = t1(:);
