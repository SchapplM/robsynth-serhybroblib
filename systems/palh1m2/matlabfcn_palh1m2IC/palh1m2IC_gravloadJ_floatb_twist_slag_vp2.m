% Calculate Gravitation load on the joints for
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:47:29
% EndTime: 2020-05-02 23:47:30
% DurationCPUTime: 1.03s
% Computational Cost: add. (967->147), mult. (905->163), div. (19->10), fcn. (601->45), ass. (0->84)
t63 = sin(qJ(5));
t71 = cos(qJ(5));
t93 = mrSges(6,1) * t71 - t63 * mrSges(6,2);
t24 = m(6) * pkin(9) + mrSges(5,1) + t93;
t52 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t64 = sin(qJ(4));
t72 = cos(qJ(4));
t101 = t24 * t72 + t52 * t64;
t128 = m(5) + m(6);
t10 = pkin(5) * t128 + mrSges(4,1) + t101;
t16 = -t24 * t64 + t52 * t72;
t14 = -mrSges(4,2) + t16;
t65 = sin(qJ(3));
t73 = cos(qJ(3));
t6 = t10 * t73 + t14 * t65;
t107 = pkin(19) - qJ(7);
t53 = -qJ(2) + t107;
t42 = pkin(4) * sin(t53);
t66 = sin(qJ(2));
t98 = -qJ(10) + t107;
t49 = -qJ(2) + t98;
t43 = sin(t49);
t44 = cos(t49);
t99 = -t43 * mrSges(11,1) + t44 * mrSges(11,2);
t137 = m(11) * (-pkin(1) * t66 + t42) + t99;
t5 = -t10 * t65 + t14 * t73;
t136 = mrSges(11,1) * t44 + mrSges(11,2) * t43;
t109 = qJ(6) - qJ(2);
t106 = qJ(4) + pkin(18);
t94 = pkin(19) + qJ(3) + t106;
t87 = t94 + t109;
t82 = -(2 * qJ(7)) - pkin(20) + t87;
t88 = t94 - t109;
t83 = pkin(20) + t88;
t135 = cos(qJ(10) - t82) + cos(qJ(10) - t83);
t67 = sin(qJ(1));
t75 = cos(qJ(1));
t96 = g(1) * t75 + g(2) * t67;
t134 = m(8) + m(4) + t128;
t34 = mrSges(6,1) * t63 + mrSges(6,2) * t71;
t133 = -mrSges(11,3) - mrSges(10,3) - t34 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t62 = sin(qJ(6));
t70 = cos(qJ(6));
t92 = -mrSges(7,1) * t70 + mrSges(7,2) * t62;
t132 = -pkin(14) * m(7) + mrSges(2,1) - t92 + (m(10) + m(9) + m(3) + t134 + m(11)) * pkin(15) + t137;
t130 = pkin(2) * m(10);
t126 = pkin(4) * cos(t53);
t7 = t101 * t73 + t16 * t65;
t74 = cos(qJ(2));
t8 = -t101 * t65 + t16 * t73;
t122 = ((-g(3) * t8 + t7 * t96) * t66 + (-g(3) * t7 - t8 * t96) * t74) / pkin(10);
t113 = t136 * t75;
t114 = t136 * t67;
t121 = (-g(1) * t113 - g(2) * t114 - g(3) * t99) / pkin(8);
t120 = m(11) * (-pkin(1) * t74 - t126);
t59 = qJ(2) + qJ(8);
t108 = sin(t59) * t130;
t115 = t96 * cos(t59) * t130 + g(3) * t108;
t105 = qJ(7) + pkin(20);
t104 = m(11) * t126;
t55 = qJ(9) + t59;
t50 = sin(t55);
t51 = cos(t55);
t102 = -(g(3) * mrSges(10,1) - mrSges(10,2) * t96) * t50 + (-mrSges(10,1) * t96 - g(3) * mrSges(10,2)) * t51;
t100 = -qJ(8) + pkin(17) + qJ(3);
t97 = t105 - t109;
t61 = sin(qJ(7));
t69 = cos(qJ(7));
t36 = -t69 * mrSges(8,1) + mrSges(8,2) * t61;
t33 = -mrSges(8,1) * t61 - mrSges(8,2) * t69;
t60 = sin(qJ(8));
t68 = cos(qJ(8));
t35 = -t68 * mrSges(9,1) + t60 * mrSges(9,2);
t32 = -mrSges(9,1) * t60 - mrSges(9,2) * t68;
t11 = (-t64 * t73 - t65 * t72) * t34;
t12 = (t64 * t65 - t72 * t73) * t34;
t90 = t11 * t66 - t12 * t74;
t85 = -qJ(7) + t88;
t84 = -qJ(7) + t87;
t79 = 0.1e1 / pkin(3);
t46 = cos(t97);
t3 = -mrSges(3,2) + t32 + t33 + t6;
t2 = t134 * pkin(1) + mrSges(3,1) - t35 - t36 - t5;
t1 = [(t133 * g(1) - t132 * g(2)) * t75 + (t132 * g(1) + t133 * g(2)) * t67 + (mrSges(10,1) * t50 + mrSges(10,2) * t51 - t2 * t66 + t3 * t74 - t108) * (g(1) * t67 - g(2) * t75); (-t3 * g(3) + t96 * t2) * t74 - (-t2 * g(3) - t96 * t3) * t66 - g(1) * (t75 * t120 + t113) - g(2) * (t67 * t120 + t114) - g(3) * t137 + t102 + t115 + (pkin(1) * sin(t105) / pkin(7) * (g(3) * t92 + t96 * (mrSges(7,1) * t62 + mrSges(7,2) * t70)) + (-pkin(3) * t46 - pkin(1) * cos(-t109)) * t79 * ((-g(3) * t36 + t33 * t96) * t66 + (-g(3) * t33 - t36 * t96) * t74 - g(1) * (-t104 * t75 + t113) - g(2) * (-t104 * t67 + t114) - g(3) * (m(11) * t42 + t99))) / t46 + (-(pkin(1) * (sin(qJ(10) - t109) + sin(qJ(10) + t109)) + (-sin(-qJ(10) + t97) + sin(qJ(10) + t97)) * pkin(3)) * pkin(4) * t122 + ((-(cos(t85) + cos(t84)) * pkin(1) + (-cos(t82) - cos(t83)) * pkin(3)) * pkin(4) + ((cos(-qJ(10) + t85) + cos(-qJ(10) + t84)) * pkin(1) + t135 * pkin(3)) * pkin(8)) * t121) / t135 * t79; (-t5 * g(3) + t96 * t6) * t66 + (-t6 * g(3) - t96 * t5) * t74 - pkin(10) * t122 + (pkin(2) * cos(t100) / pkin(12) * t102 + ((-g(3) * t35 + t32 * t96) * t66 + (-g(3) * t32 - t35 * t96) * t74 + t115) * cos(-qJ(9) + t100)) / pkin(2) / sin(qJ(9)) * pkin(6) + (-cos(qJ(3) + t98) * t122 + sin(t106) * t121) / cos(-qJ(7) - qJ(10) + t94) * pkin(5); (-t11 * t74 - t12 * t66) * g(3) + (t67 * t90 + t75 * t93) * g(2) + (-t67 * t93 + t75 * t90) * g(1);];
taug = t1(:);
