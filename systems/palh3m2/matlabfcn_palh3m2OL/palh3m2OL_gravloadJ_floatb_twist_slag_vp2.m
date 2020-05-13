% Calculate Gravitation load on the joints for
% palh3m2OL
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
% taug [10x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:42
% EndTime: 2020-05-07 04:33:49
% DurationCPUTime: 0.85s
% Computational Cost: add. (532->109), mult. (509->123), div. (0->0), fcn. (412->18), ass. (0->70)
t43 = qJ(2) + qJ(7);
t36 = sin(t43);
t38 = cos(t43);
t40 = pkin(15) - t43;
t35 = -qJ(8) + t40;
t26 = sin(t35);
t27 = cos(t35);
t68 = t27 * mrSges(9,1) + t26 * mrSges(9,2);
t118 = -m(9) * pkin(3) * cos(t40) - t38 * mrSges(8,1) + t36 * mrSges(8,2) + t68;
t47 = sin(qJ(2));
t51 = cos(qJ(2));
t136 = -t51 * mrSges(3,1) + t47 * mrSges(3,2) + t118;
t44 = qJ(2) + qJ(3);
t41 = qJ(4) + t44;
t33 = sin(t41);
t34 = cos(t41);
t50 = cos(qJ(5));
t91 = t50 * mrSges(6,1);
t135 = -mrSges(5,2) * t34 + (-pkin(8) * m(6) - mrSges(5,1) - t91) * t33;
t133 = -m(4) - m(8) - m(9);
t130 = t34 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t33;
t37 = sin(t44);
t39 = cos(t44);
t129 = mrSges(4,1) * t37 + mrSges(4,2) * t39;
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t122 = g(1) * t52 + g(2) * t48;
t78 = -t39 * mrSges(4,1) + t37 * mrSges(4,2);
t110 = pkin(10) * t34;
t112 = pkin(4) * t37;
t46 = sin(qJ(5));
t102 = mrSges(6,2) * t46;
t82 = t33 * t102;
t66 = -t34 * mrSges(6,3) - t82;
t127 = -m(6) * (-t110 + t112) - t66 - m(5) * t112;
t126 = -(t102 - t91) * t34 + t130;
t124 = t135 * t52;
t123 = t135 * t48;
t121 = -t129 * t52 + t124;
t120 = -t129 * t48 + t123;
t105 = mrSges(9,1) * t26;
t113 = pkin(3) * sin(t40);
t114 = pkin(1) * t47;
t12 = t112 - t114;
t119 = -m(5) * t12 - m(6) * (t12 - t110) - t66 - m(9) * (t113 - t114) + t105 + m(4) * t114;
t111 = pkin(4) * t39;
t75 = pkin(8) * t34 + pkin(10) * t33;
t117 = -m(6) * (-t75 - t111) + t126 - t78;
t116 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t42 = t51 * pkin(1);
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t71 = t49 * mrSges(7,1) - t45 * mrSges(7,2);
t79 = t42 - t111;
t9 = pkin(12) + t79;
t115 = pkin(6) * m(7) - m(3) * pkin(12) + m(6) * (t75 - t9) - mrSges(2,1) - m(5) * t9 - t71 - t78 + t130 + t133 * (t42 + pkin(12)) + t136;
t101 = mrSges(9,2) * t27;
t94 = t46 * t52;
t93 = t48 * t46;
t92 = t48 * t50;
t90 = t50 * t52;
t69 = -mrSges(8,1) * t36 - mrSges(8,2) * t38;
t56 = -t82 + (-m(6) * pkin(10) - mrSges(6,3)) * t34;
t11 = t52 * t101;
t10 = t48 * t101;
t4 = t34 * t90 - t93;
t3 = t34 * t94 + t92;
t2 = t34 * t92 + t94;
t1 = t34 * t93 - t90;
t5 = [(t4 * mrSges(6,1) - t3 * mrSges(6,2) + t115 * t52 + t116 * t48) * g(2) + (-t2 * mrSges(6,1) + t1 * mrSges(6,2) - t115 * t48 + t116 * t52) * g(1), (t119 * t48 - t10 + t120) * g(2) + (t119 * t52 - t11 + t121) * g(1) + (-m(5) * t79 + (-m(6) + t133) * t42 + t117 + t136) * g(3) + t122 * (m(8) * t114 + mrSges(3,1) * t47 + mrSges(3,2) * t51 - t69), (m(5) * t111 + t117) * g(3) + (t127 * t48 + t120) * g(2) + (t127 * t52 + t121) * g(1), (m(6) * t75 + t126) * g(3) + (-t56 * t48 + t123) * g(2) + (-t56 * t52 + t124) * g(1), -g(1) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(2) * (mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (mrSges(6,1) * t46 + mrSges(6,2) * t50) * t33, -g(3) * t71 + t122 * (mrSges(7,1) * t45 + mrSges(7,2) * t49), -g(1) * t11 - g(2) * t10 + t118 * g(3) + t122 * (-m(9) * t113 + t105 - t69), g(3) * t68 - g(1) * (-t52 * t105 + t11) - g(2) * (-t48 * t105 + t10), 0, 0];
taug = t5(:);
