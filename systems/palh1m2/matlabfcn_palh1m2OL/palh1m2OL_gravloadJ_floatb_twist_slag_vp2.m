% Calculate Gravitation load on the joints for
% palh1m2OL
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
% taug [13x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:42
% EndTime: 2020-05-02 21:17:46
% DurationCPUTime: 0.49s
% Computational Cost: add. (704->100), mult. (756->136), div. (0->0), fcn. (494->24), ass. (0->65)
t52 = sin(qJ(5));
t60 = cos(qJ(5));
t71 = t60 * mrSges(6,1) - t52 * mrSges(6,2);
t19 = m(6) * pkin(9) + mrSges(5,1) + t71;
t42 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t53 = sin(qJ(4));
t61 = cos(qJ(4));
t13 = -t19 * t53 + t42 * t61;
t11 = -mrSges(4,2) + t13;
t54 = sin(qJ(3));
t62 = cos(qJ(3));
t100 = -t19 * t61 - t42 * t53;
t95 = m(5) + m(6);
t8 = t95 * pkin(5) + mrSges(4,1) - t100;
t4 = t11 * t54 + t8 * t62;
t43 = -qJ(7) + pkin(19) - qJ(2);
t34 = pkin(4) * sin(t43);
t55 = sin(qJ(2));
t39 = -qJ(10) + t43;
t35 = sin(t39);
t36 = cos(t39);
t75 = -t35 * mrSges(11,1) + t36 * mrSges(11,2);
t102 = m(11) * (-t55 * pkin(1) + t34) + t75;
t3 = t11 * t62 - t8 * t54;
t101 = mrSges(11,1) * t36 + mrSges(11,2) * t35;
t56 = sin(qJ(1));
t64 = cos(qJ(1));
t74 = t64 * g(1) + t56 * g(2);
t99 = m(8) + m(4) + t95;
t51 = sin(qJ(6));
t59 = cos(qJ(6));
t70 = -mrSges(7,1) * t59 + mrSges(7,2) * t51;
t98 = pkin(14) * m(7) - mrSges(2,1) + t70 + (-m(9) - m(3) - t99 - m(11) - m(10)) * pkin(15) - t102;
t27 = t52 * mrSges(6,1) + t60 * mrSges(6,2);
t97 = -mrSges(11,3) - mrSges(10,3) - t27 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t96 = pkin(2) * m(10);
t93 = pkin(4) * cos(t43);
t63 = cos(qJ(2));
t88 = m(11) * (-t63 * pkin(1) - t93);
t85 = t101 * t56;
t84 = t101 * t64;
t48 = qJ(2) + qJ(8);
t80 = sin(t48) * t96;
t79 = m(11) * t93;
t45 = qJ(9) + t48;
t40 = sin(t45);
t41 = cos(t45);
t77 = -(g(3) * mrSges(10,1) - t74 * mrSges(10,2)) * t40 + (-t74 * mrSges(10,1) - g(3) * mrSges(10,2)) * t41;
t10 = (t53 * t54 - t61 * t62) * t27;
t9 = (-t53 * t62 - t54 * t61) * t27;
t72 = -t10 * t63 + t9 * t55;
t50 = sin(qJ(7));
t58 = cos(qJ(7));
t29 = -t58 * mrSges(8,1) + mrSges(8,2) * t50;
t26 = -t50 * mrSges(8,1) - t58 * mrSges(8,2);
t49 = sin(qJ(8));
t57 = cos(qJ(8));
t28 = -t57 * mrSges(9,1) + t49 * mrSges(9,2);
t25 = -t49 * mrSges(9,1) - t57 * mrSges(9,2);
t68 = t74 * cos(t48) * t96 + g(3) * t80 + t77;
t6 = t100 * t54 + t13 * t62;
t5 = -t100 * t62 + t13 * t54;
t2 = -mrSges(3,2) + t25 + t26 + t4;
t1 = t99 * pkin(1) + mrSges(3,1) - t28 - t29 - t3;
t7 = [(t97 * g(1) + t98 * g(2)) * t64 + (-t98 * g(1) + t97 * g(2)) * t56 + (mrSges(10,1) * t40 + mrSges(10,2) * t41 - t1 * t55 + t2 * t63 - t80) * (t56 * g(1) - t64 * g(2)), (-t2 * g(3) + t74 * t1) * t63 - (-t1 * g(3) - t74 * t2) * t55 - g(1) * (t64 * t88 + t84) - g(2) * (t56 * t88 + t85) - g(3) * t102 + t68, (-t3 * g(3) + t74 * t4) * t55 + (-t4 * g(3) - t74 * t3) * t63, (-t6 * g(3) + t74 * t5) * t55 + (-t5 * g(3) - t74 * t6) * t63, (-t10 * t55 - t9 * t63) * g(3) + (t72 * t56 + t71 * t64) * g(2) + (-t71 * t56 + t72 * t64) * g(1), t70 * g(3) + t74 * (mrSges(7,1) * t51 + mrSges(7,2) * t59), (-t29 * g(3) + t74 * t26) * t55 + (-t26 * g(3) - t74 * t29) * t63 - g(1) * (-t64 * t79 + t84) - g(2) * (-t56 * t79 + t85) - g(3) * (m(11) * t34 + t75), (-t28 * g(3) + t74 * t25) * t55 + (-t25 * g(3) - t74 * t28) * t63 + t68, t77, -g(1) * t84 - g(2) * t85 - g(3) * t75, 0, 0, 0];
taug = t7(:);
