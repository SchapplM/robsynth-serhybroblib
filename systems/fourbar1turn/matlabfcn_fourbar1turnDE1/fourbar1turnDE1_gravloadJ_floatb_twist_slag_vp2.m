% Calculate Gravitation load on the joints for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:32
% EndTime: 2020-04-12 19:25:37
% DurationCPUTime: 1.26s
% Computational Cost: add. (8857->101), mult. (12891->182), div. (612->9), fcn. (3622->10), ass. (0->87)
t51 = pkin(2) ^ 2;
t43 = cos(qJ(2));
t101 = pkin(2) * t43;
t52 = pkin(1) ^ 2;
t87 = -0.2e1 * pkin(1) * t101 + t52;
t36 = t51 + t87;
t34 = 0.1e1 / t36 ^ 2;
t119 = pkin(4) ^ 2;
t47 = 0.1e1 / t119;
t118 = pkin(3) ^ 2;
t86 = t118 - t119;
t32 = t36 - t86;
t37 = pkin(1) - t101;
t41 = sin(qJ(2));
t107 = -pkin(3) - pkin(4);
t29 = (pkin(2) - t107) * (pkin(2) + t107) + t87;
t106 = pkin(4) - pkin(3);
t30 = (pkin(2) - t106) * (pkin(2) + t106) + t87;
t55 = sqrt(-t29 * t30);
t93 = t41 * t55;
t21 = -pkin(2) * t93 + t32 * t37;
t102 = pkin(2) * t41;
t27 = t32 * t102;
t22 = t37 * t55 + t27;
t88 = t21 ^ 2 + t22 ^ 2;
t14 = t88 * t47 * t34;
t10 = t14 ^ (-0.1e1 / 0.2e1);
t46 = 0.1e1 / pkin(4);
t33 = 0.1e1 / t36;
t82 = pkin(1) * t102;
t115 = t34 * t82;
t74 = 0.2e1 * t115;
t40 = t41 ^ 2;
t109 = 0.2e1 * t40;
t96 = 0.1e1 / t55 * (-t29 - t30) * t82;
t9 = t37 * t96 + t51 * pkin(1) * t109 + (t32 * t43 + t93) * pkin(2);
t64 = t22 * t74 - t33 * t9;
t90 = t43 * t55;
t7 = t27 + (-t90 + (0.2e1 * t37 * pkin(1) - t96) * t41) * pkin(2);
t65 = t21 * t74 - t33 * t7;
t117 = t10 * t33;
t69 = 0.4e1 * t33 * t115;
t81 = 0.2e1 * t34;
t80 = ((t21 * t7 + t22 * t9) * t81 - t88 * t69) * t47 / t14 * t117;
t97 = t22 * mrSges(5,2);
t98 = t21 * mrSges(5,1);
t120 = m(4) * t102 + mrSges(3,1) * t41 + mrSges(3,2) * t43 - t46 * ((t98 / 0.2e1 + t97 / 0.2e1) * t80 + (t65 * mrSges(5,1) + t64 * mrSges(5,2)) * t10);
t50 = 0.1e1 / t118;
t31 = t36 + t86;
t38 = pkin(1) * t43 - pkin(2);
t20 = -pkin(1) * t93 - t31 * t38;
t28 = pkin(1) * t41 * t31;
t23 = -t38 * t55 + t28;
t89 = t20 ^ 2 + t23 ^ 2;
t15 = t89 * t50 * t34;
t12 = t15 ^ (-0.1e1 / 0.2e1);
t116 = t12 * t33;
t113 = -m(4) * t101 - mrSges(3,1) * t43 + mrSges(3,2) * t41;
t112 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t8 = -t38 * t96 + t52 * pkin(2) * t109 + (t31 * t43 + t93) * pkin(1);
t100 = t20 + t8;
t6 = t28 + (-t90 + (-0.2e1 * t38 * pkin(2) - t96) * t41) * pkin(1);
t99 = t23 - t6;
t111 = t100 * mrSges(4,1) - t99 * mrSges(4,2);
t110 = -(t97 + t98) * t46 * t117 + mrSges(2,1) + m(5) * pkin(1) - t113;
t44 = cos(qJ(1));
t49 = 0.1e1 / pkin(3);
t77 = t49 * t116;
t75 = t43 * t77;
t95 = t41 * t20;
t105 = (t23 * t75 + t77 * t95) * t44;
t94 = t41 * t23;
t92 = t43 * t20;
t91 = t43 * t23;
t83 = pkin(1) * pkin(2) * t34;
t79 = 0.1e1 / t15 * ((t20 * t6 + t23 * t8) * t81 - t89 * t69) * t50 * t116;
t76 = t99 * mrSges(4,1);
t70 = -t92 + t94;
t68 = t20 * t40 + t41 * t91;
t67 = -t94 / 0.2e1 + t92 / 0.2e1;
t66 = -t91 / 0.2e1 - t95 / 0.2e1;
t63 = -0.2e1 * t23 * t40 + 0.2e1 * t41 * t92;
t62 = (t67 * mrSges(4,1) + t66 * mrSges(4,2)) * t79;
t61 = (mrSges(4,1) * t63 - 0.2e1 * mrSges(4,2) * t68) * t83;
t42 = sin(qJ(1));
t3 = t42 * t20 * t75;
t1 = [(-t105 * mrSges(4,2) + (-t70 * mrSges(4,1) * t77 - t110) * t44 + t112 * t42) * g(2) + (-t3 * mrSges(4,1) + t112 * t44 + (-(-mrSges(4,1) * t94 + (-t91 - t95) * mrSges(4,2)) * t77 + t110) * t42) * g(1), (-((-t66 * mrSges(4,1) + t67 * mrSges(4,2)) * t79 + ((0.2e1 * mrSges(4,1) * t68 + mrSges(4,2) * t63) * t83 + (-t111 * t43 + (t100 * mrSges(4,2) + t76) * t41) * t33) * t12) * t49 - ((-t22 * mrSges(5,1) / 0.2e1 + t21 * mrSges(5,2) / 0.2e1) * t80 + (-t64 * mrSges(5,1) + t65 * mrSges(5,2)) * t10) * t46 + t113) * g(3) + (-t3 * mrSges(4,2) + (-(t62 + (t61 + ((t8 * mrSges(4,2) + t76) * t43 + t111 * t41) * t33) * t12) * t49 + t120) * t42) * g(2) + (-t105 * mrSges(4,1) + (-(t62 + (t61 + ((t41 * t8 - t43 * t6) * mrSges(4,1) + (t41 * t6 + t43 * t8 - t70) * mrSges(4,2)) * t33) * t12) * t49 + t120) * t44) * g(1)];
taug = t1(:);
