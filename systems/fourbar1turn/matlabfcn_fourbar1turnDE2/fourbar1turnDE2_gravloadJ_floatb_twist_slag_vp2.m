% Calculate Gravitation load on the joints for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:29
% EndTime: 2020-04-12 19:33:32
% DurationCPUTime: 0.63s
% Computational Cost: add. (4847->72), mult. (6831->128), div. (304->10), fcn. (1934->11), ass. (0->61)
t45 = pkin(2) ^ 2;
t46 = pkin(1) ^ 2;
t38 = cos(qJ(2));
t74 = pkin(2) * t38;
t67 = -0.2e1 * pkin(1) * t74 + t46;
t31 = t45 + t67;
t28 = 0.1e1 / t31;
t41 = 0.1e1 / pkin(4);
t89 = pkin(4) ^ 2;
t66 = pkin(3) ^ 2 - t89;
t27 = t31 - t66;
t32 = pkin(1) - t74;
t36 = sin(qJ(2));
t79 = -pkin(3) - pkin(4);
t24 = (pkin(2) - t79) * (pkin(2) + t79) + t67;
t78 = pkin(4) - pkin(3);
t25 = (pkin(2) - t78) * (pkin(2) + t78) + t67;
t48 = sqrt(-t24 * t25);
t70 = t36 * t48;
t16 = -pkin(2) * t70 + t27 * t32;
t75 = pkin(2) * t36;
t22 = t27 * t75;
t17 = t32 * t48 + t22;
t68 = t16 ^ 2 + t17 ^ 2;
t29 = 0.1e1 / t31 ^ 2;
t86 = t29 / t89;
t10 = t68 * t86;
t8 = t10 ^ (-0.1e1 / 0.2e1);
t55 = (mrSges(5,1) * t16 + mrSges(5,2) * t17) * t8;
t26 = t31 + t66;
t33 = pkin(1) * t38 - pkin(2);
t15 = -pkin(1) * t70 - t26 * t33;
t23 = pkin(1) * t36 * t26;
t18 = -t33 * t48 + t23;
t44 = 0.1e1 / pkin(3);
t72 = t28 * t44;
t7 = qJ(2) + atan2(t18 * t72, t15 * t72);
t5 = sin(t7);
t6 = cos(t7);
t60 = -mrSges(4,1) * t6 + mrSges(4,2) * t5;
t71 = t36 * mrSges(3,2);
t90 = t41 * t28 * t55 - m(5) * pkin(1) - (m(4) * pkin(2) + mrSges(3,1)) * t38 - mrSges(2,1) - t60 + t71;
t88 = 0.2e1 * pkin(1);
t64 = t29 * t75;
t85 = t64 * t88;
t82 = -0.2e1 * pkin(2);
t81 = 0.2e1 * t36 ^ 2;
t69 = t48 * t38;
t65 = pkin(1) * t75;
t73 = 0.1e1 / t48 * (-t24 - t25) * t65;
t3 = t22 + (-t69 + (t32 * t88 - t73) * t36) * pkin(2);
t4 = t32 * t73 + t45 * pkin(1) * t81 + (t27 * t38 + t70) * pkin(2);
t80 = 0.2e1 * (-0.2e1 * t28 * t65 * t68 + t16 * t3 + t17 * t4) * t8 / t10 * t86;
t61 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t57 = -t3 * t8 + t16 * t80 / 0.2e1;
t56 = t4 * t8 - t17 * t80 / 0.2e1;
t39 = cos(qJ(1));
t37 = sin(qJ(1));
t12 = 0.1e1 / t15 ^ 2;
t1 = 0.1e1 + (((t46 * pkin(2) * t81 - t33 * t73) * t28 + ((t26 * t38 + t70) * t28 - 0.2e1 * t18 * t64) * pkin(1)) / t15 - (t23 * t28 + (-t28 * t69 + ((t33 * t82 - t73) * t28 + t15 * t29 * t82) * t36) * pkin(1)) * t18 * t12) * pkin(3) / (t12 * t18 ^ 2 + 0.1e1) * t31 * t44;
t2 = [(t61 * g(1) + t90 * g(2)) * t39 + (-t90 * g(1) + t61 * g(2)) * t37, (-t38 * mrSges(3,1) + t71 - m(4) * t74 - t60 * t1 - ((-mrSges(5,1) * t17 + mrSges(5,2) * t16) * t8 * t85 + (t56 * mrSges(5,1) + t57 * mrSges(5,2)) * t28) * t41) * g(3) + (g(1) * t39 + g(2) * t37) * (mrSges(3,1) * t36 + mrSges(3,2) * t38 + m(4) * t75 - (mrSges(4,1) * t5 + mrSges(4,2) * t6) * t1 - t41 * (t55 * t85 + (t57 * mrSges(5,1) - t56 * mrSges(5,2)) * t28))];
taug = t2(:);
