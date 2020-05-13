% Calculate Gravitation load on the joints for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:34
% EndTime: 2020-04-12 19:18:37
% DurationCPUTime: 0.66s
% Computational Cost: add. (1849->85), mult. (2673->148), div. (120->5), fcn. (808->6), ass. (0->67)
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t71 = pkin(2) * t31;
t88 = -m(4) * t71 - t31 * mrSges(3,1) + t29 * mrSges(3,2);
t86 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t37 = pkin(2) ^ 2;
t38 = pkin(1) ^ 2;
t63 = -0.2e1 * pkin(1) * t71 + t38;
t24 = t37 + t63;
t62 = pkin(3) ^ 2 - pkin(4) ^ 2;
t21 = t24 - t62;
t25 = pkin(1) - t71;
t78 = -pkin(3) - pkin(4);
t18 = (pkin(2) - t78) * (pkin(2) + t78) + t63;
t77 = -pkin(3) + pkin(4);
t19 = (pkin(2) - t77) * (pkin(2) + t77) + t63;
t39 = sqrt(-t18 * t19);
t66 = t29 * t39;
t10 = -pkin(2) * t66 + t21 * t25;
t72 = pkin(2) * t29;
t16 = t21 * t72;
t11 = t25 * t39 + t16;
t22 = 0.1e1 / t24;
t34 = 0.1e1 / pkin(4);
t85 = (t10 * mrSges(5,1) / 0.2e1 + t11 * mrSges(5,2) / 0.2e1) * t34 * t22 - mrSges(2,1) - m(5) * pkin(1) + t88;
t28 = t29 ^ 2;
t84 = 0.2e1 * t28;
t20 = t24 + t62;
t17 = pkin(1) * t29 * t20;
t26 = pkin(1) * t31 - pkin(2);
t64 = t31 * t39;
t58 = pkin(1) * t72;
t69 = 0.1e1 / t39 * (-t18 - t19) * t58;
t1 = t17 + (-t64 + (-0.2e1 * t26 * pkin(2) - t69) * t29) * pkin(1);
t83 = -t1 / 0.2e1;
t9 = -pkin(1) * t66 - t20 * t26;
t82 = -t9 / 0.2e1;
t80 = t29 / 0.2e1;
t79 = t31 / 0.2e1;
t12 = -t26 * t39 + t17;
t30 = sin(qJ(1));
t36 = 0.1e1 / pkin(3);
t67 = t22 * t36;
t53 = t67 * t79;
t57 = t29 * t67;
t76 = (t9 * t53 - t12 * t57 / 0.2e1) * t30;
t32 = cos(qJ(1));
t75 = (t12 * t53 + t9 * t57 / 0.2e1) * t32;
t70 = t31 * t9;
t65 = t31 * t12;
t60 = m(4) * t72;
t23 = 0.1e1 / t24 ^ 2;
t59 = pkin(1) * pkin(2) * t23;
t3 = -t26 * t69 + t38 * pkin(2) * t84 + (t20 * t31 + t66) * pkin(1);
t56 = t82 - t3 / 0.2e1;
t55 = t83 + t12 / 0.2e1;
t54 = t23 * t58;
t50 = -t12 * t28 + t29 * t70;
t49 = t28 * t9 + t29 * t65;
t48 = t3 * t80 + t31 * t83;
t47 = t1 * t80 + t3 * t79;
t46 = -t65 / 0.2e1 + t29 * t82;
t45 = -t70 / 0.2e1 + t12 * t80;
t43 = -(t16 + (-t64 + (0.2e1 * t25 * pkin(1) - t69) * t29) * pkin(2)) * t22 / 0.2e1 + t10 * t54;
t42 = (t25 * t69 + t37 * pkin(1) * t84 + (t21 * t31 + t66) * pkin(2)) * t22 / 0.2e1 - t11 * t54;
t41 = (t50 * mrSges(4,1) - t49 * mrSges(4,2)) * t59;
t2 = [(-t75 * mrSges(4,2) + t86 * t30 + (-mrSges(4,1) * t45 * t67 + t85) * t32) * g(2) + (-t76 * mrSges(4,1) + t86 * t32 + (-mrSges(4,2) * t46 * t67 - t85) * t30) * g(1), -g(1) * (t75 * mrSges(4,1) + (-t60 + (t41 + (t48 * mrSges(4,1) + (-t45 + t47) * mrSges(4,2)) * t22) * t36) * t32) - g(2) * (t76 * mrSges(4,2) + (-t60 + (t41 + ((-t46 + t48) * mrSges(4,1) + t47 * mrSges(4,2)) * t22) * t36) * t30) + (-((t49 * mrSges(4,1) + t50 * mrSges(4,2)) * t59 + ((t56 * mrSges(4,1) + t55 * mrSges(4,2)) * t31 + (t55 * mrSges(4,1) - t56 * mrSges(4,2)) * t29) * t22) * t36 - (t42 * mrSges(5,1) + t43 * mrSges(5,2)) * t34 + t88) * g(3) + (mrSges(3,1) * t29 + mrSges(3,2) * t31 - t34 * (t43 * mrSges(5,1) - t42 * mrSges(5,2))) * (g(1) * t32 + g(2) * t30)];
taug = t2(:);
