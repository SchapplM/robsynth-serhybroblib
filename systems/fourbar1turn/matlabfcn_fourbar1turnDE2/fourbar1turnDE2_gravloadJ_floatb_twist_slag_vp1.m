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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:29
% EndTime: 2020-04-12 19:33:32
% DurationCPUTime: 0.69s
% Computational Cost: add. (4847->80), mult. (6831->154), div. (304->10), fcn. (1934->11), ass. (0->60)
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t81 = g(1) * t39 + g(2) * t37;
t83 = pkin(4) ^ 2;
t45 = pkin(2) ^ 2;
t46 = pkin(1) ^ 2;
t38 = cos(qJ(2));
t71 = pkin(2) * t38;
t65 = -0.2e1 * pkin(1) * t71 + t46;
t31 = t45 + t65;
t29 = 0.1e1 / t31 ^ 2;
t82 = t29 / t83;
t64 = pkin(3) ^ 2 - t83;
t27 = t31 - t64;
t32 = pkin(1) - t71;
t36 = sin(qJ(2));
t76 = -pkin(3) - pkin(4);
t24 = (pkin(2) - t76) * (pkin(2) + t76) + t65;
t75 = -pkin(3) + pkin(4);
t25 = (pkin(2) - t75) * (pkin(2) + t75) + t65;
t48 = sqrt(-t24 * t25);
t68 = t36 * t48;
t16 = -pkin(2) * t68 + t27 * t32;
t72 = pkin(2) * t36;
t22 = t27 * t72;
t17 = t32 * t48 + t22;
t66 = t16 ^ 2 + t17 ^ 2;
t10 = t66 * t82;
t8 = t10 ^ (-0.1e1 / 0.2e1);
t80 = 0.2e1 * t8;
t79 = -0.2e1 * pkin(2);
t28 = 0.1e1 / t31;
t78 = 0.2e1 * t36 ^ 2;
t67 = t48 * t38;
t63 = pkin(1) * t72;
t70 = 0.1e1 / t48 * (-t24 - t25) * t63;
t3 = t22 + (-t67 + (0.2e1 * t32 * pkin(1) - t70) * t36) * pkin(2);
t4 = t32 * t70 + t45 * pkin(1) * t78 + (t38 * t27 + t68) * pkin(2);
t77 = (-0.2e1 * t28 * t63 * t66 + t16 * t3 + t17 * t4) / t10 * t80 * t82;
t44 = 0.1e1 / pkin(3);
t69 = t28 * t44;
t62 = t29 * t72;
t26 = t31 + t64;
t33 = pkin(1) * t38 - pkin(2);
t15 = -pkin(1) * t68 - t26 * t33;
t23 = pkin(1) * t36 * t26;
t18 = -t33 * t48 + t23;
t7 = qJ(2) + atan2(t18 * t69, t15 * t69);
t5 = sin(t7);
t6 = cos(t7);
t60 = -rSges(4,1) * t6 + rSges(4,2) * t5;
t59 = rSges(3,1) * t38 - rSges(3,2) * t36;
t58 = rSges(5,1) * t16 + rSges(5,2) * t17;
t57 = -t3 * t8 + t16 * t77 / 0.2e1;
t56 = t4 * t8 - t17 * t77 / 0.2e1;
t55 = -t60 - t71;
t41 = 0.1e1 / pkin(4);
t12 = 0.1e1 / t15 ^ 2;
t1 = 0.1e1 + (((t46 * pkin(2) * t78 - t33 * t70) * t28 + ((t38 * t26 + t68) * t28 - 0.2e1 * t18 * t62) * pkin(1)) / t15 - (t23 * t28 + (-t28 * t67 + ((t33 * t79 - t70) * t28 + t15 * t29 * t79) * t36) * pkin(1)) * t18 * t12) * pkin(3) / (t12 * t18 ^ 2 + 0.1e1) * t31 * t44;
t2 = [-m(2) * (g(1) * (-rSges(2,1) * t37 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 - rSges(2,2) * t37)) - m(3) * (g(1) * (rSges(3,3) * t39 - t59 * t37) + g(2) * (rSges(3,3) * t37 + t59 * t39)) - m(4) * ((g(1) * rSges(4,3) - g(2) * t55) * t39 + (g(2) * rSges(4,3) + g(1) * t55) * t37) - m(5) * (g(1) * (rSges(5,3) * t39 - pkin(1) * t37) + g(2) * (rSges(5,3) * t37 + pkin(1) * t39) + (g(1) * t37 - g(2) * t39) * t8 * t41 * t28 * t58), (-m(3) * t59 - m(4) * (t60 * t1 + t71)) * g(3) - m(5) * ((g(3) * (-rSges(5,1) * t17 + rSges(5,2) * t16) + t81 * t58) * pkin(1) * t62 * t80 + (g(3) * (t56 * rSges(5,1) + t57 * rSges(5,2)) + t81 * (t57 * rSges(5,1) - t56 * rSges(5,2))) * t28) * t41 + t81 * (-m(3) * (-rSges(3,1) * t36 - rSges(3,2) * t38) - m(4) * (-t72 + (rSges(4,1) * t5 + rSges(4,2) * t6) * t1))];
taug = t2(:);
