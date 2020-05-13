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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:34
% EndTime: 2020-04-12 19:18:37
% DurationCPUTime: 0.73s
% Computational Cost: add. (1849->94), mult. (2673->179), div. (120->5), fcn. (808->6), ass. (0->69)
t30 = sin(qJ(1));
t72 = g(2) * t30;
t32 = cos(qJ(1));
t73 = g(1) * t32;
t84 = t73 + t72;
t29 = sin(qJ(2));
t28 = t29 ^ 2;
t83 = 0.2e1 * t28;
t37 = pkin(2) ^ 2;
t38 = pkin(1) ^ 2;
t31 = cos(qJ(2));
t70 = pkin(2) * t31;
t60 = -0.2e1 * pkin(1) * t70 + t38;
t24 = t37 + t60;
t59 = pkin(3) ^ 2 - pkin(4) ^ 2;
t20 = t24 + t59;
t17 = pkin(1) * t29 * t20;
t26 = pkin(1) * t31 - pkin(2);
t78 = -pkin(3) - pkin(4);
t18 = (pkin(2) - t78) * (pkin(2) + t78) + t60;
t77 = -pkin(3) + pkin(4);
t19 = (pkin(2) - t77) * (pkin(2) + t77) + t60;
t39 = sqrt(-t18 * t19);
t61 = t31 * t39;
t69 = t29 * pkin(2);
t58 = pkin(1) * t69;
t65 = 0.1e1 / t39 * (-t18 - t19) * t58;
t1 = t17 + (-t61 + (-0.2e1 * t26 * pkin(2) - t65) * t29) * pkin(1);
t82 = -t1 / 0.2e1;
t63 = t29 * t39;
t9 = -pkin(1) * t63 - t20 * t26;
t81 = -t9 / 0.2e1;
t80 = t29 / 0.2e1;
t79 = t31 / 0.2e1;
t12 = -t26 * t39 + t17;
t22 = 0.1e1 / t24;
t36 = 0.1e1 / pkin(3);
t64 = t22 * t36;
t53 = t64 * t79;
t57 = t29 * t64;
t76 = (t9 * t53 - t12 * t57 / 0.2e1) * t30;
t75 = (t12 * t53 + t9 * t57 / 0.2e1) * t32;
t74 = g(1) * t30;
t71 = g(2) * t32;
t68 = t30 * pkin(2);
t67 = t31 * t9;
t66 = t32 * pkin(2);
t62 = t31 * t12;
t3 = -t26 * t65 + t38 * pkin(2) * t83 + (t31 * t20 + t63) * pkin(1);
t56 = t81 - t3 / 0.2e1;
t55 = t82 + t12 / 0.2e1;
t23 = 0.1e1 / t24 ^ 2;
t54 = t23 * t58;
t52 = rSges(3,1) * t31 - rSges(3,2) * t29;
t50 = -t12 * t28 + t29 * t67;
t49 = t28 * t9 + t29 * t62;
t48 = t3 * t80 + t31 * t82;
t47 = t1 * t80 + t3 * t79;
t46 = -t62 / 0.2e1 + t29 * t81;
t45 = -t67 / 0.2e1 + t12 * t80;
t21 = t24 - t59;
t25 = pkin(1) - t70;
t10 = -pkin(2) * t63 + t21 * t25;
t16 = t21 * t69;
t43 = -(t16 + (-t61 + (0.2e1 * t25 * pkin(1) - t65) * t29) * pkin(2)) * t22 / 0.2e1 + t10 * t54;
t11 = t25 * t39 + t16;
t42 = (t25 * t65 + t37 * pkin(1) * t83 + (t31 * t21 + t63) * pkin(2)) * t22 / 0.2e1 - t11 * t54;
t34 = 0.1e1 / pkin(4);
t2 = [-m(2) * (g(1) * (-rSges(2,1) * t30 - rSges(2,2) * t32) + g(2) * (rSges(2,1) * t32 - rSges(2,2) * t30)) - m(3) * (g(1) * (rSges(3,3) * t32 - t52 * t30) + g(2) * (rSges(3,3) * t30 + t52 * t32)) - m(4) * (g(1) * (t76 * rSges(4,1) + t32 * rSges(4,3) - t31 * t68) + g(2) * (t75 * rSges(4,2) + t30 * rSges(4,3) + t31 * t66) + (t45 * rSges(4,1) * t71 + t46 * rSges(4,2) * t74) * t64) - m(5) * (g(1) * (t32 * rSges(5,3) - t30 * pkin(1)) + g(2) * (t30 * rSges(5,3) + t32 * pkin(1)) + (-t71 + t74) * t34 * t22 * (t10 * rSges(5,1) / 0.2e1 + t11 * rSges(5,2) / 0.2e1)), -m(3) * (g(3) * t52 + t84 * (-rSges(3,1) * t29 - rSges(3,2) * t31)) - m(5) * (g(3) * (t42 * rSges(5,1) + t43 * rSges(5,2)) + t84 * (t43 * rSges(5,1) - t42 * rSges(5,2))) * t34 + (-g(1) * (t75 * rSges(4,1) - t29 * t66) - g(2) * (t76 * rSges(4,2) - t29 * t68) - g(3) * t70 - ((g(3) * (t49 * rSges(4,1) + t50 * rSges(4,2)) + t84 * (t50 * rSges(4,1) - t49 * rSges(4,2))) * t23 * pkin(2) * pkin(1) + (g(3) * ((t56 * rSges(4,1) + t55 * rSges(4,2)) * t31 + (t55 * rSges(4,1) - t56 * rSges(4,2)) * t29) + (t48 * rSges(4,1) + (-t45 + t47) * rSges(4,2)) * t73 + ((-t46 + t48) * rSges(4,1) + t47 * rSges(4,2)) * t72) * t22) * t36) * m(4)];
taug = t2(:);
