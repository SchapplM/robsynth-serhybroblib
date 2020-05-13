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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:32
% EndTime: 2020-04-12 19:25:37
% DurationCPUTime: 1.36s
% Computational Cost: add. (8857->110), mult. (12891->211), div. (612->9), fcn. (3622->10), ass. (0->82)
t114 = pkin(4) ^ 2;
t113 = pkin(3) ^ 2;
t51 = pkin(2) ^ 2;
t52 = pkin(1) ^ 2;
t43 = cos(qJ(2));
t99 = pkin(2) * t43;
t82 = -0.2e1 * pkin(1) * t99 + t52;
t36 = t51 + t82;
t34 = 0.1e1 / t36 ^ 2;
t47 = 0.1e1 / t114;
t81 = t113 - t114;
t32 = t36 - t81;
t37 = pkin(1) - t99;
t41 = sin(qJ(2));
t106 = -pkin(3) - pkin(4);
t29 = (pkin(2) - t106) * (pkin(2) + t106) + t82;
t105 = pkin(4) - pkin(3);
t30 = (pkin(2) - t105) * (pkin(2) + t105) + t82;
t55 = sqrt(-t29 * t30);
t89 = t41 * t55;
t21 = -pkin(2) * t89 + t32 * t37;
t98 = t41 * pkin(2);
t27 = t32 * t98;
t22 = t37 * t55 + t27;
t83 = t21 ^ 2 + t22 ^ 2;
t14 = t83 * t47 * t34;
t112 = t14 ^ (-0.1e1 / 0.2e1) / pkin(4);
t50 = 0.1e1 / t113;
t31 = t36 + t81;
t38 = pkin(1) * t43 - pkin(2);
t20 = -pkin(1) * t89 - t31 * t38;
t28 = pkin(1) * t41 * t31;
t23 = -t38 * t55 + t28;
t84 = t20 ^ 2 + t23 ^ 2;
t15 = t84 * t50 * t34;
t111 = t15 ^ (-0.1e1 / 0.2e1) / pkin(3);
t80 = pkin(1) * t98;
t110 = t34 * t80;
t42 = sin(qJ(1));
t101 = g(2) * t42;
t44 = cos(qJ(1));
t102 = g(1) * t44;
t109 = t102 + t101;
t86 = t43 * t55;
t92 = 0.1e1 / t55 * (-t29 - t30) * t80;
t6 = t28 + (-t86 + (-0.2e1 * t38 * pkin(2) - t92) * t41) * pkin(1);
t95 = t23 - t6;
t40 = t41 ^ 2;
t107 = 0.2e1 * t40;
t8 = -t38 * t92 + t52 * pkin(2) * t107 + (t43 * t31 + t89) * pkin(1);
t96 = t20 + t8;
t108 = t96 * rSges(4,1) - t95 * rSges(4,2);
t33 = 0.1e1 / t36;
t78 = t33 * t111;
t76 = t43 * t78;
t85 = t44 * t41;
t104 = t20 * t78 * t85 + t44 * t23 * t76;
t103 = g(1) * t42;
t100 = g(2) * t44;
t97 = t42 * pkin(2);
t94 = t21 * rSges(5,1);
t93 = t22 * rSges(5,2);
t91 = t41 * t20;
t90 = t41 * t23;
t88 = t43 * t20;
t87 = t43 * t23;
t79 = 0.2e1 * t34;
t77 = t95 * rSges(4,1);
t75 = 0.2e1 * t110;
t74 = rSges(3,1) * t43 - rSges(3,2) * t41;
t71 = -t88 + t90;
t70 = 0.4e1 * t33 * t110;
t69 = t20 * t40 + t41 * t87;
t67 = -t90 / 0.2e1 + t88 / 0.2e1;
t66 = -t87 / 0.2e1 - t91 / 0.2e1;
t7 = t27 + (-t86 + (0.2e1 * t37 * pkin(1) - t92) * t41) * pkin(2);
t65 = t21 * t75 - t33 * t7;
t9 = t37 * t92 + t51 * pkin(1) * t107 + (t43 * t32 + t89) * pkin(2);
t64 = t22 * t75 - t33 * t9;
t63 = -0.2e1 * t23 * t40 + 0.2e1 * t41 * t88;
t3 = t42 * t20 * t76;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t42 - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - rSges(2,2) * t42)) - m(3) * (g(1) * (rSges(3,3) * t44 - t74 * t42) + g(2) * (rSges(3,3) * t42 + t74 * t44)) - m(4) * (g(1) * (t3 * rSges(4,1) + t44 * rSges(4,3) - t43 * t97) + g(2) * (t104 * rSges(4,2) + t42 * rSges(4,3) + t44 * t99) + (t71 * rSges(4,1) * t100 + (-rSges(4,1) * t90 + (-t87 - t91) * rSges(4,2)) * t103) * t78) - m(5) * (g(1) * (rSges(5,3) * t44 - pkin(1) * t42) + g(2) * (rSges(5,3) * t42 + pkin(1) * t44) + (-t100 + t103) * t33 * (t93 + t94) * t112), -m(3) * (g(3) * t74 + t109 * (-rSges(3,1) * t41 - rSges(3,2) * t43)) - m(5) * ((g(3) * (-t22 * rSges(5,1) / 0.2e1 + t21 * rSges(5,2) / 0.2e1) + t109 * (t94 / 0.2e1 + t93 / 0.2e1)) * t33 / t14 * ((t21 * t7 + t22 * t9) * t79 - t83 * t70) * t47 + g(3) * (-t64 * rSges(5,1) + t65 * rSges(5,2)) + t109 * (t65 * rSges(5,1) + t64 * rSges(5,2))) * t112 + (-g(1) * (t104 * rSges(4,1) - pkin(2) * t85) - g(2) * (t3 * rSges(4,2) - t41 * t97) - g(3) * t99 - ((g(3) * (0.2e1 * rSges(4,1) * t69 + rSges(4,2) * t63) + t109 * (rSges(4,1) * t63 - 0.2e1 * rSges(4,2) * t69)) * t34 * pkin(2) * pkin(1) + ((g(3) * (-t66 * rSges(4,1) + t67 * rSges(4,2)) + t109 * (t67 * rSges(4,1) + t66 * rSges(4,2))) * ((t20 * t6 + t23 * t8) * t79 - t84 * t70) * t50 / t15 + g(3) * (-t108 * t43 + (t96 * rSges(4,2) + t77) * t41) + ((t41 * t8 - t43 * t6) * rSges(4,1) + (t41 * t6 + t43 * t8 - t71) * rSges(4,2)) * t102 + ((t8 * rSges(4,2) + t77) * t43 + t108 * t41) * t101) * t33) * t111) * m(4)];
taug = t1(:);
