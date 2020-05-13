% Calculate Gravitation load on the joints for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [12x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:57
% EndTime: 2020-05-11 05:44:59
% DurationCPUTime: 0.25s
% Computational Cost: add. (362->97), mult. (210->119), div. (0->0), fcn. (138->20), ass. (0->62)
t43 = qJ(1) + qJ(2);
t37 = qJ(6) + t43;
t23 = sin(t37);
t26 = cos(t37);
t58 = -rSges(7,1) * t23 - rSges(7,2) * t26;
t67 = t26 * rSges(7,1) - rSges(7,2) * t23;
t77 = m(7) * (g(1) * t58 + g(2) * t67);
t39 = qJ(3) + t43;
t30 = qJ(9) + t39;
t15 = sin(t30);
t16 = cos(t30);
t76 = t15 * rSges(10,1) + t16 * rSges(10,2);
t38 = qJ(4) + t43;
t24 = sin(t38);
t27 = cos(t38);
t75 = t24 * rSges(5,1) + t27 * rSges(5,2);
t36 = cos(t43);
t74 = pkin(2) * t36;
t73 = pkin(3) * t36;
t47 = cos(qJ(1));
t72 = t47 * pkin(1);
t34 = sin(t43);
t71 = t34 * rSges(3,1) + t36 * rSges(3,2);
t21 = pkin(3) * t34;
t70 = t21 + t75;
t69 = -rSges(5,1) * t27 + t24 * rSges(5,2);
t68 = -rSges(10,1) * t16 + t15 * rSges(10,2);
t66 = -rSges(3,1) * t36 + t34 * rSges(3,2);
t29 = qJ(10) + t38;
t13 = sin(t29);
t14 = cos(t29);
t65 = t14 * rSges(11,1) - rSges(11,2) * t13;
t25 = sin(t39);
t28 = cos(t39);
t64 = t28 * rSges(4,1) - rSges(4,2) * t25;
t42 = qJ(1) + qJ(8);
t33 = sin(t42);
t35 = cos(t42);
t63 = t35 * rSges(9,1) - rSges(9,2) * t33;
t62 = -pkin(6) * t25 + t76;
t61 = pkin(6) * t28 + t68;
t22 = pkin(2) * t34;
t60 = t22 + t62;
t59 = -rSges(4,1) * t25 - rSges(4,2) * t28;
t57 = -rSges(9,1) * t33 - rSges(9,2) * t35;
t56 = -rSges(11,1) * t13 - rSges(11,2) * t14;
t55 = t69 - t73;
t54 = t64 - t74;
t53 = -pkin(4) * t27 + t65;
t52 = t22 + t59;
t51 = pkin(4) * t24 + t56;
t50 = t61 - t74;
t49 = t21 + t51;
t48 = t53 - t73;
t46 = cos(qJ(7));
t45 = sin(qJ(1));
t44 = sin(qJ(7));
t41 = pkin(8) + qJ(5);
t40 = t45 * pkin(1);
t32 = cos(t41);
t31 = sin(t41);
t1 = [-m(2) * (g(1) * (rSges(2,1) * t45 + rSges(2,2) * t47) + g(2) * (-rSges(2,1) * t47 + rSges(2,2) * t45)) - m(3) * (g(1) * (t40 + t71) + g(2) * (t66 - t72)) - m(4) * (g(1) * (t40 + t52) + g(2) * (t54 - t72)) - m(5) * (g(1) * (t40 + t70) + g(2) * (t55 - t72)) - m(7) * (g(1) * (t40 + t58) + g(2) * (t67 - t72)) - m(9) * (g(1) * (t40 + t57) + g(2) * (t63 - t72)) - m(10) * (g(1) * (t40 + t60) + g(2) * (t50 - t72)) - m(11) * (g(1) * (t40 + t49) + g(2) * (t48 - t72)), -m(3) * (g(1) * t71 + g(2) * t66) - m(4) * (g(1) * t52 + g(2) * t54) - m(5) * (g(1) * t70 + g(2) * t55) - t77 - m(10) * (g(1) * t60 + g(2) * t50) - m(11) * (g(1) * t49 + g(2) * t48), -m(4) * (g(1) * t59 + g(2) * t64) - m(10) * (g(1) * t62 + g(2) * t61), -m(5) * (g(1) * t75 + g(2) * t69) - m(11) * (g(1) * t51 + g(2) * t53), -m(6) * (g(1) * (-rSges(6,1) * t31 - rSges(6,2) * t32) + g(2) * (rSges(6,1) * t32 - rSges(6,2) * t31)), -t77, -m(8) * (g(1) * (rSges(8,1) * t46 - rSges(8,2) * t44) + g(2) * (rSges(8,1) * t44 + rSges(8,2) * t46)), -m(9) * (g(1) * t57 + g(2) * t63), -m(10) * (g(1) * t76 + g(2) * t68), -m(11) * (g(1) * t56 + g(2) * t65), 0, 0];
taug = t1(:);
