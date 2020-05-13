% Calculate Gravitation load on the joints for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2TE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2TE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:14:48
% EndTime: 2020-05-01 20:14:52
% DurationCPUTime: 0.71s
% Computational Cost: add. (408->141), mult. (656->206), div. (0->0), fcn. (349->34), ass. (0->77)
t73 = m(5) + m(6);
t28 = pkin(5) * t73 + m(4) * rSges(4,1);
t95 = -m(11) * rSges(11,1) - t28;
t96 = t95 * g(1);
t86 = m(4) + t73;
t57 = sin(qJ(4));
t63 = cos(qJ(4));
t80 = t57 * rSges(6,1) + t63 * rSges(6,2);
t94 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + m(6) * t80 + m(7) * rSges(7,3);
t93 = m(9) * rSges(9,2);
t92 = rSges(3,2) * m(3);
t70 = rSges(4,2) * m(4);
t91 = g(2) * t95;
t36 = pkin(2) * m(10) + m(9) * rSges(9,1);
t90 = g(2) * t36;
t54 = cos(pkin(19));
t89 = rSges(9,2) * t54;
t51 = sin(pkin(19));
t88 = t51 * rSges(9,2);
t32 = m(11) * rSges(11,2) + t70;
t47 = qJ(3) + qJ(2);
t85 = g(2) * t93;
t59 = sin(qJ(2));
t84 = -pkin(1) * t59 + pkin(15);
t39 = pkin(19) + t47;
t83 = (-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2);
t60 = sin(qJ(1));
t66 = cos(qJ(1));
t27 = t66 * g(1) + t60 * g(2);
t82 = t60 * g(1) - t66 * g(2);
t81 = t63 * rSges(6,1) - t57 * rSges(6,2);
t61 = sin(pkin(18));
t67 = cos(pkin(18));
t79 = t67 * rSges(7,1) + t61 * rSges(7,2);
t78 = t61 * rSges(7,1) - t67 * rSges(7,2);
t40 = sin(t47);
t41 = cos(t47);
t48 = sin(pkin(22));
t52 = cos(pkin(22));
t77 = rSges(11,1) * t41 - rSges(11,2) * t40 + (-(t61 * t48 + t52 * t67) * cos(pkin(21)) + (t67 * t48 - t61 * t52) * sin(pkin(21))) * pkin(4) + t84;
t76 = m(6) * t81;
t75 = m(9) * t89 + t36 * t51 + t32;
t62 = sin(pkin(17));
t68 = cos(pkin(17));
t25 = rSges(7,1) * t62 + rSges(7,2) * t68;
t26 = rSges(7,1) * t68 - rSges(7,2) * t62;
t58 = sin(qJ(3));
t64 = cos(qJ(3));
t65 = cos(qJ(2));
t71 = rSges(3,1) * m(3);
t74 = (t64 * t70 + t28 * t58 + t71 + (t25 * t61 + t26 * t67) * m(7) + t86 * pkin(1)) * t59 - (-t58 * t70 - t92 + t28 * t64 + (-t25 * t67 + t26 * t61) * m(7)) * t65 - m(2) * rSges(2,1) + m(7) * pkin(14) + (-m(3) - t86) * pkin(15);
t53 = cos(pkin(20));
t50 = sin(pkin(20));
t46 = g(1) * t93;
t45 = pkin(22) + pkin(21);
t43 = qJ(1) - t47;
t42 = qJ(1) + t47;
t38 = cos(t45);
t37 = sin(t45);
t35 = pkin(9) * m(6) + m(5) * rSges(5,1);
t34 = -qJ(1) + t39;
t33 = qJ(1) + t39;
t31 = t36 * g(1);
t30 = g(2) * t32;
t29 = t32 * g(1);
t16 = g(3) * t61 + t27 * t67;
t15 = t67 * g(3) - t61 * t27;
t14 = t50 * t35 - t53 * t83;
t13 = t35 * t53 + t50 * t83;
t12 = -m(9) * t88 + t36 * t54 - t95;
t11 = t82 * (t51 * rSges(9,1) + t89);
t10 = t82 * (rSges(9,1) * t54 - t88);
t9 = (-rSges(8,1) * t67 + rSges(8,2) * t61) * t52 + (-rSges(8,1) * t61 - rSges(8,2) * t67) * t48 + t84;
t8 = -rSges(10,1) * t59 - rSges(10,2) * t65 + pkin(15) + ((-t51 * t58 + t54 * t64) * t65 + (-t51 * t64 - t54 * t58) * t59) * pkin(2);
t2 = t12 * t64 - t75 * t58 - rSges(10,2) * m(10) - t92 + (-t62 * t79 + t68 * t78) * m(7);
t1 = t75 * t64 + t12 * t58 + rSges(10,1) * m(10) + t71 + (m(11) + m(8) + t86) * pkin(1) + (t62 * t78 + t79 * t68) * m(7);
t3 = [-(t74 * g(1) + t94 * g(2)) * t60 + (-t94 * g(1) + t74 * g(2)) * t66 - m(8) * (g(1) * (t66 * rSges(8,3) - t9 * t60) + g(2) * (t60 * rSges(8,3) + t9 * t66)) - ((-t10 * t64 + t11 * t58) * t65 - (-rSges(9,3) * g(2) + pkin(15) * g(1)) * t60 + (rSges(9,3) * g(1) + pkin(15) * g(2)) * t66 + (t10 * t58 + t11 * t64) * t59) * m(9) - m(10) * (g(1) * (t66 * rSges(10,3) - t8 * t60) + g(2) * (t60 * rSges(10,3) + t8 * t66)) - m(11) * ((g(1) * rSges(11,3) + g(2) * t77) * t66 + (g(2) * rSges(11,3) - g(1) * t77) * t60) + (-t38 * (t13 * t67 + t14 * t61 + (t61 * t50 + t67 * t53) * t76) + t37 * (-t13 * t61 + t14 * t67 + (t67 * t50 - t61 * t53) * t76)) * t82, (-g(3) * t2 + t1 * t27) * t65 + (g(3) * t1 + t2 * t27) * t59, (t29 - t91) * cos(t43) / 0.2e1 + (t30 + t96) * sin(t43) / 0.2e1 + (t46 + t90) * cos(t34) / 0.2e1 + (t31 - t85) * sin(t34) / 0.2e1 + (t46 - t90) * cos(t33) / 0.2e1 + (t31 + t85) * sin(t33) / 0.2e1 + (t29 + t91) * cos(t42) / 0.2e1 + (t30 - t96) * sin(t42) / 0.2e1 - g(3) * (t36 * cos(t39) - sin(t39) * t93 - t95 * t41 - t40 * t32), -(t81 * t82 + ((-t15 * t50 + t16 * t53) * t38 - (t15 * t53 + t50 * t16) * t37) * t80) * m(6)];
taug = t3(:);
