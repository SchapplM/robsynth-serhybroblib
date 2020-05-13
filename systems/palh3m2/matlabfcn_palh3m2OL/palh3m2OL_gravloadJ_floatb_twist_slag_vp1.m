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
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [10x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:42
% EndTime: 2020-05-07 04:33:49
% DurationCPUTime: 0.88s
% Computational Cost: add. (535->151), mult. (492->205), div. (0->0), fcn. (412->18), ass. (0->82)
t51 = sin(qJ(5));
t101 = rSges(6,2) * t51;
t55 = cos(qJ(5));
t105 = rSges(6,1) * t55;
t49 = qJ(2) + qJ(3);
t46 = qJ(4) + t49;
t38 = sin(t46);
t120 = (-t101 + t105) * t38;
t106 = pkin(10) + rSges(6,3);
t119 = t106 * t38;
t39 = cos(t46);
t118 = -t39 * rSges(5,1) + t38 * rSges(5,2);
t48 = qJ(2) + qJ(7);
t45 = pkin(15) - t48;
t40 = -qJ(8) + t45;
t32 = sin(t40);
t33 = cos(t40);
t69 = t33 * rSges(9,1) + t32 * rSges(9,2);
t63 = pkin(3) * cos(t45) - t69;
t42 = sin(t49);
t44 = cos(t49);
t117 = -t44 * rSges(4,1) + t42 * rSges(4,2);
t41 = sin(t48);
t43 = cos(t48);
t116 = t43 * rSges(8,1) - t41 * rSges(8,2);
t104 = rSges(9,1) * t32;
t78 = -t104 + pkin(3) * sin(t45);
t53 = sin(qJ(1));
t57 = cos(qJ(1));
t115 = -g(1) * t57 - g(2) * t53;
t114 = -t39 * pkin(8) - t119;
t112 = pkin(4) * t44;
t56 = cos(qJ(2));
t47 = t56 * pkin(1);
t37 = t47 + pkin(12);
t110 = g(2) * t57 * t37;
t52 = sin(qJ(2));
t107 = t52 * pkin(1);
t103 = rSges(4,2) * t44;
t102 = rSges(5,2) * t39;
t100 = rSges(9,2) * t33;
t99 = t38 * t53;
t98 = t38 * t57;
t95 = t42 * t53;
t94 = t42 * t57;
t92 = t53 * t51;
t91 = t53 * t55;
t90 = t57 * t51;
t89 = t57 * t55;
t88 = rSges(5,1) * t99 + t53 * t102;
t87 = rSges(5,1) * t98 + t57 * t102;
t86 = rSges(4,1) * t95 + t53 * t103;
t85 = rSges(4,1) * t94 + t57 * t103;
t80 = -t107 + t78;
t77 = t56 * rSges(3,1) - t52 * rSges(3,2);
t50 = sin(qJ(6));
t54 = cos(qJ(6));
t73 = t54 * rSges(7,1) - t50 * rSges(7,2);
t70 = -rSges(8,1) * t41 - rSges(8,2) * t43;
t68 = t39 * t101 - t119;
t67 = t118 - t112;
t66 = -t37 - t63;
t65 = -pkin(6) + t73;
t64 = pkin(12) + t77;
t62 = pkin(8) * t99 + t120 * t53;
t61 = pkin(8) * t98 + t120 * t57;
t59 = t68 - t112;
t58 = (g(3) * (-pkin(8) - t105) + t115 * t106) * t39;
t29 = pkin(4) * t94;
t28 = pkin(4) * t95;
t15 = pkin(4) * t42 - t107;
t14 = t57 * t100;
t13 = t53 * t100;
t12 = t37 - t112;
t7 = t57 * t15;
t6 = t53 * t15;
t5 = t57 * t12;
t4 = t39 * t89 - t92;
t3 = t39 * t90 + t91;
t2 = t39 * t91 + t90;
t1 = t39 * t92 - t89;
t8 = [-m(2) * (g(1) * (-t53 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) - t53 * rSges(2,2))) - m(3) * ((g(1) * rSges(3,3) + g(2) * t64) * t57 + (g(2) * rSges(3,3) - g(1) * t64) * t53) - m(4) * (t110 + (g(1) * rSges(4,3) + g(2) * t117) * t57 + (g(1) * (-t37 - t117) + g(2) * rSges(4,3)) * t53) - m(5) * (g(2) * t5 + (g(1) * rSges(5,3) + g(2) * t118) * t57 + (g(1) * (-t12 - t118) + g(2) * rSges(5,3)) * t53) - m(6) * ((-t4 * rSges(6,1) + t3 * rSges(6,2) + t114 * t57 + t5) * g(2) + (t2 * rSges(6,1) - t1 * rSges(6,2) + (-t114 - t12) * t53) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t65) * t57 + (g(2) * rSges(7,3) - g(1) * t65) * t53) - m(8) * (t110 + (g(1) * rSges(8,3) + g(2) * t116) * t57 + (g(1) * (-t37 - t116) + g(2) * rSges(8,3)) * t53) - m(9) * ((g(1) * rSges(9,3) - g(2) * t66) * t57 + (g(2) * rSges(9,3) + g(1) * t66) * t53), -m(3) * (g(3) * t77 - t115 * (-rSges(3,1) * t52 - rSges(3,2) * t56)) - m(4) * (g(1) * (-t57 * t107 + t85) + g(2) * (-t53 * t107 + t86) + g(3) * (t47 + t117)) - m(5) * (g(1) * (t7 + t87) + g(2) * (t6 + t88) + g(3) * (t47 + t67)) - m(6) * (g(1) * (t61 + t7) + g(2) * (t6 + t62) + g(3) * (t47 + t59) + t58) - m(8) * (g(3) * (t47 + t116) - t115 * (t70 - t107)) - m(9) * (g(1) * (t80 * t57 + t14) + g(2) * (t80 * t53 + t13) + g(3) * (t47 + t63)), -m(4) * (g(1) * t85 + g(2) * t86 + g(3) * t117) - m(5) * (g(1) * (t29 + t87) + g(2) * (t28 + t88) + g(3) * t67) - m(6) * (g(1) * (t29 + t61) + g(2) * (t28 + t62) + g(3) * t59 + t58), -m(5) * (g(1) * t87 + g(2) * t88 + g(3) * t118) - m(6) * (g(1) * t61 + g(2) * t62 + g(3) * t68 + t58), -m(6) * (g(1) * (t3 * rSges(6,1) + t4 * rSges(6,2)) + g(2) * (t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (rSges(6,1) * t51 + rSges(6,2) * t55) * t38), -m(7) * (g(3) * t73 - t115 * (-rSges(7,1) * t50 - rSges(7,2) * t54)), -m(8) * (g(3) * t116 - t115 * t70) - m(9) * (g(1) * (t78 * t57 + t14) + g(2) * (t78 * t53 + t13) + g(3) * t63), -m(9) * (g(1) * (-t57 * t104 + t14) + g(2) * (-t53 * t104 + t13) - g(3) * t69), 0, 0];
taug = t8(:);
