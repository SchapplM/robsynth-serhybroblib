% Calculate Gravitation load on the joints for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2DE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:34
% EndTime: 2020-05-07 02:14:34
% DurationCPUTime: 0.42s
% Computational Cost: add. (563->112), mult. (422->142), div. (0->0), fcn. (280->47), ass. (0->71)
t92 = m(5) + m(6);
t98 = m(9) + m(4) + m(8) + t92;
t55 = sin(qJ(4));
t60 = cos(qJ(4));
t97 = t55 * rSges(6,1) + t60 * rSges(6,2);
t57 = sin(qJ(2));
t62 = cos(qJ(2));
t27 = rSges(9,1) * t57 + t62 * rSges(9,2);
t32 = pkin(4) * t92 + m(4) * rSges(4,1);
t52 = pkin(15) + pkin(18);
t46 = qJ(2) + t52;
t33 = pkin(17) + qJ(3) + t46;
t17 = pkin(17) - atan2(-sin(t33), -cos(t33));
t34 = sin(t46);
t35 = cos(t46);
t5 = -atan2(t34, -t35) + t17;
t4 = cos(t5);
t53 = qJ(2) + qJ(3);
t44 = sin(t53);
t45 = cos(t53);
t89 = m(4) * rSges(4,2);
t96 = -m(9) * t27 * t4 - t32 * t44 - t45 * t89;
t95 = -m(7) / 0.2e1;
t94 = m(7) / 0.2e1;
t28 = t62 * rSges(9,1) - t57 * rSges(9,2);
t3 = sin(t5);
t58 = sin(qJ(1));
t63 = cos(qJ(1));
t30 = t63 * g(1) + t58 * g(2);
t91 = (t28 * t30 * t3 + (-t27 * t35 + t34 * t28) * g(3) * sin(t17)) * m(9);
t66 = rSges(7,1) * g(1);
t86 = rSges(7,2) * g(2);
t36 = t66 + t86;
t65 = rSges(7,2) * g(1);
t88 = rSges(7,1) * g(2);
t39 = t65 - t88;
t77 = pkin(14) - qJ(2);
t48 = -qJ(1) + t77;
t59 = sin(pkin(15));
t64 = cos(pkin(15));
t90 = (t59 * t36 + t39 * t64) * cos(t48) * t94 + (t36 * t64 - t59 * t39) * sin(t48) * t95;
t87 = rSges(3,2) * m(3);
t85 = rSges(6,1) * t57;
t84 = rSges(6,2) * t57;
t37 = t66 - t86;
t38 = t65 + t88;
t47 = qJ(1) + t77;
t83 = (t37 * t64 - t59 * t38) * sin(t47);
t82 = (t59 * t37 + t38 * t64) * cos(t47);
t56 = sin(qJ(3));
t61 = cos(qJ(3));
t78 = t32 * t56 + t61 * t89;
t31 = pkin(16) + t33;
t16 = atan2(-sin(t31), cos(t31));
t14 = qJ(3) + t16;
t75 = (t34 * t27 + t28 * t35) * m(9) * cos(t17);
t23 = t98 * pkin(1) + rSges(3,1) * m(3);
t12 = -qJ(4) + t14;
t11 = qJ(4) + t14;
t71 = t32 * t61 - t56 * t89;
t70 = -t23 * t62 + t57 * t87 + pkin(6) * m(7) - m(2) * rSges(2,1) - (m(3) + t98) * pkin(12);
t69 = m(2) * rSges(2,2) - m(3) * rSges(3,3) - m(4) * rSges(4,3) - m(5) * rSges(5,3) - t97 * m(6) - m(7) * rSges(7,3) - m(8) * rSges(8,3) - m(9) * rSges(9,3);
t54 = -pkin(15) + pkin(14);
t51 = pkin(17) - qJ(2);
t50 = rSges(6,1) * t62;
t49 = rSges(6,2) * t62;
t29 = t58 * g(1) - t63 * g(2);
t13 = t16 + t53;
t10 = qJ(2) + t12;
t9 = qJ(2) + t11;
t1 = [(-t82 / 0.2e1 + t83 / 0.2e1) * m(7) + (t69 * t58 + t63 * t70) * g(2) + (-t70 * t58 + t63 * t69) * g(1) + (-((pkin(10) + rSges(6,3)) * m(6) - m(5) * rSges(5,2)) * sin(t13) - (pkin(8) * m(6) + m(5) * rSges(5,1)) * cos(t13) - t32 * t45 + t44 * t89 + (-sin(t10) / 0.2e1 + sin(t9) / 0.2e1) * m(6) * rSges(6,2) + (-cos(t10) / 0.2e1 - cos(t9) / 0.2e1) * m(6) * rSges(6,1) + (-rSges(8,1) * cos(t52) - rSges(8,2) * sin(t52)) * m(8) + (-t27 * t3 - t28 * t4 + (-cos(t51) * t35 + sin(t51) * t34) * pkin(3)) * m(9)) * t29 + t90, t82 * t94 + t83 * t95 + t90 + t91 + (t23 * t57 + t87 * t62 + t96) * t30 + ((-(rSges(7,1) * t62 - rSges(7,2) * t57) * cos(t54) - (rSges(7,1) * t57 + rSges(7,2) * t62) * sin(t54)) * m(7) - t75 - (t78 - t87) * t57 - (-t71 + t23) * t62) * g(3), (-t78 * t57 + t71 * t62 - t75) * g(3) + t96 * t30 + t91, -(0.2e1 * (t60 * rSges(6,1) - t55 * rSges(6,2)) * t29 + 0.2e1 * (cos(t14) * t57 + sin(t14) * t62) * t97 * g(3) + ((t49 - t85) * cos(t12) - (t50 + t84) * sin(t12) + (t49 + t85) * cos(t11) + (t50 - t84) * sin(t11)) * t30) * m(6) / 0.2e1];
taug = t1(:);
