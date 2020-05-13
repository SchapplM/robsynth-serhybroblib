% Calculate potential energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:25
% EndTime: 2020-05-07 02:14:27
% DurationCPUTime: 1.31s
% Computational Cost: add. (695->146), mult. (356->148), div. (0->0), fcn. (182->68), ass. (0->81)
t92 = 2 * g(1);
t95 = m(5) + m(6);
t64 = m(9) + m(4) + m(8) + t95;
t94 = m(3) + t64;
t23 = (pkin(1) * t64 + mrSges(3,1));
t93 = (m(7) * pkin(6) - pkin(12) * t94 - mrSges(2,1));
t91 = pkin(3) * m(9);
t84 = mrSges(6,1) * g(1);
t83 = mrSges(6,1) * g(2);
t82 = mrSges(7,1) * g(1);
t81 = mrSges(7,1) * g(2);
t80 = mrSges(8,1) * g(1);
t79 = mrSges(8,1) * g(2);
t78 = mrSges(9,1) * g(2);
t77 = mrSges(3,2) * g(1);
t76 = mrSges(3,2) * g(2);
t75 = mrSges(4,2) * g(1);
t74 = mrSges(4,2) * g(2);
t73 = mrSges(6,2) * g(1);
t71 = mrSges(7,2) * g(1);
t70 = mrSges(7,2) * g(2);
t69 = mrSges(8,2) * g(2);
t68 = mrSges(9,2) * g(2);
t67 = g(2) * t23;
t38 = (pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3));
t66 = g(2) * t38;
t43 = (pkin(8) * m(6) + mrSges(5,1));
t65 = g(2) * t43;
t63 = pkin(15) + qJ(2);
t39 = pkin(18) + t63;
t20 = atan2(sin(t39), -cos(t39));
t19 = qJ(2) + t20;
t50 = qJ(1) - qJ(2);
t49 = qJ(1) + qJ(2);
t46 = qJ(2) + qJ(3);
t30 = pkin(17) + qJ(3) + t39;
t24 = pkin(16) + t30;
t13 = atan2(-sin(t24), cos(t24));
t12 = t13 + t46;
t59 = pkin(17) - t19;
t42 = qJ(1) - t46;
t41 = qJ(1) + t46;
t40 = pkin(14) - t63;
t9 = -qJ(4) + t12;
t8 = qJ(4) + t12;
t16 = -qJ(1) + t59;
t15 = qJ(1) + t59;
t55 = m(7) + m(2) + t94;
t54 = mrSges(9,1) * g(1);
t53 = mrSges(6,2) * g(2);
t52 = mrSges(8,2) * g(1);
t51 = mrSges(9,2) * g(1);
t48 = qJ(1) - qJ(4);
t47 = qJ(1) + qJ(4);
t37 = -qJ(1) + t40;
t36 = qJ(1) + t40;
t35 = -t73 + t83;
t34 = t73 + t83;
t33 = -t53 + t84;
t32 = t53 + t84;
t31 = (t95 * pkin(4) + mrSges(4,1));
t29 = t43 * t92;
t28 = (mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3));
t27 = g(2) * t31;
t26 = m(1) + t55;
t25 = t38 * t92;
t22 = t31 * t92;
t21 = t23 * t92;
t18 = -t20 + t50;
t17 = t20 + t49;
t14 = atan2(-sin(t30), -cos(t30));
t11 = -t13 + t42;
t10 = t13 + t41;
t7 = qJ(1) - t8;
t6 = qJ(1) + t9;
t5 = qJ(1) - t9;
t4 = qJ(1) + t8;
t3 = -t14 + t59;
t2 = -t14 + t16;
t1 = -t14 + t15;
t44 = -((cos(t16) + cos(t15)) * g(1) + g(2) * sin(t15)) * t91 / 0.2e1 + (-t21 + 2 * t76) * cos(t50) / 0.4e1 + (-t21 - 2 * t76) * cos(t49) / 0.4e1 + (g(1) * t93 - g(2) * t28) * cos(qJ(1)) + t35 * sin(t47) / 0.2e1 + t32 * cos(t47) / 0.2e1 - t33 * cos(t48) / 0.2e1 - t34 * sin(t48) / 0.2e1 + ((cos(t8) / 0.2e1 - cos(t9) / 0.2e1) * mrSges(6,2) + (sin(t8) + sin(t9)) * mrSges(6,1) / 0.2e1 - pkin(13) * m(7) - t55 * pkin(11) - mrSges(1,3) - mrSges(2,3) + sin(t59) * t91 - t38 * cos(t12) - mrSges(7,2) * cos(t40) + mrSges(7,1) * sin(t40) + t43 * sin(t12) + t31 * sin(t46) + mrSges(4,2) * cos(t46) - mrSges(8,1) * sin(t19) - mrSges(8,2) * cos(t19) - t23 * sin(qJ(2)) - mrSges(9,1) * sin(t3) + mrSges(9,2) * cos(t3) - mrSges(3,2) * cos(qJ(2)) - t26 * r_base(3)) * g(3) + (t22 + 2 * t74) * cos(t41) / 0.4e1 + (t22 - 2 * t74) * cos(t42) / 0.4e1 + (t27 - t75) * sin(t41) / 0.2e1 + (t27 + t75) * sin(t42) / 0.2e1 + (-t67 - t77) * sin(t50) / 0.2e1 + (-t67 + t77) * sin(t49) / 0.2e1 + (t51 - t78) * sin(t2) / 0.2e1 + (t51 + t78) * sin(t1) / 0.2e1 + (t52 - t79) * sin(t17) / 0.2e1 - (t52 + t79) * sin(t18) / 0.2e1 + (-t69 - t80) * cos(t17) / 0.2e1 + (t69 - t80) * cos(t18) / 0.2e1 + (-t71 - t81) * sin(t36) / 0.2e1 + (-t71 + t81) * sin(t37) / 0.2e1 + (-t70 - t82) * cos(t37) / 0.2e1 + (t70 - t82) * cos(t36) / 0.2e1 + (t25 + 2 * t65) * sin(t10) / 0.4e1 - (t25 - 2 * t65) * sin(t11) / 0.4e1 + (t29 + 2 * t66) * cos(t11) / 0.4e1 + (t29 - 2 * t66) * cos(t10) / 0.4e1 + (t54 - t68) * cos(t1) / 0.2e1 + (t54 + t68) * cos(t2) / 0.2e1 + (sin(t5) + sin(t4)) * t35 / 0.4e1 + (sin(t7) + sin(t6)) * t34 / 0.4e1 + (cos(t7) + cos(t6)) * t33 / 0.4e1 + (cos(t5) + cos(t4)) * t32 / 0.4e1 + g(2) * sin(t16) * t91 / 0.2e1 + (-g(1) * r_base(1) - g(2) * r_base(2)) * t26 + (t28 * g(1) + t93 * g(2)) * sin(qJ(1)) - (mrSges(1,1) * g(1)) - (mrSges(1,2) * g(2));
U = t44;
