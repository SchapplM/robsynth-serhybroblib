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
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:01
% EndTime: 2020-05-07 02:14:03
% DurationCPUTime: 1.24s
% Computational Cost: add. (449->122), mult. (381->160), div. (0->0), fcn. (220->52), ass. (0->67)
t81 = (m(5) + m(6));
t89 = (m(4) + m(8) + m(9) + t81);
t87 = -(2 * rSges(2,2) * m(2)) + (2 * rSges(3,3) * m(3)) + (2 * m(4) * rSges(4,3)) + (2 * m(5) * rSges(5,3)) + (2 * rSges(7,3) * m(7)) + (2 * m(8) * rSges(8,3)) + (2 * m(9) * rSges(9,3)) + 0.2e1 * (rSges(6,1) * sin(qJ(4)) + rSges(6,2) * cos(qJ(4))) * m(6);
t82 = m(9) * g(3);
t80 = m(4) * rSges(4,2);
t46 = sin(qJ(1));
t51 = cos(qJ(1));
t19 = g(1) * t51 + g(2) * t46;
t79 = m(6) * t19;
t77 = m(9) * t19;
t76 = rSges(7,1) * g(2);
t75 = rSges(3,2) * m(3);
t74 = rSges(7,2) * g(2);
t45 = sin(qJ(2));
t73 = rSges(6,2) * t45;
t50 = cos(qJ(2));
t72 = rSges(6,2) * t50;
t38 = pkin(15) + pkin(18);
t32 = qJ(2) + t38;
t24 = pkin(17) + qJ(3) + t32;
t20 = pkin(16) + t24;
t9 = atan2(-sin(t20), cos(t20));
t8 = qJ(3) + t9;
t71 = pkin(14) - qJ(2);
t10 = pkin(17) - atan2(-sin(t24), -cos(t24));
t39 = qJ(2) + qJ(3);
t66 = t79 / 0.2e1;
t6 = -qJ(4) + t8;
t5 = qJ(4) + t8;
t64 = -m(3) - t89;
t14 = t89 * pkin(1) + rSges(3,1) * m(3);
t62 = -m(2) - m(7) + t64;
t61 = (2 * m(7) * pkin(6)) - (2 * rSges(2,1) * m(2)) + (2 * t64 * pkin(12)) - 0.2e1 * t14 * t50 + 0.2e1 * t45 * t75;
t54 = rSges(7,1) * g(1);
t53 = rSges(7,2) * g(1);
t52 = cos(pkin(15));
t49 = cos(qJ(3));
t47 = sin(pkin(15));
t44 = sin(qJ(3));
t42 = cos(pkin(18));
t41 = sin(pkin(18));
t40 = -pkin(15) + pkin(14);
t37 = pkin(17) - qJ(2);
t36 = rSges(6,1) * t50;
t35 = rSges(6,1) * t45;
t34 = -qJ(1) + t71;
t33 = qJ(1) + t71;
t31 = pkin(8) * m(6) + (m(5) * rSges(5,1));
t30 = t53 - t76;
t29 = t53 + t76;
t28 = t54 - t74;
t27 = t54 + t74;
t26 = cos(t32);
t25 = sin(t32);
t23 = (-rSges(6,3) - pkin(10)) * m(6) + (m(5) * rSges(5,2));
t22 = t81 * pkin(4) + m(4) * rSges(4,1);
t21 = m(1) - t62;
t18 = rSges(9,1) * t50 - rSges(9,2) * t45;
t17 = rSges(9,1) * t45 + rSges(9,2) * t50;
t16 = -t41 * t47 + t42 * t52;
t15 = t41 * t52 + t42 * t47;
t7 = t9 + t39;
t4 = qJ(2) + t6;
t3 = qJ(2) + t5;
t2 = -atan2(t25, -t26) + t10;
t1 = pkin(17) - atan2(t15 * t50 + t16 * t45, t15 * t45 - t16 * t50);
t11 = (-g(1) * t87 + g(2) * t61) * t46 / 0.2e1 + (g(1) * t61 + g(2) * t87) * t51 / 0.2e1 - m(1) * (rSges(1,1) * g(1) + rSges(1,2) * g(2)) + ((-t17 * t26 + t18 * t25) * cos(t10) + (t17 * t25 + t18 * t26) * sin(t10)) * t82 + (t17 * sin(t2) + t18 * cos(t2)) * t77 + (cos(t4) + cos(t3)) * rSges(6,1) * t66 + (-g(1) * r_base(1) - g(2) * r_base(2)) * t21 + (-sin(t3) * t79 / 0.2e1 + sin(t4) * t66) * rSges(6,2) - ((t28 * t52 - t29 * t47) * cos(t33) + (t27 * t47 + t30 * t52) * sin(t34) + (t28 * t47 + t29 * t52) * sin(t33) + (t27 * t52 - t30 * t47) * cos(t34)) * m(7) / 0.2e1 + ((rSges(8,2) * sin(t38) + rSges(8,1) * cos(t38)) * m(8) + t22 * cos(t39) - t23 * sin(t7) + t31 * cos(t7) - sin(t39) * t80) * t19 + ((-cos(t1) * t45 + sin(t1) * t50) * t82 + (-sin(t37) * t25 + cos(t37) * t26) * t77) * pkin(3) + (-(rSges(1,3) * m(1)) - (m(2) * rSges(2,3)) + (pkin(11) * t62) + ((t35 - t72) * cos(t6) + (t35 + t72) * cos(t5) + (t36 - t73) * sin(t5) + (t36 + t73) * sin(t6)) * m(6) / 0.2e1 + ((rSges(8,1) * t45 + rSges(8,2) * t50) * t26 - (rSges(8,1) * t50 - rSges(8,2) * t45) * t25) * m(8) + (t22 * t44 + t49 * t80 - t75) * t50 - (-t22 * t49 + t44 * t80 + t14) * t45 + (t23 * t50 + t31 * t45) * cos(t8) + (-t23 * t45 + t31 * t50) * sin(t8) - (t21 * r_base(3)) + (-pkin(13) + (rSges(7,1) * t50 - rSges(7,2) * t45) * sin(t40) - (rSges(7,1) * t45 + rSges(7,2) * t50) * cos(t40)) * m(7)) * g(3);
U = t11;
