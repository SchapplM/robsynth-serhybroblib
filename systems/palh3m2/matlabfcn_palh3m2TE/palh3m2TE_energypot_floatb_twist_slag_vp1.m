% Calculate potential energy for
% palh3m2TE
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:47
% EndTime: 2020-05-07 01:41:48
% DurationCPUTime: 0.52s
% Computational Cost: add. (264->97), mult. (408->130), div. (0->0), fcn. (237->22), ass. (0->49)
t60 = g(1) * r_base(1) + g(2) * r_base(2);
t30 = sin(qJ(4));
t36 = cos(qJ(4));
t59 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t30 + rSges(6,2) * t36) * m(6) + rSges(7,3) * m(7);
t58 = m(5) + m(6);
t57 = m(4) * rSges(4,2);
t38 = cos(qJ(2));
t56 = pkin(1) * t38 + pkin(12);
t55 = -m(3) - m(4) - m(5);
t32 = sin(qJ(2));
t52 = t32 * pkin(1) + pkin(11) + r_base(3);
t51 = m(2) + m(7) - t55;
t50 = (-rSges(6,3) - pkin(10)) * m(6) + m(5) * rSges(5,2);
t33 = sin(qJ(1));
t39 = cos(qJ(1));
t49 = g(1) * t39 + g(2) * t33;
t26 = sin(pkin(18));
t28 = cos(pkin(18));
t34 = sin(pkin(15));
t40 = cos(pkin(15));
t47 = t26 * t40 + t28 * t34;
t46 = -t26 * t34 + t28 * t40;
t23 = qJ(3) + qJ(2);
t18 = sin(t23);
t19 = cos(t23);
t27 = sin(pkin(17));
t29 = cos(pkin(17));
t45 = -rSges(9,1) * t19 + rSges(9,2) * t18 + (t27 * t47 - t29 * t46) * pkin(3) + t56;
t44 = m(6) * (rSges(6,1) * t36 - rSges(6,2) * t30);
t35 = sin(pkin(14));
t41 = cos(pkin(14));
t13 = rSges(7,1) * t35 - rSges(7,2) * t41;
t14 = rSges(7,1) * t41 + rSges(7,2) * t35;
t15 = pkin(4) * t58 + m(4) * rSges(4,1);
t31 = sin(qJ(3));
t37 = cos(qJ(3));
t3 = t31 * t57 + rSges(3,1) * m(3) - t15 * t37 + (t13 * t34 + t14 * t40) * m(7) + (m(4) + t58) * pkin(1);
t4 = t37 * t57 - rSges(3,2) * m(3) + t15 * t31 + (t13 * t40 - t14 * t34) * m(7);
t43 = m(7) * pkin(6) - m(2) * rSges(2,1) - t3 * t38 - t32 * t4 + (-m(6) + t55) * pkin(12);
t25 = cos(pkin(16));
t24 = sin(pkin(16));
t22 = pkin(17) + pkin(18);
t17 = pkin(8) * m(6) + m(5) * rSges(5,1);
t8 = t24 * t17 - t25 * t50;
t7 = t17 * t25 + t24 * t50;
t6 = (-rSges(8,1) * t40 - rSges(8,2) * t34) * t28 + t26 * (rSges(8,1) * t34 - rSges(8,2) * t40) + t56;
t2 = -t8 * t34 + t7 * t40 + (-t24 * t34 + t25 * t40) * t44;
t1 = t7 * t34 + t8 * t40 + (t24 * t40 + t25 * t34) * t44;
t5 = (-g(3) * t1 + t2 * t49) * cos(t22) + (-t2 * g(3) - t49 * t1) * sin(t22) + (t43 * g(1) + t59 * g(2)) * t39 + (-t59 * g(1) + t43 * g(2)) * t33 + g(3) * t4 * t38 - g(3) * t3 * t32 - m(6) * pkin(11) * g(3) - (m(2) * rSges(2,3) + pkin(13) * m(7) + t51 * pkin(11)) * g(3) - m(8) * (g(1) * (rSges(8,3) * t33 + t39 * t6 + r_base(1)) + g(2) * (-rSges(8,3) * t39 + t33 * t6 + r_base(2)) + g(3) * ((rSges(8,1) * t26 - rSges(8,2) * t28) * t40 + t34 * (rSges(8,1) * t28 + rSges(8,2) * t26) + t52)) - m(9) * ((-g(2) * rSges(9,3) + g(1) * t45) * t39 + (g(1) * rSges(9,3) + g(2) * t45) * t33 + (-t18 * rSges(9,1) - t19 * rSges(9,2) + t52 + (t27 * t46 + t29 * t47) * pkin(3)) * g(3) + t60) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1) + (-g(3) * r_base(3) - t60) * (m(6) + m(1) + t51);
U = t5;
