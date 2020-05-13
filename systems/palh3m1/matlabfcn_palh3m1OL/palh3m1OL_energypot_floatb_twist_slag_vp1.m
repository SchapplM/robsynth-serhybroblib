% Calculate potential energy for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1OL_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:03:43
% EndTime: 2020-04-20 17:03:44
% DurationCPUTime: 0.54s
% Computational Cost: add. (232->111), mult. (185->123), div. (0->0), fcn. (153->18), ass. (0->44)
t51 = pkin(10) + rSges(6,3);
t22 = qJ(2) + qJ(3);
t18 = qJ(4) + t22;
t11 = cos(t18);
t50 = t11 * pkin(8);
t24 = sin(qJ(5));
t30 = cos(qJ(1));
t48 = t24 * t30;
t26 = sin(qJ(1));
t47 = t26 * t24;
t28 = cos(qJ(5));
t46 = t26 * t28;
t45 = t28 * t30;
t29 = cos(qJ(2));
t9 = t29 * pkin(1) + pkin(12);
t21 = qJ(2) + qJ(7);
t16 = cos(t22);
t4 = -pkin(4) * t16 + t9;
t44 = t26 * t4 + r_base(2);
t43 = t30 * t4 + r_base(1);
t42 = t26 * t9 + r_base(2);
t41 = t30 * t9 + r_base(1);
t40 = pkin(11) + r_base(3);
t17 = pkin(15) - t21;
t25 = sin(qJ(2));
t39 = t25 * pkin(1) + t40;
t14 = sin(t22);
t38 = -rSges(4,1) * t16 + rSges(4,2) * t14;
t10 = sin(t18);
t37 = -rSges(5,1) * t11 + rSges(5,2) * t10;
t13 = sin(t21);
t15 = cos(t21);
t36 = rSges(8,1) * t15 - rSges(8,2) * t13;
t12 = -qJ(8) + t17;
t7 = sin(t12);
t8 = cos(t12);
t35 = -rSges(9,1) * t8 - rSges(9,2) * t7 + pkin(3) * cos(t17) + t9;
t23 = sin(qJ(6));
t27 = cos(qJ(6));
t34 = rSges(7,1) * t27 - rSges(7,2) * t23 - pkin(6);
t33 = rSges(3,1) * t29 - rSges(3,2) * t25 + pkin(12);
t32 = -pkin(4) * t14 + t39;
t31 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(3) * (rSges(3,1) * t25 + rSges(3,2) * t29 + t40) + (-g(2) * rSges(3,3) + g(1) * t33) * t30 + (g(1) * rSges(3,3) + g(2) * t33) * t26 + t31) - m(4) * (g(1) * (t26 * rSges(4,3) + t38 * t30 + t41) + g(2) * (-rSges(4,3) * t30 + t38 * t26 + t42) + g(3) * (-rSges(4,1) * t14 - rSges(4,2) * t16 + t39)) - m(5) * (g(1) * (t26 * rSges(5,3) + t37 * t30 + t43) + g(2) * (-rSges(5,3) * t30 + t37 * t26 + t44) + g(3) * (-rSges(5,1) * t10 - rSges(5,2) * t11 + t32)) - m(6) * (g(1) * (-t30 * t50 + (-t11 * t45 + t47) * rSges(6,1) + (t11 * t48 + t46) * rSges(6,2) + t43) + g(2) * (-t26 * t50 + (-t11 * t46 - t48) * rSges(6,1) + (t11 * t47 - t45) * rSges(6,2) + t44) + g(3) * (t51 * t11 + t32) + (g(3) * (-rSges(6,1) * t28 + rSges(6,2) * t24 - pkin(8)) - (g(1) * t30 + g(2) * t26) * t51) * t10) - m(7) * (g(3) * (rSges(7,1) * t23 + rSges(7,2) * t27 + pkin(13) + t40) + (-g(2) * rSges(7,3) + g(1) * t34) * t30 + (g(1) * rSges(7,3) + g(2) * t34) * t26 + t31) - m(8) * (g(1) * (t26 * rSges(8,3) + t36 * t30 + t41) + g(2) * (-rSges(8,3) * t30 + t36 * t26 + t42) + g(3) * (rSges(8,1) * t13 + rSges(8,2) * t15 + t39)) - m(9) * (g(3) * (-pkin(3) * sin(t17) + t7 * rSges(9,1) - t8 * rSges(9,2) + t39) + (-g(2) * rSges(9,3) + g(1) * t35) * t30 + (g(1) * rSges(9,3) + g(2) * t35) * t26 + t31);
U = t1;
