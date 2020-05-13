% Calculate potential energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:20
% EndTime: 2020-05-02 21:17:21
% DurationCPUTime: 0.97s
% Computational Cost: add. (267->121), mult. (259->140), div. (0->0), fcn. (148->24), ass. (0->50)
t62 = -m(3) - m(4);
t44 = r_base(1) * g(1) + r_base(2) * g(2);
t36 = sin(qJ(1));
t41 = cos(qJ(1));
t48 = g(1) * t41 + t36 * g(2);
t33 = sin(qJ(5));
t38 = cos(qJ(5));
t61 = rSges(6,1) * t38 - rSges(6,2) * t33 + pkin(9);
t60 = pkin(15) * g(1);
t59 = pkin(15) * g(2);
t56 = -rSges(6,3) - pkin(11);
t31 = qJ(2) + qJ(3);
t30 = qJ(2) + qJ(7);
t29 = qJ(2) + qJ(8);
t35 = sin(qJ(2));
t16 = -t35 * pkin(1) + pkin(15);
t6 = pkin(5) * cos(t31) + t16;
t53 = t36 * t6 + r_base(2);
t52 = t41 * t6 + r_base(1);
t28 = pkin(13) + r_base(3);
t51 = -rSges(2,2) * m(2) + rSges(3,3) * m(3) + rSges(4,3) * m(4);
t40 = cos(qJ(2));
t50 = t40 * pkin(1) + t28;
t19 = pkin(19) - t30;
t49 = pkin(5) * sin(t31) + t50;
t25 = qJ(4) + t31;
t17 = sin(t25);
t18 = cos(t25);
t47 = rSges(5,1) * t18 - rSges(5,2) * t17;
t21 = sin(t30);
t23 = cos(t30);
t46 = -rSges(8,1) * t21 - rSges(8,2) * t23 + t16;
t12 = -qJ(10) + t19;
t10 = sin(t12);
t11 = cos(t12);
t45 = -rSges(11,1) * t10 + rSges(11,2) * t11 + pkin(4) * sin(t19) + t16;
t34 = sin(qJ(3));
t39 = cos(qJ(3));
t1 = rSges(3,1) * m(3) + (rSges(4,1) * t34 + rSges(4,2) * t39 + pkin(1)) * m(4);
t4 = -rSges(3,2) * m(3) + (rSges(4,1) * t39 - rSges(4,2) * t34) * m(4);
t43 = -rSges(2,1) * m(2) + pkin(15) * t62 + t1 * t35 - t4 * t40;
t42 = g(3) * t28 + t44;
t37 = cos(qJ(6));
t32 = sin(qJ(6));
t24 = qJ(9) + t29;
t22 = cos(t29);
t20 = sin(t29);
t8 = rSges(6,1) * t33 + rSges(6,2) * t38;
t7 = rSges(7,1) * t37 - rSges(7,2) * t32 - pkin(14);
t2 = (t43 * g(1) + t51 * g(2)) * t41 + (-t51 * g(1) + t43 * g(2)) * t36 - g(3) * t4 * t35 - g(3) * t1 * t40 - m(4) * pkin(13) * g(3) - (pkin(13) * m(3) + (pkin(13) + rSges(2,3)) * m(2)) * g(3) - m(5) * (g(1) * (t36 * rSges(5,3) + t47 * t41 + t52) + g(2) * (-rSges(5,3) * t41 + t47 * t36 + t53) + g(3) * (rSges(5,1) * t17 + rSges(5,2) * t18 + t49)) - m(6) * (g(1) * (t36 * t8 + t52) + g(2) * (-t41 * t8 + t53) + g(3) * t49 + (g(3) * t56 + t48 * t61) * t18 + (g(3) * t61 - t48 * t56) * t17) - m(7) * (g(1) * (t36 * rSges(7,3) + t41 * t7 + r_base(1)) + g(2) * (-rSges(7,3) * t41 + t7 * t36 + r_base(2)) + g(3) * (rSges(7,1) * t32 + rSges(7,2) * t37 - pkin(16) + t28)) - m(8) * (g(3) * (rSges(8,1) * t23 - rSges(8,2) * t21 + t50) + (-g(2) * rSges(8,3) + g(1) * t46) * t41 + (g(1) * rSges(8,3) + g(2) * t46) * t36 + t44) - ((g(3) * rSges(9,1) - t48 * rSges(9,2)) * t22 + (-t48 * rSges(9,1) - g(3) * rSges(9,2)) * t20 + (-rSges(9,3) * g(2) + t60) * t41 + (rSges(9,3) * g(1) + t59) * t36 + t42) * m(9) + m(10) * ((g(3) * rSges(10,1) - t48 * rSges(10,2)) * cos(t24) + (-t48 * rSges(10,1) - g(3) * rSges(10,2)) * sin(t24) + (rSges(10,3) * g(2) - t60) * t41 + (-rSges(10,3) * g(1) - t59) * t36 + (-g(3) * t22 + t48 * t20) * pkin(2) - t42) - m(11) * (g(3) * (pkin(4) * cos(t19) - t11 * rSges(11,1) - t10 * rSges(11,2) + t50) + (-g(2) * rSges(11,3) + g(1) * t45) * t41 + (g(1) * rSges(11,3) + g(2) * t45) * t36 + t44) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1) + (-g(3) * r_base(3) - t44) * (m(1) + m(2) - t62);
U = t2;
