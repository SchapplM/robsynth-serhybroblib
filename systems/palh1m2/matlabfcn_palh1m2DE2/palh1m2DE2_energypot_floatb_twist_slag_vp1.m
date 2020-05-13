% Calculate potential energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:49
% EndTime: 2020-05-02 20:57:51
% DurationCPUTime: 0.88s
% Computational Cost: add. (325->98), mult. (372->123), div. (0->0), fcn. (238->40), ass. (0->51)
t59 = m(5) + m(6);
t50 = m(4) + m(8) + m(11) + t59;
t45 = m(9) + m(3) + m(10) + t50;
t30 = sin(qJ(4));
t36 = cos(qJ(4));
t61 = -m(2) * rSges(2,2) + m(11) * rSges(11,3) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t30 + rSges(6,2) * t36) * m(6) + rSges(7,3) * m(7) + m(8) * rSges(8,3) + rSges(9,3) * m(9) + rSges(10,3) * m(10);
t19 = t59 * pkin(5) + m(4) * rSges(4,1);
t31 = sin(qJ(3));
t35 = sin(pkin(17));
t37 = cos(qJ(3));
t41 = cos(pkin(17));
t34 = sin(pkin(18));
t40 = cos(pkin(18));
t46 = rSges(7,1) * t34 - rSges(7,2) * t40;
t47 = rSges(7,1) * t40 + rSges(7,2) * t34;
t57 = m(4) * rSges(4,2);
t2 = t37 * t57 + rSges(3,1) * m(3) + t19 * t31 + (t46 * t35 + t47 * t41) * m(7) + t50 * pkin(1);
t32 = sin(qJ(2));
t38 = cos(qJ(2));
t4 = t31 * t57 + rSges(3,2) * m(3) - t19 * t37 + (t47 * t35 - t46 * t41) * m(7);
t60 = -m(2) * rSges(2,1) + pkin(14) * m(7) - t45 * pkin(15) + t2 * t32 + t38 * t4;
t56 = m(9) * rSges(9,2);
t27 = sin(pkin(19));
t29 = cos(pkin(19));
t13 = t27 * t37 + t29 * t31;
t14 = t27 * t31 - t29 * t37;
t7 = qJ(2) + atan2(t14, t13);
t25 = pkin(18) - pkin(22);
t24 = -qJ(2) + t25;
t54 = pkin(21) - atan2(cos(t24), -sin(t24));
t52 = -qJ(2) + t54;
t51 = -pkin(21) + t25;
t33 = sin(qJ(1));
t39 = cos(qJ(1));
t18 = g(1) * t39 + g(2) * t33;
t44 = m(7) + m(2) + t45;
t28 = cos(pkin(22));
t26 = sin(pkin(22));
t23 = pkin(2) * m(10) + m(9) * rSges(9,1);
t22 = -pkin(20) + t51;
t21 = (-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2);
t20 = -qJ(2) - qJ(3) + t51;
t16 = t26 * t40 - t28 * t34;
t15 = t26 * t34 + t28 * t40;
t11 = m(5) * rSges(5,1) + (rSges(6,1) * t36 - rSges(6,2) * t30 + pkin(9)) * m(6);
t10 = atan2(-sin(t20), cos(t20));
t6 = -t10 + t54;
t5 = -t10 + t52;
t3 = atan2(t14, -t13) + t7;
t1 = pkin(21) - atan2(t15 * t38 - t16 * t32, t15 * t32 + t16 * t38);
t8 = (t60 * g(1) + t61 * g(2)) * t39 + (-t61 * g(1) + t60 * g(2)) * t33 + (-g(3) * t23 + t18 * t56) * cos(t7) + (g(3) * t56 + t18 * t23) * sin(t7) + (g(3) * t21 + t18 * t11) * cos(t22) + (g(3) * t11 - t18 * t21) * sin(t22) - g(3) * t2 * t38 + g(3) * t4 * t32 - g(3) * (m(2) * rSges(2,3) - pkin(16) * m(7) + t44 * pkin(13)) + (-(-g(3) * rSges(10,1) + t18 * rSges(10,2)) * cos(t3) - (t18 * rSges(10,1) + g(3) * rSges(10,2)) * sin(t3)) * m(10) + ((t18 * rSges(8,1) + g(3) * rSges(8,2)) * cos(t25) - (-g(3) * rSges(8,1) + t18 * rSges(8,2)) * sin(t25)) * m(8) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - g(3) * rSges(1,3)) * m(1) + (-g(1) * r_base(1) - g(2) * r_base(2) - g(3) * r_base(3)) * (m(1) + t44) + ((-pkin(4) * sin(t52) - rSges(11,2) * cos(t5) + rSges(11,1) * sin(t5)) * t18 + ((rSges(11,1) * t38 - rSges(11,2) * t32) * cos(t6) + (rSges(11,1) * t32 + rSges(11,2) * t38) * sin(t6) + (-sin(t1) * t32 - cos(t1) * t38) * pkin(4)) * g(3)) * m(11);
U = t8;
