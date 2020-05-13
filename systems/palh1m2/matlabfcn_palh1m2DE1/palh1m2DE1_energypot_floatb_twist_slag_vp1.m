% Calculate potential energy for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:44
% EndTime: 2020-05-01 20:55:45
% DurationCPUTime: 0.81s
% Computational Cost: add. (322->126), mult. (507->170), div. (0->0), fcn. (312->24), ass. (0->59)
t53 = g(1) * r_base(1) + g(2) * r_base(2);
t38 = sin(qJ(4));
t44 = cos(qJ(4));
t69 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t38 + rSges(6,2) * t44) * m(6) + rSges(7,3) * m(7);
t68 = m(5) + m(6);
t67 = m(4) * rSges(4,2);
t66 = m(3) + m(4) + m(5);
t63 = pkin(13) + r_base(3);
t40 = sin(qJ(2));
t62 = -pkin(1) * t40 + pkin(15);
t46 = cos(qJ(2));
t61 = t46 * pkin(1) + t63;
t60 = m(7) + m(2) + t66;
t59 = (-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2);
t41 = sin(qJ(1));
t47 = cos(qJ(1));
t58 = g(1) * t47 + g(2) * t41;
t30 = sin(pkin(22));
t34 = cos(pkin(22));
t42 = sin(pkin(18));
t48 = cos(pkin(18));
t56 = t48 * t30 - t42 * t34;
t55 = t42 * t30 + t34 * t48;
t29 = qJ(3) + qJ(2);
t25 = sin(t29);
t26 = cos(t29);
t31 = sin(pkin(21));
t35 = cos(pkin(21));
t54 = rSges(11,1) * t26 - rSges(11,2) * t25 + (t31 * t56 - t35 * t55) * pkin(4) + t62;
t52 = m(6) * (rSges(6,1) * t44 - rSges(6,2) * t38);
t43 = sin(pkin(17));
t49 = cos(pkin(17));
t20 = rSges(7,1) * t43 + rSges(7,2) * t49;
t21 = rSges(7,1) * t49 - rSges(7,2) * t43;
t22 = pkin(5) * t68 + m(4) * rSges(4,1);
t39 = sin(qJ(3));
t45 = cos(qJ(3));
t3 = t45 * t67 + rSges(3,1) * m(3) + t22 * t39 + (t20 * t42 + t21 * t48) * m(7) + (m(4) + t68) * pkin(1);
t4 = -t39 * t67 - rSges(3,2) * m(3) + t22 * t45 + (-t20 * t48 + t21 * t42) * m(7);
t51 = -m(2) * rSges(2,1) + pkin(14) * m(7) + t3 * t40 - t4 * t46 + (-m(6) - t66) * pkin(15);
t37 = cos(pkin(19));
t36 = cos(pkin(20));
t33 = sin(pkin(19));
t32 = sin(pkin(20));
t28 = pkin(22) + pkin(21);
t24 = pkin(9) * m(6) + m(5) * rSges(5,1);
t19 = t33 * rSges(9,1) + rSges(9,2) * t37;
t18 = rSges(9,1) * t37 - t33 * rSges(9,2);
t13 = -rSges(10,2) + (-t33 * t39 + t37 * t45) * pkin(2);
t12 = rSges(10,1) + (t33 * t45 + t37 * t39) * pkin(2);
t11 = t32 * t24 - t36 * t59;
t10 = t24 * t36 + t32 * t59;
t9 = (-rSges(8,1) * t48 + rSges(8,2) * t42) * t34 + (-rSges(8,1) * t42 - rSges(8,2) * t48) * t30 + t62;
t8 = -t12 * t40 + t13 * t46 + pkin(15);
t7 = -g(3) * t18 + t19 * t58;
t6 = g(3) * t19 + t18 * t58;
t2 = t10 * t48 + t11 * t42 + (t42 * t32 + t48 * t36) * t52;
t1 = -t10 * t42 + t11 * t48 + (t48 * t32 - t42 * t36) * t52;
t5 = (-t1 * g(3) + t2 * t58) * cos(t28) + (-t2 * g(3) - t1 * t58) * sin(t28) + (t51 * g(1) + t69 * g(2)) * t47 + (-t69 * g(1) + t51 * g(2)) * t41 - t4 * g(3) * t40 - t3 * g(3) * t46 - m(6) * pkin(13) * g(3) - (m(2) * rSges(2,3) - pkin(16) * m(7) + t60 * pkin(13)) * g(3) - m(8) * (g(1) * (t41 * rSges(8,3) + t9 * t47 + r_base(1)) + g(2) * (-t47 * rSges(8,3) + t9 * t41 + r_base(2)) + g(3) * (-(rSges(8,1) * t34 + rSges(8,2) * t30) * t42 + (t30 * rSges(8,1) - t34 * rSges(8,2)) * t48 + t61)) - ((-t39 * t7 + t6 * t45) * t46 + (-rSges(9,3) * g(2) + pkin(15) * g(1)) * t47 + (rSges(9,3) * g(1) + g(2) * pkin(15)) * t41 + g(3) * t63 + (-t39 * t6 - t7 * t45) * t40 + t53) * m(9) - m(10) * (g(1) * (t41 * rSges(10,3) + t8 * t47 + r_base(1)) + g(2) * (-t47 * rSges(10,3) + t8 * t41 + r_base(2)) + g(3) * (t12 * t46 + t13 * t40 + t63)) - m(11) * ((-g(2) * rSges(11,3) + g(1) * t54) * t47 + (g(1) * rSges(11,3) + g(2) * t54) * t41 + t53 + (t25 * rSges(11,1) + t26 * rSges(11,2) + t61 + (t31 * t55 + t35 * t56) * pkin(4)) * g(3)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1) + (-g(3) * r_base(3) - t53) * (m(6) + m(1) + t60);
U = t5;
