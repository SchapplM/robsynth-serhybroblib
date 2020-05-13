% Calculate potential energy for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:17
% EndTime: 2020-04-15 19:28:18
% DurationCPUTime: 0.71s
% Computational Cost: add. (284->137), mult. (224->151), div. (0->0), fcn. (184->22), ass. (0->54)
t64 = pkin(11) + rSges(6,3);
t31 = qJ(2) + qJ(3);
t25 = qJ(4) + t31;
t16 = cos(t25);
t63 = t16 * pkin(9);
t33 = sin(qJ(5));
t39 = cos(qJ(1));
t61 = t33 * t39;
t35 = sin(qJ(1));
t60 = t35 * t33;
t37 = cos(qJ(5));
t59 = t35 * t37;
t58 = t37 * t39;
t30 = qJ(2) + qJ(7);
t29 = qJ(2) + qJ(8);
t34 = sin(qJ(2));
t12 = -t34 * pkin(1) + pkin(15);
t23 = cos(t31);
t4 = pkin(5) * t23 + t12;
t57 = t35 * t4 + r_base(2);
t56 = t39 * t4 + r_base(1);
t55 = t35 * t12 + r_base(2);
t54 = t39 * t12 + r_base(1);
t53 = pkin(13) + r_base(3);
t52 = t35 * pkin(15) + r_base(2);
t51 = t39 * pkin(15) + r_base(1);
t38 = cos(qJ(2));
t50 = t38 * pkin(1) + t53;
t17 = pkin(19) - t30;
t20 = sin(t31);
t49 = pkin(5) * t20 + t50;
t48 = -rSges(3,1) * t34 - rSges(3,2) * t38;
t47 = rSges(4,1) * t23 - rSges(4,2) * t20;
t14 = sin(t25);
t46 = rSges(5,1) * t16 - rSges(5,2) * t14;
t19 = sin(t30);
t22 = cos(t30);
t45 = -rSges(8,1) * t19 - rSges(8,2) * t22;
t18 = sin(t29);
t21 = cos(t29);
t44 = -rSges(9,1) * t18 - rSges(9,2) * t21;
t10 = -qJ(10) + t17;
t8 = sin(t10);
t9 = cos(t10);
t43 = -rSges(11,1) * t8 + rSges(11,2) * t9 + pkin(4) * sin(t17) + t12;
t24 = qJ(9) + t29;
t13 = sin(t24);
t15 = cos(t24);
t42 = -pkin(2) * t18 + rSges(10,1) * t13 + rSges(10,2) * t15 + pkin(15);
t32 = sin(qJ(6));
t36 = cos(qJ(6));
t41 = rSges(7,1) * t36 - rSges(7,2) * t32 - pkin(14);
t40 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t39 - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + rSges(2,2) * t39 + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t35 * rSges(3,3) + t39 * t48 + t51) + g(2) * (-rSges(3,3) * t39 + t35 * t48 + t52) + g(3) * (rSges(3,1) * t38 - rSges(3,2) * t34 + t53)) - m(4) * (g(1) * (t35 * rSges(4,3) + t39 * t47 + t54) + g(2) * (-rSges(4,3) * t39 + t35 * t47 + t55) + g(3) * (rSges(4,1) * t20 + rSges(4,2) * t23 + t50)) - m(5) * (g(1) * (t35 * rSges(5,3) + t39 * t46 + t56) + g(2) * (-rSges(5,3) * t39 + t35 * t46 + t57) + g(3) * (rSges(5,1) * t14 + rSges(5,2) * t16 + t49)) - m(6) * (g(1) * (t39 * t63 + (t16 * t58 + t60) * rSges(6,1) + (-t16 * t61 + t59) * rSges(6,2) + t56) + g(2) * (t35 * t63 + (t16 * t59 - t61) * rSges(6,1) + (-t16 * t60 - t58) * rSges(6,2) + t57) + g(3) * (-t64 * t16 + t49) + (g(3) * (rSges(6,1) * t37 - rSges(6,2) * t33 + pkin(9)) + (g(1) * t39 + g(2) * t35) * t64) * t14) - m(7) * (g(3) * (rSges(7,1) * t32 + rSges(7,2) * t36 - pkin(16) + t53) + (-g(2) * rSges(7,3) + g(1) * t41) * t39 + (g(1) * rSges(7,3) + g(2) * t41) * t35 + t40) - m(8) * (g(1) * (t35 * rSges(8,3) + t39 * t45 + t54) + g(2) * (-rSges(8,3) * t39 + t35 * t45 + t55) + g(3) * (rSges(8,1) * t22 - rSges(8,2) * t19 + t50)) - m(9) * (g(1) * (t35 * rSges(9,3) + t39 * t44 + t51) + g(2) * (-rSges(9,3) * t39 + t35 * t44 + t52) + g(3) * (t21 * rSges(9,1) - rSges(9,2) * t18 + t53)) - m(10) * (g(3) * (pkin(2) * t21 - rSges(10,1) * t15 + rSges(10,2) * t13 + t53) + (-g(2) * rSges(10,3) + g(1) * t42) * t39 + (g(1) * rSges(10,3) + g(2) * t42) * t35 + t40) - m(11) * (g(3) * (pkin(4) * cos(t17) - t9 * rSges(11,1) - t8 * rSges(11,2) + t50) + (-g(2) * rSges(11,3) + g(1) * t43) * t39 + (g(1) * rSges(11,3) + g(2) * t43) * t35 + t40);
U = t1;
