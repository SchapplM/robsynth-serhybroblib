% Calculate Gravitation load on the joints for
% palh3m2TE
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2TE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:42:09
% EndTime: 2020-05-07 01:42:10
% DurationCPUTime: 0.34s
% Computational Cost: add. (295->95), mult. (493->143), div. (0->0), fcn. (270->22), ass. (0->53)
t66 = m(5) + m(6);
t62 = m(4) + t66;
t39 = sin(qJ(4));
t45 = cos(qJ(4));
t58 = t39 * rSges(6,1) + t45 * rSges(6,2);
t67 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + t58 * m(6) + m(7) * rSges(7,3);
t65 = m(4) * rSges(4,2);
t64 = rSges(3,2) * m(3);
t47 = cos(qJ(2));
t63 = pkin(1) * t47 + pkin(12);
t61 = (-rSges(6,3) - pkin(10)) * m(6) + m(5) * rSges(5,2);
t42 = sin(qJ(1));
t48 = cos(qJ(1));
t23 = t48 * g(1) + t42 * g(2);
t60 = t42 * g(1) - t48 * g(2);
t59 = t45 * rSges(6,1) - t39 * rSges(6,2);
t43 = sin(pkin(15));
t49 = cos(pkin(15));
t57 = t49 * rSges(7,1) - t43 * rSges(7,2);
t56 = t43 * rSges(7,1) + t49 * rSges(7,2);
t33 = qJ(3) + qJ(2);
t36 = sin(pkin(18));
t38 = cos(pkin(18));
t55 = rSges(9,1) * cos(t33) - rSges(9,2) * sin(t33) - (-(-t36 * t43 + t38 * t49) * cos(pkin(17)) + (t36 * t49 + t43 * t38) * sin(pkin(17))) * pkin(3) - t63;
t24 = pkin(4) * t66 + m(4) * rSges(4,1);
t54 = m(6) * t59;
t44 = sin(pkin(14));
t50 = cos(pkin(14));
t21 = rSges(7,1) * t44 - rSges(7,2) * t50;
t22 = rSges(7,1) * t50 + rSges(7,2) * t44;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t46 = cos(qJ(3));
t52 = rSges(3,1) * m(3);
t53 = -t41 * (t46 * t65 - t64 + t24 * t40 + (t21 * t49 - t22 * t43) * m(7)) - t47 * (t40 * t65 - t24 * t46 + t52 + (t21 * t43 + t22 * t49) * m(7) + t62 * pkin(1)) + pkin(6) * m(7) - m(2) * rSges(2,1) + (-m(3) - t62) * pkin(12);
t35 = cos(pkin(16));
t34 = sin(pkin(16));
t32 = pkin(17) + pkin(18);
t28 = cos(t32);
t27 = sin(t32);
t26 = pkin(8) * m(6) + m(5) * rSges(5,1);
t25 = m(9) * rSges(9,2) + t65;
t19 = m(9) * rSges(9,1) + t24;
t14 = -t43 * g(3) + t23 * t49;
t13 = t49 * g(3) + t43 * t23;
t12 = t34 * t26 - t61 * t35;
t11 = t26 * t35 + t61 * t34;
t10 = (-rSges(8,1) * t49 - rSges(8,2) * t43) * t38 + t36 * (rSges(8,1) * t43 - rSges(8,2) * t49) + t63;
t9 = t19 * g(3) - t23 * t25;
t8 = g(3) * t25 + t23 * t19;
t4 = -t64 + t19 * t40 + t25 * t46 + (t57 * t44 - t56 * t50) * m(7);
t3 = -t19 * t46 + t25 * t40 + (m(8) + m(9) + t62) * pkin(1) + t52 + (t56 * t44 + t57 * t50) * m(7);
t1 = [-(t53 * g(1) + t67 * g(2)) * t42 + (-t67 * g(1) + t53 * g(2)) * t48 - m(8) * (g(1) * (t48 * rSges(8,3) - t10 * t42) + g(2) * (t42 * rSges(8,3) + t10 * t48)) - m(9) * ((g(1) * rSges(9,3) - g(2) * t55) * t48 + (g(2) * rSges(9,3) + g(1) * t55) * t42) + (-t28 * (t11 * t49 - t12 * t43 + (-t34 * t43 + t35 * t49) * t54) + t27 * (t11 * t43 + t12 * t49 + (t34 * t49 + t35 * t43) * t54)) * t60, (-t4 * g(3) + t23 * t3) * t41 - t47 * (g(3) * t3 + t23 * t4), (-t8 * t40 + t9 * t46) * t47 - t41 * (t40 * t9 + t8 * t46), -m(6) * (t59 * t60 + ((-t34 * t13 + t14 * t35) * t28 - (t13 * t35 + t34 * t14) * t27) * t58)];
taug = t1(:);
