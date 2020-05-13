% Calculate Gravitation load on the joints for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:43:56
% EndTime: 2020-05-03 01:43:59
% DurationCPUTime: 0.46s
% Computational Cost: add. (448->90), mult. (474->115), div. (0->0), fcn. (290->14), ass. (0->55)
t36 = qJ(2) + qJ(3);
t22 = pkin(2) * cos(t36);
t43 = cos(qJ(2));
t35 = t43 * pkin(4);
t10 = t22 + t35 + pkin(1);
t33 = qJ(4) + t36;
t28 = cos(t33);
t19 = pkin(5) * t28;
t87 = t19 + t10;
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t86 = t41 * rSges(7,1) - rSges(7,2) * t37;
t39 = sin(qJ(2));
t42 = cos(qJ(3));
t30 = rSges(4,1) * m(4) * t42;
t38 = sin(qJ(3));
t72 = rSges(4,2) * t38;
t7 = m(3) * rSges(3,1) + t30 + (pkin(4) - t72) * m(4);
t12 = (rSges(4,1) * t38 + rSges(4,2) * t42) * m(4);
t9 = rSges(3,2) * m(3) + t12;
t85 = t9 * t39 - t7 * t43;
t40 = sin(qJ(1));
t44 = cos(qJ(1));
t62 = t44 * g(1) + t40 * g(2);
t14 = pkin(3) + t86;
t29 = qJ(5) + t33;
t20 = sin(t29);
t21 = cos(t29);
t84 = -t20 * rSges(7,3) + t14 * t21;
t27 = sin(t33);
t66 = pkin(2) * sin(t36);
t8 = -pkin(5) * t27 - t66;
t83 = -rSges(6,1) * t20 - rSges(6,2) * t21 + t8;
t82 = t21 * rSges(6,1) - t20 * rSges(6,2);
t81 = t28 * rSges(5,1) - t27 * rSges(5,2);
t80 = -rSges(5,1) * t27 - rSges(5,2) * t28 - t66;
t25 = m(6) * rSges(6,2) + m(7) * rSges(7,3);
t6 = m(6) * rSges(6,1) + t14 * m(7);
t79 = (g(3) * t25 + t62 * t6) * t20 + (-g(3) * t6 + t62 * t25) * t21;
t78 = m(5) * rSges(5,2);
t77 = t39 * pkin(4);
t67 = t19 + t22;
t65 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + rSges(4,3) * m(4);
t63 = t22 + t81;
t15 = t37 * rSges(7,1) + t41 * rSges(7,2);
t56 = t67 + t84;
t55 = t67 + t82;
t54 = t82 + t87;
t51 = t10 + t81;
t48 = -m(2) * rSges(2,1) + (-m(3) - m(4)) * pkin(1) + t85;
t47 = -rSges(7,3) * t21 - t14 * t20 + t8;
t46 = t84 + t87;
t17 = rSges(5,1) * m(5) + (m(6) + m(7)) * pkin(5);
t13 = -m(4) * t72 + t30;
t1 = [-(t48 * g(1) + t65 * g(2)) * t40 + (-t65 * g(1) + t48 * g(2)) * t44 - m(5) * ((g(1) * rSges(5,3) + g(2) * t51) * t44 + (g(2) * rSges(5,3) - g(1) * t51) * t40) - m(6) * ((g(1) * rSges(6,3) + g(2) * t54) * t44 + (g(2) * rSges(6,3) - g(1) * t54) * t40) - m(7) * ((-g(1) * t46 - g(2) * t15) * t40 + (-t15 * g(1) + g(2) * t46) * t44), (-m(5) * (t35 + t63) - m(6) * (t35 + t55) - m(7) * (t35 + t56) + t85) * g(3) + (t39 * t7 + t43 * t9 - m(5) * (-t77 + t80) - m(6) * (-t77 + t83) - m(7) * (t47 - t77)) * t62, (-m(5) * t63 - m(6) * t55 - m(7) * t56 + t12 * t39 - t13 * t43) * g(3) + (-m(5) * t80 - m(6) * t83 - m(7) * t47 + t12 * t43 + t13 * t39) * t62, (g(3) * t78 + t62 * t17) * t27 - (g(3) * t17 - t62 * t78) * t28 + t79, t79, (-t86 * (-t40 * g(1) + t44 * g(2)) + (g(3) * t20 + t62 * t21) * t15) * m(7)];
taug = t1(:);
