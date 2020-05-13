% Calculate Gravitation load on the joints for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [12x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:54
% EndTime: 2020-05-09 23:18:56
% DurationCPUTime: 0.30s
% Computational Cost: add. (356->79), mult. (212->82), div. (0->0), fcn. (138->20), ass. (0->53)
t70 = pkin(6) * m(10) + mrSges(4,1);
t44 = qJ(1) + qJ(2);
t38 = qJ(6) + t44;
t24 = sin(t38);
t27 = cos(t38);
t2 = t27 * mrSges(7,1) - t24 * mrSges(7,2);
t35 = sin(t44);
t37 = cos(t44);
t39 = qJ(4) + t44;
t25 = sin(t39);
t28 = cos(t39);
t30 = qJ(10) + t39;
t14 = sin(t30);
t15 = cos(t30);
t58 = t15 * mrSges(11,1) - t14 * mrSges(11,2);
t51 = -t25 * mrSges(5,2) + (pkin(4) * m(11) + mrSges(5,1)) * t28 - t58;
t63 = m(11) + m(5);
t65 = m(4) + m(10);
t40 = qJ(3) + t44;
t26 = sin(t40);
t29 = cos(t40);
t31 = qJ(9) + t40;
t16 = sin(t31);
t17 = cos(t31);
t59 = -t17 * mrSges(10,1) + t16 * mrSges(10,2);
t66 = t26 * mrSges(4,2) - t70 * t29 - t59;
t68 = -t2 + (t65 * pkin(2) + t63 * pkin(3) + mrSges(3,1)) * t37 - t35 * mrSges(3,2) + t51 + t66;
t1 = -t24 * mrSges(7,1) - t27 * mrSges(7,2);
t54 = -t14 * mrSges(11,1) - t15 * mrSges(11,2);
t52 = -t25 * mrSges(5,1) - t28 * mrSges(5,2) - t54;
t64 = t16 * mrSges(10,1) + t17 * mrSges(10,2);
t56 = t29 * mrSges(4,2) - t64;
t67 = -t35 * mrSges(3,1) + t26 * mrSges(4,1) - t37 * mrSges(3,2) - t1 + t52 + t56;
t13 = pkin(4) * t25;
t22 = pkin(3) * t35;
t46 = sin(qJ(1));
t41 = t46 * pkin(1);
t62 = t22 + t41;
t61 = -m(3) - m(7) - m(9);
t23 = pkin(2) * t35;
t60 = -pkin(6) * t26 + t23;
t43 = qJ(1) + qJ(8);
t34 = sin(t43);
t36 = cos(t43);
t57 = t36 * mrSges(9,1) - t34 * mrSges(9,2);
t55 = -t34 * mrSges(9,1) - t36 * mrSges(9,2);
t48 = cos(qJ(1));
t47 = cos(qJ(7));
t45 = sin(qJ(7));
t42 = pkin(8) + qJ(5);
t33 = cos(t42);
t32 = sin(t42);
t3 = [((mrSges(2,1) + (-t61 + t63 + t65) * pkin(1)) * t48 - t46 * mrSges(2,2) - t57 + t68) * g(2) + ((t61 * pkin(1) - mrSges(2,1)) * t46 - m(11) * (t13 + t62) - m(4) * (t23 + t41) - m(5) * t62 - t48 * mrSges(2,2) - m(10) * (t41 + t60) - t55 + t67) * g(1), t68 * g(2) + ((-pkin(2) * m(4) - pkin(3) * m(5)) * t35 - m(11) * (t22 + t13) - m(10) * t60 + t67) * g(1), t66 * g(2) + (t70 * t26 + t56) * g(1), t51 * g(2) + (-m(11) * t13 + t52) * g(1), -g(1) * (-t32 * mrSges(6,1) - t33 * mrSges(6,2)) - g(2) * (t33 * mrSges(6,1) - t32 * mrSges(6,2)), -g(1) * t1 - g(2) * t2, -g(1) * (t47 * mrSges(8,1) - t45 * mrSges(8,2)) - g(2) * (t45 * mrSges(8,1) + t47 * mrSges(8,2)), -g(1) * t55 - g(2) * t57, -g(1) * t64 - g(2) * t59, -g(1) * t54 - g(2) * t58, 0, 0];
taug = t3(:);
