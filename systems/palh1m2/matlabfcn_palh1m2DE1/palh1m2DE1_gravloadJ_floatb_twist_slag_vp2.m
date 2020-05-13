% Calculate Gravitation load on the joints for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2DE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:37
% EndTime: 2020-05-01 20:57:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (495->77), mult. (527->99), div. (0->0), fcn. (376->20), ass. (0->52)
t74 = m(5) + m(6);
t39 = sin(qJ(4));
t45 = cos(qJ(4));
t73 = mrSges(6,1) * t39 + mrSges(6,2) * t45;
t21 = pkin(9) * m(6) + mrSges(6,1) * t45 - mrSges(6,2) * t39 + mrSges(5,1);
t29 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t33 = sin(pkin(20));
t37 = cos(pkin(20));
t12 = m(11) * pkin(4) + t21 * t37 - t29 * t33;
t32 = sin(pkin(21));
t36 = cos(pkin(21));
t55 = t33 * t21 + t29 * t37;
t72 = t12 * t32 + t55 * t36 - mrSges(8,2);
t71 = m(11) + m(4) + m(8) + t74;
t30 = pkin(2) * m(10) + mrSges(9,1);
t34 = sin(pkin(19));
t38 = cos(pkin(19));
t20 = mrSges(9,2) * t38 + t30 * t34 + mrSges(11,2) + mrSges(4,2);
t40 = sin(qJ(3));
t46 = cos(qJ(3));
t52 = t74 * pkin(5) - t34 * mrSges(9,2) + t30 * t38 + mrSges(11,1) + mrSges(4,1);
t68 = t71 * pkin(1) + t20 * t46 + t40 * t52 + mrSges(3,1) + mrSges(10,1);
t67 = -t20 * t40 + t46 * t52 - mrSges(3,2) - mrSges(10,2);
t43 = sin(pkin(18));
t49 = cos(pkin(18));
t42 = sin(qJ(1));
t48 = cos(qJ(1));
t59 = g(1) * t48 + t42 * g(2);
t18 = t43 * g(3) + t59 * t49;
t19 = g(3) * t49 - t43 * t59;
t31 = sin(pkin(22));
t35 = cos(pkin(22));
t11 = t31 * t18 + t19 * t35;
t56 = t18 * t35 - t31 * t19;
t66 = (-t11 * t32 + t56 * t36) * t37;
t58 = mrSges(7,1) * t49 + mrSges(7,2) * t43;
t57 = mrSges(7,1) * t43 - mrSges(7,2) * t49;
t54 = -mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3) + mrSges(11,3) + t73;
t44 = sin(pkin(17));
t50 = cos(pkin(17));
t24 = mrSges(7,1) * t44 + mrSges(7,2) * t50;
t25 = mrSges(7,1) * t50 - mrSges(7,2) * t44;
t41 = sin(qJ(2));
t47 = cos(qJ(2));
t7 = t12 * t36 - t32 * t55 + mrSges(8,1);
t51 = pkin(14) * m(7) + (t47 * t24 + t25 * t41 - t72 * t31 + t7 * t35) * t49 - t67 * t47 + (t24 * t41 - t47 * t25 + t7 * t31 + t72 * t35) * t43 + t68 * t41 - mrSges(2,1) + (-m(9) - m(10) - m(3) - t71) * pkin(15);
t26 = g(1) * t42 - g(2) * t48;
t8 = -g(3) * t52 + t59 * t20;
t6 = g(3) * t20 + t52 * t59;
t5 = -t58 * t44 + t57 * t50 + t67;
t4 = t57 * t44 + t58 * t50 + t68;
t1 = [(-t42 * t54 + t51 * t48) * g(2) + (-t42 * t51 - t54 * t48) * g(1), (-g(3) * t5 + t59 * t4) * t47 + (t4 * g(3) + t59 * t5) * t41, (t6 * t40 + t8 * t46) * t47 + (-t8 * t40 + t6 * t46) * t41, (-mrSges(6,1) * t26 - mrSges(6,2) * t66) * t45 - t39 * (mrSges(6,1) * t66 - mrSges(6,2) * t26) + t73 * t33 * (t11 * t36 + t56 * t32)];
taug = t1(:);
