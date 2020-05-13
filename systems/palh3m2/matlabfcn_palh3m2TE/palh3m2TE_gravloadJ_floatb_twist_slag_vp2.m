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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2TE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:46:25
% EndTime: 2020-05-07 01:46:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (294->84), mult. (362->109), div. (0->0), fcn. (272->22), ass. (0->46)
t53 = m(5) + m(6);
t23 = t53 * pkin(4) + mrSges(4,1);
t49 = cos(qJ(2));
t25 = pkin(1) * t49 + pkin(12);
t60 = m(4) + t53;
t22 = mrSges(9,1) + t23;
t44 = sin(qJ(1));
t50 = cos(qJ(1));
t21 = g(1) * t50 + t44 * g(2);
t41 = sin(qJ(4));
t47 = cos(qJ(4));
t59 = t47 * mrSges(6,1) - t41 * mrSges(6,2);
t58 = t41 * mrSges(6,1) + t47 * mrSges(6,2);
t45 = sin(pkin(15));
t51 = cos(pkin(15));
t57 = t51 * mrSges(7,1) - mrSges(7,2) * t45;
t56 = mrSges(7,1) * t45 + mrSges(7,2) * t51;
t55 = -mrSges(8,3) - mrSges(9,3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - t58;
t26 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t31 = pkin(8) * m(6) + mrSges(5,1);
t35 = sin(pkin(16));
t36 = cos(pkin(16));
t11 = t26 * t36 + t35 * t31;
t12 = -t35 * t26 + t31 * t36;
t46 = sin(pkin(14));
t52 = cos(pkin(14));
t19 = mrSges(7,1) * t46 - mrSges(7,2) * t52;
t20 = mrSges(7,1) * t52 + mrSges(7,2) * t46;
t33 = pkin(17) + pkin(18);
t27 = sin(t33);
t28 = cos(t33);
t34 = qJ(3) + qJ(2);
t37 = sin(pkin(18));
t39 = cos(pkin(18));
t42 = sin(qJ(3));
t43 = sin(qJ(2));
t48 = cos(qJ(3));
t54 = m(8) * t25 + m(9) * ((-(-t37 * t45 + t39 * t51) * cos(pkin(17)) + (t37 * t51 + t45 * t39) * sin(pkin(17))) * pkin(3) + t25) - mrSges(9,1) * cos(t34) + mrSges(9,2) * sin(t34) + (t11 * t51 + t12 * t45 + t59 * (t35 * t51 + t36 * t45)) * t27 - (-t11 * t45 + t12 * t51 + t59 * (-t35 * t45 + t36 * t51)) * t28 + (t60 * pkin(1) + t42 * mrSges(4,2) + t19 * t45 + t20 * t51 - t23 * t48 + mrSges(3,1)) * t49 + t43 * (t48 * mrSges(4,2) + t19 * t51 - t20 * t45 + t23 * t42 - mrSges(3,2)) - pkin(6) * m(7) + mrSges(2,1) - (-m(3) - t60) * pkin(12) + (-mrSges(8,1) * t51 - mrSges(8,2) * t45) * t39 + t37 * (mrSges(8,1) * t45 - mrSges(8,2) * t51);
t40 = mrSges(4,2) + mrSges(9,2);
t14 = g(3) * t51 + t45 * t21;
t13 = -t45 * g(3) + t21 * t51;
t10 = g(3) * t22 - t21 * t40;
t8 = g(3) * t40 + t21 * t22;
t5 = t22 * t42 + t40 * t48 + t57 * t46 - t56 * t52 - mrSges(3,2);
t3 = -t22 * t48 + t40 * t42 + t57 * t52 + (m(8) + m(9) + t60) * pkin(1) + mrSges(3,1) + t56 * t46;
t1 = [(t55 * t44 - t54 * t50) * g(2) + (t54 * t44 + t55 * t50) * g(1), (-t5 * g(3) + t21 * t3) * t43 - (t3 * g(3) + t21 * t5) * t49, (t10 * t48 - t42 * t8) * t49 - (t42 * t10 + t8 * t48) * t43, -t59 * (g(1) * t44 - g(2) * t50) + (-(t13 * t36 - t14 * t35) * t28 + (t13 * t35 + t14 * t36) * t27) * t58];
taug = t1(:);
