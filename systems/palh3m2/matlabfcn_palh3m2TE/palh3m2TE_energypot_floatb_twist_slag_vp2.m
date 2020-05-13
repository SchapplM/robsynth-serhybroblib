% Calculate potential energy for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:42:04
% EndTime: 2020-05-07 01:42:04
% DurationCPUTime: 0.38s
% Computational Cost: add. (264->88), mult. (325->97), div. (0->0), fcn. (239->22), ass. (0->48)
t55 = m(5) + m(6);
t54 = m(4) + t55;
t53 = -m(3) - t54;
t51 = m(2) + m(7) - t53;
t14 = m(1) + t51;
t58 = -t14 - m(8) - m(9);
t32 = sin(qJ(4));
t38 = cos(qJ(4));
t57 = mrSges(6,1) * t32 + mrSges(6,2) * t38 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t40 = cos(qJ(2));
t17 = pkin(1) * t40 + pkin(12);
t25 = qJ(3) + qJ(2);
t19 = sin(t25);
t20 = cos(t25);
t28 = sin(pkin(18));
t29 = sin(pkin(17));
t37 = sin(pkin(14));
t43 = cos(pkin(14));
t12 = mrSges(7,1) * t37 - mrSges(7,2) * t43;
t13 = mrSges(7,1) * t43 + mrSges(7,2) * t37;
t15 = t55 * pkin(4) + mrSges(4,1);
t33 = sin(qJ(3));
t36 = sin(pkin(15));
t39 = cos(qJ(3));
t42 = cos(pkin(15));
t3 = t54 * pkin(1) + t33 * mrSges(4,2) + t12 * t36 + t13 * t42 - t15 * t39 + mrSges(3,1);
t30 = cos(pkin(18));
t31 = cos(pkin(17));
t34 = sin(qJ(2));
t47 = -t28 * t36 + t30 * t42;
t48 = t28 * t42 + t30 * t36;
t5 = t39 * mrSges(4,2) + t12 * t42 - t13 * t36 + t15 * t33 - mrSges(3,2);
t56 = -m(8) * t17 - m(9) * ((t48 * t29 - t47 * t31) * pkin(3) + t17) + t20 * mrSges(9,1) - t19 * mrSges(9,2) + m(7) * pkin(6) + t53 * pkin(12) - t3 * t40 - t34 * t5 - mrSges(2,1) - (-mrSges(8,1) * t42 - mrSges(8,2) * t36) * t30 - t28 * (mrSges(8,1) * t36 - mrSges(8,2) * t42);
t52 = t34 * pkin(1) + pkin(11) + r_base(3);
t35 = sin(qJ(1));
t41 = cos(qJ(1));
t50 = g(1) * t41 + g(2) * t35;
t49 = mrSges(6,1) * t38 - mrSges(6,2) * t32;
t27 = cos(pkin(16));
t26 = sin(pkin(16));
t24 = pkin(17) + pkin(18);
t21 = pkin(8) * m(6) + mrSges(5,1);
t18 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t8 = -t26 * t18 + t21 * t27;
t7 = t18 * t27 + t26 * t21;
t2 = -t7 * t36 + t8 * t42 + t49 * (-t26 * t36 + t27 * t42);
t1 = t8 * t36 + t7 * t42 + t49 * (t26 * t42 + t27 * t36);
t4 = (-t1 * g(3) + t50 * t2) * cos(t24) + (-g(3) * t2 - t50 * t1) * sin(t24) + g(3) * t5 * t40 - g(3) * t3 * t34 - g(3) * t14 * r_base(3) - (m(7) * pkin(13) + t51 * pkin(11) + mrSges(1,3) + mrSges(2,3)) * g(3) - g(3) * (m(8) * t52 + (t28 * mrSges(8,1) - mrSges(8,2) * t30) * t42 + t36 * (mrSges(8,1) * t30 + mrSges(8,2) * t28)) - g(3) * (-t19 * mrSges(9,1) - t20 * mrSges(9,2) + (t52 + (t47 * t29 + t48 * t31) * pkin(3)) * m(9)) + (t56 * t35 + t57 * t41 + t58 * r_base(2) - mrSges(1,2)) * g(2) + (-t57 * t35 + t56 * t41 + t58 * r_base(1) - mrSges(1,1)) * g(1);
U = t4;
