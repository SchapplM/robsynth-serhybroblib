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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:19
% EndTime: 2020-05-02 21:17:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (296->87), mult. (248->77), div. (0->0), fcn. (134->24), ass. (0->38)
t49 = g(1) * r_base(1) + r_base(2) * g(2);
t44 = m(5) + m(6);
t48 = m(8) + m(4) + t44;
t12 = -qJ(7) + pkin(19) - qJ(2);
t18 = sin(qJ(6));
t22 = sin(qJ(2));
t26 = cos(qJ(6));
t39 = -m(9) - m(3) - t48;
t10 = -qJ(10) + t12;
t8 = sin(t10);
t9 = cos(t10);
t47 = pkin(14) * m(7) - t26 * mrSges(7,1) + mrSges(7,2) * t18 - mrSges(2,1) - m(11) * (-t22 * pkin(1) + pkin(4) * sin(t12)) + t8 * mrSges(11,1) - t9 * mrSges(11,2) + (t39 - m(11) - m(10)) * pkin(15);
t19 = sin(qJ(5));
t27 = cos(qJ(5));
t46 = mrSges(6,1) * t19 + mrSges(6,2) * t27 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t15 = qJ(2) + qJ(8);
t38 = pkin(13) + r_base(3);
t23 = sin(qJ(1));
t31 = cos(qJ(1));
t37 = g(1) * t31 + t23 * g(2);
t11 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t20 = sin(qJ(4));
t28 = cos(qJ(4));
t5 = m(6) * pkin(9) + mrSges(6,1) * t27 - mrSges(6,2) * t19 + mrSges(5,1);
t35 = t11 * t28 - t20 * t5 - mrSges(4,2);
t34 = m(7) + m(2) - t39;
t30 = cos(qJ(2));
t29 = cos(qJ(3));
t25 = cos(qJ(7));
t24 = cos(qJ(8));
t21 = sin(qJ(3));
t17 = sin(qJ(7));
t16 = sin(qJ(8));
t13 = qJ(9) + t15;
t3 = t44 * pkin(5) + t11 * t20 + t28 * t5 + mrSges(4,1);
t2 = -t17 * mrSges(8,1) - t16 * mrSges(9,1) - t25 * mrSges(8,2) - t24 * mrSges(9,2) + t35 * t21 + t3 * t29 - mrSges(3,2);
t1 = t48 * pkin(1) + t25 * mrSges(8,1) + t24 * mrSges(9,1) - mrSges(8,2) * t17 - t16 * mrSges(9,2) + t3 * t21 - t35 * t29 + mrSges(3,1);
t4 = (-t2 * g(3) + t37 * t1) * t22 + (-t1 * g(3) - t37 * t2) * t30 - mrSges(7,2) * t26 * g(3) - mrSges(7,1) * t18 * g(3) - (-m(7) * pkin(16) + t34 * pkin(13) + mrSges(1,3) + mrSges(2,3)) * g(3) + (g(3) * mrSges(10,1) - t37 * mrSges(10,2)) * cos(t13) + (-t37 * mrSges(10,1) - g(3) * mrSges(10,2)) * sin(t13) - g(3) * (m(11) * (pkin(4) * cos(t12) + t30 * pkin(1) + t38) - t9 * mrSges(11,1) - t8 * mrSges(11,2)) + (-m(11) * r_base(2) + t47 * t23 + t46 * t31 - mrSges(1,2)) * g(2) + (-m(11) * r_base(1) - t46 * t23 + t47 * t31 - mrSges(1,1)) * g(1) + (-g(3) * r_base(3) - t49) * (m(1) + t34) + (-g(3) * t38 - t49 + (t37 * sin(t15) - cos(t15) * g(3)) * pkin(2)) * m(10);
U = t4;
