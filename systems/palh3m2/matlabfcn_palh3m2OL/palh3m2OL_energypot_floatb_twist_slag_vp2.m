% Calculate potential energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:22
% EndTime: 2020-05-07 04:33:23
% DurationCPUTime: 0.54s
% Computational Cost: add. (232->79), mult. (203->62), div. (0->0), fcn. (153->18), ass. (0->34)
t53 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t26 = sin(qJ(5));
t30 = cos(qJ(5));
t52 = pkin(8) * m(6) + mrSges(6,1) * t30 - mrSges(6,2) * t26 + mrSges(5,1);
t51 = -m(4) - m(8);
t50 = -m(5) - m(6);
t49 = -m(7) - m(2) - m(3);
t47 = -m(1) - m(9) + t49;
t31 = cos(qJ(2));
t10 = t31 * pkin(1) + pkin(12);
t24 = qJ(2) + qJ(3);
t19 = qJ(4) + t24;
t11 = sin(t19);
t12 = cos(t19);
t23 = qJ(2) + qJ(7);
t14 = sin(t23);
t15 = sin(t24);
t16 = cos(t23);
t17 = cos(t24);
t18 = pkin(15) - t23;
t25 = sin(qJ(6));
t27 = sin(qJ(2));
t29 = cos(qJ(6));
t13 = -qJ(8) + t18;
t8 = sin(t13);
t9 = cos(t13);
t46 = m(7) * pkin(6) - m(3) * pkin(12) - t31 * mrSges(3,1) - t29 * mrSges(7,1) + t27 * mrSges(3,2) + t25 * mrSges(7,2) - mrSges(2,1) - m(9) * (pkin(3) * cos(t18) + t10) + t9 * mrSges(9,1) + t8 * mrSges(9,2) - t16 * mrSges(8,1) + t14 * mrSges(8,2) + t17 * mrSges(4,1) - t15 * mrSges(4,2) + t52 * t12 - t53 * t11;
t45 = -t26 * mrSges(6,1) - t30 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t22 = pkin(11) + r_base(3);
t7 = t27 * pkin(1) + t22;
t32 = cos(qJ(1));
t28 = sin(qJ(1));
t4 = -pkin(4) * t17 + t10;
t1 = (m(9) * pkin(3) * sin(t18) - m(1) * r_base(3) - m(7) * pkin(13) - t27 * mrSges(3,1) + t15 * mrSges(4,1) - t25 * mrSges(7,1) - t14 * mrSges(8,1) - t8 * mrSges(9,1) - t31 * mrSges(3,2) + t17 * mrSges(4,2) - t29 * mrSges(7,2) - t16 * mrSges(8,2) + t9 * mrSges(9,2) - mrSges(1,3) - mrSges(2,3) + (-m(9) + t51) * t7 + t50 * (-pkin(4) * t15 + t7) + t49 * t22 + t53 * t12 + t52 * t11) * g(3) + (-mrSges(1,2) + t51 * (t10 * t28 + r_base(2)) + t50 * (t28 * t4 + r_base(2)) + t47 * r_base(2) - t45 * t32 + t46 * t28) * g(2) + (-mrSges(1,1) + t51 * (t10 * t32 + r_base(1)) + t50 * (t32 * t4 + r_base(1)) + t47 * r_base(1) + t46 * t32 + t45 * t28) * g(1);
U = t1;
