% Calculate potential energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:48
% EndTime: 2020-05-11 05:44:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (190->89), mult. (134->70), div. (0->0), fcn. (68->22), ass. (0->38)
t39 = pkin(5) * m(6);
t28 = sin(qJ(1));
t38 = t28 * pkin(1);
t30 = cos(qJ(1));
t37 = t30 * pkin(1);
t26 = qJ(1) + qJ(2);
t36 = -m(7) - m(9) - m(3);
t23 = qJ(3) + t26;
t22 = qJ(4) + t26;
t18 = sin(t26);
t35 = -pkin(2) * t18 - t38;
t20 = cos(t26);
t34 = -pkin(2) * t20 - t37;
t33 = -pkin(3) * t18 - t38;
t32 = -pkin(3) * t20 - t37;
t31 = -m(6) - m(11) - m(10) - m(5) - m(4) - m(8) - m(2) - m(1);
t29 = cos(qJ(7));
t27 = sin(qJ(7));
t25 = qJ(1) + qJ(8);
t24 = pkin(8) + qJ(5);
t21 = qJ(6) + t26;
t19 = cos(t25);
t17 = sin(t25);
t16 = cos(t24);
t15 = sin(t24);
t14 = qJ(9) + t23;
t13 = qJ(10) + t22;
t12 = cos(t23);
t11 = cos(t22);
t10 = cos(t21);
t9 = sin(t23);
t8 = sin(t22);
t7 = sin(t21);
t4 = cos(t14);
t3 = sin(t14);
t2 = cos(t13);
t1 = sin(t13);
t5 = (-mrSges(1,3) - mrSges(2,3) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) + (t31 + t36) * r_base(3)) * g(3) + (t36 * (r_base(2) - t38) - sin(pkin(8)) * t39 + t31 * r_base(2) - m(11) * (-pkin(4) * t8 + t33) - m(5) * t33 - m(10) * (pkin(6) * t9 + t35) - m(4) * t35 + t30 * mrSges(2,2) - t15 * mrSges(6,1) - t16 * mrSges(6,2) - t17 * mrSges(9,1) + t18 * mrSges(3,1) - t19 * mrSges(9,2) + t20 * mrSges(3,2) - t27 * mrSges(8,2) + t28 * mrSges(2,1) + t29 * mrSges(8,1) - t2 * mrSges(11,2) + t3 * mrSges(10,1) + t4 * mrSges(10,2) - t7 * mrSges(7,1) + t8 * mrSges(5,1) - t9 * mrSges(4,1) - t10 * mrSges(7,2) + t11 * mrSges(5,2) - t12 * mrSges(4,2) - t1 * mrSges(11,1) - mrSges(1,2)) * g(2) + (-cos(pkin(8)) * t39 + t36 * (r_base(1) - t37) + t31 * r_base(1) - m(11) * (-pkin(4) * t11 + t32) - m(5) * t32 - m(10) * (pkin(6) * t12 + t34) - m(4) * t34 - m(8) * pkin(7) - t29 * mrSges(8,2) + t30 * mrSges(2,1) + t15 * mrSges(6,2) - t16 * mrSges(6,1) + t17 * mrSges(9,2) - t18 * mrSges(3,2) - t19 * mrSges(9,1) + t20 * mrSges(3,1) - t27 * mrSges(8,1) - t28 * mrSges(2,2) - t2 * mrSges(11,1) - t3 * mrSges(10,2) + t4 * mrSges(10,1) + t7 * mrSges(7,2) - t8 * mrSges(5,2) + t9 * mrSges(4,2) - t10 * mrSges(7,1) + t11 * mrSges(5,1) - t12 * mrSges(4,1) + t1 * mrSges(11,2) - mrSges(1,1)) * g(1);
U = t5;
