% Calculate potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:12:46
% EndTime: 2020-05-03 01:12:46
% DurationCPUTime: 0.25s
% Computational Cost: add. (200->58), mult. (146->46), div. (0->0), fcn. (102->14), ass. (0->29)
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t47 = -pkin(3) * m(7) - t28 * mrSges(7,1) + t24 * mrSges(7,2) - mrSges(6,1);
t46 = mrSges(6,2) + mrSges(7,3);
t45 = -m(7) - m(6);
t41 = m(3) + m(4);
t19 = m(1) + m(2) + t41;
t44 = -m(5) - t19;
t43 = t24 * mrSges(7,1) + t28 * mrSges(7,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t23 = qJ(2) + qJ(3);
t18 = qJ(4) + t23;
t17 = qJ(5) + t18;
t11 = sin(t17);
t12 = cos(t17);
t15 = sin(t18);
t16 = cos(t18);
t26 = sin(qJ(2));
t30 = cos(qJ(2));
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t4 = pkin(4) * m(4) + t29 * mrSges(4,1) - mrSges(4,2) * t25 + mrSges(3,1);
t5 = pkin(2) * cos(t23) + t30 * pkin(4) + pkin(1);
t6 = t25 * mrSges(4,1) + t29 * mrSges(4,2) + mrSges(3,2);
t42 = -t41 * pkin(1) - m(5) * t5 - t16 * mrSges(5,1) + t15 * mrSges(5,2) + t46 * t11 + t47 * t12 + t26 * t6 - t30 * t4 - mrSges(2,1);
t37 = pkin(2) * sin(t23) + t26 * pkin(4) + r_base(3);
t31 = cos(qJ(1));
t27 = sin(qJ(1));
t3 = pkin(5) * t16 + t5;
t1 = (-mrSges(1,2) + t44 * r_base(2) + t45 * (t27 * t3 + r_base(2)) - t43 * t31 + t42 * t27) * g(2) + (-mrSges(1,1) + t44 * r_base(1) + t45 * (t31 * t3 + r_base(1)) + t42 * t31 + t43 * t27) * g(1) + (-m(5) * t37 - t15 * mrSges(5,1) - t16 * mrSges(5,2) - t19 * r_base(3) - t4 * t26 - t6 * t30 - mrSges(1,3) - mrSges(2,3) + t45 * (pkin(5) * t15 + t37) - t46 * t12 + t47 * t11) * g(3);
U = t1;
