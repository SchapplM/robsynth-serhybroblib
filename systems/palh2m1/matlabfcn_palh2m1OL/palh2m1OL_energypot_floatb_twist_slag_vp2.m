% Calculate potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:58:40
% EndTime: 2020-05-02 23:58:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (171->44), mult. (165->44), div. (0->0), fcn. (72->10), ass. (0->23)
t24 = m(5) + m(6);
t23 = m(4) + t24;
t22 = m(3) + t23;
t21 = m(2) + t22;
t14 = sin(qJ(1));
t18 = cos(qJ(1));
t20 = g(1) * t18 + g(2) * t14;
t11 = sin(qJ(5));
t15 = cos(qJ(5));
t19 = mrSges(6,1) * t11 + mrSges(6,2) * t15 + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t17 = cos(qJ(3));
t16 = cos(qJ(4));
t13 = sin(qJ(3));
t12 = sin(qJ(4));
t9 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t7 = t23 * pkin(2) + mrSges(3,1);
t6 = t22 * pkin(1) + mrSges(2,1);
t5 = pkin(4) * m(6) + mrSges(6,1) * t15 - mrSges(6,2) * t11 + mrSges(5,1);
t4 = t12 * t5 - t16 * t9 + mrSges(4,2);
t3 = t24 * pkin(3) + t12 * t9 + t16 * t5 + mrSges(4,1);
t2 = g(3) * t3 + t20 * t4;
t1 = -t4 * g(3) + t20 * t3;
t8 = (mrSges(3,2) * g(3) - t1 * t17 + t2 * t13 - t20 * t7) * cos(qJ(2)) + (t20 * mrSges(3,2) + g(3) * t7 + t1 * t13 + t2 * t17) * sin(qJ(2)) + (-g(1) * t6 - t19 * g(2)) * t18 + (t19 * g(1) - g(2) * t6) * t14 - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) - g(3) * (mrSges(1,3) + mrSges(2,3)) + (-g(1) * r_base(1) - g(2) * r_base(2) - g(3) * r_base(3)) * (m(1) + t21) - g(3) * t21 * pkin(5);
U = t8;
