% Calculate Gravitation load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:07:26
% EndTime: 2020-05-03 00:07:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (442->56), mult. (505->75), div. (0->0), fcn. (288->10), ass. (0->33)
t44 = m(5) + m(6);
t45 = m(4) + t44;
t43 = cos(qJ(5));
t23 = sin(qJ(5));
t15 = -pkin(4) * m(6) - t43 * mrSges(6,1) + t23 * mrSges(6,2) - mrSges(5,1);
t20 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t24 = sin(qJ(4));
t28 = cos(qJ(4));
t11 = t44 * pkin(3) - t15 * t28 + t20 * t24 + mrSges(4,1);
t27 = sin(qJ(1));
t31 = cos(qJ(1));
t19 = g(1) * t31 + t27 * g(2);
t38 = t15 * t24 + t20 * t28 - mrSges(4,2);
t1 = g(3) * t38 + t19 * t11;
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t3 = g(3) * t11 - t19 * t38;
t42 = t1 * t25 + t3 * t29;
t17 = t45 * pkin(2) + mrSges(3,1);
t41 = t1 * t29 - t3 * t25;
t40 = t24 * g(3) - t19 * t28;
t26 = sin(qJ(2));
t30 = cos(qJ(2));
t39 = t26 * (t11 * t25 - t29 * t38 + mrSges(3,2)) - t30 * (t11 * t29 + t38 * t25 + t17) - mrSges(2,1) - (m(3) + t45) * pkin(1);
t37 = mrSges(6,1) * t23 + mrSges(6,2) * t43 + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t14 = g(3) * t28 + t24 * t19;
t36 = t30 * (-t14 * t25 - t40 * t29) + (-t14 * t29 + t40 * t25) * t26;
t18 = g(1) * t27 - g(2) * t31;
t10 = g(3) * t15 + t19 * t20;
t9 = g(3) * t20 - t19 * t15;
t5 = t24 * t10 + t9 * t28;
t4 = -t10 * t28 + t24 * t9;
t2 = [(t37 * t27 + t39 * t31) * g(2) + (-t39 * t27 + t37 * t31) * g(1), -(mrSges(3,2) * g(3) - t19 * t17 - t41) * t26 + (t19 * mrSges(3,2) + g(3) * t17 + t42) * t30, t41 * t26 + t42 * t30, (t5 * t25 + t4 * t29) * t30 + t26 * (-t25 * t4 + t5 * t29), (mrSges(6,1) * t18 + t36 * mrSges(6,2)) * t43 + (t36 * mrSges(6,1) - mrSges(6,2) * t18) * t23];
taug = t2(:);
