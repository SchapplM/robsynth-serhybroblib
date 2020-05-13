% Calculate Gravitation load on the joints for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
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
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:43:56
% EndTime: 2020-05-03 01:43:59
% DurationCPUTime: 0.44s
% Computational Cost: add. (478->77), mult. (402->81), div. (0->0), fcn. (295->14), ass. (0->46)
t37 = sin(qJ(6));
t74 = cos(qJ(6));
t12 = t74 * mrSges(7,1) - t37 * mrSges(7,2);
t93 = mrSges(6,1) + t12;
t36 = mrSges(6,2) + mrSges(7,3);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t13 = t41 * mrSges(4,1) - mrSges(4,2) * t38;
t5 = pkin(4) * m(4) + mrSges(3,1) + t13;
t34 = qJ(2) + qJ(3);
t28 = qJ(4) + t34;
t25 = qJ(5) + t28;
t19 = sin(t25);
t20 = cos(t25);
t22 = sin(t28);
t23 = cos(t28);
t86 = -t23 * mrSges(5,1) + t22 * mrSges(5,2) + t36 * t19 - t93 * t20;
t10 = -t38 * mrSges(4,1) - t41 * mrSges(4,2);
t9 = mrSges(3,2) - t10;
t92 = t39 * t9 - t42 * t5 + t86;
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t61 = g(1) * t43 + t40 * g(2);
t90 = mrSges(5,1) * t22 + mrSges(5,2) * t23 + t93 * t19 + t36 * t20;
t4 = -pkin(3) * m(7) - t93;
t88 = (g(3) * t4 + t61 * t36) * t20 + (g(3) * t36 - t61 * t4) * t19;
t11 = t37 * mrSges(7,1) + t74 * mrSges(7,2);
t84 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) + t11;
t18 = pkin(5) * t23;
t21 = pkin(2) * cos(t34);
t68 = t42 * pkin(4) + t21;
t6 = pkin(1) + t68;
t1 = t18 + t6;
t14 = t20 * pkin(3);
t81 = m(5) * t6 + m(6) * t1 + m(7) * (t1 + t14) + mrSges(2,1) + (m(3) + m(4)) * pkin(1) - t92;
t78 = t39 * pkin(4);
t69 = t21 + t18;
t67 = pkin(2) * sin(t34);
t66 = pkin(3) * t19;
t65 = t18 + t68;
t3 = -pkin(5) * t22 - t67;
t24 = (m(6) + m(7)) * pkin(5) + mrSges(5,1);
t2 = t3 - t78;
t7 = [(t84 * t40 - t81 * t43) * g(2) + (t81 * t40 + t84 * t43) * g(1), (-m(5) * t68 - m(6) * t65 - m(7) * (t14 + t65) + t92) * g(3) + t61 * (t39 * t5 + t42 * t9 - m(5) * (-t67 - t78) - m(6) * t2 - m(7) * (t2 - t66) + t90), (-t13 * t42 - t10 * t39 - m(5) * t21 - m(6) * t69 - m(7) * (t14 + t69) + t86) * g(3) + t61 * (-t10 * t42 + t13 * t39 + m(5) * t67 - m(6) * t3 - m(7) * (t3 - t66) + t90), (mrSges(5,2) * g(3) + t61 * t24) * t22 - (-t61 * mrSges(5,2) + t24 * g(3)) * t23 + t88, t88, t12 * (g(1) * t40 - g(2) * t43) + (g(3) * t19 + t61 * t20) * t11];
taug = t7(:);
