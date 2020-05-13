% Calculate Gravitation load on the joints for
% palh2m2IC
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:35
% EndTime: 2020-05-03 06:51:36
% DurationCPUTime: 0.44s
% Computational Cost: add. (480->49), mult. (402->51), div. (0->0), fcn. (295->13), ass. (0->32)
t86 = m(6) + m(7);
t88 = m(5) + t86;
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t60 = g(1) * t43 + t40 * g(2);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t13 = t41 * mrSges(4,1) - mrSges(4,2) * t38;
t5 = pkin(4) * m(4) + mrSges(3,1) + t13;
t10 = -t38 * mrSges(4,1) - t41 * mrSges(4,2);
t9 = mrSges(3,2) - t10;
t83 = t39 * t9 - t42 * t5;
t37 = sin(qJ(6));
t74 = cos(qJ(6));
t11 = t37 * mrSges(7,1) + t74 * mrSges(7,2);
t82 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) + t11;
t12 = t74 * mrSges(7,1) - t37 * mrSges(7,2);
t34 = qJ(2) + qJ(3);
t28 = qJ(4) + t34;
t25 = qJ(5) + t28;
t19 = sin(t25);
t20 = cos(t25);
t22 = sin(t28);
t23 = cos(t28);
t21 = pkin(2) * cos(t34);
t68 = t42 * pkin(4) + t21;
t6 = pkin(1) + t68;
t81 = m(5) * t6 + t23 * mrSges(5,1) - t22 * mrSges(5,2) + mrSges(2,1) + (m(3) + m(4)) * pkin(1) - t83 + t86 * (pkin(5) * t23 + t6) + (pkin(3) * m(7) + mrSges(6,1) + t12) * t20 + (-mrSges(6,2) - mrSges(7,3)) * t19;
t24 = t86 * pkin(5) + mrSges(5,1);
t1 = [(t82 * t40 - t81 * t43) * g(2) + (t81 * t40 + t82 * t43) * g(1); (t10 * t39 + t13 * t42 + t83 + t88 * (t21 - t68)) * g(3) + ((t9 + t10) * t42 + (t88 * pkin(4) - t13 + t5) * t39) * t60; (mrSges(5,2) * g(3) + t60 * t24) * t22 - (-t60 * mrSges(5,2) + t24 * g(3)) * t23; t12 * (g(1) * t40 - g(2) * t43) + (g(3) * t19 + t60 * t20) * t11;];
taug = t1(:);
