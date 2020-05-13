% Calculate Gravitation load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2DE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:27
% EndTime: 2020-05-03 01:06:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->32), mult. (94->39), div. (0->0), fcn. (32->8), ass. (0->20)
t11 = sin(qJ(4));
t14 = sin(qJ(1));
t15 = cos(qJ(4));
t18 = cos(qJ(1));
t7 = mrSges(7,1) * g(1) + mrSges(7,2) * g(2);
t8 = mrSges(7,1) * g(2) - mrSges(7,2) * g(1);
t29 = (t11 * t8 + t15 * t7) * t14 + (t7 * t11 - t8 * t15) * t18;
t28 = m(6) + m(7);
t27 = m(5) + t28;
t25 = m(4) + t27;
t24 = t28 * pkin(5) + mrSges(5,1);
t23 = g(1) * t18 + g(2) * t14;
t21 = t25 * pkin(4) + mrSges(3,1);
t12 = sin(qJ(3));
t13 = sin(qJ(2));
t16 = cos(qJ(3));
t17 = cos(qJ(2));
t20 = mrSges(3,2) * t13 + mrSges(5,2) * t12 - t16 * t24 - t17 * t21 - (m(3) + t25) * pkin(1) - pkin(3) * m(7) - mrSges(2,1) - mrSges(4,1) - mrSges(6,1) - t27 * pkin(2);
t9 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t1 = [(t9 * t14 + t20 * t18) * g(2) + (-t20 * t14 + t9 * t18) * g(1) + t29, (mrSges(3,2) * g(3) + t21 * t23) * t13 - (-t23 * mrSges(3,2) + g(3) * t21) * t17, (g(3) * mrSges(5,2) + t23 * t24) * t12 - t16 * (-t23 * mrSges(5,2) + g(3) * t24), t29];
taug = t1(:);
