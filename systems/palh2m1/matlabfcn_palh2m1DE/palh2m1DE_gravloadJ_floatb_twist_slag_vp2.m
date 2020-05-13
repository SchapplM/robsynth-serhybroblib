% Calculate Gravitation load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1DE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (123->33), mult. (139->44), div. (0->0), fcn. (56->8), ass. (0->22)
t35 = m(5) + m(6);
t37 = m(4) + t35;
t17 = sin(qJ(3));
t21 = cos(qJ(3));
t32 = t35 * pkin(3) + mrSges(4,1);
t6 = mrSges(4,2) * t21 + t17 * t32 + mrSges(3,2);
t11 = mrSges(6,1) * g(1) + mrSges(6,2) * g(2);
t12 = mrSges(6,1) * g(2) - mrSges(6,2) * g(1);
t16 = sin(qJ(4));
t19 = sin(qJ(1));
t20 = cos(qJ(4));
t23 = cos(qJ(1));
t36 = (t11 * t20 + t12 * t16) * t19 + (t11 * t16 - t12 * t20) * t23;
t3 = -mrSges(4,2) * t17 + t37 * pkin(2) + t21 * t32 + mrSges(3,1);
t31 = g(1) * t23 + g(2) * t19;
t18 = sin(qJ(2));
t22 = cos(qJ(2));
t29 = t22 * t3 - t6 * t18 + pkin(4) * m(6) + mrSges(2,1) + mrSges(5,1) + (m(3) + t37) * pkin(1);
t14 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t4 = t31 * mrSges(4,2) + g(3) * t32;
t2 = -g(3) * mrSges(4,2) + t31 * t32;
t1 = [(t14 * t19 - t29 * t23) * g(2) + (t14 * t23 + t29 * t19) * g(1) + t36, (-g(3) * t6 + t31 * t3) * t18 + t22 * (g(3) * t3 + t31 * t6), (t17 * t2 + t21 * t4) * t22 + t18 * (-t17 * t4 + t2 * t21), t36];
taug = t1(:);
