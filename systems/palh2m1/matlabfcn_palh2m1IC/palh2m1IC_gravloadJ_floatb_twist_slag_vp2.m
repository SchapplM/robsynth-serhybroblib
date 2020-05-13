% Calculate Gravitation load on the joints for
% palh2m1IC
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:54
% EndTime: 2020-05-03 01:03:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (527->55), mult. (607->69), div. (0->0), fcn. (342->10), ass. (0->31)
t46 = m(5) + m(6);
t49 = m(4) + t46;
t25 = sin(qJ(5));
t45 = cos(qJ(5));
t17 = -pkin(4) * m(6) - t45 * mrSges(6,1) + t25 * mrSges(6,2) - mrSges(5,1);
t29 = sin(qJ(1));
t33 = cos(qJ(1));
t21 = g(1) * t33 + t29 * g(2);
t22 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t11 = g(3) * t22 - t21 * t17;
t12 = g(3) * t17 + t21 * t22;
t26 = sin(qJ(4));
t30 = cos(qJ(4));
t13 = t46 * pkin(3) - t17 * t30 + t22 * t26 + mrSges(4,1);
t40 = t17 * t26 + t22 * t30 - mrSges(4,2);
t48 = g(3) * t40 - t11 * t30 - t26 * t12 + t21 * t13;
t47 = g(3) * t13 - t26 * t11 + t12 * t30 - t21 * t40;
t19 = t49 * pkin(2) + mrSges(3,1);
t27 = sin(qJ(3));
t31 = cos(qJ(3));
t44 = t48 * t27 + t47 * t31;
t43 = t26 * g(3) - t21 * t30;
t42 = -t47 * t27 + t48 * t31;
t28 = sin(qJ(2));
t32 = cos(qJ(2));
t41 = t28 * (t13 * t27 - t31 * t40 + mrSges(3,2)) - t32 * (t13 * t31 + t40 * t27 + t19) - mrSges(2,1) - (m(3) + t49) * pkin(1);
t39 = mrSges(6,1) * t25 + mrSges(6,2) * t45 + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t16 = g(3) * t30 + t26 * t21;
t38 = (-t16 * t27 - t43 * t31) * t32 + (-t16 * t31 + t43 * t27) * t28;
t20 = g(1) * t29 - g(2) * t33;
t1 = [(t39 * t29 + t41 * t33) * g(2) + (-t41 * t29 + t39 * t33) * g(1); (t21 * mrSges(3,2) + g(3) * t19 + t44) * t32 + (-mrSges(3,2) * g(3) + t21 * t19 + t42) * t28; t42 * t28 + t44 * t32; (mrSges(6,1) * t20 + t38 * mrSges(6,2)) * t45 + (t38 * mrSges(6,1) - mrSges(6,2) * t20) * t25;];
taug = t1(:);
