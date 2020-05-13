% Calculate Gravitation load on the joints for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:35
% EndTime: 2020-05-07 11:33:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (75->30), mult. (117->41), div. (4->3), fcn. (90->11), ass. (0->23)
t10 = qJ(2) + qJ(3);
t8 = sin(t10);
t9 = cos(t10);
t39 = mrSges(4,1) * t8 + mrSges(4,2) * t9;
t12 = sin(qJ(2));
t15 = cos(qJ(2));
t22 = -t9 * mrSges(4,1) + t8 * mrSges(4,2);
t34 = m(4) * pkin(2);
t38 = t12 * mrSges(3,2) - t22 - (mrSges(3,1) + t34) * t15;
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t36 = g(1) * t16 + g(2) * t13;
t11 = sin(qJ(4));
t14 = cos(qJ(4));
t19 = t14 * mrSges(5,1) - t11 * mrSges(5,2);
t35 = m(5) * pkin(1) + mrSges(2,1) + t19 - t38;
t33 = t39 * t13;
t32 = t39 * t16;
t25 = -qJ(4) + qJ(2);
t24 = t12 * t34;
t21 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t7 = sin(qJ(3) + t25);
t1 = [(t21 * g(1) - t35 * g(2)) * t16 + (t35 * g(1) + t21 * g(2)) * t13; -g(1) * (-t16 * t24 + t32) - g(2) * (-t13 * t24 + t33) + t38 * g(3) + ((-pkin(3) * t7 + pkin(2) * sin(t25)) / pkin(3) * (-g(1) * t32 - g(2) * t33 - g(3) * t22) + sin(qJ(3)) * pkin(2) / pkin(4) * (-g(3) * t19 + t36 * (mrSges(5,1) * t11 + mrSges(5,2) * t14))) / t7 + t36 * (mrSges(3,1) * t12 + mrSges(3,2) * t15);];
taug = t1(:);
