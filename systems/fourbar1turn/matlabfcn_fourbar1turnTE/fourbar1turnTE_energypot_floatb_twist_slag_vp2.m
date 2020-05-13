% Calculate potential energy for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnTE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:28
% EndTime: 2020-04-12 19:18:28
% DurationCPUTime: 0.34s
% Computational Cost: add. (439->52), mult. (615->62), div. (36->3), fcn. (195->6), ass. (0->33)
t46 = -m(2) - m(3) - m(4) - m(5);
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t36 = pkin(2) * t15;
t31 = (-0.2e1 * t36 + pkin(1)) * pkin(1);
t37 = -pkin(3) + pkin(4);
t38 = -pkin(3) - pkin(4);
t5 = sqrt(-((pkin(2) - t38) * (pkin(2) + t38) + t31) * ((pkin(2) - t37) * (pkin(2) + t37) + t31));
t34 = t13 * t5;
t29 = pkin(2) ^ 2 + t31;
t30 = pkin(3) ^ 2 - pkin(4) ^ 2;
t7 = t29 - t30;
t9 = pkin(1) - t36;
t2 = -pkin(2) * t34 + t9 * t7;
t10 = pkin(1) * t15 - pkin(2);
t6 = t29 + t30;
t1 = -pkin(1) * t34 - t10 * t6;
t8 = 0.1e1 / t29;
t32 = 0.1e1 / pkin(3) * t8;
t39 = t13 / 0.2e1;
t4 = pkin(1) * t13 * t6 - t10 * t5;
t26 = (t4 * t39 - t15 * t1 / 0.2e1) * t32;
t35 = t13 * pkin(2);
t3 = t35 * t7 + t9 * t5;
t33 = 0.1e1 / pkin(4) * t8;
t40 = -mrSges(5,2) / 0.2e1;
t42 = (t15 * t4 / 0.2e1 + t1 * t39) * t32;
t45 = -mrSges(2,1) - m(5) * pkin(1) - (-t2 * mrSges(5,1) / 0.2e1 + t3 * t40) * t33 - mrSges(3,1) * t15 + mrSges(3,2) * t13 - m(4) * t36 - mrSges(4,1) * t26 - mrSges(4,2) * t42;
t44 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t43 = -m(1) + t46;
t16 = cos(qJ(1));
t14 = sin(qJ(1));
t11 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t13 * mrSges(3,1) - t15 * mrSges(3,2) - m(4) * t35 + mrSges(4,1) * t42 - mrSges(4,2) * t26 - (t3 * mrSges(5,1) / 0.2e1 + t2 * t40) * t33 + t46 * (pkin(5) + r_base(3))) * g(3) + (t45 * t14 - t44 * t16 + t43 * r_base(2) - mrSges(1,2)) * g(2) + (t44 * t14 + t45 * t16 + t43 * r_base(1) - mrSges(1,1)) * g(1);
U = t11;
