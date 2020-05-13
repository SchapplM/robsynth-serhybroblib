% Calculate potential energy for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:23
% EndTime: 2020-04-12 19:33:23
% DurationCPUTime: 0.31s
% Computational Cost: add. (697->52), mult. (965->57), div. (60->5), fcn. (299->11), ass. (0->31)
t44 = pkin(4) ^ 2;
t43 = -m(2) - m(3) - m(4) - m(5);
t42 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t41 = -m(1) + t43;
t19 = cos(qJ(2));
t14 = pkin(1) * t19 - pkin(2);
t17 = sin(qJ(2));
t37 = pkin(2) * t19;
t33 = (-0.2e1 * t37 + pkin(1)) * pkin(1);
t12 = pkin(2) ^ 2 + t33;
t10 = 0.1e1 / t12;
t34 = t10 / pkin(3);
t38 = -pkin(3) + pkin(4);
t39 = -pkin(3) - pkin(4);
t7 = sqrt(-((pkin(2) - t39) * (pkin(2) + t39) + t33) * ((pkin(2) - t38) * (pkin(2) + t38) + t33));
t35 = t17 * t7;
t32 = pkin(3) ^ 2 - t44;
t8 = t12 + t32;
t3 = qJ(2) + atan2((pkin(1) * t17 * t8 - t14 * t7) * t34, (-pkin(1) * t35 - t14 * t8) * t34);
t1 = sin(t3);
t2 = cos(t3);
t13 = pkin(1) - t37;
t9 = t12 - t32;
t5 = -pkin(2) * t35 + t13 * t9;
t36 = t17 * pkin(2);
t6 = t13 * t7 + t36 * t9;
t31 = t10 / pkin(4) * ((t5 ^ 2 + t6 ^ 2) / t44 / t12 ^ 2) ^ (-0.1e1 / 0.2e1);
t40 = -m(5) * pkin(1) - mrSges(2,1) - (-mrSges(5,1) * t5 - mrSges(5,2) * t6) * t31 - m(4) * t37 + t2 * mrSges(4,1) - t1 * mrSges(4,2) - mrSges(3,1) * t19 + mrSges(3,2) * t17;
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t4 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t17 * mrSges(3,1) - t19 * mrSges(3,2) - m(4) * t36 + t1 * mrSges(4,1) + t2 * mrSges(4,2) - (mrSges(5,1) * t6 - mrSges(5,2) * t5) * t31 + t43 * (pkin(5) + r_base(3))) * g(3) + (t40 * t18 - t42 * t20 + t41 * r_base(2) - mrSges(1,2)) * g(2) + (t42 * t18 + t40 * t20 + t41 * r_base(1) - mrSges(1,1)) * g(1);
U = t4;
