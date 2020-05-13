% Calculate potential energy for
% fourbar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:02:33
% EndTime: 2020-04-24 20:02:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (124->46), mult. (177->46), div. (8->3), fcn. (44->4), ass. (0->23)
t29 = -0.1e1 / pkin(4) / 0.2e1;
t28 = -0.1e1 / pkin(3) / 0.2e1;
t27 = -pkin(3) - pkin(4);
t26 = -pkin(3) + pkin(4);
t11 = cos(qJ(1));
t25 = pkin(2) * t11;
t10 = sin(qJ(1));
t24 = t10 * pkin(2);
t23 = (-0.2e1 * t25 + pkin(1)) * pkin(1);
t22 = pkin(3) ^ 2 - pkin(4) ^ 2;
t12 = m(1) + m(2);
t21 = -m(3) - m(4) - t12;
t20 = pkin(2) ^ 2 + t23;
t19 = -pkin(1) + t25;
t8 = 0.1e1 / t20;
t7 = t20 - t22;
t6 = t20 + t22;
t5 = -t19 * mrSges(3,1) + mrSges(3,2) * t24;
t4 = mrSges(3,1) * t24 + t19 * mrSges(3,2);
t3 = mrSges(4,1) * t24 + t19 * mrSges(4,2);
t2 = -t19 * mrSges(4,1) + mrSges(4,2) * t24;
t1 = sqrt(-((pkin(2) - t27) * (pkin(2) + t27) + t23) * ((pkin(2) - t26) * (pkin(2) + t26) + t23));
t9 = (t21 * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3)) * g(3) + (-mrSges(2,2) * t11 - mrSges(1,2) + (-m(3) * pkin(2) - mrSges(2,1)) * t10 + t21 * r_base(2) + ((t5 * t1 - t6 * t4) * t28 + (t2 * t1 + t7 * t3) * t29) * t8) * g(2) + (-mrSges(2,1) * t11 + mrSges(2,2) * t10 - t12 * r_base(1) - mrSges(1,1) - m(3) * (r_base(1) + t25) - m(4) * (pkin(1) + r_base(1)) + ((t4 * t1 + t6 * t5) * t28 + (t3 * t1 - t7 * t2) * t29) * t8) * g(1);
U = t9;
