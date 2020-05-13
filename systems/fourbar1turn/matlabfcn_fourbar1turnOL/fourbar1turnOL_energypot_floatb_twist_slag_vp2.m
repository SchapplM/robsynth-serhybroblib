% Calculate potential energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
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
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnOL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (61->30), mult. (77->25), div. (0->0), fcn. (47->8), ass. (0->15)
t19 = -m(4) * pkin(2) - mrSges(3,1);
t17 = -m(2) - m(3) - m(4) - m(5);
t4 = qJ(2) + qJ(3);
t1 = sin(t4);
t2 = cos(t4);
t5 = sin(qJ(4));
t6 = sin(qJ(2));
t8 = cos(qJ(4));
t9 = cos(qJ(2));
t16 = -m(5) * pkin(1) + t2 * mrSges(4,1) - t8 * mrSges(5,1) + mrSges(3,2) * t6 - t1 * mrSges(4,2) + t5 * mrSges(5,2) + t19 * t9 - mrSges(2,1);
t15 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t14 = -m(1) + t17;
t10 = cos(qJ(1));
t7 = sin(qJ(1));
t3 = (-m(1) * r_base(3) + t1 * mrSges(4,1) - mrSges(5,1) * t5 - mrSges(3,2) * t9 + t2 * mrSges(4,2) - mrSges(5,2) * t8 - mrSges(1,3) - mrSges(2,3) + t19 * t6 + t17 * (pkin(5) + r_base(3))) * g(3) + (-t15 * t10 + t14 * r_base(2) + t16 * t7 - mrSges(1,2)) * g(2) + (t16 * t10 + t14 * r_base(1) + t15 * t7 - mrSges(1,1)) * g(1);
U = t3;
