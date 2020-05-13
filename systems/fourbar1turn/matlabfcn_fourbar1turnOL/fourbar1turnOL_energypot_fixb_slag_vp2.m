% Calculate potential energy for
% fourbar1turnOL
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnOL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (46->27), mult. (70->22), div. (0->0), fcn. (47->8), ass. (0->13)
t31 = -m(4) * pkin(2) - mrSges(3,1);
t30 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t22 = qJ(2) + qJ(3);
t20 = sin(t22);
t21 = cos(t22);
t23 = sin(qJ(4));
t24 = sin(qJ(2));
t26 = cos(qJ(4));
t27 = cos(qJ(2));
t29 = -m(5) * pkin(1) + mrSges(4,1) * t21 - mrSges(5,1) * t26 + mrSges(3,2) * t24 - mrSges(4,2) * t20 + mrSges(5,2) * t23 + t31 * t27 - mrSges(2,1);
t28 = cos(qJ(1));
t25 = sin(qJ(1));
t1 = (t20 * mrSges(4,1) - t23 * mrSges(5,1) - t27 * mrSges(3,2) + t21 * mrSges(4,2) - t26 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) + t31 * t24 + (-m(2) - m(3) - m(4) - m(5)) * pkin(5)) * g(3) + (t29 * t25 - t30 * t28 - mrSges(1,2)) * g(2) + (t30 * t25 + t29 * t28 - mrSges(1,1)) * g(1);
U = t1;
