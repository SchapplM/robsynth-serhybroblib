% Calculate potential energy for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar2OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_energypot_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:15
% EndTime: 2020-04-24 20:32:15
% DurationCPUTime: 0.02s
% Computational Cost: add. (24->21), mult. (30->25), div. (0->0), fcn. (8->6), ass. (0->6)
t11 = cos(qJ(2));
t10 = sin(qJ(2));
t9 = m(3) * pkin(2) + mrSges(2,1);
t8 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t7 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t1 = (-mrSges(2,2) * g(2) - t9 * g(1) + t8 * t10 + t7 * t11) * cos(qJ(1)) + (mrSges(2,2) * g(1) - g(2) * t9 - t7 * t10 + t8 * t11) * sin(qJ(1)) + (-mrSges(4,1) * g(1) - mrSges(4,2) * g(2)) * cos(qJ(3)) + (-mrSges(4,1) * g(2) + mrSges(4,2) * g(1)) * sin(qJ(3)) + (-m(4) * pkin(1) - mrSges(1,1)) * g(1) - mrSges(1,2) * g(2) - g(3) * (mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3));
U = t1;
