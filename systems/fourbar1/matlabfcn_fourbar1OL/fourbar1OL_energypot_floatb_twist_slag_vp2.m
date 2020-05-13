% Calculate potential energy for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1OL_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:17
% EndTime: 2020-04-24 20:10:17
% DurationCPUTime: 0.09s
% Computational Cost: add. (36->27), mult. (36->29), div. (0->0), fcn. (8->6), ass. (0->6)
t6 = cos(qJ(2));
t5 = sin(qJ(2));
t3 = m(3) * pkin(2) + mrSges(2,1);
t2 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t1 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t4 = (-mrSges(2,2) * g(2) - t3 * g(1) + t1 * t6 + t2 * t5) * cos(qJ(1)) + (mrSges(2,2) * g(1) - g(2) * t3 - t1 * t5 + t2 * t6) * sin(qJ(1)) + (-mrSges(4,1) * g(1) - mrSges(4,2) * g(2)) * cos(qJ(3)) + (-mrSges(4,1) * g(2) + mrSges(4,2) * g(1)) * sin(qJ(3)) + (-m(4) * pkin(1) - mrSges(1,1)) * g(1) - mrSges(1,2) * g(2) - g(3) * (mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3)) + (-g(1) * r_base(1) - g(2) * r_base(2) - g(3) * r_base(3)) * (m(1) + m(2) + m(3) + m(4));
U = t4;
