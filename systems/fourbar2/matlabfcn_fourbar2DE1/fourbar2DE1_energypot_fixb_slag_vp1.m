% Calculate potential energy for
% fourbar2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:22
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar2DE1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE1_energypot_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2DE1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE1_energypot_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2DE1_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2DE1_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:22:52
% EndTime: 2020-04-24 20:22:53
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->15), mult. (28->23), div. (0->0), fcn. (2->2), ass. (0->3)
t5 = m(2) * rSges(2,2) + m(4) * rSges(4,2);
t4 = m(2) * rSges(2,1) + m(3) * pkin(2) + m(4) * rSges(4,1);
t1 = (-t4 * g(1) - g(2) * t5) * cos(qJ(1)) + (t5 * g(1) - g(2) * t4) * sin(qJ(1)) + (-m(1) * rSges(1,1) - m(3) * rSges(3,1) - m(4) * pkin(1)) * g(1) + (-m(1) * rSges(1,2) - m(3) * rSges(3,2)) * g(2) - g(3) * (m(1) * rSges(1,3) + m(2) * rSges(2,3) + m(3) * rSges(3,3) + m(4) * rSges(4,3));
U = t1;
