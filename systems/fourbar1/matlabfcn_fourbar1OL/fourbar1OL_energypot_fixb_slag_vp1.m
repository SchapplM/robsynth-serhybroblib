% Calculate potential energy for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1OL_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:18
% EndTime: 2020-04-24 20:10:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->21), mult. (33->32), div. (0->0), fcn. (8->6), ass. (0->4)
t9 = cos(qJ(1));
t8 = sin(qJ(1));
t7 = qJ(1) + qJ(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t9 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t8 + g(3) * rSges(2,3)) - ((-rSges(3,1) * g(1) - rSges(3,2) * g(2)) * cos(t7) + (-rSges(3,1) * g(2) + rSges(3,2) * g(1)) * sin(t7) + g(3) * rSges(3,3) + (g(1) * t9 + g(2) * t8) * pkin(2)) * m(3) - m(4) * ((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * cos(qJ(3)) + (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * sin(qJ(3)) + g(1) * pkin(1) + g(3) * rSges(4,3));
U = t1;
