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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnOL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.08s
% Computational Cost: add. (46->36), mult. (67->50), div. (0->0), fcn. (47->8), ass. (0->13)
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t26 = rSges(3,1) * t22 - rSges(3,2) * t19;
t18 = sin(qJ(4));
t21 = cos(qJ(4));
t25 = rSges(5,1) * t21 - rSges(5,2) * t18 + pkin(1);
t17 = qJ(2) + qJ(3);
t15 = sin(t17);
t16 = cos(t17);
t24 = -rSges(4,1) * t16 + rSges(4,2) * t15 + pkin(2) * t22;
t23 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t20 * rSges(3,3) + t26 * t23) + g(2) * (-t23 * rSges(3,3) + t26 * t20) + g(3) * (t19 * rSges(3,1) + t22 * rSges(3,2) + pkin(5))) - m(4) * (g(3) * (-t15 * rSges(4,1) - t16 * rSges(4,2) + t19 * pkin(2) + pkin(5)) + (-g(2) * rSges(4,3) + g(1) * t24) * t23 + (g(1) * rSges(4,3) + g(2) * t24) * t20) - m(5) * (g(3) * (t18 * rSges(5,1) + t21 * rSges(5,2) + pkin(5)) + (-g(2) * rSges(5,3) + g(1) * t25) * t23 + (g(1) * rSges(5,3) + g(2) * t25) * t20);
U = t1;
