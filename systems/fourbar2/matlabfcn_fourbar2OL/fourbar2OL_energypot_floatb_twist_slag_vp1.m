% Calculate potential energy for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar2OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2OL_energypot_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:16
% EndTime: 2020-04-24 20:32:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->31), mult. (40->34), div. (0->0), fcn. (8->6), ass. (0->5)
t6 = g(1) * r_base(1) + g(2) * r_base(2);
t5 = cos(qJ(1));
t4 = sin(qJ(1));
t3 = qJ(1) + qJ(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t5 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t4 + g(3) * (r_base(3) + rSges(2,3)) + t6) * m(2) - ((-rSges(3,1) * g(1) - rSges(3,2) * g(2)) * cos(t3) + (-rSges(3,1) * g(2) + rSges(3,2) * g(1)) * sin(t3) + g(3) * (r_base(3) + rSges(3,3)) + (g(1) * t5 + g(2) * t4) * pkin(2) + t6) * m(3) - ((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * cos(qJ(3)) + (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * sin(qJ(3)) + g(1) * pkin(1) + (r_base(3) + rSges(4,3)) * g(3) + t6) * m(4);
U = t1;
