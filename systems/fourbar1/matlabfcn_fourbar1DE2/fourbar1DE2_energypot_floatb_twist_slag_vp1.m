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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE2_energypot_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:04:42
% EndTime: 2020-04-24 20:04:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (124->51), mult. (172->53), div. (8->3), fcn. (44->4), ass. (0->21)
t11 = cos(qJ(1));
t25 = pkin(2) * t11;
t23 = (-0.2e1 * t25 + pkin(1)) * pkin(1);
t21 = pkin(2) ^ 2 + t23;
t28 = 0.1e1 / t21 / 0.2e1;
t27 = -pkin(3) - pkin(4);
t26 = -pkin(3) + pkin(4);
t10 = sin(qJ(1));
t24 = t10 * pkin(2);
t22 = pkin(3) ^ 2 - pkin(4) ^ 2;
t20 = 0.1e1 / pkin(4) * t28;
t19 = 0.1e1 / pkin(3) * t28;
t18 = -pkin(1) + t25;
t7 = t21 - t22;
t6 = t21 + t22;
t5 = -t18 * rSges(3,1) + rSges(3,2) * t24;
t4 = rSges(3,1) * t24 + t18 * rSges(3,2);
t3 = rSges(4,1) * t24 + t18 * rSges(4,2);
t2 = -t18 * rSges(4,1) + rSges(4,2) * t24;
t1 = sqrt(-((pkin(2) - t27) * (pkin(2) + t27) + t23) * ((pkin(2) - t26) * (pkin(2) + t26) + t23));
t8 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t11 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t10 + g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (t25 + r_base(1) + (t4 * t1 + t6 * t5) * t19) + g(2) * (t24 + r_base(2) + (t5 * t1 - t6 * t4) * t19) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (pkin(1) + r_base(1) + (t3 * t1 - t7 * t2) * t20) + g(2) * (r_base(2) + (t2 * t1 + t7 * t3) * t20) + g(3) * (r_base(3) + rSges(4,3)));
U = t8;
