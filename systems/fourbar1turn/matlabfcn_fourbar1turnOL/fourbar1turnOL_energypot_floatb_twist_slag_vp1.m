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

function U = fourbar1turnOL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (61->47), mult. (67->52), div. (0->0), fcn. (47->8), ass. (0->15)
t14 = pkin(5) + r_base(3);
t5 = sin(qJ(2));
t8 = cos(qJ(2));
t13 = rSges(3,1) * t8 - rSges(3,2) * t5;
t4 = sin(qJ(4));
t7 = cos(qJ(4));
t12 = rSges(5,1) * t7 - rSges(5,2) * t4 + pkin(1);
t11 = g(1) * r_base(1) + g(2) * r_base(2);
t3 = qJ(2) + qJ(3);
t1 = sin(t3);
t2 = cos(t3);
t10 = -rSges(4,1) * t2 + rSges(4,2) * t1 + pkin(2) * t8;
t9 = cos(qJ(1));
t6 = sin(qJ(1));
t15 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t9 - t6 * rSges(2,2) + r_base(1)) + g(2) * (t6 * rSges(2,1) + rSges(2,2) * t9 + r_base(2)) + g(3) * (rSges(2,3) + t14)) - m(3) * (g(1) * (t6 * rSges(3,3) + t13 * t9 + r_base(1)) + g(2) * (-rSges(3,3) * t9 + t13 * t6 + r_base(2)) + g(3) * (rSges(3,1) * t5 + rSges(3,2) * t8 + t14)) - m(4) * (g(3) * (-rSges(4,1) * t1 - rSges(4,2) * t2 + pkin(2) * t5 + t14) + (-g(2) * rSges(4,3) + g(1) * t10) * t9 + (g(1) * rSges(4,3) + g(2) * t10) * t6 + t11) - m(5) * (g(3) * (rSges(5,1) * t4 + rSges(5,2) * t7 + t14) + (-g(2) * rSges(5,3) + g(1) * t12) * t9 + (g(1) * rSges(5,3) + g(2) * t12) * t6 + t11);
U = t15;
