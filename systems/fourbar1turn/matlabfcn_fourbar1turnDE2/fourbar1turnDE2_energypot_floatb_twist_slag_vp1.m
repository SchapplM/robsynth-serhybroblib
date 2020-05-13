% Calculate potential energy for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:23
% EndTime: 2020-04-12 19:33:23
% DurationCPUTime: 0.29s
% Computational Cost: add. (697->70), mult. (955->84), div. (60->5), fcn. (299->11), ass. (0->29)
t39 = pkin(4) ^ 2;
t38 = -pkin(3) - pkin(4);
t37 = -pkin(3) + pkin(4);
t18 = cos(qJ(2));
t36 = pkin(2) * t18;
t16 = sin(qJ(2));
t35 = t16 * pkin(2);
t32 = (-0.2e1 * t36 + pkin(1)) * pkin(1);
t7 = sqrt(-((pkin(2) - t38) * (pkin(2) + t38) + t32) * ((pkin(2) - t37) * (pkin(2) + t37) + t32));
t34 = t16 * t7;
t12 = pkin(2) ^ 2 + t32;
t10 = 0.1e1 / t12;
t33 = t10 / pkin(3);
t31 = pkin(3) ^ 2 - t39;
t30 = pkin(5) + r_base(3);
t28 = rSges(3,1) * t18 - rSges(3,2) * t16;
t14 = pkin(1) * t18 - pkin(2);
t8 = t12 + t31;
t3 = qJ(2) + atan2((pkin(1) * t16 * t8 - t14 * t7) * t33, (-pkin(1) * t34 - t14 * t8) * t33);
t1 = sin(t3);
t2 = cos(t3);
t27 = -rSges(4,1) * t2 + rSges(4,2) * t1 + t36;
t19 = cos(qJ(1));
t17 = sin(qJ(1));
t13 = pkin(1) - t36;
t9 = t12 - t31;
t6 = t13 * t7 + t9 * t35;
t5 = -pkin(2) * t34 + t13 * t9;
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t17 * rSges(3,3) + t28 * t19 + r_base(1)) + g(2) * (-t19 * rSges(3,3) + t28 * t17 + r_base(2)) + g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t30)) - m(4) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (-t1 * rSges(4,1) - t2 * rSges(4,2) + t30 + t35) + (-g(2) * rSges(4,3) + g(1) * t27) * t19 + (g(1) * rSges(4,3) + g(2) * t27) * t17) - m(5) * (g(1) * (t17 * rSges(5,3) + t19 * pkin(1) + r_base(1)) + g(2) * (-t19 * rSges(5,3) + t17 * pkin(1) + r_base(2)) + g(3) * t30 + (g(3) * (rSges(5,1) * t6 - rSges(5,2) * t5) + (g(1) * t19 + g(2) * t17) * (-rSges(5,1) * t5 - rSges(5,2) * t6)) * ((t5 ^ 2 + t6 ^ 2) / t39 / t12 ^ 2) ^ (-0.1e1 / 0.2e1) / pkin(4) * t10);
U = t4;
