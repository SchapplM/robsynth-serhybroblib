% Calculate potential energy for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnTE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:28
% EndTime: 2020-04-12 19:18:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (439->70), mult. (605->89), div. (36->3), fcn. (195->6), ass. (0->31)
t14 = cos(qJ(2));
t10 = pkin(1) * t14 - pkin(2);
t12 = sin(qJ(2));
t35 = pkin(2) * t14;
t31 = (-0.2e1 * t35 + pkin(1)) * pkin(1);
t36 = -pkin(3) + pkin(4);
t37 = -pkin(3) - pkin(4);
t5 = sqrt(-((pkin(2) - t37) * (pkin(2) + t37) + t31) * ((pkin(2) - t36) * (pkin(2) + t36) + t31));
t33 = t12 * t5;
t28 = pkin(2) ^ 2 + t31;
t30 = pkin(3) ^ 2 - pkin(4) ^ 2;
t6 = t28 + t30;
t1 = -pkin(1) * t33 - t10 * t6;
t8 = 0.1e1 / t28;
t32 = 0.1e1 / pkin(3) * t8;
t38 = t12 / 0.2e1;
t4 = pkin(1) * t12 * t6 - t10 * t5;
t41 = (t14 * t4 / 0.2e1 + t1 * t38) * t32;
t24 = (t4 * t38 - t14 * t1 / 0.2e1) * t32;
t40 = rSges(4,1) * t24 + rSges(4,2) * t41 + t35;
t39 = -rSges(5,2) / 0.2e1;
t34 = t12 * pkin(2);
t29 = pkin(5) + r_base(3);
t27 = rSges(3,1) * t14 - rSges(3,2) * t12;
t15 = cos(qJ(1));
t13 = sin(qJ(1));
t9 = pkin(1) - t35;
t7 = t28 - t30;
t3 = t34 * t7 + t9 * t5;
t2 = -pkin(2) * t33 + t9 * t7;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t15 * rSges(2,1) - t13 * rSges(2,2) + r_base(1)) + g(2) * (t13 * rSges(2,1) + t15 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t29)) - m(3) * (g(1) * (t13 * rSges(3,3) + t15 * t27 + r_base(1)) + g(2) * (-t15 * rSges(3,3) + t13 * t27 + r_base(2)) + g(3) * (t12 * rSges(3,1) + t14 * rSges(3,2) + t29)) - m(4) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (-rSges(4,1) * t41 + rSges(4,2) * t24 + t29 + t34) + (-g(2) * rSges(4,3) + g(1) * t40) * t15 + (g(1) * rSges(4,3) + g(2) * t40) * t13) - m(5) * (g(1) * (t13 * rSges(5,3) + t15 * pkin(1) + r_base(1)) + g(2) * (-t15 * rSges(5,3) + t13 * pkin(1) + r_base(2)) + g(3) * t29 + (g(3) * (t3 * rSges(5,1) / 0.2e1 + t2 * t39) + (g(1) * t15 + g(2) * t13) * (-t2 * rSges(5,1) / 0.2e1 + t3 * t39)) * t8 / pkin(4));
U = t11;
