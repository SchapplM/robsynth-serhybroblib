% Calculate potential energy for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:25
% EndTime: 2020-04-12 19:25:25
% DurationCPUTime: 0.40s
% Computational Cost: add. (1213->72), mult. (1739->93), div. (108->6), fcn. (519->10), ass. (0->31)
t45 = pkin(4) ^ 2;
t44 = pkin(3) ^ 2;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t43 = g(1) * t19 + g(2) * t17;
t42 = -pkin(3) - pkin(4);
t41 = -pkin(3) + pkin(4);
t18 = cos(qJ(2));
t38 = pkin(2) * t18;
t16 = sin(qJ(2));
t37 = t16 * pkin(2);
t35 = (-0.2e1 * t38 + pkin(1)) * pkin(1);
t7 = sqrt(-((pkin(2) - t42) * (pkin(2) + t42) + t35) * ((pkin(2) - t41) * (pkin(2) + t41) + t35));
t36 = t16 * t7;
t34 = t44 - t45;
t33 = pkin(5) + r_base(3);
t12 = pkin(2) ^ 2 + t35;
t14 = pkin(1) * t18 - pkin(2);
t8 = t12 + t34;
t3 = -pkin(1) * t36 - t14 * t8;
t6 = pkin(1) * t16 * t8 - t14 * t7;
t31 = t16 * t6 - t18 * t3;
t30 = t16 * t3 + t18 * t6;
t29 = rSges(3,1) * t18 - rSges(3,2) * t16;
t13 = pkin(1) - t38;
t11 = 0.1e1 / t12 ^ 2;
t10 = 0.1e1 / t12;
t9 = t12 - t34;
t5 = t13 * t7 + t37 * t9;
t4 = -pkin(2) * t36 + t13 * t9;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t17 * rSges(3,3) + t19 * t29 + r_base(1)) + g(2) * (-t19 * rSges(3,3) + t17 * t29 + r_base(2)) + g(3) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t33)) - m(4) * (g(1) * (t17 * rSges(4,3) + t19 * t38 + r_base(1)) + g(2) * (-t19 * rSges(4,3) + t17 * t38 + r_base(2)) + g(3) * (t33 + t37) + (g(3) * (-rSges(4,1) * t30 + rSges(4,2) * t31) + t43 * (rSges(4,1) * t31 + rSges(4,2) * t30)) / pkin(3) * ((t3 ^ 2 + t6 ^ 2) / t44 * t11) ^ (-0.1e1 / 0.2e1) * t10) - m(5) * (g(1) * (t17 * rSges(5,3) + t19 * pkin(1) + r_base(1)) + g(2) * (-t19 * rSges(5,3) + t17 * pkin(1) + r_base(2)) + g(3) * t33 + (g(3) * (rSges(5,1) * t5 - rSges(5,2) * t4) + t43 * (-rSges(5,1) * t4 - rSges(5,2) * t5)) / pkin(4) * t10 * ((t4 ^ 2 + t5 ^ 2) / t45 * t11) ^ (-0.1e1 / 0.2e1));
U = t1;
