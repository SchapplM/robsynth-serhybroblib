% Calculate potential energy for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbarprisDE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_energypot_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:44:52
% EndTime: 2020-05-07 09:44:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (101->48), mult. (163->61), div. (4->3), fcn. (4->2), ass. (0->17)
t5 = (rSges(2,1) * m(2) + rSges(3,3) * m(3));
t31 = (m(3) * qJ(1) + t5);
t21 = (qJ(1) + pkin(3));
t29 = -g(1) * r_base(1) - g(2) * r_base(2);
t28 = 2 * pkin(3);
t27 = m(4) / 0.2e1;
t26 = pkin(3) * m(3);
t4 = -rSges(2,2) * m(2) + rSges(3,1) * m(3);
t25 = g(2) * t4;
t24 = -m(2) - m(3);
t23 = (pkin(2) - pkin(3)) * (pkin(2) + pkin(3));
t18 = -pkin(2) - t21;
t17 = -pkin(2) + t21;
t12 = (pkin(1) ^ 2);
t11 = (qJ(1) ^ 2);
t1 = sqrt(-((pkin(1) + t18) * (pkin(1) + t17) * (pkin(1) - t17) * (pkin(1) - t18)));
t2 = (((t4 * g(1) + t31 * g(2)) * t1 + ((-2 * pkin(3) * t25 + (-2 * m(2) * t12 + (2 * t5 + t26) * pkin(3) + (-pkin(2) ^ 2 - t12) * m(3)) * g(1)) * qJ(1)) + (((t24 * t28 + t5) * g(1) - t25) * t12) - ((t5 * g(1) - t25) * t23) + ((-t25 + (2 * t26 + t31) * g(1)) * t11)) / t21 / 0.2e1 + (0.2e1 * ((-rSges(4,3) - r_base(3)) * g(3) + t29) * t27 + (-rSges(2,3) * m(2) + rSges(3,2) * m(3)) * g(3) + (-(rSges(1,1) * g(1)) - (rSges(1,2) * g(2)) - rSges(1,3) * g(3)) * m(1) + (m(1) - t24) * (-g(3) * r_base(3) + t29)) * pkin(1) + ((rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t1 + ((qJ(1) * t28 + t11 - t12 - t23) * (rSges(4,1) * g(1) + rSges(4,2) * g(2)))) * t27 / pkin(2)) / pkin(1);
U = t2;
