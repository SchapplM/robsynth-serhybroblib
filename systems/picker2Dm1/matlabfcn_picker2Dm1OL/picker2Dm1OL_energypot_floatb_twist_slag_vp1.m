% Calculate potential energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1OL_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:48
% EndTime: 2020-05-11 05:44:49
% DurationCPUTime: 0.30s
% Computational Cost: add. (190->106), mult. (112->96), div. (0->0), fcn. (68->22), ass. (0->35)
t24 = qJ(1) + qJ(2);
t21 = qJ(3) + t24;
t20 = qJ(4) + t24;
t26 = sin(qJ(1));
t34 = -t26 * pkin(1) + r_base(2);
t28 = cos(qJ(1));
t33 = -t28 * pkin(1) + r_base(1);
t16 = sin(t24);
t32 = -pkin(2) * t16 + t34;
t18 = cos(t24);
t31 = -pkin(2) * t18 + t33;
t30 = -pkin(3) * t16 + t34;
t29 = -pkin(3) * t18 + t33;
t27 = cos(qJ(7));
t25 = sin(qJ(7));
t23 = qJ(1) + qJ(8);
t22 = pkin(8) + qJ(5);
t19 = qJ(6) + t24;
t17 = cos(t23);
t15 = sin(t23);
t14 = cos(t22);
t13 = sin(t22);
t12 = qJ(9) + t21;
t11 = qJ(10) + t20;
t10 = cos(t21);
t9 = cos(t20);
t8 = cos(t19);
t7 = sin(t21);
t6 = sin(t20);
t5 = sin(t19);
t4 = cos(t12);
t3 = sin(t12);
t2 = cos(t11);
t1 = sin(t11);
t35 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (-t28 * rSges(2,1) + t26 * rSges(2,2) + r_base(1)) + g(2) * (-t26 * rSges(2,1) - t28 * rSges(2,2) + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (-t18 * rSges(3,1) + t16 * rSges(3,2) + t33) + g(2) * (-t16 * rSges(3,1) - t18 * rSges(3,2) + t34) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (t10 * rSges(4,1) - t7 * rSges(4,2) + t31) + g(2) * (t7 * rSges(4,1) + t10 * rSges(4,2) + t32) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (-t9 * rSges(5,1) + t6 * rSges(5,2) + t29) + g(2) * (-t6 * rSges(5,1) - t9 * rSges(5,2) + t30) + g(3) * (r_base(3) + rSges(5,3))) - m(6) * (g(1) * (cos(pkin(8)) * pkin(5) + r_base(1) + t14 * rSges(6,1) - t13 * rSges(6,2)) + g(2) * (sin(pkin(8)) * pkin(5) + r_base(2) + t13 * rSges(6,1) + t14 * rSges(6,2)) + g(3) * (r_base(3) + rSges(6,3))) - m(7) * (g(1) * (t8 * rSges(7,1) - t5 * rSges(7,2) + t33) + g(2) * (t5 * rSges(7,1) + t8 * rSges(7,2) + t34) + g(3) * (r_base(3) + rSges(7,3))) - m(8) * (g(1) * (t25 * rSges(8,1) + t27 * rSges(8,2) + pkin(7) + r_base(1)) + g(2) * (-t27 * rSges(8,1) + t25 * rSges(8,2) + r_base(2)) + g(3) * (r_base(3) + rSges(8,3))) - m(9) * (g(1) * (t17 * rSges(9,1) - t15 * rSges(9,2) + t33) + g(2) * (t15 * rSges(9,1) + t17 * rSges(9,2) + t34) + g(3) * (r_base(3) + rSges(9,3))) - m(10) * (g(1) * (pkin(6) * t10 - t4 * rSges(10,1) + t3 * rSges(10,2) + t31) + g(2) * (pkin(6) * t7 - t3 * rSges(10,1) - t4 * rSges(10,2) + t32) + g(3) * (r_base(3) + rSges(10,3))) - m(11) * (g(1) * (-pkin(4) * t9 + t2 * rSges(11,1) - t1 * rSges(11,2) + t29) + g(2) * (-pkin(4) * t6 + t1 * rSges(11,1) + t2 * rSges(11,2) + t30) + g(3) * (r_base(3) + rSges(11,3)));
U = t35;
