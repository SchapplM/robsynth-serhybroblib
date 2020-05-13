% Calculate potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:00:28
% EndTime: 2020-05-03 00:00:28
% DurationCPUTime: 0.28s
% Computational Cost: add. (146->59), mult. (150->69), div. (0->0), fcn. (81->12), ass. (0->28)
t36 = -m(3) - m(4);
t34 = -m(2) + t36;
t16 = sin(qJ(5));
t20 = cos(qJ(5));
t33 = t20 * rSges(6,1) - t16 * rSges(6,2) + pkin(4);
t30 = pkin(6) + rSges(6,3);
t15 = qJ(2) + qJ(3);
t19 = sin(qJ(1));
t22 = cos(qJ(2));
t5 = pkin(3) * cos(t15) + t22 * pkin(2) + pkin(1);
t29 = t19 * t5 + r_base(2);
t23 = cos(qJ(1));
t28 = t23 * t5 + r_base(1);
t27 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t13 = qJ(4) + t15;
t10 = cos(t13);
t9 = sin(t13);
t26 = rSges(5,1) * t10 - rSges(5,2) * t9;
t18 = sin(qJ(2));
t25 = -t18 * pkin(2) - pkin(3) * sin(t15) + pkin(5) + r_base(3);
t17 = sin(qJ(3));
t21 = cos(qJ(3));
t1 = rSges(3,1) * m(3) + (rSges(4,1) * t21 - rSges(4,2) * t17 + pkin(2)) * m(4);
t4 = rSges(3,2) * m(3) + (rSges(4,1) * t17 + rSges(4,2) * t21) * m(4);
t24 = -m(2) * rSges(2,1) + t36 * pkin(1) - t1 * t22 + t4 * t18;
t14 = m(1) - t34;
t6 = t16 * rSges(6,1) + t20 * rSges(6,2);
t2 = (t24 * g(1) - t27 * g(2)) * t23 + (t27 * g(1) + t24 * g(2)) * t19 - m(5) * (g(1) * (-t19 * rSges(5,3) + t26 * t23 + t28) + g(2) * (t23 * rSges(5,3) + t26 * t19 + t29)) - m(6) * (g(1) * (-t19 * t6 + t28) + g(2) * (t23 * t6 + t29) + (t33 * t10 + t30 * t9) * (g(1) * t23 + g(2) * t19)) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t14 + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (t4 * t22 + t1 * t18 - t14 * r_base(3) - m(2) * rSges(2,3) - m(1) * rSges(1,3) - m(5) * (-t9 * rSges(5,1) - t10 * rSges(5,2) + t25) - m(6) * (t30 * t10 - t33 * t9 + t25) + t34 * pkin(5)) * g(3);
U = t2;
