% Calculate potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:18:27
% EndTime: 2020-05-03 01:18:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (188->66), mult. (171->79), div. (0->0), fcn. (101->14), ass. (0->32)
t45 = -m(3) - m(4);
t22 = qJ(2) + qJ(3);
t18 = qJ(4) + t22;
t15 = cos(t18);
t29 = cos(qJ(2));
t4 = pkin(2) * cos(t22) + t29 * pkin(4) + pkin(1);
t44 = pkin(5) * t15 + t4;
t43 = -g(1) * r_base(1) - g(2) * r_base(2);
t39 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t25 = sin(qJ(2));
t38 = t25 * pkin(4) + pkin(2) * sin(t22) + r_base(3);
t14 = sin(t18);
t37 = pkin(5) * t14 + t38;
t23 = sin(qJ(6));
t27 = cos(qJ(6));
t36 = t23 * rSges(7,1) + t27 * rSges(7,2);
t16 = qJ(5) + t18;
t7 = sin(t16);
t8 = cos(t16);
t35 = rSges(6,1) * t8 - rSges(6,2) * t7 + t44;
t34 = rSges(5,1) * t15 - rSges(5,2) * t14 + t4;
t24 = sin(qJ(3));
t28 = cos(qJ(3));
t2 = m(3) * rSges(3,1) + (rSges(4,1) * t28 - rSges(4,2) * t24 + pkin(4)) * m(4);
t3 = rSges(3,2) * m(3) + (rSges(4,1) * t24 + rSges(4,2) * t28) * m(4);
t32 = t45 * pkin(1) - m(2) * rSges(2,1) - t2 * t29 + t3 * t25;
t5 = t27 * rSges(7,1) - rSges(7,2) * t23 + pkin(3);
t31 = -t7 * rSges(7,3) + t5 * t8 + t44;
t30 = cos(qJ(1));
t26 = sin(qJ(1));
t19 = m(1) + m(2) - t45;
t1 = (t32 * g(1) + t39 * g(2)) * t30 + (-t39 * g(1) + t32 * g(2)) * t26 - m(5) * ((-g(2) * rSges(5,3) + g(1) * t34) * t30 + (g(1) * rSges(5,3) + g(2) * t34) * t26 - t43) - m(6) * ((-g(2) * rSges(6,3) + g(1) * t35) * t30 + (g(1) * rSges(6,3) + g(2) * t35) * t26 - t43) - m(7) * ((g(1) * t31 + g(2) * t36) * t30 + (-g(1) * t36 + g(2) * t31) * t26 - t43) + t43 * t19 + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (-t3 * t29 - t2 * t25 - t19 * r_base(3) - m(1) * rSges(1,3) - m(2) * rSges(2,3) - m(5) * (t14 * rSges(5,1) + t15 * rSges(5,2) + t38) - m(6) * (t7 * rSges(6,1) + t8 * rSges(6,2) + t37) - m(7) * (t8 * rSges(7,3) + t5 * t7 + t37)) * g(3);
U = t1;
