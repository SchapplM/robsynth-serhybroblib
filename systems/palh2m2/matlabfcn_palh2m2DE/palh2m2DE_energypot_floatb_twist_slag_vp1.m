% Calculate potential energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2DE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:21
% EndTime: 2020-05-03 01:06:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (99->45), mult. (114->54), div. (0->0), fcn. (18->8), ass. (0->17)
t19 = m(6) + m(7);
t9 = pkin(1) + pkin(2);
t18 = m(4) + m(5) + t19;
t12 = sin(qJ(3));
t13 = sin(qJ(2));
t15 = cos(qJ(3));
t16 = cos(qJ(2));
t2 = t18 * pkin(4) + m(3) * rSges(3,1);
t3 = t19 * pkin(5) + m(5) * rSges(5,1);
t17 = -m(2) * rSges(2,1) - (rSges(4,1) + pkin(1)) * m(4) - (rSges(6,1) + t9) * m(6) - (pkin(3) + t9) * m(7) - t3 * t15 - t2 * t16 + (rSges(5,2) * t12 - t9) * m(5) + (rSges(3,2) * t13 - pkin(1)) * m(3);
t14 = cos(qJ(4));
t11 = sin(qJ(4));
t6 = -rSges(7,1) * g(2) + rSges(7,2) * g(1);
t5 = rSges(7,1) * g(1) + rSges(7,2) * g(2);
t4 = m(1) + m(2) + m(3) + t18;
t1 = m(2) * rSges(2,2) - m(3) * rSges(3,3) - m(4) * rSges(4,3) - m(5) * rSges(5,3) - m(6) * rSges(6,3);
t7 = (-g(2) * t1 + (t6 * t11 - t5 * t14) * m(7) + t17 * g(1)) * cos(qJ(1)) + (t1 * g(1) + (t5 * t11 + t6 * t14) * m(7) + t17 * g(2)) * sin(qJ(1)) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t4 + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (-m(3) * t16 * rSges(3,2) - m(5) * t15 * rSges(5,2) - m(1) * rSges(1,3) - m(2) * rSges(2,3) - m(4) * rSges(4,2) - m(6) * rSges(6,2) - m(7) * rSges(7,3) - t3 * t12 - t2 * t13 - t4 * r_base(3)) * g(3);
U = t7;
