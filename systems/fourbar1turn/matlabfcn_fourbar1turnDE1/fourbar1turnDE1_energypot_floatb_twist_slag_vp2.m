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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1turnDE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:25
% EndTime: 2020-04-12 19:25:25
% DurationCPUTime: 0.36s
% Computational Cost: add. (1213->54), mult. (1749->67), div. (108->6), fcn. (519->10), ass. (0->34)
t48 = pkin(4) ^ 2;
t47 = pkin(3) ^ 2;
t46 = -m(2) - m(3) - m(4) - m(5);
t45 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t44 = -m(1) + t46;
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t14 = pkin(1) * t19 - pkin(2);
t39 = pkin(2) * t19;
t37 = (-0.2e1 * t39 + pkin(1)) * pkin(1);
t41 = -pkin(3) + pkin(4);
t42 = -pkin(3) - pkin(4);
t7 = sqrt(-((pkin(2) - t42) * (pkin(2) + t42) + t37) * ((pkin(2) - t41) * (pkin(2) + t41) + t37));
t38 = t17 * t7;
t12 = pkin(2) ^ 2 + t37;
t36 = t47 - t48;
t8 = t12 + t36;
t3 = -pkin(1) * t38 - t14 * t8;
t6 = pkin(1) * t17 * t8 - t14 * t7;
t32 = t17 * t3 + t19 * t6;
t33 = t17 * t6 - t19 * t3;
t10 = 0.1e1 / t12;
t11 = 0.1e1 / t12 ^ 2;
t34 = t10 * ((t3 ^ 2 + t6 ^ 2) / t47 * t11) ^ (-0.1e1 / 0.2e1) / pkin(3);
t13 = pkin(1) - t39;
t9 = t12 - t36;
t4 = -pkin(2) * t38 + t13 * t9;
t40 = pkin(2) * t17;
t5 = t13 * t7 + t9 * t40;
t35 = ((t4 ^ 2 + t5 ^ 2) / t48 * t11) ^ (-0.1e1 / 0.2e1) * t10 / pkin(4);
t43 = -m(4) * t39 - m(5) * pkin(1) - mrSges(2,1) - (t33 * mrSges(4,1) + t32 * mrSges(4,2)) * t34 - (-mrSges(5,1) * t4 - mrSges(5,2) * t5) * t35 - mrSges(3,1) * t19 + mrSges(3,2) * t17;
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t17 * mrSges(3,1) - t19 * mrSges(3,2) - m(4) * t40 - (-t32 * mrSges(4,1) + t33 * mrSges(4,2)) * t34 - (mrSges(5,1) * t5 - mrSges(5,2) * t4) * t35 + t46 * (pkin(5) + r_base(3))) * g(3) + (t43 * t18 - t45 * t20 + t44 * r_base(2) - mrSges(1,2)) * g(2) + (t45 * t18 + t43 * t20 + t44 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
