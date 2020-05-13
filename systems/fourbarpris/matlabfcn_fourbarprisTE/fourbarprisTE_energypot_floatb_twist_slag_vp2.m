% Calculate potential energy for
% fourbarprisTE
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
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbarprisTE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisTE_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 08:59:59
% EndTime: 2020-05-07 09:00:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (101->55), mult. (137->67), div. (3->4), fcn. (2->2), ass. (0->26)
t8 = (mrSges(3,1) - mrSges(2,2));
t9 = (mrSges(2,1) + mrSges(3,3));
t29 = (t9 * g(1) - g(2) * t8);
t11 = (pkin(3) ^ 2);
t36 = -2 * pkin(3);
t35 = (0.1e1 / pkin(1) / pkin(2));
t34 = (pkin(2) * m(3));
t28 = (mrSges(4,2) * g(2));
t3 = (mrSges(4,1) * g(1) + t28);
t33 = pkin(3) * t3;
t32 = (m(3) * g(1));
t30 = (m(2) + m(3));
t13 = (pkin(2) ^ 2);
t27 = t3 * t13;
t26 = qJ(1) + pkin(3);
t25 = -2 * (mrSges(1,1) * g(1) + mrSges(1,2) * g(2) - g(3) * (mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - mrSges(4,3))) * pkin(1);
t24 = -pkin(2) - t26;
t23 = -pkin(2) + t26;
t22 = -2 * (m(1) + m(4) + t30) / t35;
t21 = mrSges(4,1) * g(2) - mrSges(4,2) * g(1);
t20 = g(1) * r_base(1) * t22;
t19 = g(2) * r_base(2) * t22;
t18 = g(3) * r_base(3) * t22;
t14 = pkin(1) ^ 2;
t12 = pkin(2) * t13;
t1 = (((g(2) * t34 + t21) * qJ(1) + (g(1) * t8 + g(2) * t9) * pkin(2) + pkin(3) * t21) * sqrt(-((pkin(1) + t24) * (pkin(1) + t23) * (pkin(1) - t23) * (pkin(1) - t24))) + ((t20 + t19 + t18 - t12 * t32 - t27 + (t25 - t29 * t36 + (-2 * m(2) * t14 + (-t14 + t11) * m(3)) * g(1)) * pkin(2) - (-3 * t11 + t14) * t3) * qJ(1)) + (pkin(3) * t20) + (pkin(3) * t19) + (pkin(3) * t18) - (t29 * t12) - (pkin(3) * t27) + (((g(1) * t30 * t36 + t29) * t14 + pkin(3) * t25 + t29 * t11) * pkin(2)) - ((pkin(1) - pkin(3)) * (pkin(1) + pkin(3)) * t33) + (((t28 + (mrSges(4,1) + t34) * g(1)) * qJ(1) + (2 * pkin(3) * t32 + t29) * pkin(2) + 3 * t33) * qJ(1) ^ 2)) / t26 * t35 / 0.2e1;
U = t1;
