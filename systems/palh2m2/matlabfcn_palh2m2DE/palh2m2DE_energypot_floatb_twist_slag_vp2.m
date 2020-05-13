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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2DE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:14
% EndTime: 2020-05-03 01:06:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (99->43), mult. (74->35), div. (0->0), fcn. (18->8), ass. (0->18)
t19 = m(5) + m(4);
t18 = m(7) + m(6);
t8 = pkin(1) + pkin(2);
t17 = t18 + t19;
t1 = t17 * pkin(4) + mrSges(3,1);
t11 = sin(qJ(3));
t12 = sin(qJ(2));
t14 = cos(qJ(3));
t15 = cos(qJ(2));
t3 = t18 * pkin(5) + mrSges(5,1);
t16 = -(m(3) + t19) * pkin(1) - m(5) * pkin(2) - t8 * m(6) - (pkin(3) + t8) * m(7) + t12 * mrSges(3,2) + mrSges(5,2) * t11 - t1 * t15 - t3 * t14 - mrSges(2,1) - mrSges(4,1) - mrSges(6,1);
t13 = cos(qJ(4));
t10 = sin(qJ(4));
t6 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t5 = -mrSges(7,1) * g(2) + mrSges(7,2) * g(1);
t4 = mrSges(7,1) * g(1) + mrSges(7,2) * g(2);
t2 = m(1) + m(2) + m(3) + t17;
t7 = (t16 * g(1) - g(2) * t6 + t5 * t10 - t4 * t13) * cos(qJ(1)) + (t6 * g(1) + t16 * g(2) + t4 * t10 + t5 * t13) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t2 + (-mrSges(3,2) * t15 - mrSges(5,2) * t14 - t1 * t12 - t3 * t11 - t2 * r_base(3) - mrSges(4,2) - mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(7,3)) * g(3);
U = t7;
