% Calculate potential energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1DE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:40
% EndTime: 2020-05-02 23:51:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (98->41), mult. (85->35), div. (0->0), fcn. (24->8), ass. (0->18)
t18 = m(5) + m(6);
t17 = m(4) + t18;
t16 = m(2) + m(3) + t17;
t12 = cos(qJ(3));
t3 = t18 * pkin(3) + mrSges(4,1);
t9 = sin(qJ(3));
t1 = -mrSges(4,2) * t9 + t17 * pkin(2) + t3 * t12 + mrSges(3,1);
t10 = sin(qJ(2));
t13 = cos(qJ(2));
t2 = mrSges(4,2) * t12 + t3 * t9 + mrSges(3,2);
t15 = (-pkin(4) - pkin(1)) * m(6) - t1 * t13 + t2 * t10 - mrSges(2,1) - mrSges(5,1) + (-m(3) - m(4) - m(5)) * pkin(1);
t11 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t6 = m(1) + t16;
t5 = -mrSges(6,1) * g(2) + mrSges(6,2) * g(1);
t4 = mrSges(6,1) * g(1) + mrSges(6,2) * g(2);
t14 = (t15 * g(1) - g(2) * t7 - t4 * t11 + t5 * t8) * cos(qJ(1)) + (t7 * g(1) + t15 * g(2) + t5 * t11 + t4 * t8) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t6 + (-m(6) * pkin(6) - t16 * pkin(5) + t1 * t10 + t2 * t13 - t6 * r_base(3) + mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(6,3)) * g(3);
U = t14;
