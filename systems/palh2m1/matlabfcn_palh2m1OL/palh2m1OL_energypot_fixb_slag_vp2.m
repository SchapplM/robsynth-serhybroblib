% Calculate potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:58:40
% EndTime: 2020-05-02 23:58:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (153->40), mult. (159->40), div. (0->0), fcn. (72->10), ass. (0->22)
t46 = m(5) + m(6);
t45 = m(4) + t46;
t44 = m(3) + t45;
t37 = sin(qJ(1));
t41 = cos(qJ(1));
t43 = g(1) * t41 + g(2) * t37;
t34 = sin(qJ(5));
t38 = cos(qJ(5));
t42 = mrSges(6,1) * t34 + mrSges(6,2) * t38 + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t40 = cos(qJ(3));
t39 = cos(qJ(4));
t36 = sin(qJ(3));
t35 = sin(qJ(4));
t32 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t31 = t45 * pkin(2) + mrSges(3,1);
t30 = t44 * pkin(1) + mrSges(2,1);
t29 = pkin(4) * m(6) + t38 * mrSges(6,1) - t34 * mrSges(6,2) + mrSges(5,1);
t28 = t29 * t35 - t32 * t39 + mrSges(4,2);
t27 = t46 * pkin(3) + t29 * t39 + t32 * t35 + mrSges(4,1);
t26 = g(3) * t27 + t43 * t28;
t25 = -t28 * g(3) + t43 * t27;
t1 = (mrSges(3,2) * g(3) - t25 * t40 + t26 * t36 - t43 * t31) * cos(qJ(2)) + (t43 * mrSges(3,2) + g(3) * t31 + t25 * t36 + t26 * t40) * sin(qJ(2)) + (-t30 * g(1) - t42 * g(2)) * t41 + (t42 * g(1) - g(2) * t30) * t37 - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) - g(3) * (mrSges(1,3) + mrSges(2,3)) - g(3) * (m(2) + t44) * pkin(5);
U = t1;
