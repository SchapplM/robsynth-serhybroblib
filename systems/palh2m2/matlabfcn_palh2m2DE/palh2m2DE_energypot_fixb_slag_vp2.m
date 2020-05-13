% Calculate potential energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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

function U = palh2m2DE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:14
% EndTime: 2020-05-03 01:06:14
% DurationCPUTime: 0.14s
% Computational Cost: add. (78->39), mult. (68->31), div. (0->0), fcn. (18->8), ass. (0->16)
t36 = m(5) + m(4);
t35 = m(7) + m(6);
t26 = pkin(1) + pkin(2);
t20 = mrSges(3,1) + (t35 + t36) * pkin(4);
t21 = t35 * pkin(5) + mrSges(5,1);
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t32 = cos(qJ(3));
t33 = cos(qJ(2));
t34 = -pkin(1) * (m(3) + t36) - m(5) * pkin(2) - m(6) * t26 - m(7) * (pkin(3) + t26) + t30 * mrSges(3,2) + mrSges(5,2) * t29 - t20 * t33 - t21 * t32 - mrSges(2,1) - mrSges(4,1) - mrSges(6,1);
t31 = cos(qJ(4));
t28 = sin(qJ(4));
t24 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t23 = -mrSges(7,1) * g(2) + mrSges(7,2) * g(1);
t22 = mrSges(7,1) * g(1) + mrSges(7,2) * g(2);
t1 = (t34 * g(1) - g(2) * t24 - t22 * t31 + t23 * t28) * cos(qJ(1)) + (t24 * g(1) + t34 * g(2) + t22 * t28 + t23 * t31) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (-mrSges(3,2) * t33 - mrSges(5,2) * t32 - t20 * t30 - t21 * t29 - mrSges(4,2) - mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(7,3)) * g(3);
U = t1;
