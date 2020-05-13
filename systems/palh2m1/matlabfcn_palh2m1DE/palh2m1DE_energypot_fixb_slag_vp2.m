% Calculate potential energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1DE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:40
% EndTime: 2020-05-02 23:51:40
% DurationCPUTime: 0.17s
% Computational Cost: add. (80->38), mult. (79->31), div. (0->0), fcn. (24->8), ass. (0->16)
t34 = m(5) + m(6);
t33 = m(4) + t34;
t21 = t34 * pkin(3) + mrSges(4,1);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t19 = -mrSges(4,2) * t26 + t33 * pkin(2) + t21 * t29 + mrSges(3,1);
t20 = mrSges(4,2) * t29 + t21 * t26 + mrSges(3,2);
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t32 = (-pkin(4) - pkin(1)) * m(6) - t19 * t30 + t20 * t27 - mrSges(2,1) - mrSges(5,1) + (-m(3) - m(4) - m(5)) * pkin(1);
t28 = cos(qJ(4));
t25 = sin(qJ(4));
t24 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t23 = -mrSges(6,1) * g(2) + mrSges(6,2) * g(1);
t22 = mrSges(6,1) * g(1) + mrSges(6,2) * g(2);
t1 = (t32 * g(1) - g(2) * t24 - t22 * t28 + t23 * t25) * cos(qJ(1)) + (t24 * g(1) + t32 * g(2) + t22 * t25 + t23 * t28) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (t20 * t30 + t19 * t27 - m(6) * pkin(6) + mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(6,3) - (m(2) + m(3) + t33) * pkin(5)) * g(3);
U = t1;
