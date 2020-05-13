% Calculate potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
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
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:12:46
% EndTime: 2020-05-03 01:12:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (179->49), mult. (140->40), div. (0->0), fcn. (102->14), ass. (0->25)
t86 = -m(7) - m(6);
t85 = mrSges(6,2) + mrSges(7,3);
t68 = sin(qJ(6));
t72 = cos(qJ(6));
t84 = -pkin(3) * m(7) - t72 * mrSges(7,1) + t68 * mrSges(7,2) - mrSges(6,1);
t83 = t68 * mrSges(7,1) + t72 * mrSges(7,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t49 = pkin(4) * m(4) + t73 * mrSges(4,1) - mrSges(4,2) * t69 + mrSges(3,1);
t67 = qJ(2) + qJ(3);
t74 = cos(qJ(2));
t50 = pkin(2) * cos(t67) + t74 * pkin(4) + pkin(1);
t51 = t69 * mrSges(4,1) + t73 * mrSges(4,2) + mrSges(3,2);
t63 = qJ(4) + t67;
t62 = qJ(5) + t63;
t56 = sin(t62);
t57 = cos(t62);
t60 = sin(t63);
t61 = cos(t63);
t70 = sin(qJ(2));
t82 = -m(5) * t50 - mrSges(5,1) * t61 + mrSges(5,2) * t60 - t49 * t74 + t51 * t70 - mrSges(2,1) - (m(3) + m(4)) * pkin(1) + t86 * (pkin(5) * t61 + t50) + t84 * t57 + t85 * t56;
t81 = pkin(2) * sin(t67) + t70 * pkin(4);
t75 = cos(qJ(1));
t71 = sin(qJ(1));
t1 = (t82 * t71 - t83 * t75 - mrSges(1,2)) * g(2) + (t83 * t71 + t82 * t75 - mrSges(1,1)) * g(1) + (-m(5) * t81 - t60 * mrSges(5,1) - t61 * mrSges(5,2) - t49 * t70 - t51 * t74 - mrSges(1,3) - mrSges(2,3) + t86 * (pkin(5) * t60 + t81) - t85 * t57 + t84 * t56) * g(3);
U = t1;
