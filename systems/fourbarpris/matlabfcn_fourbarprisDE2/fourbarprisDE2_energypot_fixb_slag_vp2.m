% Calculate potential energy for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
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
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbarprisDE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energypot_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energypot_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE2_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:43:59
% EndTime: 2020-05-07 09:43:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (77->47), mult. (104->60), div. (3->3), fcn. (2->2), ass. (0->20)
t43 = (mrSges(3,1) - mrSges(2,2));
t44 = (mrSges(2,1) + mrSges(3,3));
t57 = (t44 * g(1) - g(2) * t43);
t46 = (pkin(3) ^ 2);
t65 = -2 * pkin(3);
t64 = (pkin(2) * m(3));
t63 = (m(3) * g(1));
t61 = (mrSges(4,2) * g(2));
t39 = (mrSges(4,1) * g(1) + t61);
t62 = pkin(3) * t39;
t48 = (pkin(2) ^ 2);
t59 = t39 * t48;
t58 = qJ(1) + pkin(3);
t56 = -2 * (mrSges(1,1) * g(1) + mrSges(1,2) * g(2) - g(3) * (mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - mrSges(4,3))) * pkin(1);
t55 = -pkin(2) - t58;
t54 = -pkin(2) + t58;
t53 = mrSges(4,1) * g(2) - mrSges(4,2) * g(1);
t49 = pkin(1) ^ 2;
t47 = pkin(2) * t48;
t1 = (((g(2) * t64 + t53) * qJ(1) + (t43 * g(1) + g(2) * t44) * pkin(2) + pkin(3) * t53) * sqrt(-((pkin(1) + t55) * (pkin(1) + t54) * (pkin(1) - t54) * (pkin(1) - t55))) + ((-t47 * t63 - t59 + (t56 - t57 * t65 + (-2 * m(2) * t49 + (-t49 + t46) * m(3)) * g(1)) * pkin(2) - (-3 * t46 + t49) * t39) * qJ(1)) - (t57 * t47) - (pkin(3) * t59) + (((g(1) * (m(2) + m(3)) * t65 + t57) * t49 + pkin(3) * t56 + t57 * t46) * pkin(2)) - ((pkin(1) - pkin(3)) * (pkin(1) + pkin(3)) * t62) + (((t61 + (mrSges(4,1) + t64) * g(1)) * qJ(1) + (2 * pkin(3) * t63 + t57) * pkin(2) + 3 * t62) * qJ(1) ^ 2)) / pkin(1) / t58 / pkin(2) / 0.2e1;
U = t1;
