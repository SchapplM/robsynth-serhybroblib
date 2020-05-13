% Calculate potential energy for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbar1DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energypot_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE1_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:54:40
% EndTime: 2020-04-24 19:54:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (112->38), mult. (168->42), div. (8->3), fcn. (44->4), ass. (0->22)
t57 = -0.1e1 / pkin(4) / 0.2e1;
t56 = -0.1e1 / pkin(3) / 0.2e1;
t55 = -pkin(3) - pkin(4);
t54 = -pkin(3) + pkin(4);
t40 = cos(qJ(1));
t53 = pkin(2) * t40;
t39 = sin(qJ(1));
t52 = t39 * pkin(2);
t51 = (-0.2e1 * t53 + pkin(1)) * pkin(1);
t50 = pkin(3) ^ 2 - pkin(4) ^ 2;
t49 = pkin(2) ^ 2 + t51;
t48 = -m(3) * pkin(2) - mrSges(2,1);
t47 = -pkin(1) + t53;
t37 = 0.1e1 / t49;
t36 = t49 - t50;
t35 = t49 + t50;
t34 = -t47 * mrSges(3,1) + mrSges(3,2) * t52;
t33 = mrSges(3,1) * t52 + t47 * mrSges(3,2);
t32 = mrSges(4,1) * t52 + t47 * mrSges(4,2);
t31 = -t47 * mrSges(4,1) + mrSges(4,2) * t52;
t30 = sqrt(-((pkin(2) - t55) * (pkin(2) + t55) + t51) * ((pkin(2) - t54) * (pkin(2) + t54) + t51));
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3)) * g(3) + (-mrSges(2,2) * t40 - mrSges(1,2) + t48 * t39 + ((t34 * t30 - t35 * t33) * t56 + (t31 * t30 + t36 * t32) * t57) * t37) * g(2) + (-m(4) * pkin(1) + mrSges(2,2) * t39 - mrSges(1,1) + t48 * t40 + ((t33 * t30 + t35 * t34) * t56 + (t32 * t30 - t36 * t31) * t57) * t37) * g(1);
U = t1;
