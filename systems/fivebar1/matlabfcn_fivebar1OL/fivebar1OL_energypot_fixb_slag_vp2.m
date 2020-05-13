% Calculate potential energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:53
% EndTime: 2020-04-27 06:12:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->30), mult. (45->34), div. (0->0), fcn. (12->8), ass. (0->11)
t22 = cos(qJ(2));
t21 = cos(qJ(4));
t20 = sin(qJ(2));
t19 = sin(qJ(4));
t18 = pkin(2) * m(3) + mrSges(2,1);
t17 = pkin(3) * m(5) + mrSges(4,1);
t16 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t15 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t14 = -mrSges(5,1) * g(2) + mrSges(5,2) * g(1);
t13 = mrSges(5,1) * g(1) + mrSges(5,2) * g(2);
t1 = (-mrSges(2,2) * g(2) - t18 * g(1) + t15 * t22 + t16 * t20) * cos(qJ(1)) + (-mrSges(4,2) * g(2) - t17 * g(1) - t13 * t21 + t14 * t19) * cos(qJ(3)) + (mrSges(2,2) * g(1) - g(2) * t18 - t15 * t20 + t16 * t22) * sin(qJ(1)) + (mrSges(4,2) * g(1) - g(2) * t17 + t13 * t19 + t14 * t21) * sin(qJ(3)) + (-mrSges(1,1) + (-m(4) - m(5)) * pkin(1)) * g(1) - mrSges(1,2) * g(2) - g(3) * (mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3));
U = t1;
