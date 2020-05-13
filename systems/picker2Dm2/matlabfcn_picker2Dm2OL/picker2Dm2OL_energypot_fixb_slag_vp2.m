% Calculate potential energy for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_energypot_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:48
% EndTime: 2020-05-09 23:18:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (157->74), mult. (118->52), div. (0->0), fcn. (68->22), ass. (0->36)
t74 = pkin(5) * m(6);
t73 = m(4) + m(10);
t72 = m(11) + m(5);
t63 = qJ(1) + qJ(2);
t71 = -m(10) * pkin(6) - mrSges(4,1);
t60 = qJ(3) + t63;
t59 = qJ(4) + t63;
t70 = m(11) * pkin(4) + mrSges(5,1);
t69 = t73 * pkin(2) + t72 * pkin(3) + mrSges(3,1);
t68 = mrSges(2,1) + (m(3) + m(7) + m(9) + t72 + t73) * pkin(1);
t67 = cos(qJ(1));
t66 = cos(qJ(7));
t65 = sin(qJ(1));
t64 = sin(qJ(7));
t62 = qJ(1) + qJ(8);
t61 = pkin(8) + qJ(5);
t58 = qJ(6) + t63;
t57 = cos(t63);
t56 = cos(t62);
t55 = sin(t63);
t54 = sin(t62);
t53 = cos(t61);
t52 = sin(t61);
t51 = qJ(9) + t60;
t50 = qJ(10) + t59;
t49 = cos(t60);
t48 = cos(t59);
t47 = cos(t58);
t46 = sin(t60);
t45 = sin(t59);
t44 = sin(t58);
t43 = cos(t51);
t42 = sin(t51);
t41 = cos(t50);
t40 = sin(t50);
t1 = (-mrSges(10,3) - mrSges(9,3) - mrSges(8,3) - mrSges(7,3) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) - mrSges(3,3) - mrSges(11,3) - mrSges(2,3) - mrSges(1,3)) * g(3) + (t67 * mrSges(2,2) - t56 * mrSges(9,2) + t57 * mrSges(3,2) - t64 * mrSges(8,2) + t66 * mrSges(8,1) - t40 * mrSges(11,1) - t41 * mrSges(11,2) + t42 * mrSges(10,1) + t43 * mrSges(10,2) - t44 * mrSges(7,1) - t47 * mrSges(7,2) + t48 * mrSges(5,2) - t49 * mrSges(4,2) - t52 * mrSges(6,1) - t53 * mrSges(6,2) - t54 * mrSges(9,1) - sin(pkin(8)) * t74 - mrSges(1,2) + t71 * t46 + t70 * t45 + t69 * t55 + t68 * t65) * g(2) + (t54 * mrSges(9,2) - t55 * mrSges(3,2) - t56 * mrSges(9,1) - t64 * mrSges(8,1) - t65 * mrSges(2,2) - t66 * mrSges(8,2) + t40 * mrSges(11,2) - t41 * mrSges(11,1) - t42 * mrSges(10,2) + t43 * mrSges(10,1) + t44 * mrSges(7,2) - t45 * mrSges(5,2) + t46 * mrSges(4,2) - t47 * mrSges(7,1) + t52 * mrSges(6,2) - t53 * mrSges(6,1) - m(8) * pkin(7) - cos(pkin(8)) * t74 - mrSges(1,1) + t71 * t49 + t70 * t48 + t69 * t57 + t68 * t67) * g(1);
U = t1;
