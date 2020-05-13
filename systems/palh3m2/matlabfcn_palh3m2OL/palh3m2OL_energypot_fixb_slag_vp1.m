% Calculate potential energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2OL_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:23
% EndTime: 2020-05-07 04:33:24
% DurationCPUTime: 0.49s
% Computational Cost: add. (205->97), mult. (185->121), div. (0->0), fcn. (153->18), ass. (0->42)
t96 = pkin(10) + rSges(6,3);
t73 = qJ(2) + qJ(3);
t69 = qJ(4) + t73;
t62 = cos(t69);
t95 = t62 * pkin(8);
t75 = sin(qJ(5));
t81 = cos(qJ(1));
t93 = t75 * t81;
t77 = sin(qJ(1));
t92 = t77 * t75;
t79 = cos(qJ(5));
t91 = t77 * t79;
t90 = t79 * t81;
t76 = sin(qJ(2));
t89 = t76 * pkin(1) + pkin(11);
t80 = cos(qJ(2));
t60 = t80 * pkin(1) + pkin(12);
t72 = qJ(2) + qJ(7);
t68 = pkin(15) - t72;
t65 = sin(t73);
t88 = -pkin(4) * t65 + t89;
t67 = cos(t73);
t87 = -rSges(4,1) * t67 + rSges(4,2) * t65;
t61 = sin(t69);
t86 = -rSges(5,1) * t62 + rSges(5,2) * t61;
t64 = sin(t72);
t66 = cos(t72);
t85 = rSges(8,1) * t66 - rSges(8,2) * t64;
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t84 = rSges(7,1) * t78 - rSges(7,2) * t74 - pkin(6);
t83 = rSges(3,1) * t80 - rSges(3,2) * t76 + pkin(12);
t63 = -qJ(8) + t68;
t58 = sin(t63);
t59 = cos(t63);
t82 = -rSges(9,1) * t59 - rSges(9,2) * t58 + pkin(3) * cos(t68) + t60;
t57 = t81 * t60;
t56 = t77 * t60;
t55 = -pkin(4) * t67 + t60;
t53 = t81 * t55;
t52 = t77 * t55;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + rSges(2,2) * t81) + g(3) * (pkin(11) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t76 + rSges(3,2) * t80 + pkin(11)) + (-g(2) * rSges(3,3) + g(1) * t83) * t81 + (g(1) * rSges(3,3) + g(2) * t83) * t77) - m(4) * (g(1) * (t77 * rSges(4,3) + t81 * t87 + t57) + g(2) * (-rSges(4,3) * t81 + t77 * t87 + t56) + g(3) * (-rSges(4,1) * t65 - rSges(4,2) * t67 + t89)) - m(5) * (g(1) * (t77 * rSges(5,3) + t81 * t86 + t53) + g(2) * (-rSges(5,3) * t81 + t77 * t86 + t52) + g(3) * (-rSges(5,1) * t61 - rSges(5,2) * t62 + t88)) - m(6) * (g(1) * (-t81 * t95 + t53 + (-t62 * t90 + t92) * rSges(6,1) + (t62 * t93 + t91) * rSges(6,2)) + g(2) * (-t77 * t95 + t52 + (-t62 * t91 - t93) * rSges(6,1) + (t62 * t92 - t90) * rSges(6,2)) + g(3) * (t96 * t62 + t88) + (g(3) * (-rSges(6,1) * t79 + rSges(6,2) * t75 - pkin(8)) - (g(1) * t81 + g(2) * t77) * t96) * t61) - m(7) * (g(3) * (rSges(7,1) * t74 + rSges(7,2) * t78 + pkin(11) + pkin(13)) + (-g(2) * rSges(7,3) + g(1) * t84) * t81 + (g(1) * rSges(7,3) + g(2) * t84) * t77) - m(8) * (g(1) * (t77 * rSges(8,3) + t81 * t85 + t57) + g(2) * (-rSges(8,3) * t81 + t77 * t85 + t56) + g(3) * (rSges(8,1) * t64 + rSges(8,2) * t66 + t89)) - m(9) * (g(3) * (-pkin(3) * sin(t68) + t58 * rSges(9,1) - t59 * rSges(9,2) + t89) + (-g(2) * rSges(9,3) + g(1) * t82) * t81 + (g(1) * rSges(9,3) + g(2) * t82) * t77);
U = t1;
