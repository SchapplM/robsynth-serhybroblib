% Calculate potential energy for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m2TE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energypot_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:48
% EndTime: 2020-05-07 01:41:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (237->89), mult. (402->126), div. (0->0), fcn. (237->22), ass. (0->47)
t89 = sin(qJ(4));
t95 = cos(qJ(4));
t115 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t89 + rSges(6,2) * t95) * m(6) + rSges(7,3) * m(7);
t114 = m(5) + m(6);
t113 = m(4) * rSges(4,2);
t91 = sin(qJ(2));
t112 = t91 * pkin(1) + pkin(11);
t97 = cos(qJ(2));
t111 = pkin(1) * t97 + pkin(12);
t110 = -m(3) - m(4) - m(5);
t92 = sin(qJ(1));
t98 = cos(qJ(1));
t109 = g(1) * t98 + g(2) * t92;
t108 = (-rSges(6,3) - pkin(10)) * m(6) + m(5) * rSges(5,2);
t85 = sin(pkin(18));
t87 = cos(pkin(18));
t93 = sin(pkin(15));
t99 = cos(pkin(15));
t106 = t85 * t99 + t87 * t93;
t105 = -t85 * t93 + t87 * t99;
t82 = qJ(3) + qJ(2);
t77 = sin(t82);
t78 = cos(t82);
t86 = sin(pkin(17));
t88 = cos(pkin(17));
t104 = -rSges(9,1) * t78 + rSges(9,2) * t77 + (-t105 * t88 + t106 * t86) * pkin(3) + t111;
t103 = m(6) * (rSges(6,1) * t95 - rSges(6,2) * t89);
t100 = cos(pkin(14));
t94 = sin(pkin(14));
t73 = rSges(7,1) * t94 - rSges(7,2) * t100;
t74 = rSges(7,1) * t100 + rSges(7,2) * t94;
t75 = pkin(4) * t114 + m(4) * rSges(4,1);
t90 = sin(qJ(3));
t96 = cos(qJ(3));
t63 = t90 * t113 + rSges(3,1) * m(3) - t75 * t96 + (t73 * t93 + t74 * t99) * m(7) + (m(4) + t114) * pkin(1);
t64 = t96 * t113 - rSges(3,2) * m(3) + t75 * t90 + (t73 * t99 - t74 * t93) * m(7);
t102 = m(7) * pkin(6) - m(2) * rSges(2,1) - t63 * t97 - t64 * t91 + (-m(6) + t110) * pkin(12);
t84 = cos(pkin(16));
t83 = sin(pkin(16));
t81 = pkin(17) + pkin(18);
t76 = pkin(8) * m(6) + m(5) * rSges(5,1);
t68 = -t108 * t84 + t83 * t76;
t67 = t108 * t83 + t76 * t84;
t66 = (-rSges(8,1) * t99 - rSges(8,2) * t93) * t87 + t85 * (rSges(8,1) * t93 - rSges(8,2) * t99) + t111;
t62 = t67 * t99 - t68 * t93 + (-t83 * t93 + t84 * t99) * t103;
t61 = t67 * t93 + t68 * t99 + (t83 * t99 + t84 * t93) * t103;
t1 = (-g(3) * t61 + t109 * t62) * cos(t81) + (-t62 * g(3) - t109 * t61) * sin(t81) + (t102 * g(1) + t115 * g(2)) * t98 + (-t115 * g(1) + t102 * g(2)) * t92 + g(3) * t64 * t97 - g(3) * t63 * t91 - m(6) * pkin(11) * g(3) - ((m(2) + m(7) - t110) * pkin(11) + m(2) * rSges(2,3) + pkin(13) * m(7)) * g(3) - m(8) * (g(1) * (t92 * rSges(8,3) + t66 * t98) + g(2) * (-t98 * rSges(8,3) + t66 * t92) + g(3) * ((t85 * rSges(8,1) - rSges(8,2) * t87) * t99 + t93 * (rSges(8,1) * t87 + rSges(8,2) * t85) + t112)) - m(9) * ((-g(2) * rSges(9,3) + g(1) * t104) * t98 + (g(1) * rSges(9,3) + g(2) * t104) * t92 + (-t77 * rSges(9,1) - t78 * rSges(9,2) + t112 + (t105 * t86 + t106 * t88) * pkin(3)) * g(3)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1);
U = t1;
