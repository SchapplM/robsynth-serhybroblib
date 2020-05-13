% Calculate potential energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1OL_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:18
% EndTime: 2020-04-15 19:28:19
% DurationCPUTime: 0.64s
% Computational Cost: add. (251->121), mult. (224->149), div. (0->0), fcn. (184->22), ass. (0->52)
t120 = pkin(11) + rSges(6,3);
t95 = qJ(2) + qJ(3);
t89 = qJ(4) + t95;
t80 = cos(t89);
t119 = t80 * pkin(9);
t97 = sin(qJ(5));
t99 = sin(qJ(1));
t117 = t99 * t97;
t102 = cos(qJ(2));
t116 = t102 * pkin(1) + pkin(13);
t103 = cos(qJ(1));
t115 = t103 * t97;
t101 = cos(qJ(5));
t114 = t99 * t101;
t113 = t101 * t103;
t94 = qJ(2) + qJ(7);
t93 = qJ(2) + qJ(8);
t84 = sin(t95);
t112 = pkin(5) * t84 + t116;
t98 = sin(qJ(2));
t76 = -t98 * pkin(1) + pkin(15);
t81 = pkin(19) - t94;
t87 = cos(t95);
t111 = rSges(4,1) * t87 - rSges(4,2) * t84;
t78 = sin(t89);
t110 = rSges(5,1) * t80 - rSges(5,2) * t78;
t83 = sin(t94);
t86 = cos(t94);
t109 = -rSges(8,1) * t83 - rSges(8,2) * t86;
t82 = sin(t93);
t85 = cos(t93);
t108 = -rSges(9,1) * t82 - rSges(9,2) * t85;
t107 = -rSges(3,1) * t98 - rSges(3,2) * t102;
t88 = qJ(9) + t93;
t77 = sin(t88);
t79 = cos(t88);
t106 = -pkin(2) * t82 + rSges(10,1) * t77 + rSges(10,2) * t79 + pkin(15);
t100 = cos(qJ(6));
t96 = sin(qJ(6));
t105 = rSges(7,1) * t100 - rSges(7,2) * t96 - pkin(14);
t74 = -qJ(10) + t81;
t72 = sin(t74);
t73 = cos(t74);
t104 = -rSges(11,1) * t72 + rSges(11,2) * t73 + pkin(4) * sin(t81) + t76;
t91 = t103 * pkin(15);
t90 = t99 * pkin(15);
t70 = t103 * t76;
t69 = t99 * t76;
t68 = pkin(5) * t87 + t76;
t66 = t103 * t68;
t65 = t99 * t68;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t103 - t99 * rSges(2,2)) + g(2) * (t99 * rSges(2,1) + rSges(2,2) * t103) + g(3) * (pkin(13) + rSges(2,3))) - m(3) * (g(1) * (t99 * rSges(3,3) + t103 * t107 + t91) + g(2) * (-rSges(3,3) * t103 + t107 * t99 + t90) + g(3) * (rSges(3,1) * t102 - rSges(3,2) * t98 + pkin(13))) - m(4) * (g(1) * (t99 * rSges(4,3) + t103 * t111 + t70) + g(2) * (-rSges(4,3) * t103 + t111 * t99 + t69) + g(3) * (rSges(4,1) * t84 + rSges(4,2) * t87 + t116)) - m(5) * (g(1) * (t99 * rSges(5,3) + t103 * t110 + t66) + g(2) * (-rSges(5,3) * t103 + t110 * t99 + t65) + g(3) * (rSges(5,1) * t78 + rSges(5,2) * t80 + t112)) - m(6) * (g(1) * (t103 * t119 + t66 + (t113 * t80 + t117) * rSges(6,1) + (-t115 * t80 + t114) * rSges(6,2)) + g(2) * (t99 * t119 + t65 + (t114 * t80 - t115) * rSges(6,1) + (-t117 * t80 - t113) * rSges(6,2)) + g(3) * (-t120 * t80 + t112) + (g(3) * (rSges(6,1) * t101 - rSges(6,2) * t97 + pkin(9)) + (g(1) * t103 + g(2) * t99) * t120) * t78) - m(7) * (g(3) * (rSges(7,1) * t96 + rSges(7,2) * t100 + pkin(13) - pkin(16)) + (g(1) * rSges(7,3) + g(2) * t105) * t99 + (-g(2) * rSges(7,3) + g(1) * t105) * t103) - m(8) * (g(1) * (t99 * rSges(8,3) + t103 * t109 + t70) + g(2) * (-rSges(8,3) * t103 + t109 * t99 + t69) + g(3) * (rSges(8,1) * t86 - rSges(8,2) * t83 + t116)) - m(9) * (g(1) * (t99 * rSges(9,3) + t103 * t108 + t91) + g(2) * (-rSges(9,3) * t103 + t108 * t99 + t90) + g(3) * (rSges(9,1) * t85 - rSges(9,2) * t82 + pkin(13))) - m(10) * (g(3) * (pkin(2) * t85 - rSges(10,1) * t79 + rSges(10,2) * t77 + pkin(13)) + (g(1) * rSges(10,3) + g(2) * t106) * t99 + (-g(2) * rSges(10,3) + g(1) * t106) * t103) - m(11) * (g(3) * (pkin(4) * cos(t81) - t73 * rSges(11,1) - t72 * rSges(11,2) + t116) + (g(1) * rSges(11,3) + g(2) * t104) * t99 + (-g(2) * rSges(11,3) + g(1) * t104) * t103);
U = t1;
