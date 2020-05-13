% Calculate potential energy for
% palh1m2OL
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:21
% EndTime: 2020-05-02 21:17:22
% DurationCPUTime: 0.64s
% Computational Cost: add. (234->108), mult. (249->135), div. (0->0), fcn. (148->24), ass. (0->47)
t101 = cos(qJ(1));
t96 = sin(qJ(1));
t106 = t101 * g(1) + t96 * g(2);
t116 = pkin(15) * g(1);
t115 = pkin(15) * g(2);
t114 = g(3) * pkin(13);
t93 = sin(qJ(5));
t98 = cos(qJ(5));
t112 = t98 * rSges(6,1) - t93 * rSges(6,2) + pkin(9);
t111 = -rSges(6,3) - pkin(11);
t100 = cos(qJ(2));
t109 = t100 * pkin(1) + pkin(13);
t91 = qJ(2) + qJ(3);
t90 = qJ(2) + qJ(7);
t89 = qJ(2) + qJ(8);
t108 = pkin(5) * sin(t91) + t109;
t95 = sin(qJ(2));
t78 = -t95 * pkin(1) + pkin(15);
t107 = -rSges(2,2) * m(2) + rSges(3,3) * m(3) + rSges(4,3) * m(4);
t81 = pkin(19) - t90;
t87 = qJ(4) + t91;
t79 = sin(t87);
t80 = cos(t87);
t105 = rSges(5,1) * t80 - rSges(5,2) * t79;
t83 = sin(t90);
t85 = cos(t90);
t104 = -rSges(8,1) * t83 - rSges(8,2) * t85 + t78;
t74 = -qJ(10) + t81;
t72 = sin(t74);
t73 = cos(t74);
t103 = -rSges(11,1) * t72 + rSges(11,2) * t73 + pkin(4) * sin(t81) + t78;
t94 = sin(qJ(3));
t99 = cos(qJ(3));
t63 = rSges(3,1) * m(3) + (rSges(4,1) * t94 + rSges(4,2) * t99 + pkin(1)) * m(4);
t66 = -rSges(3,2) * m(3) + (rSges(4,1) * t99 - rSges(4,2) * t94) * m(4);
t102 = -rSges(2,1) * m(2) - t66 * t100 + t63 * t95 + (-m(3) - m(4)) * pkin(15);
t97 = cos(qJ(6));
t92 = sin(qJ(6));
t86 = qJ(9) + t89;
t84 = cos(t89);
t82 = sin(t89);
t70 = t93 * rSges(6,1) + t98 * rSges(6,2);
t69 = t97 * rSges(7,1) - t92 * rSges(7,2) - pkin(14);
t68 = pkin(5) * cos(t91) + t78;
t65 = t101 * t68;
t64 = t96 * t68;
t1 = (t102 * g(1) + t107 * g(2)) * t101 + (-g(1) * t107 + t102 * g(2)) * t96 - g(3) * t66 * t95 - g(3) * t63 * t100 - m(4) * t114 - (pkin(13) * m(3) + (pkin(13) + rSges(2,3)) * m(2)) * g(3) - m(5) * (g(1) * (t96 * rSges(5,3) + t101 * t105 + t65) + g(2) * (-t101 * rSges(5,3) + t105 * t96 + t64) + g(3) * (t79 * rSges(5,1) + t80 * rSges(5,2) + t108)) - m(6) * (g(1) * (t96 * t70 + t65) + g(2) * (-t101 * t70 + t64) + g(3) * t108 + (g(3) * t111 + t106 * t112) * t80 + (g(3) * t112 - t106 * t111) * t79) - m(7) * (g(1) * (t96 * rSges(7,3) + t69 * t101) + g(2) * (-t101 * rSges(7,3) + t69 * t96) + g(3) * (t92 * rSges(7,1) + t97 * rSges(7,2) + pkin(13) - pkin(16))) - m(8) * (g(3) * (t85 * rSges(8,1) - t83 * rSges(8,2) + t109) + (g(1) * rSges(8,3) + g(2) * t104) * t96 + (-g(2) * rSges(8,3) + g(1) * t104) * t101) - ((g(3) * rSges(9,1) - rSges(9,2) * t106) * t84 + (-rSges(9,1) * t106 - g(3) * rSges(9,2)) * t82 + (-rSges(9,3) * g(2) + t116) * t101 + (rSges(9,3) * g(1) + t115) * t96 + t114) * m(9) + m(10) * ((g(3) * rSges(10,1) - rSges(10,2) * t106) * cos(t86) + (-t106 * rSges(10,1) - g(3) * rSges(10,2)) * sin(t86) + (rSges(10,3) * g(2) - t116) * t101 + (-rSges(10,3) * g(1) - t115) * t96 - t114 + (-g(3) * t84 + t106 * t82) * pkin(2)) - m(11) * (g(3) * (pkin(4) * cos(t81) - t73 * rSges(11,1) - t72 * rSges(11,2) + t109) + (g(1) * rSges(11,3) + g(2) * t103) * t96 + (-g(2) * rSges(11,3) + g(1) * t103) * t101) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1);
U = t1;
