% Calculate potential energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:51
% EndTime: 2020-05-02 20:57:53
% DurationCPUTime: 0.85s
% Computational Cost: add. (292->95), mult. (366->119), div. (0->0), fcn. (238->40), ass. (0->50)
t118 = m(5) + m(6);
t109 = m(4) + m(8) + m(11) + t118;
t104 = m(9) + m(3) + m(10) + t109;
t90 = sin(qJ(4));
t96 = cos(qJ(4));
t120 = -m(2) * rSges(2,2) + m(11) * rSges(11,3) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + (rSges(6,1) * t90 + rSges(6,2) * t96) * m(6) + rSges(7,3) * m(7) + m(8) * rSges(8,3) + rSges(9,3) * m(9) + rSges(10,3) * m(10);
t101 = cos(pkin(17));
t100 = cos(pkin(18));
t94 = sin(pkin(18));
t105 = rSges(7,1) * t100 + rSges(7,2) * t94;
t106 = rSges(7,1) * t94 - rSges(7,2) * t100;
t116 = m(4) * rSges(4,2);
t79 = t118 * pkin(5) + m(4) * rSges(4,1);
t91 = sin(qJ(3));
t95 = sin(pkin(17));
t97 = cos(qJ(3));
t63 = t97 * t116 + rSges(3,1) * m(3) + t79 * t91 + (t105 * t101 + t106 * t95) * m(7) + t109 * pkin(1);
t65 = t91 * t116 + rSges(3,2) * m(3) - t79 * t97 + (-t106 * t101 + t105 * t95) * m(7);
t92 = sin(qJ(2));
t98 = cos(qJ(2));
t119 = -m(2) * rSges(2,1) + pkin(14) * m(7) - t104 * pkin(15) + t63 * t92 + t65 * t98;
t115 = m(9) * rSges(9,2);
t87 = sin(pkin(19));
t89 = cos(pkin(19));
t74 = t87 * t97 + t89 * t91;
t75 = t87 * t91 - t89 * t97;
t68 = qJ(2) + atan2(t75, t74);
t85 = pkin(18) - pkin(22);
t84 = -qJ(2) + t85;
t113 = pkin(21) - atan2(cos(t84), -sin(t84));
t111 = -qJ(2) + t113;
t110 = -pkin(21) + t85;
t93 = sin(qJ(1));
t99 = cos(qJ(1));
t78 = g(1) * t99 + g(2) * t93;
t88 = cos(pkin(22));
t86 = sin(pkin(22));
t83 = pkin(2) * m(10) + m(9) * rSges(9,1);
t82 = -pkin(20) + t110;
t81 = (-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2);
t80 = -qJ(2) - qJ(3) + t110;
t77 = t100 * t86 - t88 * t94;
t76 = t100 * t88 + t86 * t94;
t72 = m(5) * rSges(5,1) + (rSges(6,1) * t96 - rSges(6,2) * t90 + pkin(9)) * m(6);
t71 = atan2(-sin(t80), cos(t80));
t67 = -t71 + t113;
t66 = -t71 + t111;
t64 = atan2(t75, -t74) + t68;
t62 = pkin(21) - atan2(t76 * t98 - t77 * t92, t76 * t92 + t77 * t98);
t1 = (g(3) * t115 + t78 * t83) * sin(t68) + (-g(3) * t83 + t78 * t115) * cos(t68) + (g(3) * t81 + t78 * t72) * cos(t82) + (g(3) * t72 - t78 * t81) * sin(t82) + (-t120 * g(1) + t119 * g(2)) * t93 + (t119 * g(1) + t120 * g(2)) * t99 - g(3) * (m(2) * rSges(2,3) - pkin(16) * m(7) + (m(7) + m(2) + t104) * pkin(13)) + g(3) * t65 * t92 - g(3) * t63 * t98 + (-(-g(3) * rSges(10,1) + t78 * rSges(10,2)) * cos(t64) - (t78 * rSges(10,1) + g(3) * rSges(10,2)) * sin(t64)) * m(10) + ((t78 * rSges(8,1) + g(3) * rSges(8,2)) * cos(t85) - (-g(3) * rSges(8,1) + t78 * rSges(8,2)) * sin(t85)) * m(8) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - g(3) * rSges(1,3)) * m(1) + ((-pkin(4) * sin(t111) - rSges(11,2) * cos(t66) + rSges(11,1) * sin(t66)) * t78 + ((rSges(11,1) * t98 - rSges(11,2) * t92) * cos(t67) + (rSges(11,1) * t92 + rSges(11,2) * t98) * sin(t67) + (-cos(t62) * t98 - sin(t62) * t92) * pkin(4)) * g(3)) * m(11);
U = t1;
