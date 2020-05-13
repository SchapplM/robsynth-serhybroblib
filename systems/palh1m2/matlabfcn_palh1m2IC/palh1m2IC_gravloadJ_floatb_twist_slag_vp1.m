% Calculate Gravitation load on the joints for
% palh1m2IC
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:46:48
% EndTime: 2020-05-02 23:46:50
% DurationCPUTime: 1.93s
% Computational Cost: add. (830->203), mult. (733->242), div. (19->10), fcn. (506->45), ass. (0->106)
t76 = sin(qJ(5));
t81 = cos(qJ(5));
t35 = rSges(6,1) * t81 - rSges(6,2) * t76;
t172 = pkin(9) + t35;
t169 = pkin(11) + rSges(6,3);
t77 = sin(qJ(3));
t82 = cos(qJ(3));
t23 = (rSges(4,1) * t82 - rSges(4,2) * t77) * m(4);
t13 = -rSges(3,2) * m(3) + t23;
t110 = -rSges(4,1) * t77 - rSges(4,2) * t82;
t7 = rSges(3,1) * m(3) + (pkin(1) - t110) * m(4);
t78 = sin(qJ(2));
t83 = cos(qJ(2));
t171 = -t13 * t83 + t7 * t78;
t79 = sin(qJ(1));
t84 = cos(qJ(1));
t115 = g(1) * t84 + g(2) * t79;
t74 = qJ(2) + qJ(3);
t67 = qJ(4) + t74;
t56 = sin(t67);
t58 = cos(t67);
t170 = t169 * t56 + t172 * t58;
t130 = -qJ(7) + pkin(19);
t119 = -qJ(10) + t130;
t50 = -qJ(2) + t119;
t41 = sin(t50);
t42 = cos(t50);
t167 = rSges(11,1) * t42 + rSges(11,2) * t41;
t166 = g(2) * t84;
t131 = qJ(4) + pkin(18);
t111 = pkin(19) + qJ(3) + t131;
t132 = qJ(6) - qJ(2);
t100 = t111 + t132;
t92 = -(2 * qJ(7)) - pkin(20) + t100;
t101 = t111 - t132;
t93 = pkin(20) + t101;
t165 = cos(qJ(10) - t92) + cos(qJ(10) - t93);
t120 = -rSges(11,1) * t41 + rSges(11,2) * t42;
t60 = -qJ(2) + t130;
t117 = t120 + pkin(4) * sin(t60);
t164 = t58 * rSges(5,1) - rSges(5,2) * t56;
t51 = pkin(5) * cos(t74);
t159 = pkin(1) * t78;
t54 = pkin(15) - t159;
t17 = t51 + t54;
t160 = t17 * t166;
t158 = pkin(1) * t83;
t157 = pkin(4) * cos(t60);
t156 = pkin(5) * sin(t74);
t155 = pkin(15) * g(1);
t154 = pkin(15) * g(2);
t108 = -rSges(5,1) * t56 - rSges(5,2) * t58;
t102 = t108 * t84;
t103 = t108 * t79;
t140 = t58 * t84;
t136 = t169 * t140;
t141 = t58 * t79;
t137 = t169 * t141;
t90 = t115 * t172 * t56;
t151 = (-m(5) * (g(1) * t102 + g(2) * t103 + g(3) * t164) - m(6) * (g(1) * t136 + g(2) * t137 + g(3) * t170 - t90)) / pkin(10);
t73 = qJ(2) + qJ(7);
t65 = cos(t73);
t146 = rSges(8,1) * t65;
t62 = sin(t73);
t143 = rSges(8,2) * t62;
t139 = t167 * t79;
t138 = t167 * t84;
t72 = qJ(2) + qJ(8);
t129 = qJ(7) + pkin(20);
t127 = m(11) * (g(1) * t138 + g(2) * t139 + g(3) * t120) / pkin(8);
t66 = qJ(9) + t72;
t55 = sin(t66);
t57 = cos(t66);
t124 = -(g(3) * rSges(10,1) - t115 * rSges(10,2)) * t55 + (-t115 * rSges(10,1) - g(3) * rSges(10,2)) * t57;
t61 = sin(t72);
t64 = cos(t72);
t126 = -(-(g(3) * rSges(9,1) - t115 * rSges(9,2)) * t61 + (-t115 * rSges(9,1) - g(3) * rSges(9,2)) * t64) * m(9) + m(10) * ((g(3) * t61 + t115 * t64) * pkin(2) + t124);
t125 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t75 = sin(qJ(6));
t80 = cos(qJ(6));
t122 = t80 * rSges(7,1) - rSges(7,2) * t75;
t121 = -qJ(8) + pkin(17) + qJ(3);
t118 = t164 + t51;
t116 = t129 - t132;
t114 = g(1) * t79 - t166;
t113 = -t146 - t158;
t34 = -t76 * rSges(6,1) - t81 * rSges(6,2);
t107 = rSges(8,1) * t62 + rSges(8,2) * t65;
t105 = t107 - t54;
t104 = -t54 - t117;
t97 = -m(2) * rSges(2,1) + (-m(3) - m(4)) * pkin(15) + t171;
t96 = t51 + t170;
t95 = -qJ(7) + t101;
t94 = -qJ(7) + t100;
t88 = 0.1e1 / pkin(3);
t46 = cos(t116);
t45 = cos(qJ(9) - t121);
t38 = t84 * t143;
t37 = t79 * t143;
t26 = -pkin(14) + t122;
t24 = -t156 - t158;
t22 = t110 * m(4);
t16 = -t157 - t158;
t11 = t84 * t24;
t10 = t79 * t24;
t1 = [-(t97 * g(1) + t125 * g(2)) * t79 + (-t125 * g(1) + t97 * g(2)) * t84 - m(5) * (t160 + (g(1) * rSges(5,3) + g(2) * t164) * t84 + (g(1) * (-t164 - t17) + g(2) * rSges(5,3)) * t79) - m(6) * (t160 + (-g(1) * t34 + g(2) * t170) * t84 + (g(1) * (-t17 - t170) - g(2) * t34) * t79) - m(7) * (g(1) * (rSges(7,3) * t84 - t26 * t79) + g(2) * (rSges(7,3) * t79 + t26 * t84)) - m(8) * ((g(1) * rSges(8,3) - g(2) * t105) * t84 + (g(2) * rSges(8,3) + g(1) * t105) * t79) - (-(-rSges(9,3) * g(2) + t155) * t79 + (rSges(9,3) * g(1) + t154) * t84 + (rSges(9,1) * t61 + rSges(9,2) * t64) * t114) * m(9) + m(10) * (-(rSges(10,3) * g(2) - t155) * t79 + (-rSges(10,3) * g(1) - t154) * t84 + (-pkin(2) * t61 + t55 * rSges(10,1) + t57 * rSges(10,2)) * t114) - m(11) * ((g(1) * rSges(11,3) - g(2) * t104) * t84 + (g(2) * rSges(11,3) + g(1) * t104) * t79); -m(5) * (g(1) * (t11 + t102) + g(2) * (t10 + t103)) - m(6) * (g(1) * (t11 + t136) + g(2) * (t10 + t137) - t90) - m(8) * (g(1) * (t113 * t84 + t38) + g(2) * (t113 * t79 + t37)) - m(11) * (g(1) * (t16 * t84 + t138) + g(2) * (t16 * t79 + t139)) + t126 + (-m(5) * (t118 - t159) - m(6) * (t96 - t159) - m(8) * (-t107 - t159) - m(11) * (t117 - t159) + t171) * g(3) + (-pkin(1) * sin(t129) / pkin(7) * m(7) * (g(3) * t122 + t115 * (-rSges(7,1) * t75 - rSges(7,2) * t80)) + (-pkin(3) * t46 - pkin(1) * cos(-t132)) * t88 * (-m(8) * (g(1) * (-t84 * t146 + t38) + g(2) * (-t79 * t146 + t37) - g(3) * t107) - m(11) * (g(1) * (-t84 * t157 + t138) + g(2) * (-t79 * t157 + t139) + g(3) * t117))) / t46 + (-(pkin(1) * (sin(qJ(10) - t132) + sin(qJ(10) + t132)) + (-sin(-qJ(10) + t116) + sin(qJ(10) + t116)) * pkin(3)) * pkin(4) * t151 - ((-(cos(t95) + cos(t94)) * pkin(1) + (-cos(t92) - cos(t93)) * pkin(3)) * pkin(4) + ((cos(-qJ(10) + t95) + cos(-qJ(10) + t94)) * pkin(1) + t165 * pkin(3)) * pkin(8)) * t127) / t165 * t88 + t115 * (t13 * t78 + t7 * t83); -m(6) * (g(1) * (-t84 * t156 + t136) + g(2) * (-t79 * t156 + t137) - t90) - pkin(10) * t151 + (t45 * t126 + (-pkin(12) * t45 + pkin(2) * cos(t121)) / pkin(12) * m(10) * t124) * pkin(6) / sin(qJ(9)) / pkin(2) + (-m(5) * t118 - m(6) * t96 - t22 * t78 - t23 * t83) * g(3) + (-m(5) * (t108 - t156) - t22 * t83 + t23 * t78) * t115 + (-cos(qJ(3) + t119) * t151 - sin(t131) * t127) / cos(-qJ(7) - qJ(10) + t111) * pkin(5); -m(6) * (g(1) * (t34 * t140 + t35 * t79) + g(2) * (t34 * t141 - t35 * t84) + g(3) * t34 * t56);];
taug = t1(:);
