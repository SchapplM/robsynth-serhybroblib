% Calculate Gravitation load on the joints for
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2DE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:29
% EndTime: 2020-05-02 20:58:34
% DurationCPUTime: 0.84s
% Computational Cost: add. (761->154), mult. (430->154), div. (0->0), fcn. (352->72), ass. (0->100)
t100 = m(5) + m(6);
t113 = m(8) + m(4) + t100;
t147 = m(11) + t113;
t140 = pkin(5) * t100 + mrSges(4,1);
t146 = g(1) * t140;
t78 = pkin(18) - pkin(22);
t65 = -qJ(2) + t78;
t112 = pkin(21) - qJ(2) - atan2(cos(t65), -sin(t65));
t109 = -pkin(21) + t78;
t79 = qJ(2) + qJ(3);
t45 = t109 - t79;
t21 = atan2(-sin(t45), cos(t45));
t15 = -t21 + t112;
t77 = qJ(3) + pkin(19);
t64 = qJ(2) + t77;
t76 = pkin(2) * m(10) + mrSges(9,1);
t145 = -t140 * cos(t79) - mrSges(11,2) * cos(t15) + mrSges(4,2) * sin(t79) + mrSges(9,2) * sin(t64) - t76 * cos(t64);
t118 = mrSges(11,2) * g(2);
t119 = mrSges(11,2) * g(1);
t120 = mrSges(11,1) * g(2);
t121 = mrSges(11,1) * g(1);
t23 = -qJ(1) + t112;
t14 = -t21 + t23;
t144 = -(t118 + t121) * cos(t14) / 0.2e1 + (-t119 + t120) * sin(t14) / 0.2e1;
t84 = sin(qJ(4));
t90 = cos(qJ(4));
t107 = mrSges(6,1) * t84 + mrSges(6,2) * t90;
t143 = -mrSges(11,3) - t107 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3);
t88 = sin(pkin(18));
t94 = cos(pkin(18));
t105 = mrSges(7,1) * t88 - mrSges(7,2) * t94;
t106 = mrSges(7,1) * t94 + mrSges(7,2) * t88;
t85 = sin(qJ(3));
t86 = sin(qJ(2));
t89 = sin(pkin(17));
t91 = cos(qJ(3));
t92 = cos(qJ(2));
t95 = cos(pkin(17));
t142 = pkin(14) * m(7) + (pkin(1) * t113 + mrSges(4,2) * t91 + t105 * t89 + t106 * t95 + t140 * t85 + mrSges(3,1)) * t86 - (-t85 * mrSges(4,2) + t105 * t95 - t106 * t89 + t140 * t91 - mrSges(3,2)) * t92 - mrSges(2,1) + (-m(9) - m(10) - m(3) - t147) * pkin(15);
t22 = qJ(1) + t112;
t13 = -t21 + t22;
t141 = (-t118 + t121) * cos(t13) / 0.2e1 - (-t119 - t120) * sin(t13) / 0.2e1;
t80 = qJ(1) + qJ(2);
t139 = sin(t80) / 0.2e1;
t137 = m(11) * g(1);
t136 = m(11) * g(2);
t135 = mrSges(7,1) * g(2);
t134 = mrSges(10,1) * g(1);
t133 = mrSges(10,1) * g(2);
t132 = mrSges(3,2) * g(1);
t131 = mrSges(3,2) * g(2);
t130 = mrSges(7,2) * g(1);
t129 = mrSges(7,2) * g(2);
t128 = mrSges(9,2) * g(2);
t127 = mrSges(10,2) * g(1);
t126 = mrSges(10,2) * g(2);
t125 = g(2) * t140;
t124 = g(2) * t76;
t82 = sin(pkin(19));
t83 = cos(pkin(19));
t33 = t82 * t91 + t83 * t85;
t34 = t82 * t85 - t83 * t91;
t16 = qJ(2) + atan2(t34, t33);
t116 = -t137 / 0.2e1;
t115 = t137 / 0.2e1;
t114 = -t136 / 0.2e1;
t62 = sin(t77);
t63 = cos(t77);
t20 = qJ(2) + atan2(-t63, t62) + atan2(-t63, -t62);
t72 = qJ(2) + pkin(17) - pkin(18);
t111 = mrSges(6,1) * t90 - mrSges(6,2) * t84;
t87 = sin(qJ(1));
t93 = cos(qJ(1));
t108 = g(1) * t87 - g(2) * t93;
t53 = qJ(1) + t64;
t54 = -qJ(1) + t64;
t57 = t76 * g(1);
t73 = qJ(1) + t79;
t74 = qJ(1) - t79;
t96 = mrSges(9,2) * g(1);
t97 = mrSges(4,2) * g(2);
t98 = mrSges(4,2) * g(1);
t102 = (t57 + t128) * sin(t53) / 0.2e1 + (t57 - t128) * sin(t54) / 0.2e1 + (t97 + t146) * sin(t73) / 0.2e1 + (t97 - t146) * sin(t74) / 0.2e1 + (t96 + t124) * cos(t54) / 0.2e1 + (t96 - t124) * cos(t53) / 0.2e1 + (t98 + t125) * cos(t74) / 0.2e1 + (t98 - t125) * cos(t73) / 0.2e1 + mrSges(11,1) * sin(t15) * g(3) - t141 + t144;
t99 = mrSges(7,1) * g(1);
t81 = qJ(1) - qJ(2);
t71 = cos(t81);
t70 = cos(t80);
t68 = sin(t81);
t59 = -qJ(1) + t72;
t58 = qJ(1) + t72;
t47 = -pkin(20) + t109;
t44 = cos(t47);
t43 = sin(t47);
t42 = pkin(1) * t147 + mrSges(3,1);
t39 = g(2) * t42;
t38 = t42 * g(1);
t19 = qJ(1) - t20;
t18 = qJ(1) + t20;
t1 = atan2(t34, -t33) + t16;
t2 = [(g(1) * t143 + g(2) * t142) * t93 + (-g(1) * t142 + g(2) * t143) * t87 + (t114 * t68 + t115 * t70 + t116 * t71 + t136 * t139) * pkin(1) + (-cos(t16) * mrSges(9,2) - sin(t16) * t76 + sin(t1) * mrSges(10,1) + cos(t1) * mrSges(10,2) - t44 * (pkin(9) * m(6) + mrSges(5,1) + t111) - t43 * (pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3)) - cos(t78) * mrSges(8,1) + sin(t78) * mrSges(8,2)) * t108 + t141 + t144 + ((sin(t23) + sin(t22)) * t114 + cos(t22) * t116 + cos(t23) * t115) * pkin(4), t102 + (mrSges(7,1) * sin(t72) - mrSges(10,1) * sin(t20) + mrSges(3,2) * t92 + mrSges(7,2) * cos(t72) - mrSges(10,2) * cos(t20) + t42 * t86 + t145) * g(3) + (t99 + t129) * cos(t58) / 0.2e1 + (t99 - t129) * cos(t59) / 0.2e1 + (-t130 + t135) * sin(t58) / 0.2e1 + (-t130 - t135) * sin(t59) / 0.2e1 + (t39 - t132) * t139 + (t39 + t132) * t68 / 0.2e1 + (t38 + t131) * t70 / 0.2e1 + (t38 - t131) * t71 / 0.2e1 + (t127 - t133) * sin(t18) / 0.2e1 + (-t126 - t134) * cos(t18) / 0.2e1 + (-t127 - t133) * sin(t19) / 0.2e1 + (t126 - t134) * cos(t19) / 0.2e1, g(3) * t145 + t102, -t111 * t108 + (-(g(1) * t93 + g(2) * t87) * t44 - g(3) * t43) * t107];
taug = t2(:);
