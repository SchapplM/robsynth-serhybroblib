% Calculate Gravitation load on the joints for
% palh3m1IC
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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:32:02
% EndTime: 2020-04-20 17:32:04
% DurationCPUTime: 1.32s
% Computational Cost: add. (841->136), mult. (758->161), div. (16->6), fcn. (587->28), ass. (0->94)
t65 = qJ(2) + qJ(7);
t55 = sin(t65);
t57 = cos(t65);
t64 = -qJ(7) + pkin(15);
t61 = -qJ(2) + t64;
t59 = -qJ(8) + t64;
t53 = -qJ(2) + t59;
t38 = sin(t53);
t39 = cos(t53);
t92 = t39 * mrSges(9,1) + t38 * mrSges(9,2);
t168 = -t57 * mrSges(8,1) + t55 * mrSges(8,2) - m(9) * pkin(3) * cos(t61) + t92;
t72 = cos(qJ(5));
t130 = mrSges(6,1) * t72;
t110 = qJ(3) + qJ(4);
t62 = qJ(2) + t110;
t51 = sin(t62);
t52 = cos(t62);
t167 = -mrSges(5,2) * t52 + (-pkin(8) * m(6) - mrSges(5,1) - t130) * t51;
t101 = -t52 * mrSges(5,1) + t51 * mrSges(5,2);
t66 = qJ(2) + qJ(3);
t56 = sin(t66);
t58 = cos(t66);
t102 = -t58 * mrSges(4,1) + t56 * mrSges(4,2);
t69 = sin(qJ(2));
t73 = cos(qJ(2));
t165 = -t73 * mrSges(3,1) + t69 * mrSges(3,2) - t101 - t102 + t168;
t74 = cos(qJ(1));
t145 = t167 * t74;
t155 = mrSges(4,1) * t56 + mrSges(4,2) * t58;
t163 = -t155 * t74 + t145;
t70 = sin(qJ(1));
t144 = t167 * t70;
t162 = -t155 * t70 + t144;
t160 = -m(9) - m(8) - m(4);
t136 = pkin(4) * t58;
t99 = pkin(8) * t52 + pkin(10) * t51;
t159 = m(6) * (-t99 - t136);
t129 = mrSges(9,1) * t38;
t135 = pkin(10) * t52;
t138 = pkin(3) * sin(t61);
t139 = pkin(1) * t69;
t137 = pkin(4) * t56;
t24 = t137 - t139;
t68 = sin(qJ(5));
t126 = mrSges(6,2) * t68;
t106 = t51 * t126;
t89 = -t52 * mrSges(6,3) - t106;
t154 = m(4) * t139 - m(5) * t24 - m(6) * (t24 - t135) - t89 - m(9) * (t138 - t139) + t129;
t60 = pkin(16) + t65;
t45 = sin(t60);
t19 = -pkin(2) * t45 - t139;
t48 = cos(t60);
t63 = t73 * pkin(1);
t20 = pkin(2) * t48 + t63;
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t6 = 0.1e1 / (-t45 * t71 + t48 * t67) / pkin(5) / pkin(2);
t152 = pkin(2) * t6 * (-t19 * t48 - t20 * t45);
t124 = t51 * mrSges(6,3);
t88 = -t124 + (t126 - t130) * t52;
t147 = -t101 - t88;
t78 = -t106 + (-pkin(10) * m(6) - mrSges(6,3)) * t52;
t151 = pkin(7) * ((m(6) * t99 + t147) * g(3) + (-t70 * t78 + t144) * g(2) + (-t74 * t78 + t145) * g(1));
t125 = mrSges(9,2) * t39;
t22 = t70 * t125;
t23 = t74 * t125;
t149 = pkin(9) * (g(3) * t92 - g(2) * (-t129 * t70 + t22) - g(1) * (-t129 * t74 + t23));
t148 = pkin(5) * t6 * (t19 * t71 + t20 * t67);
t146 = -m(5) * t137 - m(6) * (-t135 + t137) - t89;
t142 = -mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(3,3) - mrSges(4,3) + mrSges(2,2);
t105 = t63 - t136;
t21 = pkin(12) + t105;
t95 = t71 * mrSges(7,1) - t67 * mrSges(7,2);
t141 = m(6) * (-t21 + t99) - mrSges(2,1) + t124 - m(5) * t21 - m(3) * pkin(12) + pkin(6) * m(7) - t95 + t160 * (t63 + pkin(12)) + t165;
t118 = t68 * t74;
t117 = t70 * t68;
t116 = t70 * t72;
t115 = t72 * t74;
t93 = -mrSges(8,1) * t55 - mrSges(8,2) * t57;
t54 = pkin(14) + t110;
t47 = cos(t59);
t44 = sin(t59);
t43 = cos(t54);
t42 = sin(t54);
t17 = pkin(9) * t43 + cos(qJ(3)) * pkin(4);
t16 = -pkin(9) * t42 - sin(qJ(3)) * pkin(4);
t12 = -t47 * pkin(7) + pkin(3) * cos(t64);
t11 = -t44 * pkin(7) + pkin(3) * sin(t64);
t10 = t115 * t52 - t117;
t9 = t118 * t52 + t116;
t8 = t116 * t52 + t118;
t7 = t117 * t52 - t115;
t5 = 0.1e1 / (t42 * t47 + t43 * t44) / pkin(9) / pkin(7);
t1 = [(t10 * mrSges(6,1) - t9 * mrSges(6,2) + t141 * t74 + t142 * t70) * g(2) + (-t8 * mrSges(6,1) + t7 * mrSges(6,2) - t141 * t70 + t142 * t74) * g(1); (-m(5) * t105 - t159 - t88 - t95 * t152 - t168 * t148 + (-m(6) + t160) * t63 + t165) * g(3) + ((-t11 * t47 + t12 * t44) * t151 + (-t11 * t43 - t12 * t42) * t149) * t5 * t148 + (t148 * t22 + t154 * t70 + t162 - t22) * g(2) + (t148 * t23 + t154 * t74 + t163 - t23) * g(1) + ((mrSges(7,1) * t67 + mrSges(7,2) * t71) * t152 - (-m(9) * t138 + t129 - t93) * t148 + m(8) * t139 + mrSges(3,1) * t69 + mrSges(3,2) * t73 - t93) * (g(1) * t74 + g(2) * t70); ((t16 * t47 - t17 * t44) * t151 + (t16 * t43 + t17 * t42) * t149) * t5 + (m(5) * t136 - t102 + t147 - t159) * g(3) + (t146 * t70 + t162) * g(2) + (t146 * t74 + t163) * g(1); -g(1) * (mrSges(6,1) * t9 + mrSges(6,2) * t10) - g(2) * (mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(3) * (mrSges(6,1) * t68 + mrSges(6,2) * t72) * t51;];
taug = t1(:);
