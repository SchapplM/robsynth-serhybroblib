% Calculate Gravitation load on the joints for
% palh1m1IC
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:03:01
% EndTime: 2020-04-15 20:03:03
% DurationCPUTime: 1.61s
% Computational Cost: add. (1003->167), mult. (921->197), div. (20->8), fcn. (707->38), ass. (0->118)
t202 = pkin(2) * m(10) + mrSges(9,1);
t100 = cos(qJ(2));
t96 = sin(qJ(2));
t171 = pkin(1) * t96;
t88 = -qJ(7) + pkin(19);
t74 = -qJ(10) + t88;
t65 = -qJ(2) + t74;
t48 = sin(t65);
t49 = cos(t65);
t133 = -t48 * mrSges(11,1) + t49 * mrSges(11,2);
t92 = qJ(2) + qJ(7);
t80 = sin(t92);
t84 = cos(t92);
t181 = mrSges(8,1) * t80 + mrSges(8,2) * t84 - t133;
t90 = qJ(8) + qJ(9);
t86 = qJ(2) + t90;
t68 = sin(t86);
t70 = cos(t86);
t148 = t68 * mrSges(10,1) + t70 * mrSges(10,2);
t91 = qJ(2) + qJ(8);
t79 = sin(t91);
t83 = cos(t91);
t185 = mrSges(9,2) * t83 + t202 * t79 - t148;
t77 = -qJ(2) + t88;
t47 = pkin(4) * sin(t77);
t200 = t185 - m(11) * (t47 - t171) + mrSges(3,1) * t96 + mrSges(3,2) * t100 + t181;
t199 = m(8) + m(4);
t159 = mrSges(10,2) * t68;
t198 = -mrSges(9,2) * t79 + t202 * t83 + t159;
t95 = sin(qJ(5));
t154 = t95 * mrSges(6,2);
t140 = qJ(3) + qJ(4);
t87 = qJ(2) + t140;
t69 = sin(t87);
t71 = cos(t87);
t197 = t69 * t154 + t71 * (pkin(11) * m(6) + mrSges(6,3));
t195 = -t71 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t69;
t161 = mrSges(5,2) * t71;
t93 = qJ(2) + qJ(3);
t81 = sin(t93);
t85 = cos(t93);
t193 = -mrSges(4,1) * t81 - mrSges(5,1) * t69 - mrSges(4,2) * t85 - t161;
t192 = mrSges(11,1) * t49 + mrSges(11,2) * t48;
t101 = cos(qJ(1));
t99 = cos(qJ(5));
t151 = t99 * mrSges(6,1);
t102 = t161 + (m(6) * pkin(9) + mrSges(5,1) + t151) * t69;
t131 = t197 * t101;
t97 = sin(qJ(1));
t132 = t197 * t97;
t147 = t71 * pkin(9) + t69 * pkin(11);
t180 = -(t151 - t154) * t71 + t195;
t190 = pkin(8) * ((-m(6) * t147 + t180) * g(3) + (t102 * t97 - t132) * g(2) + (t101 * t102 - t131) * g(1));
t149 = t192 * t101;
t150 = t192 * t97;
t189 = pkin(10) * (-g(1) * t149 - g(2) * t150 - g(3) * t133);
t162 = mrSges(10,1) * t70;
t36 = t97 * t162;
t187 = t198 * t97 - t36;
t39 = t101 * t162;
t186 = t198 * t101 - t39;
t160 = mrSges(8,2) * t80;
t184 = -t101 * t160 - t149;
t183 = -t97 * t160 - t150;
t135 = t85 * mrSges(4,1) - t81 * mrSges(4,2);
t76 = pkin(20) + t92;
t60 = sin(t76);
t26 = pkin(3) * t60 + t171;
t165 = pkin(1) * t100;
t63 = cos(t76);
t27 = pkin(3) * t63 + t165;
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t9 = 0.1e1 / (-t60 * t94 - t63 * t98) / pkin(7) / pkin(3);
t182 = (t26 * t94 + t27 * t98) * pkin(7) * t9;
t179 = -g(1) * t101 - g(2) * t97;
t178 = -t135 + t180;
t177 = (t26 * t63 - t27 * t60) * t9 * pkin(3);
t138 = t69 * t151;
t163 = mrSges(8,1) * t84;
t167 = pkin(9) * t69;
t169 = pkin(4) * cos(t77);
t168 = pkin(5) * t81;
t35 = -t165 - t168;
t176 = (m(11) * t169 + t163) * t182 - m(6) * (t35 - t167) + t138 + m(8) * t165 + t163 - m(11) * (-t165 - t169);
t175 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) - mrSges(11,3);
t128 = mrSges(7,1) * t98 - mrSges(7,2) * t94;
t66 = pkin(5) * t85;
t136 = t66 - t171;
t28 = pkin(15) + t136;
t174 = m(6) * (t28 + t147) + mrSges(2,1) - m(7) * pkin(14) + t128 + m(5) * t28 + t135 - t195 + t199 * (pkin(15) - t171) + (m(11) + m(3) + m(10) + m(9)) * pkin(15) - t200;
t153 = t97 * t95;
t152 = t97 * t99;
t143 = t101 * t95;
t142 = t101 * t99;
t103 = m(6) * (-t167 - t168) - t138;
t89 = qJ(3) + pkin(17);
t82 = cos(t90);
t78 = sin(t90);
t75 = pkin(18) + t140;
t73 = cos(t89);
t72 = sin(t89);
t62 = cos(t75);
t59 = sin(t75);
t57 = cos(t74);
t56 = sin(t74);
t34 = t82 * pkin(12) - cos(qJ(8)) * pkin(2);
t33 = -t78 * pkin(12) + sin(qJ(8)) * pkin(2);
t24 = pkin(10) * t62 + cos(qJ(3)) * pkin(5);
t23 = pkin(10) * t59 + sin(qJ(3)) * pkin(5);
t18 = t57 * pkin(8) - pkin(4) * cos(t88);
t17 = t56 * pkin(8) - pkin(4) * sin(t88);
t16 = t142 * t71 + t153;
t15 = -t143 * t71 + t152;
t14 = -t152 * t71 + t143;
t13 = t153 * t71 + t142;
t6 = 0.1e1 / (-t56 * t59 + t57 * t62) / pkin(10) / pkin(8);
t1 = [(-t16 * mrSges(6,1) - t15 * mrSges(6,2) - t174 * t101 + t175 * t97) * g(2) + (-t14 * mrSges(6,1) - t13 * mrSges(6,2) + t175 * t101 + t174 * t97) * g(1); -((-t17 * t59 + t18 * t62) * t189 + (t17 * t57 - t18 * t56) * t190) * t6 * t182 + (-m(6) * (t136 + t147) - m(5) * t136 + (-m(11) * t47 + t181) * t182 - t128 * t177 + t199 * t171 + t178 + t200) * g(3) + (t176 * t97 + t183 * t182 - t132 + t183 + t187) * g(2) + (t176 * t101 + t184 * t182 - t131 + t184 + t186) * g(1) + (-(mrSges(7,1) * t94 + mrSges(7,2) * t98) * t177 - mrSges(3,1) * t100 + mrSges(3,2) * t96 + m(5) * t35 - m(4) * t165 + t193) * t179; -g(1) * (t103 * t101 + t131) - g(2) * (t103 * t97 + t132) + ((t23 * t56 - t24 * t57) * t190 + (-t23 * t62 + t24 * t59) * t189) * t6 + (-m(5) * t66 - m(6) * (t66 + t147) + t178) * g(3) + ((-t72 * t78 - t73 * t82) * pkin(12) * (t186 * g(1) + t187 * g(2) + t185 * g(3)) + (-t33 * t72 + t34 * t73) * (-g(3) * t148 - g(1) * (-t101 * t159 + t39) - g(2) * (-t159 * t97 + t36))) * pkin(6) / (t33 * t82 + t34 * t78) / pkin(12) + (-m(5) * t168 + t193) * t179; -g(1) * (mrSges(6,1) * t15 - mrSges(6,2) * t16) - g(2) * (-mrSges(6,1) * t13 + mrSges(6,2) * t14) - g(3) * (-mrSges(6,1) * t95 - mrSges(6,2) * t99) * t69;];
taug = t1(:);
