% Calculate Gravitation load on the joints for
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [13x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:36
% EndTime: 2020-04-15 19:28:42
% DurationCPUTime: 1.00s
% Computational Cost: add. (652->134), mult. (624->147), div. (0->0), fcn. (498->22), ass. (0->88)
t164 = pkin(2) * m(10) + mrSges(9,1);
t65 = sin(qJ(2));
t137 = pkin(1) * t65;
t61 = qJ(2) + qJ(7);
t51 = pkin(19) - t61;
t44 = -qJ(10) + t51;
t33 = sin(t44);
t34 = cos(t44);
t100 = -t33 * mrSges(11,1) + t34 * mrSges(11,2);
t53 = sin(t61);
t56 = cos(t61);
t146 = t53 * mrSges(8,1) + t56 * mrSges(8,2) - t100;
t60 = qJ(2) + qJ(8);
t58 = qJ(9) + t60;
t47 = sin(t58);
t49 = cos(t58);
t111 = t47 * mrSges(10,1) + t49 * mrSges(10,2);
t52 = sin(t60);
t55 = cos(t60);
t148 = t55 * mrSges(9,2) + t164 * t52 - t111;
t32 = pkin(4) * sin(t51);
t69 = cos(qJ(2));
t162 = t148 - m(11) * (t32 - t137) + t65 * mrSges(3,1) + t69 * mrSges(3,2) + t146;
t161 = m(8) + m(4);
t124 = mrSges(10,2) * t47;
t160 = -mrSges(9,2) * t52 + t164 * t55 + t124;
t64 = sin(qJ(5));
t119 = t64 * mrSges(6,2);
t62 = qJ(2) + qJ(3);
t59 = qJ(4) + t62;
t48 = sin(t59);
t50 = cos(t59);
t159 = t48 * t119 + t50 * (pkin(11) * m(6) + mrSges(6,3));
t126 = mrSges(5,2) * t50;
t54 = sin(t62);
t57 = cos(t62);
t157 = mrSges(4,1) * t54 + mrSges(5,1) * t48 + mrSges(4,2) * t57 + t126;
t156 = -t50 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t48;
t154 = mrSges(11,1) * t34 + mrSges(11,2) * t33;
t66 = sin(qJ(1));
t70 = cos(qJ(1));
t153 = g(1) * t70 + g(2) * t66;
t127 = mrSges(10,1) * t49;
t21 = t66 * t127;
t150 = t160 * t66 - t21;
t24 = t70 * t127;
t149 = t160 * t70 - t24;
t110 = t50 * pkin(9) + t48 * pkin(11);
t68 = cos(qJ(5));
t115 = t68 * mrSges(6,1);
t145 = -(t115 - t119) * t50 + t156;
t112 = t154 * t70;
t125 = mrSges(8,2) * t53;
t144 = -t70 * t125 - t112;
t113 = t154 * t66;
t143 = -t66 * t125 - t113;
t102 = t57 * mrSges(4,1) - t54 * mrSges(4,2);
t105 = t48 * t115;
t128 = mrSges(8,1) * t56;
t132 = pkin(9) * t48;
t134 = pkin(4) * cos(t51);
t136 = pkin(1) * t69;
t133 = pkin(5) * t54;
t20 = -t133 - t136;
t141 = -m(11) * (-t134 - t136) - m(6) * (t20 - t132) + t105 + m(8) * t136 + t128;
t140 = -t102 + t145;
t139 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) - mrSges(11,3);
t45 = pkin(5) * t57;
t103 = t45 - t137;
t15 = pkin(15) + t103;
t63 = sin(qJ(6));
t67 = cos(qJ(6));
t93 = t67 * mrSges(7,1) - t63 * mrSges(7,2);
t138 = m(6) * (t15 + t110) + mrSges(2,1) - m(7) * pkin(14) + t93 + m(5) * t15 + t102 - t156 + t161 * (pkin(15) - t137) + (m(11) + m(10) + m(9) + m(3)) * pkin(15) - t162;
t118 = t64 * t70;
t117 = t66 * t64;
t116 = t66 * t68;
t114 = t68 * t70;
t99 = t159 * t66;
t98 = t159 * t70;
t85 = m(11) * t134 + t128;
t72 = m(6) * (-t132 - t133) - t105;
t71 = t126 + (m(6) * pkin(9) + mrSges(5,1) + t115) * t48;
t9 = t50 * t114 + t117;
t8 = -t50 * t118 + t116;
t7 = -t50 * t116 + t118;
t6 = t50 * t117 + t114;
t1 = [(-t9 * mrSges(6,1) - t8 * mrSges(6,2) - t138 * t70 + t139 * t66) * g(2) + (-t7 * mrSges(6,1) - t6 * mrSges(6,2) + t138 * t66 + t139 * t70) * g(1), (t141 * t66 + t143 + t150 - t99) * g(2) + (t141 * t70 + t144 + t149 - t98) * g(1) + (-m(5) * t103 - m(6) * (t103 + t110) + t161 * t137 + t140 + t162) * g(3) + t153 * (m(4) * t136 - m(5) * t20 + mrSges(3,1) * t69 - mrSges(3,2) * t65 + t157), -g(1) * (t72 * t70 + t98) - g(2) * (t72 * t66 + t99) + (-m(5) * t45 - m(6) * (t45 + t110) + t140) * g(3) + t153 * (m(5) * t133 + t157), (-m(6) * t110 + t145) * g(3) + (t71 * t66 - t99) * g(2) + (t71 * t70 - t98) * g(1), -g(1) * (mrSges(6,1) * t8 - mrSges(6,2) * t9) - g(2) * (-mrSges(6,1) * t6 + mrSges(6,2) * t7) - g(3) * (-mrSges(6,1) * t64 - mrSges(6,2) * t68) * t48, -g(3) * t93 + t153 * (mrSges(7,1) * t63 + mrSges(7,2) * t67), (-m(11) * t32 + t146) * g(3) + (t85 * t66 + t143) * g(2) + (t85 * t70 + t144) * g(1), t149 * g(1) + t150 * g(2) + t148 * g(3), -g(3) * t111 - g(1) * (-t70 * t124 + t24) - g(2) * (-t66 * t124 + t21), -g(1) * t112 - g(2) * t113 - g(3) * t100, 0, 0, 0];
taug = t1(:);
