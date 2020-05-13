% Calculate joint inertia matrix for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE1_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE1_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:03:12
% EndTime: 2020-05-07 02:03:13
% DurationCPUTime: 1.50s
% Computational Cost: add. (574->229), mult. (600->240), div. (0->0), fcn. (240->90), ass. (0->101)
t126 = pkin(17) + pkin(18);
t114 = pkin(15) + t126;
t104 = (pkin(16) + t114);
t180 = pkin(3) * m(9);
t143 = cos(qJ(2));
t181 = 0.2e1 * t143;
t147 = m(5) + m(6);
t178 = pkin(8) * mrSges(6,2);
t177 = (pkin(10) * mrSges(6,3));
t140 = sin(pkin(15));
t176 = pkin(8) * t140;
t144 = cos(pkin(15));
t175 = pkin(8) * t144;
t137 = sin(qJ(4));
t174 = t137 * mrSges(6,2);
t141 = cos(qJ(4));
t173 = t141 * mrSges(6,1);
t119 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t142 = cos(qJ(3));
t108 = -pkin(4) * t142 + pkin(1);
t138 = sin(qJ(3));
t139 = sin(qJ(2));
t171 = t138 * t139;
t66 = pkin(4) * t171 + t108 * t143;
t111 = m(4) + m(8) + m(9) + t147;
t172 = t111 * pkin(1) ^ 2;
t170 = t138 * t143;
t125 = t147 * pkin(4) ^ 2;
t128 = qJ(3) + qJ(2);
t169 = qJ(4) - qJ(2);
t127 = pkin(18) + pkin(15);
t102 = pkin(4) * t147 + mrSges(4,1) + mrSges(9,1);
t167 = t102 * pkin(1) * t142;
t165 = t125 + Ifges(4,3) + Ifges(9,3);
t164 = -pkin(10) * m(6) - mrSges(6,3);
t120 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t163 = 2 * t104;
t162 = -qJ(2) + t114;
t161 = qJ(2) + t114;
t93 = -qJ(2) + t104;
t92 = qJ(2) + t104;
t160 = mrSges(6,1) * t176 + t120 * t144;
t89 = -qJ(4) + t93;
t88 = -qJ(4) + t92;
t87 = qJ(4) + t93;
t86 = qJ(4) + t92;
t115 = sin(t128);
t116 = cos(t128);
t121 = qJ(4) + t128;
t122 = qJ(3) - t169;
t159 = -t116 * (Ifges(4,6) + Ifges(9,6)) + (-Ifges(4,5) - Ifges(9,5)) * t115 + (mrSges(5,3) * t115 + (-cos(t121) / 0.2e1 + cos(t122) / 0.2e1) * mrSges(6,1) + (sin(t121) + sin(t122)) * mrSges(6,2) / 0.2e1) * pkin(4);
t155 = pkin(8) ^ 2;
t154 = pkin(10) ^ 2;
t152 = 0.2e1 * qJ(2);
t151 = 0.2e1 * qJ(4);
t146 = pkin(8) * mrSges(6,1);
t136 = mrSges(4,2) + mrSges(9,2);
t133 = cos(pkin(16));
t132 = sin(pkin(16));
t130 = qJ(2) + qJ(4);
t129 = t152 + qJ(3);
t124 = pkin(8) * m(6) + mrSges(5,1);
t123 = -qJ(2) + pkin(14) - pkin(15);
t118 = -qJ(2) + t127;
t117 = qJ(2) + t127;
t113 = cos(t126);
t112 = sin(t126);
t110 = mrSges(5,2) + t164;
t109 = 0.2e1 * t128;
t107 = cos(t123);
t106 = sin(t123);
t105 = 0.2e1 * t127;
t101 = 0.2e1 * t123;
t100 = pkin(1) * t136 * t138;
t99 = qJ(3) + t161;
t98 = -qJ(3) + t162;
t97 = -qJ(4) + t104;
t96 = qJ(4) + t104;
t95 = -qJ(4) + t163;
t94 = qJ(4) + t163;
t91 = -qJ(3) + t93;
t90 = qJ(3) + t92;
t85 = 2 * t104;
t83 = mrSges(6,1) * t137 + mrSges(6,2) * t141;
t78 = -qJ(3) + t89;
t77 = -qJ(3) + t87;
t76 = qJ(3) + t88;
t75 = qJ(3) + t86;
t74 = 0.2e1 * t97;
t73 = 0.2e1 * t96;
t72 = t139 * t142 + t170;
t71 = t142 * t143 - t171;
t70 = mrSges(6,2) * t176 + t119 * t144;
t69 = mrSges(6,1) * t175 - t120 * t140;
t68 = mrSges(6,2) * t175 - t119 * t140;
t65 = -pkin(4) * t170 + t108 * t139;
t64 = t140 * t71 + t144 * t72;
t63 = -t140 * t72 + t144 * t71;
t62 = -t140 * t65 + t144 * t66;
t61 = t140 * t66 + t144 * t65;
t1 = [Ifges(3,4) * sin(t152) + ((cos(t90) + cos(t91)) * t124 + (sin(t90) + sin(t91)) * t110 + (-sin(t75) / 0.2e1 + sin(t76) / 0.2e1 - sin(t77) / 0.2e1 + sin(t78) / 0.2e1) * mrSges(6,2) + (cos(t75) / 0.2e1 + cos(t76) / 0.2e1 + cos(t77) / 0.2e1 + cos(t78) / 0.2e1) * mrSges(6,1)) * pkin(4) + (-Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2) + t125) * cos(t109) / 0.2e1 + (-Ifges(3,1) + Ifges(3,2) + t172) * cos(t152) / 0.2e1 + (-Ifges(7,1) + Ifges(7,2)) * cos(t101) / 0.2e1 + (-Ifges(8,1) + Ifges(8,2)) * cos(t105) / 0.2e1 + (-t120 - t178) * sin(t94) / 0.2e1 + (-t120 + t178) * sin(t95) / 0.2e1 + t125 / 0.2e1 + (0.2e1 * (-t154 + t155) * m(6) - (4 * t177) - (2 * Ifges(5,1)) - Ifges(6,1) + (2 * Ifges(5,2)) - Ifges(6,2) + (2 * Ifges(6,3))) * cos(t85) / 0.4e1 + t177 + (pkin(12) * t111 * t181 + t136 * sin(t129) + (-cos(t92) - cos(t93)) * t124 + (-sin(t92) - sin(t93)) * t110 + (-cos(t129) - t142) * t102 + (-sin(t117) - sin(t118)) * mrSges(8,2) + (-cos(t117) - cos(t118)) * mrSges(8,1) + (-cos(t162) - cos(t161)) * t180 + (sin(t86) / 0.2e1 + sin(t87) / 0.2e1 - sin(t88) / 0.2e1 - sin(t89) / 0.2e1) * mrSges(6,2) + (-cos(t86) / 0.2e1 - cos(t87) / 0.2e1 - cos(t88) / 0.2e1 - cos(t89) / 0.2e1) * mrSges(6,1)) * pkin(1) + (-0.2e1 * cos(t114) * t180 + mrSges(3,1) * t181 - 0.2e1 * mrSges(8,1) * cos(t127) - 0.2e1 * mrSges(3,2) * t139 - 0.2e1 * mrSges(8,2) * sin(t127) - 0.2e1 * t110 * sin(t104) - 0.2e1 * t124 * cos(t104) - 0.2e1 * t102 * t116 + 0.2e1 * t136 * t115 + (m(3) + t111) * pkin(12) + (sin(t96) - sin(t97)) * mrSges(6,2) + (-cos(t96) - cos(t97)) * mrSges(6,1)) * pkin(12) + (Ifges(4,4) + Ifges(9,4)) * sin(t109) + ((-sin(t99) + sin(t98)) * mrSges(9,2) + (cos(t99) + cos(t98)) * mrSges(9,1) + (0.1e1 / 0.2e1 + cos(0.2e1 * t114) / 0.2e1) * t180) * pkin(3) + (t146 - t119) * cos(t94) / 0.2e1 + (t146 + t119) * cos(t95) / 0.2e1 - pkin(8) * t174 + (pkin(6) * m(7) - 0.2e1 * mrSges(7,1) * t107 - 0.2e1 * mrSges(7,2) * t106) * pkin(6) + t100 + pkin(8) * t173 + (-cos(t73) / 0.8e1 - cos(t74) / 0.8e1 + cos(t151) / 0.4e1) * (Ifges(6,1) - Ifges(6,2)) + (t154 + t155) * m(6) / 0.2e1 + (t164 * pkin(8) - Ifges(5,4)) * sin(t85) + t172 / 0.2e1 - Ifges(7,4) * sin(t101) - Ifges(8,4) * sin(t105) + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.2e1 + Ifges(8,1) / 0.2e1 + Ifges(9,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.2e1 + Ifges(8,2) / 0.2e1 + Ifges(9,2) / 0.2e1 + Ifges(2,3) + Ifges(6,3) / 0.2e1 + (-sin(t151) / 0.2e1 + sin(t73) / 0.4e1 - sin(t74) / 0.4e1) * Ifges(6,4); Ifges(7,6) * t107 - Ifges(7,5) * t106 + Ifges(3,5) * t139 + Ifges(3,6) * t143 + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t139 + (-sin(-t169) / 0.2e1 - sin(t130) / 0.2e1) * mrSges(6,2) + (-cos(-t169) / 0.2e1 + cos(t130) / 0.2e1) * mrSges(6,1)) * pkin(1) + t159; Ifges(3,3) + Ifges(7,3) + 0.2e1 * t100 + t165 - 0.2e1 * t167 + t172; t159; t100 + t165 - t167; t165; ((-t132 * t160 + t69 * t133) * t141 + (t132 * t70 - t133 * t68) * t137 + Ifges(6,3) * (-t132 * t140 + t133 * t144)) * t113 + ((-t132 * t69 - t160 * t133) * t141 + (t132 * t68 + t133 * t70) * t137 - Ifges(6,3) * (t132 * t144 + t133 * t140)) * t112 - (t173 - t174) * (pkin(12) + t66); -((t132 * t62 + t133 * t61) * t113 + t112 * (-t132 * t61 + t133 * t62)) * t83; t83 * pkin(4) * ((t132 * t63 + t133 * t64) * t113 + (-t132 * t64 + t133 * t63) * t112); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
