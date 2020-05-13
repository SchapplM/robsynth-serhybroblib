% Calculate joint inertia matrix for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:43:08
% EndTime: 2020-05-03 01:43:10
% DurationCPUTime: 1.49s
% Computational Cost: add. (1206->213), mult. (1841->298), div. (0->0), fcn. (1545->10), ass. (0->106)
t79 = sin(qJ(6));
t84 = cos(qJ(6));
t159 = -mrSges(7,1) * t79 - mrSges(7,2) * t84;
t100 = -mrSges(6,3) - t159;
t137 = Ifges(7,6) * t84;
t120 = mrSges(7,1) * pkin(3);
t74 = t84 ^ 2;
t160 = 0.2e1 * Ifges(7,4) * t74 - t79 * t120 - Ifges(7,4) - Ifges(6,5);
t119 = mrSges(7,2) * pkin(3);
t77 = Ifges(7,2) - Ifges(7,1);
t51 = t77 * t79 + t119;
t33 = t51 * t84 - t160;
t85 = cos(qJ(5));
t154 = t33 * t85;
t138 = Ifges(7,5) * t79;
t62 = -Ifges(6,6) + t138;
t80 = sin(qJ(5));
t158 = -t100 * pkin(5) - (t62 + t137) * t80 - Ifges(5,5) - t154;
t136 = Ifges(7,6) * t85;
t19 = t160 * t80 + (-t51 * t80 + t136) * t84 + t62 * t85 - Ifges(5,6);
t81 = sin(qJ(4));
t86 = cos(qJ(4));
t6 = -t158 * t86 + t19 * t81;
t96 = -mrSges(5,3) + t100;
t3 = t96 * pkin(2) + Ifges(4,5) + t6;
t7 = -t158 * t81 - t19 * t86;
t5 = Ifges(4,6) + t7;
t82 = sin(qJ(3));
t87 = cos(qJ(3));
t173 = t3 * t87 - t5 * t82;
t121 = pkin(2) * t86;
t122 = pkin(2) * t81;
t78 = mrSges(6,2) + mrSges(7,3);
t128 = t78 * t80;
t127 = t79 * mrSges(7,2);
t53 = pkin(3) * m(7) + mrSges(6,1) - t127;
t69 = t84 * mrSges(7,1);
t114 = t53 + t69;
t130 = t114 * t85;
t170 = -t128 + t130;
t89 = m(6) + m(7);
t66 = t89 * pkin(5) + mrSges(5,1);
t152 = -t66 - t170;
t34 = t114 * t80 + t78 * t85 + mrSges(5,2);
t167 = m(5) + t89;
t116 = pkin(5) * t128;
t139 = Ifges(7,4) * t79;
t156 = -0.2e1 * t79;
t102 = 0.2e1 * (t120 + t139) * t84 + pkin(3) ^ 2 * m(7) + t119 * t156 + t77 * t74 + Ifges(7,1) + Ifges(6,3);
t95 = t89 * pkin(5) ^ 2 + Ifges(5,3) + t102;
t93 = 0.2e1 * pkin(5) * t130 - 0.2e1 * t116 + t95;
t92 = pkin(2) ^ 2 * t167 + Ifges(4,3) + t93;
t171 = -0.2e1 * t152 * t121 - 0.2e1 * t34 * t122 + t92;
t169 = -t152 * t86 - t34 * t81;
t168 = pkin(2) * t167 + mrSges(4,1);
t161 = (t80 * t119 - t136) * t79 - Ifges(7,3) * t80;
t110 = pkin(3) * t85 + pkin(5);
t26 = (t110 * mrSges(7,1) + Ifges(7,5) * t80) * t84 + (-t110 * mrSges(7,2) - Ifges(7,6) * t80) * t79 + t85 * Ifges(7,3);
t49 = -t85 * Ifges(7,5) + t80 * t120;
t157 = (-t49 * t84 + t161) * t86 - t26 * t81;
t147 = pkin(4) * t82;
t56 = t78 * t147;
t146 = pkin(4) * t87;
t60 = pkin(2) + t146;
t24 = -t114 * t60 + t56;
t151 = t170 * pkin(5);
t148 = t3 * t82 + t5 * t87;
t145 = pkin(5) * t80;
t143 = mrSges(5,2) * t82;
t83 = sin(qJ(2));
t88 = cos(qJ(2));
t46 = t82 * t88 + t83 * t87;
t125 = t82 * t83;
t47 = t87 * t88 - t125;
t27 = -t46 * t81 + t47 * t86;
t28 = t46 * t86 + t47 * t81;
t14 = t27 * t80 + t28 * t85;
t141 = mrSges(7,3) * t14;
t131 = t114 * t82;
t129 = t66 * t82;
t126 = t81 * t82;
t117 = t78 * pkin(2);
t41 = pkin(4) * t131;
t124 = -t41 - 0.2e1 * t117;
t123 = pkin(2) * mrSges(5,2);
t118 = t66 * pkin(2);
t58 = pkin(5) * t86 + pkin(2);
t23 = (pkin(5) * t126 - t58 * t87 - pkin(4)) * t88 + (pkin(5) * t81 * t87 + t58 * t82) * t83 - pkin(1);
t109 = -0.2e1 * pkin(2) * t114 + t56;
t104 = -t137 - t138;
t59 = pkin(5) + t121;
t98 = t85 * t122 + t59 * t80;
t61 = -t88 * pkin(4) - pkin(1);
t52 = t69 - t127;
t50 = -Ifges(6,6) - t104;
t37 = (-pkin(2) * t87 - pkin(4)) * t88 + pkin(2) * t125 - pkin(1);
t30 = -t117 + (-t78 * t87 - t131) * pkin(4);
t29 = -t60 * t78 - t41;
t22 = t33 * t80 - t50 * t85;
t21 = t50 * t80 + t154;
t13 = t27 * t85 - t28 * t80;
t12 = t26 * t86 + (pkin(2) * mrSges(7,1) - t49 * t81) * t84 + t161 * t81 - pkin(2) * t127;
t11 = t21 * t81 + t22 * t86;
t10 = t21 * t86 - t22 * t81;
t9 = -pkin(3) * t13 + t23;
t1 = [0.2e1 * t37 * (-mrSges(5,1) * t27 + mrSges(5,2) * t28) + t27 * (Ifges(5,4) * t28 + Ifges(5,2) * t27) + t28 * (Ifges(5,1) * t28 + Ifges(5,4) * t27) + t46 * (Ifges(4,1) * t46 + Ifges(4,4) * t47) + t47 * (Ifges(4,4) * t46 + Ifges(4,2) * t47) + 0.2e1 * t61 * (-mrSges(4,1) * t47 + mrSges(4,2) * t46) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t88 + mrSges(3,2) * t83) + t83 * (Ifges(3,1) * t83 + Ifges(3,4) * t88) + t88 * (Ifges(3,4) * t83 + Ifges(3,2) * t88) + m(3) * pkin(1) ^ 2 + m(6) * t23 ^ 2 + m(5) * t37 ^ 2 + m(4) * t61 ^ 2 + Ifges(2,3) + (-0.2e1 * t23 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3)) * t13) * t13 + (0.2e1 * t23 * mrSges(6,2) + (Ifges(6,1) - t79 * (Ifges(7,4) * t84 - Ifges(7,2) * t79) + t84 * (Ifges(7,1) * t84 - t139)) * t14 + 0.2e1 * (Ifges(7,5) * t84 - Ifges(7,6) * t79 + Ifges(6,4)) * t13) * t14 + ((-mrSges(7,2) * t13 - t79 * t141) * t156 - 0.2e1 * t84 * (mrSges(7,1) * t13 - t84 * t141) + m(7) * (t79 ^ 2 + t74) * t9) * t9; (Ifges(3,6) + t148) * t88 + (Ifges(3,5) + (-mrSges(4,3) + t96) * pkin(4) + t173) * t83; 0.2e1 * (t168 + t169) * t146 + 0.2e1 * (t152 * t81 - t34 * t86 - mrSges(4,2)) * t147 + Ifges(3,3) + (m(4) + t167) * pkin(4) ^ 2 + t171; t148 * t88 + t83 * t173; (-t109 * t85 + t124 * t80 + 0.2e1 * t118 + (-t152 * t87 - t143) * pkin(4)) * t86 + t92 + (t124 * t85 + t109 * t80 - 0.2e1 * t123 + (-t34 * t87 - t129) * pkin(4)) * t81 + t168 * t146 - mrSges(4,2) * t147; t171; (t6 * t82 + t7 * t87) * t88 + (t6 * t87 - t7 * t82) * t83; (-t24 * t85 + t29 * t80 + t118) * t86 + (t24 * t80 + t29 * t85 - t123) * t81 + ((t66 * t87 - t143) * t86 + (-mrSges(5,2) * t87 - t129) * t81) * pkin(4) + t93; t169 * pkin(2) + 0.2e1 * t151 + t95; t93; (t10 * t82 + t11 * t87) * t88 + t83 * (t10 * t87 - t11 * t82); (pkin(5) * t114 - t24 * t86 + t30 * t81) * t85 + (-pkin(5) * t78 + t24 * t81 + t30 * t86) * t80 + t102; (pkin(5) * t53 + t59 * t69) * t85 - t116 + ((t53 * t86 - t78 * t81) * t85 + (-t114 * t81 - t78 * t86) * t80) * pkin(2) + t102; t151 + t102; t102; (pkin(4) * t52 + t12 * t87 + t157 * t82) * t88 + pkin(1) * t52 + (-t82 * t12 + t157 * t87) * t83; t104 + t159 * ((t86 * t147 + t60 * t81) * t85 + (-pkin(4) * t126 + t60 * t86 + pkin(5)) * t80); (-mrSges(7,2) * t98 - Ifges(7,6)) * t84 - t79 * (mrSges(7,1) * t98 + Ifges(7,5)); (-mrSges(7,2) * t145 - Ifges(7,6)) * t84 - t79 * (mrSges(7,1) * t145 + Ifges(7,5)); t104; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
