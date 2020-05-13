% Calculate joint inertia matrix for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:07:13
% EndTime: 2020-05-03 00:07:15
% DurationCPUTime: 2.16s
% Computational Cost: add. (799->199), mult. (1211->247), div. (0->0), fcn. (458->8), ass. (0->113)
t146 = Ifges(5,1) + Ifges(6,2);
t72 = sin(qJ(5));
t161 = mrSges(6,2) * t72;
t134 = pkin(4) * t161;
t159 = Ifges(6,4) * t72;
t76 = cos(qJ(5));
t153 = (mrSges(6,1) * pkin(4) + t159) * t76;
t68 = t76 ^ 2;
t71 = Ifges(6,1) - Ifges(6,2);
t45 = t71 * t68;
t205 = 0.2e1 * t153 - 0.2e1 * t134 - t45;
t206 = Ifges(5,2) + Ifges(6,3) - t146 + t205;
t164 = (pkin(6) * mrSges(6,3));
t108 = -(2 * t164) + t206;
t53 = mrSges(6,2) * pkin(6) - Ifges(6,6);
t54 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t83 = pkin(4) * m(6);
t192 = pkin(4) * mrSges(6,3) + pkin(6) * t83 - t53 * t72 + t54 * t76 + Ifges(5,4);
t77 = cos(qJ(4));
t69 = t77 ^ 2;
t156 = t192 * t69;
t118 = t76 * mrSges(6,1) - t161;
t27 = -mrSges(5,1) - t118 - t83;
t154 = t27 * t77;
t52 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t73 = sin(qJ(4));
t38 = t52 * t73;
t117 = t38 - t154;
t62 = 2 * t164;
t115 = Ifges(6,1) + Ifges(5,3) + t62 + t205;
t95 = pkin(6) ^ 2;
t96 = pkin(4) ^ 2;
t51 = (t95 + t96) * m(6);
t204 = t51 + t115;
t199 = 2 * pkin(3);
t97 = pkin(3) ^ 2;
t186 = (t97 * m(5)) + t38 * t199;
t202 = Ifges(4,2) + t186;
t107 = 0.2e1 * t108;
t143 = t95 + t97;
t50 = -t96 + t143;
t163 = t50 * m(6);
t178 = -Ifges(4,1) + t202;
t201 = -0.2e1 * t163 - 0.2e1 * t178 + t107;
t170 = m(6) * (-t95 + t96);
t200 = t107 + 0.2e1 * t170;
t24 = t27 * t73;
t187 = t24 - mrSges(4,2);
t74 = sin(qJ(3));
t196 = (-t52 * t77 - t187) * t74;
t81 = m(5) + m(6);
t195 = -pkin(3) * t81 - mrSges(4,1);
t194 = t117 * pkin(3) + t204;
t162 = mrSges(6,1) * t72;
t14 = -0.2e1 * Ifges(6,4) * t68 + pkin(4) * t162 + (mrSges(6,2) * pkin(4) - t71 * t72) * t76 + Ifges(6,4) - Ifges(5,5);
t127 = -t53 * t76 - t54 * t72;
t184 = -Ifges(5,6) + t127;
t4 = t14 * t73 + t184 * t77;
t189 = -Ifges(4,6) + t4;
t181 = t62 - t170 - t206;
t180 = -mrSges(6,2) * t76 - t162;
t25 = t38 - t195;
t165 = pkin(4) * t73;
t179 = (mrSges(6,2) * t165 - t53 * t77) * t72 - Ifges(6,3) * t73;
t167 = pkin(3) * t27;
t176 = -0.8e1 * t192 * t73 + 0.4e1 * t167;
t175 = -2 * pkin(2);
t174 = 0.2e1 * t25;
t173 = 0.2e1 * t52;
t172 = 0.2e1 * t74;
t82 = m(4) + m(5);
t75 = sin(qJ(2));
t169 = pkin(1) * t75;
t168 = pkin(2) * t27;
t166 = pkin(3) * t52;
t157 = t14 * t77;
t155 = t25 * t74;
t152 = t52 * t74;
t150 = t74 * t75;
t149 = t74 * t77;
t122 = 0.4e1 * t192;
t145 = -t122 * t73 + 0.2e1 * t167;
t142 = -0.2e1 * t169;
t141 = t187 * t175;
t139 = 0.2e1 * t152;
t78 = cos(qJ(3));
t112 = mrSges(5,3) - t180;
t16 = t184 * t73;
t9 = pkin(3) * t112 - Ifges(4,5) - t16;
t137 = t14 * t149 + t189 * t78 + t9 * t74;
t133 = t74 * t156;
t131 = t96 + t143;
t130 = pkin(4) * t77 + pkin(3);
t129 = pkin(3) * t73 + pkin(6);
t125 = (mrSges(6,1) * t130 + t54 * t73) * t76 + (-mrSges(6,2) * t130 - t53 * t73) * t72 + Ifges(6,3) * t77;
t123 = 0.2e1 * t192;
t119 = t189 * t74;
t114 = mrSges(6,1) * t165 - t54 * t77;
t111 = Ifges(4,3) + t115;
t22 = pkin(3) * t24;
t110 = -Ifges(4,4) - t22 + t192;
t103 = -0.2e1 * pkin(3) * t154 + t111 + t186;
t101 = -t200 * t69 + t108 - t178;
t99 = pkin(1) ^ 2;
t98 = pkin(2) ^ 2;
t79 = cos(qJ(2));
t70 = t78 ^ 2;
t59 = t82 * t98;
t8 = t108 + t170;
t6 = (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t68 + t153 - t134 + (t96 / 0.2e1 - t95 / 0.2e1) * m(6) - t164 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t3 = -t16 + t157;
t1 = (t9 + t157) * t75 * t78;
t2 = [(t145 * t77 + t101 - t163) * t70 + (-0.4e1 * t133 + (t169 * t173 + 0.4e1 * (t6 * t73 - t166 / 0.2e1) * t74) * t77 - t187 * t142 + t110 * t172) * t78 + t8 * t69 + (t123 * t73 + 0.2e1 * (pkin(1) * t150 - pkin(3)) * t27) * t77 + (mrSges(3,2) + t155) * t142 + t45 - 0.2e1 * t76 * t159 + (t99 + t143) * m(6) + t62 + (m(3) + t82) * t99 + Ifges(3,1) + Ifges(2,3) + t146 + ((-0.8e1 * t6 * t69 * t150 + ((pkin(2) * t173 + t176 * t74) * t75 - 0.2e1 * pkin(1) * t27) * t77 + (t201 * t74 - t141) * t75 + pkin(1) * t174) * t78 + ((t168 * t172 + t200 * t73 - 0.2e1 * t166) * t75 + pkin(1) * t139) * t77 + (t155 * t175 + (2 * Ifges(3,4)) - 0.2e1 * Ifges(4,4) + t123 - 0.2e1 * t22) * t75 + 0.2e1 * (t187 * t74 + mrSges(3,1) + (m(6) + t82) * pkin(2)) * pkin(1) - 0.4e1 * (((t8 * t73 - t166) * t77 + t110 - 0.2e1 * t156) * t70 + t156) * t75 + ((0.8e1 * t133 + (0.4e1 * (t181 * t73 + t166) * t74 - 0.2e1 * t168) * t77 + (0.4e1 * Ifges(4,4) + 0.4e1 * t22 - t122) * t74 + pkin(2) * t174) * t78 + (-t176 * t77 - 0.4e1 * t181 * t69 - t201) * t70 + t101 + t59 - t74 * t141 + (pkin(2) * t139 + t145) * t77 - Ifges(3,1) + Ifges(3,2) + (t98 - t50) * m(6)) * t79) * t79 + t202; t1 + (-Ifges(3,6) + t137) * t79 + (-Ifges(3,5) - t119 + (mrSges(4,3) + t112) * pkin(2)) * t75; t103 + t59 + 0.2e1 * (-(-t117 + t195) * t78 - t196) * pkin(2) + Ifges(3,3) + (t98 + t131) * m(6); -t119 * t75 + t137 * t79 + t1; t81 * t97 + t51 + t117 * t199 + ((-t154 + t25) * t78 - t196) * pkin(2) + t111; m(6) * t131 + t103; (t3 * t74 + t4 * t78) * t79 + (t3 * t78 - t4 * t74) * t75; ((-t27 * t78 + t152) * t77 + (t27 * t74 + t52 * t78) * t73) * pkin(2) + t194; t194; t204; (t125 * t78 + (mrSges(6,1) * pkin(2) - t114 * t74) * t76 + t179 * t74 - pkin(2) * t161) * t79 + t118 * pkin(1) + (-(t114 * t76 - t179) * t78 - t125 * t74) * t75; t72 * Ifges(6,5) + Ifges(6,6) * t76 + t180 * (pkin(2) * t149 + (pkin(2) * t78 + pkin(3)) * t73 + pkin(6)); (-mrSges(6,2) * t129 + Ifges(6,6)) * t76 - t72 * (mrSges(6,1) * t129 - Ifges(6,5)); t127; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
