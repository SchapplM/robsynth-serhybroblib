% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [2x2]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnTE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:38
% EndTime: 2020-04-12 19:18:46
% DurationCPUTime: 2.92s
% Computational Cost: add. (34297->155), mult. (47230->349), div. (1646->15), fcn. (13143->4), ass. (0->137)
t79 = sin(qJ(2));
t80 = cos(qJ(2));
t165 = pkin(2) * t80;
t89 = pkin(1) ^ 2;
t138 = -0.2e1 * pkin(1) * t165 + t89;
t172 = -pkin(3) - pkin(4);
t63 = (pkin(2) - t172) * (pkin(2) + t172) + t138;
t171 = -pkin(3) + pkin(4);
t64 = (pkin(2) - t171) * (pkin(2) + t171) + t138;
t149 = t63 * t64;
t90 = sqrt(-t149);
t144 = t79 * t90;
t55 = pkin(2) * t144;
t194 = pkin(3) ^ 2;
t137 = -pkin(4) ^ 2 + t194;
t88 = pkin(2) ^ 2;
t74 = t88 + t138;
t70 = t74 - t137;
t75 = pkin(1) - t165;
t48 = t70 * t75 - t55;
t43 = 0.1e1 / t48 ^ 2;
t166 = pkin(2) * t79;
t59 = t70 * t166;
t49 = t75 * t90 + t59;
t154 = t43 * t49;
t132 = pkin(1) * t166;
t72 = 0.1e1 / t74 ^ 2;
t121 = t72 * t132;
t106 = t49 * t121;
t71 = 0.1e1 / t74;
t176 = t71 / 0.2e1;
t169 = pkin(1) * t88;
t78 = t79 ^ 2;
t128 = t78 * t169;
t139 = t70 * t165 + t55;
t57 = 0.1e1 / t90;
t151 = t57 * t75;
t113 = pkin(1) * pkin(2) * (-t63 - t64);
t53 = t79 * t113;
t28 = t151 * t53 + 0.2e1 * t128 + t139;
t100 = t176 * t28 - t106;
t83 = 0.1e1 / pkin(4);
t21 = t100 * t83;
t42 = 0.1e1 / t48;
t160 = t21 * t42;
t107 = t48 * t121;
t141 = t80 * t90;
t152 = t57 * t53;
t26 = t59 + (-t141 + (0.2e1 * t75 * pkin(1) - t152) * t79) * pkin(2);
t101 = t176 * t26 - t107;
t19 = t101 * t83;
t196 = 0.2e1 * t19 * t154 - 0.2e1 * t160;
t167 = pkin(2) * t72;
t133 = pkin(1) * t167;
t142 = t80 * t79;
t69 = t74 + t137;
t61 = pkin(1) * t79 * t69;
t170 = pkin(1) * t80;
t76 = -pkin(2) + t170;
t50 = -t76 * t90 + t61;
t178 = t50 / 0.2e1;
t54 = pkin(1) * t144;
t47 = -t69 * t76 - t54;
t181 = t47 / 0.2e1;
t145 = t78 * t89;
t127 = pkin(2) * t145;
t140 = t69 * t170 + t54;
t150 = t57 * t76;
t27 = -t150 * t53 + 0.2e1 * t127 + t140;
t183 = t27 / 0.2e1;
t187 = -0.2e1 * t76;
t25 = t61 + (-t141 + (pkin(2) * t187 - t152) * t79) * pkin(1);
t184 = -t25 / 0.2e1;
t85 = 0.1e1 / pkin(3);
t11 = ((t142 * t47 - t50 * t78) * t133 + ((t178 + t184) * t80 + (t183 + t181) * t79) * t71) * t85;
t185 = -t11 / 0.2e1;
t146 = t71 * t85;
t153 = t47 * t80;
t174 = t79 / 0.2e1;
t30 = (t50 * t174 - t153 / 0.2e1) * t146;
t182 = -t30 / 0.2e1;
t156 = t26 * t42 * t43;
t45 = t49 ^ 2;
t35 = t43 * t45 + 0.1e1;
t193 = (t154 * t28 - t156 * t45) / t35 ^ 2;
t92 = t47 ^ 2;
t40 = 0.1e1 / t92;
t155 = t40 * t50;
t39 = 0.1e1 / t47;
t158 = t25 * t39 * t40;
t46 = t50 ^ 2;
t36 = t40 * t46 + 0.1e1;
t192 = (t155 * t27 - t158 * t46) / t36 ^ 2;
t186 = 0.2e1 * t79;
t191 = t152 * t186 + t141;
t111 = -0.2e1 * t121;
t22 = (t111 * t50 + t27 * t71) * t85;
t159 = t22 * t39;
t20 = (t111 * t47 + t25 * t71) * t85;
t190 = -t20 * t155 + t159;
t189 = t71 ^ 2;
t180 = -t48 / 0.2e1;
t179 = t49 / 0.2e1;
t177 = -t71 / 0.2e1;
t175 = -t79 / 0.2e1;
t173 = -t80 / 0.2e1;
t168 = pkin(2) * t71;
t164 = pkin(3) * t74;
t163 = pkin(4) * t74;
t162 = Ifges(3,4) * t79;
t157 = t25 * t79;
t148 = t71 * t80;
t143 = t80 * t50;
t10 = ((-t157 / 0.2e1 + t27 * t173) * t71 + (t142 * t50 + t47 * t78) * t133) * t85 + t30;
t122 = t88 * t142;
t29 = (t47 * t175 - t143 / 0.2e1) * t146;
t1 = m(4) * t122 - (-mrSges(4,1) * t30 + mrSges(4,2) * t29) * t166 + (-mrSges(4,1) * t11 + mrSges(4,2) * t10) * t165 + (Ifges(3,1) * t80 - t162) * t175 + (Ifges(3,2) * t80 + t162) * t174 - t10 * t29 * Ifges(4,1) + (0.2e1 * Ifges(3,4) * t80 + (Ifges(3,1) - Ifges(3,2)) * t79) * t173 + (t11 * t182 + t30 * t185) * Ifges(4,2) + 0.2e1 * (t10 * t182 + t29 * t185) * Ifges(4,4) + (pkin(1) * (mrSges(5,1) * t101 + mrSges(5,2) * t100) + ((-t49 * (Ifges(5,1) * t100 - Ifges(5,4) * t101) / 0.4e1 + t48 * (Ifges(5,4) * t100 - Ifges(5,2) * t101) / 0.4e1) * t71 + (-t28 * t189 / 0.4e1 + t106 * t176) * (Ifges(5,1) * t179 + Ifges(5,4) * t180) + (t26 * t189 / 0.4e1 + t107 * t177) * (Ifges(5,4) * t179 + Ifges(5,2) * t180)) * t83) * t83;
t136 = t1 * qJD(1);
t131 = t79 * t169;
t33 = 0.1e1 / t36;
t130 = t33 * t164;
t31 = 0.1e1 / t35;
t129 = t31 * t163;
t124 = 0.4e1 * t57 / t149 * t53 ^ 2;
t123 = t88 * t145;
t119 = t72 * t131;
t73 = t71 * t72;
t114 = t73 * t123;
t112 = -t124 / 0.4e1;
t105 = 0.8e1 * t114;
t51 = t113 * t80 - 0.4e1 * t123;
t99 = (t79 * t112 + 0.2e1 * (t175 * t51 - t80 * t53) * t57) * t71;
t9 = t190 * t130 + 0.1e1;
t8 = t129 * t196;
t3 = 0.2e1 * t190 * pkin(3) * t33 * t132 + 0.2e1 * (-t159 * t192 + (t158 * t33 + t40 * t192) * t20 * t50) * t164 + ((-t20 * t27 - t22 * t25) * t40 + (((0.6e1 * t89 * pkin(2) * t142 + t76 * t112 - t51 * t150 - t61) * t71 + t50 * t105 + (t191 * t71 + (-0.4e1 * t27 * t79 - 0.2e1 * t143) * t167) * pkin(1)) * t39 - ((0.4e1 * t127 + t140) * t71 + t47 * t105 + (t99 + (t148 * t187 + (-0.2e1 * t153 - 0.4e1 * t157) * t72) * pkin(2)) * pkin(1)) * t155) * t85) * t130;
t2 = 0.2e1 * pkin(4) * t31 * t132 * t196 + 0.4e1 * (t160 * t193 + (-t156 * t31 - t193 * t43) * t19 * t49) * t163 + 0.2e1 * ((t19 * t28 + t21 * t26) * t43 + (-((t75 * t124 / 0.4e1 + t51 * t151 - t59 + t191 * pkin(2)) * t176 + 0.4e1 * t49 * t114 + (0.3e1 * t71 * t122 + (-0.2e1 * t28 * t79 - t49 * t80) * t167) * pkin(1)) * t42 - ((0.4e1 * t128 + t139) * t177 - 0.4e1 * t48 * t114 + (-t99 / 0.2e1 + (-t75 * t148 + (t186 * t26 + t48 * t80) * t72) * pkin(1)) * pkin(2)) * t154) * t83) * t129;
t4 = [-t1 * qJD(2), -t136 + (((Ifges(5,5) * t179 + Ifges(5,6) * t180) * t71 * t2 + (Ifges(5,5) * t100 - Ifges(5,6) * t101) * t8) * t83 + Ifges(3,5) * t80 - Ifges(3,6) * t79 + ((-t29 * t47 + t30 * t50) * t119 + (t27 * t182 + t50 * t185 + t25 * t29 / 0.2e1 + t10 * t181) * t168) * t85 * mrSges(4,3) + (t9 * t10 + t3 * t29) * Ifges(4,5) + (t11 * t9 + t3 * t30) * Ifges(4,6)) * qJD(2); t136, (t9 * Ifges(4,3) * t3 + t8 * Ifges(5,3) * t2 + m(4) * ((t25 * t47 + t27 * t50) * t88 * t72 + 0.2e1 * (-t46 - t92) * pkin(2) * t73 * t131) / t194 / 0.4e1 + ((mrSges(4,1) * t47 - mrSges(4,2) * t50) * t9 * t119 + ((t178 * t3 + t183 * t9) * mrSges(4,2) + (t9 * t184 - t47 * t3 / 0.2e1) * mrSges(4,1)) * t168) * t85) * qJD(2);];
Cq = t4;
