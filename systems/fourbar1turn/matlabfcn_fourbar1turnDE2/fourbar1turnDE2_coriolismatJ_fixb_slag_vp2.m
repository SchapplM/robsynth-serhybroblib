% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnDE2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:36
% EndTime: 2020-04-12 19:33:51
% DurationCPUTime: 5.16s
% Computational Cost: add. (74633->178), mult. (105864->409), div. (4520->21), fcn. (29265->8), ass. (0->165)
t102 = pkin(1) ^ 2;
t92 = cos(qJ(2));
t211 = pkin(2) * t92;
t163 = -0.2e1 * pkin(1) * t211 + t102;
t217 = -pkin(3) - pkin(4);
t75 = (pkin(2) - t217) * (pkin(2) + t217) + t163;
t216 = pkin(4) - pkin(3);
t76 = (pkin(2) - t216) * (pkin(2) + t216) + t163;
t179 = t75 * t76;
t105 = sqrt(-t179);
t91 = sin(qJ(2));
t165 = t105 * t91;
t67 = pkin(2) * t165;
t230 = pkin(3) ^ 2;
t231 = pkin(4) ^ 2;
t169 = t230 - t231;
t101 = pkin(2) ^ 2;
t86 = t101 + t163;
t82 = t86 - t169;
t87 = pkin(1) - t211;
t60 = t82 * t87 - t67;
t55 = 0.1e1 / t60 ^ 2;
t212 = pkin(2) * t91;
t71 = t82 * t212;
t61 = t105 * t87 + t71;
t188 = t55 * t61;
t160 = pkin(1) * t212;
t84 = 0.1e1 / t86 ^ 2;
t144 = t84 * t160;
t128 = t61 * t144;
t90 = t91 ^ 2;
t156 = t101 * t90 * pkin(1);
t170 = t82 * t211 + t67;
t69 = 0.1e1 / t105;
t181 = t69 * t87;
t137 = pkin(1) * pkin(2) * (-t75 - t76);
t65 = t91 * t137;
t34 = t181 * t65 + 0.2e1 * t156 + t170;
t83 = 0.1e1 / t86;
t190 = t34 * t83;
t95 = 0.1e1 / pkin(4);
t28 = (t190 / 0.2e1 - t128) * t95;
t54 = 0.1e1 / t60;
t196 = t28 * t54;
t129 = t60 * t144;
t164 = t92 * t105;
t182 = t69 * t65;
t32 = t71 + (-t164 + (0.2e1 * t87 * pkin(1) - t182) * t91) * pkin(2);
t191 = t32 * t83;
t26 = (-t191 / 0.2e1 + t129) * t95;
t233 = -0.2e1 * t26 * t188 - 0.2e1 * t196;
t185 = t61 * t34;
t53 = t60 ^ 2;
t192 = t32 * t54 / t53;
t57 = t61 ^ 2;
t47 = t55 * t57 + 0.1e1;
t229 = (t185 * t55 - t192 * t57) / t47 ^ 2;
t66 = pkin(1) * t165;
t81 = t86 + t169;
t215 = pkin(1) * t92;
t88 = -pkin(2) + t215;
t59 = -t81 * t88 - t66;
t51 = 0.1e1 / t59 ^ 2;
t73 = pkin(1) * t91 * t81;
t62 = -t105 * t88 + t73;
t189 = t51 * t62;
t222 = -0.2e1 * t88;
t31 = t73 + (-t164 + (pkin(2) * t222 - t182) * t91) * pkin(1);
t49 = t59 ^ 2;
t50 = 0.1e1 / t59;
t194 = t31 * t50 / t49;
t166 = t102 * t90;
t155 = pkin(2) * t166;
t171 = t81 * t215 + t66;
t180 = t69 * t88;
t33 = -t180 * t65 + 0.2e1 * t155 + t171;
t58 = t62 ^ 2;
t48 = t51 * t58 + 0.1e1;
t228 = (t189 * t33 - t194 * t58) / t48 ^ 2;
t85 = t83 * t84;
t143 = t85 * t160;
t159 = 0.2e1 * t84;
t172 = t53 + t57;
t96 = 0.1e1 / t231;
t41 = t172 * t96 * t84;
t35 = t41 ^ (-0.1e1 / 0.2e1);
t200 = ((t32 * t60 + t185) * t159 - 0.4e1 * t172 * t143) * t96 * t35 / t41;
t227 = t200 * t35;
t221 = 0.2e1 * t91;
t226 = t182 * t221 + t164;
t134 = -0.2e1 * t144;
t98 = 0.1e1 / pkin(3);
t29 = (t134 * t62 + t33 * t83) * t98;
t195 = t29 * t50;
t27 = (t134 * t59 + t31 * t83) * t98;
t225 = -t27 * t189 + t195;
t224 = t35 ^ 2;
t173 = t49 + t58;
t99 = 0.1e1 / t230;
t42 = t173 * t99 * t84;
t39 = 0.1e1 / t42;
t37 = t42 ^ (-0.1e1 / 0.2e1);
t220 = -t61 / 0.2e1;
t219 = -t91 / 0.2e1;
t214 = pkin(2) * t83;
t213 = pkin(2) * t84;
t210 = pkin(3) * t86;
t209 = pkin(4) * t86;
t208 = Ifges(5,1) * t61;
t207 = Ifges(3,4) * t91;
t205 = Ifges(5,4) * t60;
t204 = Ifges(5,4) * t61;
t203 = Ifges(5,5) * t61;
t202 = Ifges(5,2) * t60;
t201 = Ifges(5,6) * t60;
t121 = t173 * t143;
t122 = t31 * t59 + t33 * t62;
t21 = (t122 * t159 - 0.4e1 * t121) * t99;
t199 = t21 * t37 * t39;
t150 = t37 * t83 * t98;
t183 = t62 * t92;
t187 = t59 * t91;
t24 = (-t183 - t187) * t150;
t198 = t24 * t59;
t186 = t59 * t92;
t184 = t62 * t91;
t30 = t150 * t184;
t25 = -t150 * t186 + t30;
t197 = t25 * t62;
t193 = t31 * t91;
t178 = t83 * t92;
t175 = t92 * t91;
t174 = t33 + t59;
t115 = -0.2e1 * t128 + t190;
t116 = 0.2e1 * t129 - t191;
t123 = -t202 + t204;
t124 = -t205 + t208;
t167 = t101 * t91;
t147 = t92 * t167;
t148 = t83 * t200;
t145 = pkin(1) * t159;
t133 = pkin(2) * t145;
t154 = t83 * t199;
t8 = t30 + ((t187 / 0.2e1 + t183 / 0.2e1) * t154 + ((t175 * t62 + t59 * t90) * t133 + (-t174 * t92 - t193) * t83) * t37) * t98;
t9 = ((-t184 / 0.2e1 + t186 / 0.2e1) * t154 + ((t175 * t59 - t62 * t90) * t133 + ((-t31 + t62) * t92 + t174 * t91) * t83) * t37) * t98;
t1 = m(4) * t147 - (-mrSges(4,1) * t25 + mrSges(4,2) * t24) * t212 + (-mrSges(4,1) * t9 + mrSges(4,2) * t8) * t211 - t8 * t24 * Ifges(4,1) - t9 * Ifges(4,2) * t25 + (Ifges(3,1) * t92 - t207) * t219 + t91 * (Ifges(3,2) * t92 + t207) / 0.2e1 - (0.2e1 * Ifges(3,4) * t92 + (Ifges(3,1) - Ifges(3,2)) * t91) * t92 / 0.2e1 + (-t24 * t9 - t25 * t8) * Ifges(4,4) + (pkin(1) * ((mrSges(5,2) * t220 - t60 * mrSges(5,1) / 0.2e1) * t148 + (-mrSges(5,1) * t116 + mrSges(5,2) * t115) * t35) + ((t60 * ((-t204 / 0.2e1 + t202 / 0.2e1) * t148 + (Ifges(5,4) * t115 + Ifges(5,2) * t116) * t35) / 0.2e1 + ((-t208 / 0.2e1 + t205 / 0.2e1) * t148 + (Ifges(5,1) * t115 + Ifges(5,4) * t116) * t35) * t220) * t35 + (-t123 * t60 + t124 * t61) * t224 * t144 + ((-t34 * t224 / 0.2e1 + t61 * t227 / 0.4e1) * t124 + (t32 * t224 / 0.2e1 - t60 * t227 / 0.4e1) * t123) * t83) * t83 * t95) * t95;
t168 = t1 * qJD(1);
t45 = 0.1e1 / t48;
t158 = t45 * t210;
t43 = 0.1e1 / t47;
t157 = t43 * t209;
t149 = 0.4e1 * t69 / t179 * t65 ^ 2;
t146 = t101 * t166;
t136 = t85 * t146;
t135 = -t149 / 0.4e1;
t127 = 0.8e1 * t136;
t125 = mrSges(4,1) * t59 - mrSges(4,2) * t62;
t119 = t37 * t145 * t167;
t63 = t137 * t92 - 0.4e1 * t146;
t114 = (t91 * t135 + 0.2e1 * (t219 * t63 - t92 * t65) * t69) * t83;
t15 = t158 * t225 + 0.1e1;
t14 = t157 * t233;
t5 = 0.2e1 * t225 * pkin(3) * t45 * t160 + 0.2e1 * (-t195 * t228 + (t194 * t45 + t228 * t51) * t27 * t62) * t210 + ((-t27 * t33 - t29 * t31) * t51 + (((0.6e1 * t102 * pkin(2) * t175 + t88 * t135 - t63 * t180 - t73) * t83 + t62 * t127 + (t226 * t83 + (-0.4e1 * t33 * t91 - 0.2e1 * t183) * t213) * pkin(1)) * t50 - ((0.4e1 * t155 + t171) * t83 + t59 * t127 + (t114 + (t178 * t222 + (-0.2e1 * t186 - 0.4e1 * t193) * t84) * pkin(2)) * pkin(1)) * t189) * t98) * t158;
t4 = 0.2e1 * pkin(4) * t43 * t160 * t233 + 0.4e1 * (t196 * t229 + (t43 * t192 + t55 * t229) * t26 * t61) * t209 + 0.2e1 * ((-t26 * t34 + t28 * t32) * t55 + (-((t87 * t149 / 0.4e1 + t63 * t181 - t71 + t226 * pkin(2)) * t83 / 0.2e1 + 0.4e1 * t61 * t136 + (0.3e1 * t83 * t147 + (-0.2e1 * t34 * t91 - t61 * t92) * t213) * pkin(1)) * t54 - (-(0.4e1 * t156 + t170) * t83 / 0.2e1 - 0.4e1 * t60 * t136 + (-t114 / 0.2e1 + (-t87 * t178 + (t221 * t32 + t60 * t92) * t84) * pkin(1)) * pkin(2)) * t188) * t95) * t157;
t2 = [-t1 * qJD(2), -t168 + (((-t201 + t203) * t83 * t35 * t4 + ((-t203 / 0.2e1 + t201 / 0.2e1) * t148 + (Ifges(5,5) * t115 + Ifges(5,6) * t116) * t35) * t14) * t95 + Ifges(3,5) * t92 - Ifges(3,6) * t91 + ((t197 - t198) * t119 + ((-t198 / 0.2e1 + t197 / 0.2e1) * t199 + (t24 * t31 - t25 * t33 + t59 * t8 - t62 * t9) * t37) * t214) * t98 * mrSges(4,3) + (t15 * t8 + t5 * t24) * Ifges(4,5) + (t15 * t9 + t25 * t5) * Ifges(4,6)) * qJD(2); t168, (t15 * Ifges(4,3) * t5 + t14 * Ifges(5,3) * t4 + m(4) * (-0.8e1 * t39 * t121 + (0.4e1 * t122 * t39 - 0.2e1 * t173 / t42 ^ 2 * t21) * t84) * t101 * t99 / 0.4e1 + (t125 * t15 * t119 + (-t125 * t5 * t37 + ((t33 * t37 - t62 * t199 / 0.2e1) * mrSges(4,2) + (-t31 * t37 + t59 * t199 / 0.2e1) * mrSges(4,1)) * t15) * t214) * t98) * qJD(2);];
Cq = t2;
