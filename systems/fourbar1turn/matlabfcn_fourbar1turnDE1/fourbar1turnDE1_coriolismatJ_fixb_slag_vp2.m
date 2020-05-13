% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnDE1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:39
% EndTime: 2020-04-12 19:25:53
% DurationCPUTime: 4.91s
% Computational Cost: add. (74633->178), mult. (105864->404), div. (4520->21), fcn. (29265->8), ass. (0->166)
t78 = sin(qJ(2));
t194 = pkin(2) * t78;
t145 = pkin(1) * t194;
t79 = cos(qJ(2));
t193 = pkin(2) * t79;
t89 = pkin(1) ^ 2;
t149 = -0.2e1 * pkin(1) * t193 + t89;
t88 = pkin(2) ^ 2;
t72 = t88 + t149;
t70 = 0.1e1 / t72 ^ 2;
t128 = t70 * t145;
t200 = -pkin(3) - pkin(4);
t65 = (pkin(2) - t200) * (pkin(2) + t200) + t149;
t199 = pkin(4) - pkin(3);
t66 = (pkin(2) - t199) * (pkin(2) + t199) + t149;
t162 = t65 * t66;
t92 = sqrt(-t162);
t156 = t78 * t92;
t57 = pkin(2) * t156;
t212 = pkin(3) ^ 2;
t213 = pkin(4) ^ 2;
t148 = t212 - t213;
t68 = t72 - t148;
t73 = pkin(1) - t193;
t50 = t68 * t73 - t57;
t111 = t50 * t128;
t155 = t79 * t92;
t203 = pkin(1) * pkin(2);
t118 = (-t65 - t66) * t203;
t55 = t78 * t118;
t59 = 0.1e1 / t92;
t165 = t59 * t55;
t61 = t68 * t194;
t22 = t61 + (-t155 + (0.2e1 * t73 * pkin(1) - t165) * t78) * pkin(2);
t69 = 0.1e1 / t72;
t173 = t22 * t69;
t82 = 0.1e1 / pkin(4);
t16 = (-t173 / 0.2e1 + t111) * t82;
t45 = 0.1e1 / t50 ^ 2;
t51 = t73 * t92 + t61;
t168 = t45 * t51;
t110 = t51 * t128;
t77 = t78 ^ 2;
t139 = t88 * t77 * pkin(1);
t150 = t68 * t193 + t57;
t164 = t59 * t73;
t24 = t164 * t55 + 0.2e1 * t139 + t150;
t172 = t24 * t69;
t18 = (t172 / 0.2e1 - t110) * t82;
t44 = 0.1e1 / t50;
t178 = t18 * t44;
t215 = -0.2e1 * t16 * t168 - 0.2e1 * t178;
t43 = t50 ^ 2;
t174 = t22 * t44 / t43;
t47 = t51 ^ 2;
t37 = t45 * t47 + 0.1e1;
t211 = (t168 * t24 - t174 * t47) / t37 ^ 2;
t152 = t43 + t47;
t83 = 0.1e1 / t213;
t31 = t152 * t83 * t70;
t25 = t31 ^ (-0.1e1 / 0.2e1);
t171 = t25 * t69;
t204 = 0.2e1 * t78;
t210 = t165 * t204 + t155;
t56 = pkin(1) * t156;
t67 = t72 + t148;
t197 = pkin(1) * t79;
t74 = -pkin(2) + t197;
t49 = -t67 * t74 - t56;
t41 = 0.1e1 / t49 ^ 2;
t198 = pkin(1) * t78;
t63 = t67 * t198;
t52 = -t74 * t92 + t63;
t169 = t41 * t52;
t116 = -0.2e1 * t128;
t205 = -0.2e1 * t74;
t21 = t63 + (-t155 + (pkin(2) * t205 - t165) * t78) * pkin(1);
t85 = 0.1e1 / pkin(3);
t17 = (t116 * t49 + t21 * t69) * t85;
t160 = t77 * t89;
t138 = pkin(2) * t160;
t151 = t67 * t197 + t56;
t163 = t59 * t74;
t23 = -t163 * t55 + 0.2e1 * t138 + t151;
t19 = (t116 * t52 + t23 * t69) * t85;
t40 = 0.1e1 / t49;
t177 = t19 * t40;
t209 = -t17 * t169 + t177;
t71 = t69 * t70;
t127 = t71 * t145;
t144 = 0.2e1 * t70;
t140 = 0.1e1 / t31 * ((t22 * t50 + t24 * t51) * t144 - 0.4e1 * t152 * t127) * t83 * t171;
t170 = t25 * t82;
t208 = t128 * t170 + t82 * t140 / 0.4e1;
t39 = t49 ^ 2;
t48 = t52 ^ 2;
t153 = t39 + t48;
t86 = 0.1e1 / t212;
t32 = t153 * t86 * t70;
t207 = 0.1e1 / t32;
t27 = t32 ^ (-0.1e1 / 0.2e1);
t202 = -t78 / 0.2e1;
t196 = pkin(2) * t69;
t195 = pkin(2) * t70;
t192 = pkin(3) * t72;
t191 = pkin(4) * t72;
t104 = t70 * (t21 * t49 + t23 * t52);
t106 = t153 * t127;
t9 = (0.2e1 * t104 - 0.4e1 * t106) * t86;
t190 = t27 * t207 * t9;
t189 = mrSges(4,1) * t49;
t188 = mrSges(4,2) * t52;
t187 = Ifges(5,1) * t51;
t186 = Ifges(3,4) * t78;
t185 = Ifges(5,4) * t50;
t184 = Ifges(5,4) * t51;
t183 = Ifges(5,5) * t51;
t182 = Ifges(5,2) * t50;
t181 = Ifges(5,6) * t50;
t134 = t27 * t69 * t85;
t159 = t78 * t49;
t166 = t52 * t79;
t12 = (-t159 - t166) * t134;
t180 = t12 * t49;
t167 = t49 * t79;
t158 = t78 * t52;
t20 = t134 * t158;
t13 = -t134 * t167 + t20;
t179 = t13 * t52;
t176 = t21 * t40 / t39;
t175 = t21 * t78;
t161 = t69 * t79;
t157 = t78 * t79;
t154 = t23 + t49;
t102 = -0.2e1 * t110 + t172;
t103 = 0.2e1 * t111 - t173;
t129 = t171 / 0.2e1;
t130 = -t171 / 0.2e1;
t131 = t88 * t157;
t135 = t69 * t170;
t115 = t144 * t203;
t141 = t69 * t190;
t4 = t20 + ((t159 / 0.2e1 + t166 / 0.2e1) * t141 + ((t157 * t52 + t49 * t77) * t115 + (-t154 * t79 - t175) * t69) * t27) * t85;
t5 = ((-t158 / 0.2e1 + t167 / 0.2e1) * t141 + ((t157 * t49 - t52 * t77) * t115 + ((-t21 + t52) * t79 + t154 * t78) * t69) * t27) * t85;
t1 = -m(4) * t131 + (-mrSges(4,1) * t13 + mrSges(4,2) * t12) * t194 - (-mrSges(4,1) * t5 + mrSges(4,2) * t4) * t193 + t78 * (Ifges(3,1) * t79 - t186) / 0.2e1 + (Ifges(3,2) * t79 + t186) * t202 + t4 * t12 * Ifges(4,1) + t5 * Ifges(4,2) * t13 + (0.2e1 * Ifges(3,4) * t79 + (Ifges(3,1) - Ifges(3,2)) * t78) * t79 / 0.2e1 + (t12 * t5 + t4 * t13) * Ifges(4,4) + (-pkin(1) * ((-t51 * mrSges(5,2) / 0.2e1 - t50 * mrSges(5,1) / 0.2e1) * t140 + (-mrSges(5,1) * t103 + mrSges(5,2) * t102) * t25) + (t51 * ((-t187 / 0.2e1 + t185 / 0.2e1) * t140 + (Ifges(5,1) * t102 + Ifges(5,4) * t103) * t25) * t129 + t50 * ((-t184 / 0.2e1 + t182 / 0.2e1) * t140 + (Ifges(5,4) * t102 + Ifges(5,2) * t103) * t25) * t130) * t82) * t82 + ((t24 * t82 * t129 - t208 * t51) * (-t185 + t187) + (t22 * t82 * t130 + t208 * t50) * (-t182 + t184)) * t135;
t147 = t1 * qJD(1);
t38 = t41 * t48 + 0.1e1;
t35 = 0.1e1 / t38;
t143 = t35 * t192;
t33 = 0.1e1 / t37;
t142 = t33 * t191;
t133 = 0.4e1 * t59 / t162 * t55 ^ 2;
t132 = t88 * t160;
t125 = t88 * t144;
t119 = t71 * t132;
t117 = -t133 / 0.4e1;
t113 = 0.2e1 * t17 * t52 * t192;
t109 = 0.8e1 * t119;
t105 = t27 * t125 * t198;
t53 = t118 * t79 - 0.4e1 * t132;
t101 = (t78 * t117 + 0.2e1 * (t202 * t53 - t79 * t55) * t59) * t69;
t7 = t209 * t143 + 0.1e1;
t6 = t142 * t215;
t3 = (t113 * t41 - 0.2e1 * t177 * t192) / t38 ^ 2 * (t169 * t23 - t176 * t48) + (0.2e1 * pkin(3) * t145 * t209 + t113 * t176) * t35 + ((-t17 * t23 - t19 * t21) * t41 + (((0.6e1 * t89 * pkin(2) * t157 + t74 * t117 - t53 * t163 - t63) * t69 + t52 * t109 + (t210 * t69 + (-0.4e1 * t23 * t78 - 0.2e1 * t166) * t195) * pkin(1)) * t40 - ((0.4e1 * t138 + t151) * t69 + t49 * t109 + (t101 + (t161 * t205 + (-0.2e1 * t167 - 0.4e1 * t175) * t70) * pkin(2)) * pkin(1)) * t169) * t85) * t143;
t2 = 0.2e1 * pkin(4) * t33 * t145 * t215 + 0.4e1 * (t178 * t211 + (t33 * t174 + t45 * t211) * t16 * t51) * t191 + 0.2e1 * ((-t16 * t24 + t18 * t22) * t45 + (-((t73 * t133 / 0.4e1 + t53 * t164 - t61 + t210 * pkin(2)) * t69 / 0.2e1 + 0.4e1 * t51 * t119 + (0.3e1 * t69 * t131 + (-0.2e1 * t24 * t78 - t51 * t79) * t195) * pkin(1)) * t44 - (-(0.4e1 * t139 + t150) * t69 / 0.2e1 - 0.4e1 * t50 * t119 + (-t101 / 0.2e1 + (-t73 * t161 + (t204 * t22 + t50 * t79) * t70) * pkin(1)) * pkin(2)) * t168) * t82) * t142;
t8 = [t1 * qJD(2), t147 + (Ifges(3,5) * t79 - Ifges(3,6) * t78 + ((t179 - t180) * t105 + ((-t180 / 0.2e1 + t179 / 0.2e1) * t190 + (t12 * t21 - t13 * t23 + t4 * t49 - t5 * t52) * t27) * t196) * t85 * mrSges(4,3) + (-t181 + t183) * t2 * t135 + (t3 * t12 + t7 * t4) * Ifges(4,5) + ((-t183 / 0.2e1 + t181 / 0.2e1) * t140 + (Ifges(5,5) * t102 + Ifges(5,6) * t103) * t25) * t82 * t6 + (t13 * t3 + t5 * t7) * Ifges(4,6)) * qJD(2); -t147, (t7 * Ifges(4,3) * t3 + t6 * Ifges(5,3) * t2 + m(4) * (-t153 * t9 / t32 ^ 2 * t125 + (0.4e1 * t104 - 0.8e1 * t106) * t207 * t88) * t86 / 0.4e1 + ((-t188 + t189) * t7 * t105 + ((-t188 / 0.2e1 + t189 / 0.2e1) * t7 * t190 + ((t23 * t7 + t3 * t52) * mrSges(4,2) + (-t21 * t7 - t3 * t49) * mrSges(4,1)) * t27) * t196) * t85) * qJD(2);];
Cq = t8;
