% Calculate joint inertia matrix for
% palh2m1IC
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1IC_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:54
% EndTime: 2020-05-03 01:03:57
% DurationCPUTime: 2.00s
% Computational Cost: add. (1184->218), mult. (1741->258), div. (0->0), fcn. (764->8), ass. (0->121)
t137 = cos(qJ(4));
t129 = t137 ^ 2;
t113 = mrSges(6,2) * pkin(6) - Ifges(6,6);
t114 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t145 = pkin(4) * m(6);
t254 = -pkin(4) * mrSges(6,3) - pkin(6) * t145 + t113 * t132 - t114 * t136 - Ifges(5,4);
t211 = t254 * t129;
t217 = Ifges(6,4) * t132;
t206 = (mrSges(6,1) * pkin(4) + t217) * t136;
t219 = mrSges(6,2) * t132;
t198 = pkin(4) * t219;
t128 = t136 ^ 2;
t131 = Ifges(6,1) - Ifges(6,2);
t230 = pkin(6) * mrSges(6,3);
t255 = t131 * t128 + Ifges(5,1) + Ifges(6,2) + 0.2e1 * t230;
t252 = Ifges(5,2) + Ifges(6,3) - 0.2e1 * t198 - t255;
t261 = -0.2e1 * t206 - t252;
t166 = 0.2e1 * t261;
t158 = pkin(4) ^ 2;
t157 = pkin(6) ^ 2;
t159 = pkin(3) ^ 2;
t202 = t157 + t159;
t111 = -t158 + t202;
t224 = t111 * m(6);
t258 = 2 * pkin(3);
t112 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t133 = sin(qJ(4));
t98 = t112 * t133;
t207 = t159 * m(5) + t98 * t258;
t256 = Ifges(4,2) + t207;
t239 = -Ifges(4,1) + t256;
t260 = -0.2e1 * t224 - 0.2e1 * t239 - t166;
t229 = m(6) * (-t157 + t158);
t259 = t166 - 0.2e1 * t229;
t215 = t136 * mrSges(6,1);
t176 = t215 - t219;
t248 = -mrSges(5,1) - t145;
t86 = -t176 + t248;
t83 = t86 * t133;
t78 = -t83 + mrSges(4,2);
t141 = m(5) + m(6);
t84 = pkin(3) * t141 + mrSges(4,1) + t98;
t223 = (t157 + t158) * m(6);
t134 = sin(qJ(3));
t138 = cos(qJ(3));
t210 = t86 * t137;
t251 = (-(-t112 * t137 + t78) * t134 + (-t210 + t84) * t138) * pkin(2);
t250 = 0.2e1 * t86;
t220 = mrSges(6,1) * t132;
t73 = -0.2e1 * Ifges(6,4) * t128 + pkin(4) * t220 + (mrSges(6,2) * pkin(4) - t131 * t132) * t136 + Ifges(6,4) - Ifges(5,5);
t249 = t73 * t133;
t203 = t134 * t137;
t214 = t137 * t73;
t103 = t114 * t132;
t80 = -t113 * t136 - Ifges(5,6) - t103;
t75 = t80 * t133;
t59 = -t75 + t214;
t243 = -mrSges(6,2) * t136 - t220;
t172 = mrSges(5,3) - t243;
t68 = t172 * pkin(3) - Ifges(4,5) - t75;
t246 = t73 * t203 + (-t59 + t68) * t134 - Ifges(4,6) * t138;
t244 = -t229 + t261;
t225 = pkin(4) * t133;
t240 = (mrSges(6,2) * t225 - t113 * t137) * t132 - Ifges(6,3) * t133;
t231 = pkin(3) * t86;
t236 = 0.8e1 * t254 * t133 + 0.4e1 * t231;
t235 = 0.2e1 * t84;
t234 = -0.2e1 * t86;
t233 = 0.2e1 * t112;
t142 = m(4) + m(5);
t135 = sin(qJ(2));
t228 = pkin(1) * t135;
t227 = pkin(2) * t134;
t226 = pkin(3) * t112;
t179 = -0.4e1 * t254;
t221 = -t179 * t133 + 0.2e1 * t231;
t97 = -t219 - t248;
t185 = t97 + t215;
t205 = t134 * t112;
t213 = (t185 * pkin(3) + (-t138 * t86 + t205) * pkin(2)) * t137;
t212 = (t226 + (t112 * t138 - t185 * t134) * pkin(2)) * t133;
t76 = t80 * t137;
t209 = (0.2e1 * t217 + (pkin(3) * t137 + 0.2e1 * pkin(4)) * mrSges(6,1)) * t136;
t208 = Ifges(6,6) + t113;
t204 = t134 * t135;
t201 = -0.2e1 * t228;
t200 = -0.2e1 * t227;
t195 = 0.2e1 * t205;
t192 = t134 * t211;
t190 = Ifges(4,3) + t223;
t189 = t158 + t202;
t188 = pkin(4) * t137 + pkin(3);
t187 = pkin(3) * t133 + pkin(6);
t181 = t136 * (t188 * mrSges(6,1) + t114 * t133) + Ifges(6,3) * t137 + (-t188 * mrSges(6,2) - t113 * t133) * t132;
t180 = -0.2e1 * t254;
t174 = mrSges(6,1) * t225 - t114 * t137;
t81 = pkin(3) * t83;
t171 = -Ifges(4,4) - t81 - t254;
t170 = -0.2e1 * pkin(3) * t210 + t190 + t207 - 0.2e1 * t223;
t90 = 0.2e1 * t206;
t168 = t90 + t252;
t162 = t129 * t259 + t168 - t239;
t161 = pkin(1) ^ 2;
t160 = pkin(2) ^ 2;
t139 = cos(qJ(2));
t130 = t138 ^ 2;
t119 = t142 * t160;
t70 = Ifges(4,6) - t249;
t67 = t168 + t229;
t64 = (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t128 + t206 - t198 + (t158 / 0.2e1 - t157 / 0.2e1) * m(6) - t230 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t63 = -t132 * (t187 * mrSges(6,1) - Ifges(6,5)) + t103 + (-t187 * mrSges(6,2) + t208) * t136;
t57 = t132 * Ifges(6,5) + t208 * t136 + t103 + t243 * (pkin(2) * t203 + (pkin(2) * t138 + pkin(3)) * t133 + pkin(6));
t56 = (t68 + t214) * t135 * t138;
t54 = -t134 * (t76 + t249) + t138 * t59;
t53 = (t181 * t138 + (mrSges(6,1) * pkin(2) - t174 * t134) * t136 + t240 * t134 - pkin(2) * t219) * t139 + t176 * pkin(1) + (-(t174 * t136 - t240) * t138 - t181 * t134) * t135;
t52 = -t223 - t212 - t209 - t213 + t141 * t159 + t90 + (t98 + (t234 - t97) * t137) * pkin(3) + t251 + t190;
t51 = t56 + t246 * t139 + (-t54 + (t70 - t76) * t134) * t135;
t50 = t56 + (-Ifges(3,6) + t246) * t139 + (-t80 * t203 + t70 * t134 - Ifges(3,5) + (mrSges(4,3) + t172) * pkin(2) - t54) * t135;
t1 = [(t221 * t137 + t162 - t224) * t130 + (0.4e1 * t192 + (t228 * t233 + 0.4e1 * (t64 * t133 - t226 / 0.2e1) * t134) * t137 + t78 * t201 + 0.2e1 * t134 * t171) * t138 + t67 * t129 + (t180 * t133 + (pkin(1) * t204 - pkin(3)) * t250) * t137 + (t84 * t134 + mrSges(3,2)) * t201 - 0.2e1 * t136 * t217 + (t161 + t202) * m(6) + (m(3) + t142) * t161 + Ifges(3,1) + Ifges(2,3) + t255 + t256 + ((-0.8e1 * t64 * t129 * t204 + ((pkin(2) * t233 + t236 * t134) * t135 + pkin(1) * t234) * t137 + (-0.2e1 * t78 * pkin(2) + t260 * t134) * t135 + pkin(1) * t235) * t138 + ((-t133 * t259 + t227 * t250 - 0.2e1 * t226) * t135 + pkin(1) * t195) * t137 + (t84 * t200 + (2 * Ifges(3,4)) - 0.2e1 * Ifges(4,4) + t180 - 0.2e1 * t81) * t135 + 0.2e1 * (-t78 * t134 + mrSges(3,1) + (m(6) + t142) * pkin(2)) * pkin(1) - 0.4e1 * (((t133 * t67 - t226) * t137 + t171 + 0.2e1 * t211) * t130 - t211) * t135 + ((-0.8e1 * t192 + (0.4e1 * (t133 * t244 + t226) * t134 + pkin(2) * t234) * t137 + (0.4e1 * Ifges(4,4) + 0.4e1 * t81 - t179) * t134 + pkin(2) * t235) * t138 + (-0.4e1 * t129 * t244 - t236 * t137 - t260) * t130 + (t160 - t111) * m(6) + (pkin(2) * t195 + t221) * t137 + t119 + t78 * t200 - Ifges(3,1) + Ifges(3,2) + t162) * t139) * t139, t50, t51, t53; t50, -0.2e1 * t213 - 0.2e1 * t212 + Ifges(3,3) + t119 + (t160 + t189) * m(6) + 0.2e1 * t251 + t170, t52, t57; t51, t52, 0.2e1 * t90 + t189 * m(6) - 0.2e1 * t209 + (-t137 * t97 - t98) * t258 + t170, t63; t53, t57, t63, Ifges(6,3);];
Mq = t1;
