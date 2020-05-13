% Calculate joint inertia matrix for
% palh2m2IC
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2IC_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2IC_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:36
% EndTime: 2020-05-03 06:51:37
% DurationCPUTime: 1.42s
% Computational Cost: add. (1730->151), mult. (2597->224), div. (0->0), fcn. (2089->10), ass. (0->75)
t165 = sin(qJ(6));
t170 = cos(qJ(6));
t222 = mrSges(7,1) * t165 + mrSges(7,2) * t170;
t180 = -mrSges(6,3) + t222;
t236 = t180 * pkin(5) + Ifges(5,5);
t175 = m(6) + m(7);
t235 = (m(5) + t175);
t159 = t175 * pkin(5);
t155 = t159 + mrSges(5,1);
t166 = sin(qJ(5));
t171 = cos(qJ(5));
t193 = pkin(3) * t166;
t232 = (mrSges(7,2) * t193 - Ifges(7,6) * t171) * t165 - Ifges(7,3) * t166;
t164 = mrSges(6,2) + mrSges(7,3);
t208 = t165 * mrSges(7,2);
t145 = pkin(3) * m(7) + mrSges(6,1) - t208;
t157 = t170 * mrSges(7,1);
t137 = -t145 - t157;
t203 = t137 * t171;
t223 = t164 * t166 - t155 + t203;
t189 = pkin(3) * t171 + pkin(5);
t119 = (mrSges(7,1) * t189 + Ifges(7,5) * t166) * t170 + (-mrSges(7,2) * t189 - Ifges(7,6) * t166) * t165 + t171 * Ifges(7,3);
t141 = mrSges(7,1) * t193 - t171 * Ifges(7,5);
t167 = sin(qJ(4));
t172 = cos(qJ(4));
t229 = t119 * t167 + (t141 * t170 - t232) * t172;
t161 = t170 ^ 2;
t228 = -2 * pkin(2);
t219 = 2 * pkin(2);
t168 = sin(qJ(3));
t218 = pkin(4) * t168;
t173 = cos(qJ(3));
t217 = pkin(4) * t173;
t215 = mrSges(5,2) * t168;
t169 = sin(qJ(2));
t174 = cos(qJ(2));
t138 = t168 * t174 + t169 * t173;
t199 = t168 * t169;
t139 = t173 * t174 - t199;
t120 = -t138 * t167 + t139 * t172;
t121 = t138 * t172 + t139 * t167;
t107 = t120 * t166 + t121 * t171;
t213 = mrSges(7,3) * t107;
t207 = -Ifges(5,6) * t167 + t236 * t172;
t127 = -t137 * t166 + t164 * t171 + mrSges(5,2);
t205 = t127 * t167;
t204 = t137 * t168;
t202 = t155 * t168;
t200 = t167 * t168;
t198 = Ifges(5,6) * t172 + t236 * t167;
t134 = pkin(4) * t204;
t150 = pkin(2) + t217;
t192 = t164 * pkin(2);
t196 = -t150 * t164 + t192 - (-t164 * t173 + t204) * pkin(4) + t134;
t195 = t134 - 0.2e1 * t192;
t194 = pkin(2) * t172;
t190 = -0.2e1 * t157;
t148 = pkin(5) * t172 + pkin(2);
t116 = (pkin(5) * t200 - t148 * t173 - pkin(4)) * t174 + (pkin(5) * t167 * t173 + t148 * t168) * t169 - pkin(1);
t185 = t137 * t219 + t164 * t218;
t151 = -t174 * pkin(4) - pkin(1);
t149 = pkin(5) + t194;
t144 = t157 - t208;
t132 = -0.2e1 * pkin(5) * t203;
t130 = (-pkin(2) * t173 - pkin(4)) * t174 + pkin(2) * t199 - pkin(1);
t124 = t222 * t166 * pkin(5);
t106 = t120 * t171 - t121 * t166;
t105 = t119 * t172 + (pkin(2) * mrSges(7,1) - t141 * t167) * t170 + t232 * t167 - pkin(2) * t208;
t102 = -pkin(3) * t106 + t116;
t93 = t222 * ((pkin(4) * t200 - t150 * t172 - pkin(5) + t149) * t166 + (-t172 * t218 + (pkin(2) - t150) * t167) * t171);
t92 = (pkin(4) * t144 + t105 * t173 - t229 * t168) * t174 + pkin(1) * t144 + (-t168 * t105 - t229 * t173) * t169;
t91 = t132 + (t149 * t157 + (t137 - t145 + t190) * pkin(5)) * t171 + (t196 * t166 + (t155 * t173 - t215) * pkin(4) + (t137 + t145) * t171 * pkin(2)) * t172 + (t196 * t171 + (-mrSges(5,2) * t173 - t202) * pkin(4)) * t167;
t90 = (t168 * t207 + t173 * t198) * t174 + (-t168 * t198 + t173 * t207) * t169;
t89 = (Ifges(3,5) + (-mrSges(4,3) - mrSges(5,3) + t180) * pkin(4)) * t169 + Ifges(3,6) * t174;
t1 = [0.2e1 * t130 * (-mrSges(5,1) * t120 + mrSges(5,2) * t121) + t120 * (Ifges(5,4) * t121 + Ifges(5,2) * t120) + t121 * (Ifges(5,1) * t121 + Ifges(5,4) * t120) + t138 * (Ifges(4,1) * t138 + Ifges(4,4) * t139) + t139 * (Ifges(4,4) * t138 + Ifges(4,2) * t139) + 0.2e1 * t151 * (-mrSges(4,1) * t139 + mrSges(4,2) * t138) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t174 + mrSges(3,2) * t169) + t169 * (Ifges(3,1) * t169 + Ifges(3,4) * t174) + t174 * (Ifges(3,4) * t169 + Ifges(3,2) * t174) + m(3) * pkin(1) ^ 2 + m(6) * t116 ^ 2 + m(5) * t130 ^ 2 + m(4) * t151 ^ 2 + Ifges(2,3) + (-0.2e1 * t116 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3)) * t106) * t106 + (0.2e1 * t116 * mrSges(6,2) + (Ifges(7,1) * t161 + Ifges(6,1) + (-0.2e1 * Ifges(7,4) * t170 + Ifges(7,2) * t165) * t165) * t107 + 0.2e1 * (Ifges(7,5) * t170 - Ifges(7,6) * t165 + Ifges(6,4)) * t106) * t107 + 0.2e1 * (-t165 * (-mrSges(7,2) * t106 - t165 * t213) - t170 * (mrSges(7,1) * t106 - t170 * t213)) * t102 + m(7) * (t165 ^ 2 + t161) * t102 ^ 2, t89, t90, t92; t89, 0.2e1 * ((pkin(2) * t235) - t172 * t223 - t205) * t217 + 0.2e1 * (-t127 * t172 + t167 * t223) * t218 - 0.2e1 * (-t185 * t171 + t195 * t166 + t155 * t219 + (-t223 * t173 - t215) * pkin(4)) * t172 - 0.2e1 * (t195 * t171 + t185 * t166 + mrSges(5,2) * t228 + (-t127 * t173 - t202) * pkin(4)) * t167 - 0.4e1 * pkin(2) * t205 + t235 * t217 * t228 + Ifges(3,3) + (m(4) + t235) * pkin(4) ^ 2 - 0.4e1 * t223 * t194, t91, t93; t90, t91, Ifges(5,3) + t132 + (t159 + (-0.2e1 * t145 + t190) * t171) * pkin(5), -t124; t92, t93, -t124, Ifges(7,3);];
Mq = t1;
