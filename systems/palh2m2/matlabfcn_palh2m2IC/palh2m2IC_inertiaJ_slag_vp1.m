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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2IC_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2IC_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:27
% EndTime: 2020-05-03 06:51:28
% DurationCPUTime: 0.83s
% Computational Cost: add. (726->181), mult. (965->239), div. (0->0), fcn. (348->50), ass. (0->99)
t187 = rSges(7,1) ^ 2;
t227 = -t187 / 0.2e1;
t166 = sin(qJ(6));
t172 = cos(qJ(6));
t111 = rSges(7,1) * t172 - rSges(7,2) * t166;
t161 = qJ(3) + qJ(4);
t148 = qJ(2) + t161;
t136 = qJ(5) + t148;
t116 = sin(t136);
t118 = cos(t136);
t225 = -rSges(7,3) * t116 + (pkin(3) + t111) * t118;
t224 = t172 ^ 2;
t222 = pkin(5) * m(7);
t168 = sin(qJ(4));
t169 = sin(qJ(3));
t170 = sin(qJ(2));
t174 = cos(qJ(4));
t175 = cos(qJ(3));
t176 = cos(qJ(2));
t221 = pkin(5) * (t176 * (-t168 * t169 + t174 * t175) - t170 * (t168 * t175 + t169 * t174));
t220 = m(3) * rSges(3,3);
t219 = m(5) * rSges(5,3);
t218 = m(6) * rSges(6,3);
t217 = rSges(7,2) * m(7);
t167 = sin(qJ(5));
t216 = pkin(5) * t167;
t189 = pkin(5) ^ 2;
t215 = m(6) * t189;
t214 = rSges(7,1) * rSges(7,3);
t213 = rSges(7,2) * rSges(7,3);
t157 = qJ(6) + qJ(2);
t146 = qJ(3) + t157;
t132 = qJ(4) + t146;
t212 = rSges(7,1) * cos(t132);
t211 = rSges(7,2) * sin(t132);
t103 = (pkin(2) * t175 + pkin(4)) * t176 - t170 * t169 * pkin(2) + pkin(1);
t177 = cos(qJ(1));
t102 = t103 * t177;
t209 = t177 * t221 + t102;
t208 = pkin(3) * rSges(7,2);
t207 = (pkin(2) * t174 + pkin(5)) * t167;
t206 = t168 * cos(qJ(5));
t160 = qJ(5) - qJ(6);
t159 = qJ(5) + qJ(6);
t158 = -qJ(6) + qJ(2);
t147 = qJ(3) + t158;
t133 = qJ(4) + t147;
t205 = (sin(t133) * rSges(7,2) + cos(t133) * rSges(7,1)) * t222 / 0.2e1;
t204 = m(7) * t214;
t185 = rSges(7,2) ^ 2;
t149 = t185 + t187;
t202 = 0.2e1 * t118;
t145 = qJ(4) + t160;
t144 = qJ(4) + t159;
t180 = 0.2e1 * t189;
t198 = t180 / 0.2e1 + t185 / 0.2e1;
t197 = -t103 - t221;
t196 = m(7) * t149 + Icges(7,3);
t124 = sin(t148);
t125 = cos(t148);
t194 = rSges(5,1) * t125 - rSges(5,2) * t124;
t193 = rSges(6,1) * t118 - rSges(6,2) * t116;
t162 = qJ(2) + qJ(3);
t139 = sin(t162);
t143 = cos(t162);
t192 = t176 * pkin(4) + rSges(4,1) * t143 - rSges(4,2) * t139 + pkin(1);
t191 = rSges(3,1) * t176 - rSges(3,2) * t170 + pkin(1);
t190 = pkin(4) ^ 2;
t183 = 0.2e1 * qJ(6);
t171 = sin(qJ(1));
t163 = pkin(3) * rSges(7,1);
t141 = cos(t158);
t140 = cos(t157);
t138 = sin(t158);
t137 = sin(t157);
t135 = qJ(3) + t145;
t134 = qJ(3) + t144;
t120 = qJ(5) + t133;
t119 = qJ(5) + t132;
t112 = rSges(2,1) * t177 - rSges(2,2) * t171;
t110 = rSges(7,1) * t166 + rSges(7,2) * t172;
t109 = -rSges(2,1) * t171 - rSges(2,2) * t177;
t101 = (-m(5) * sin(t161) * rSges(5,2) + cos(t161) * (m(5) * rSges(5,1) + (m(6) + m(7)) * pkin(5))) * pkin(4);
t100 = rSges(3,3) * t171 + t191 * t177;
t99 = rSges(3,3) * t177 - t191 * t171;
t97 = rSges(4,3) * t171 + t192 * t177;
t96 = rSges(4,3) * t177 - t192 * t171;
t94 = rSges(5,3) * t171 + t194 * t177 + t102;
t93 = rSges(5,3) * t177 + (-t103 - t194) * t171;
t92 = -t216 * t217 * t172 + ((rSges(7,3) - t216) * rSges(7,1) * m(7) - t204) * t166;
t91 = rSges(6,3) * t171 + t193 * t177 + t209;
t90 = rSges(6,3) * t177 + (-t193 + t197) * t171;
t89 = ((-rSges(4,3) * m(4) - t218 - t219) * pkin(4) - rSges(3,1) * t220 + Icges(3,5)) * t170 - t176 * (rSges(3,2) * t220 - Icges(3,6)) + ((t137 / 0.2e1 + t138 / 0.2e1) * rSges(7,2) + (t141 / 0.2e1 - t140 / 0.2e1) * rSges(7,1)) * pkin(4) * m(7);
t88 = -t125 * (rSges(5,2) * t219 - Icges(5,6)) + (-t212 / 0.2e1 + t211 / 0.2e1) * t222 + (-pkin(5) * t218 - rSges(5,1) * t219 + Icges(5,5)) * t124 + t205;
t87 = -t110 * t171 + t225 * t177 + t209;
t86 = -t110 * t177 + (t197 - t225) * t171;
t85 = t204 * t166 + ((t172 * t207 + (-sin(t160) / 0.2e1 - sin(t159) / 0.2e1) * pkin(5) + (-sin(t135) / 0.2e1 - sin(t134) / 0.2e1) * pkin(4) + (-sin(t145) / 0.2e1 - sin(t144) / 0.2e1 + t172 * t206) * pkin(2)) * rSges(7,2) + (t166 * (-rSges(7,3) + t207) + (-cos(t160) / 0.2e1 + cos(t159) / 0.2e1) * pkin(5) + (-cos(t135) / 0.2e1 + cos(t134) / 0.2e1) * pkin(4) + (-cos(t145) / 0.2e1 + cos(t144) / 0.2e1 + t166 * t206) * pkin(2)) * rSges(7,1)) * m(7);
t84 = ((t163 + t213) * m(7) - Icges(7,6)) * cos(t120) / 0.2e1 + ((t208 - t214) * m(7) + Icges(7,5)) * sin(t120) / 0.2e1 + ((t163 - t213) * m(7) + Icges(7,6)) * cos(t119) / 0.2e1 + ((-t208 - t214) * m(7) + Icges(7,5)) * sin(t119) / 0.2e1 + t196 * t202 / 0.2e1 + t205 + (0.2e1 * pkin(1) * t111 + (-t211 + t212) * pkin(5) + ((-t137 + t138) * rSges(7,2) + (t140 + t141) * rSges(7,1)) * pkin(4) + ((sin(t147) - sin(t146)) * rSges(7,2) + (cos(t147) + cos(t146)) * rSges(7,1)) * pkin(2)) * m(7) / 0.2e1;
t1 = [Icges(3,2) * t176 ^ 2 + Icges(4,2) * t143 ^ 2 + Icges(5,2) * t125 ^ 2 + Icges(2,3) + (Icges(3,1) * t170 + 0.2e1 * Icges(3,4) * t176) * t170 + (Icges(4,1) * t139 + 0.2e1 * Icges(4,4) * t143) * t139 + (Icges(5,1) * t124 + 0.2e1 * Icges(5,4) * t125) * t124 + (Icges(7,3) + Icges(6,2)) * t118 ^ 2 + ((Icges(7,1) * t224 + Icges(6,1) + (-0.2e1 * Icges(7,4) * t172 + Icges(7,2) * t166) * t166) * t116 + (Icges(7,5) * t172 - Icges(7,6) * t166 + Icges(6,4)) * t202) * t116 + m(7) * (t86 ^ 2 + t87 ^ 2) + m(6) * (t90 ^ 2 + t91 ^ 2) + m(3) * (t100 ^ 2 + t99 ^ 2) + m(2) * (t109 ^ 2 + t112 ^ 2) + m(5) * (t93 ^ 2 + t94 ^ 2) + m(4) * (t96 ^ 2 + t97 ^ 2), t89, t88, t84; t89, Icges(3,3) + (m(4) + m(5)) * t190 + (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) - t215 + (-t180 / 0.2e1 - t149 / 0.2e1 + t190 + t198 - t227) * m(7) + (t190 + t189) * m(6), t101, t85; t88, t101, Icges(5,3) + t215 - Icges(7,2) / 0.2e1 + Icges(7,1) / 0.2e1 + (t198 + t227) * m(7) + (-cos(t183) / 0.2e1 + t224) * (m(7) * (-t185 + t187) - Icges(7,1) + Icges(7,2)) + (rSges(5,1) ^ 2 + rSges(5,2) ^ 2) * m(5) + 0.2e1 * (-sin(t183) / 0.2e1 + t166 * t172) * (-rSges(7,1) * t217 + Icges(7,4)), t92; t84, t85, t92, t196;];
Mq = t1;
