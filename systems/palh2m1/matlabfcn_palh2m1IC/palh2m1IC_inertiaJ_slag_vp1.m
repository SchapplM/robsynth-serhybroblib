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
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1IC_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1IC_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:48
% EndTime: 2020-05-03 01:03:49
% DurationCPUTime: 0.56s
% Computational Cost: add. (677->156), mult. (1008->214), div. (0->0), fcn. (292->32), ass. (0->98)
t157 = qJ(2) + qJ(3);
t153 = -qJ(5) + qJ(2);
t210 = rSges(6,2) * m(6);
t190 = t210 / 0.2e1;
t211 = m(6) * rSges(6,1);
t191 = t211 / 0.2e1;
t220 = (sin(t153) * t190 + cos(t153) * t191) * pkin(2);
t176 = rSges(6,2) ^ 2;
t177 = rSges(6,1) ^ 2;
t207 = m(6) * (-t176 + t177);
t116 = -Icges(6,1) + Icges(6,2) + t207;
t130 = -rSges(6,1) * t210 + Icges(6,4);
t168 = rSges(6,3) + pkin(6);
t152 = t168 ^ 2;
t158 = sin(qJ(5));
t163 = cos(qJ(5));
t175 = 0.2e1 * qJ(5);
t178 = pkin(4) ^ 2;
t179 = pkin(3) ^ 2;
t171 = rSges(6,1) * pkin(4);
t194 = t171 * m(6);
t214 = -0.2e1 * t163;
t215 = t163 ^ 2;
t218 = -Icges(6,1) / 0.2e1 + Icges(6,2) / 0.2e1;
t219 = (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3) + t116 * t215 + (t152 + t176 + t178) * m(6) + t194 * t214 + t179 * m(5) + 0.2e1 * (t130 * t158 + t194) * t163 - t218 - t130 * sin(t175);
t160 = sin(qJ(3));
t165 = cos(qJ(3));
t217 = (-m(4) * rSges(4,2) * t160 + (m(4) * rSges(4,1) + (m(5) + m(6)) * pkin(3)) * t165) * pkin(2);
t143 = qJ(4) + t157;
t124 = sin(t143);
t126 = cos(t143);
t109 = rSges(6,1) * t163 - rSges(6,2) * t158;
t205 = pkin(4) + t109;
t216 = t168 * t124 + t205 * t126;
t213 = m(3) * rSges(3,3);
t212 = m(5) * rSges(5,3);
t209 = rSges(6,2) * pkin(4);
t208 = rSges(4,3) * m(4);
t204 = rSges(6,1) * t158;
t203 = rSges(6,2) * t163;
t202 = rSges(6,2) * t168;
t144 = t168 * rSges(6,1);
t199 = m(6) * t202;
t156 = qJ(2) + qJ(5);
t141 = qJ(3) + t153;
t188 = pkin(3) * t190;
t197 = sin(t141) * t188 + pkin(3) * cos(t141) * t191;
t145 = t176 + t177;
t195 = pkin(3) * t211;
t193 = 0.2e1 * t126;
t142 = qJ(3) + t156;
t189 = m(6) * t145 + Icges(6,3);
t149 = 0.2e1 * t152;
t173 = 0.2e1 * t178;
t113 = t149 + t173 + t145;
t186 = rSges(5,1) * t126 - rSges(5,2) * t124;
t161 = sin(qJ(2));
t166 = cos(qJ(2));
t185 = rSges(3,1) * t166 - rSges(3,2) * t161 + pkin(1);
t159 = sin(qJ(4));
t164 = cos(qJ(4));
t184 = pkin(2) * t160 * t164 + (pkin(2) * t165 + pkin(3)) * t159;
t138 = sin(t157);
t140 = cos(t157);
t183 = rSges(4,1) * t140 - rSges(4,2) * t138 + t166 * pkin(2) + pkin(1);
t123 = sin(t142);
t125 = cos(t142);
t81 = (rSges(4,2) * t208 - Icges(4,6)) * t140 + t123 * t188 - t125 * t195 / 0.2e1 + t197 + (rSges(4,1) * t208 + pkin(3) * t212 - Icges(4,5)) * t138;
t151 = cos(t175);
t181 = -t116 * t151 / 0.2e1 + t219 + (-t113 / 0.2e1 + t179) * m(6);
t180 = pkin(2) ^ 2;
t167 = cos(qJ(1));
t162 = sin(qJ(1));
t155 = qJ(4) - qJ(5);
t154 = qJ(4) + qJ(5);
t139 = cos(t156);
t137 = sin(t156);
t134 = qJ(4) + t141;
t133 = qJ(4) + t142;
t110 = rSges(2,1) * t167 - rSges(2,2) * t162;
t108 = t203 + t204;
t107 = -rSges(2,1) * t162 - rSges(2,2) * t167;
t95 = (pkin(3) * t165 + pkin(2)) * t166 - t161 * t160 * pkin(3) + pkin(1);
t92 = t95 * t167;
t91 = -rSges(3,3) * t162 + t185 * t167;
t90 = -rSges(3,3) * t167 - t185 * t162;
t89 = -rSges(4,3) * t162 + t183 * t167;
t88 = -rSges(4,3) * t167 - t183 * t162;
t87 = -rSges(5,3) * t162 + t186 * t167 + t92;
t86 = -rSges(5,3) * t167 + (-t186 - t95) * t162;
t85 = ((-pkin(3) * t159 - t168) * t210 + t199) * t163 - t159 * t195 * t158;
t84 = t199 * t163 + (-t184 * t204 + (-t168 - t184) * t203) * m(6);
t83 = -t108 * t162 + t216 * t167 + t92;
t82 = -t108 * t167 + (-t216 - t95) * t162;
t80 = ((t208 + t212) * pkin(2) + rSges(3,1) * t213 - Icges(3,5)) * t161 + t166 * (rSges(3,2) * t213 - Icges(3,6)) + (-rSges(6,1) * t139 / 0.2e1 + rSges(6,2) * t137 / 0.2e1) * pkin(2) * m(6) + t81 + t220;
t79 = ((t171 - t202) * m(6) + Icges(6,6)) * cos(t134) / 0.2e1 + ((t144 + t209) * m(6) - Icges(6,5)) * sin(t134) / 0.2e1 + ((t171 + t202) * m(6) - Icges(6,6)) * cos(t133) / 0.2e1 + ((t144 - t209) * m(6) - Icges(6,5)) * sin(t133) / 0.2e1 + t189 * t193 / 0.2e1 + m(6) * (pkin(1) * t109 + (-pkin(3) * t123 / 0.2e1 - pkin(2) * t137 / 0.2e1) * rSges(6,2) + (pkin(3) * t125 / 0.2e1 + pkin(2) * t139 / 0.2e1) * rSges(6,1)) + t197 + t220;
t78 = t181 + ((sin(t155) / 0.2e1 - sin(t154) / 0.2e1) * rSges(6,2) + (cos(t154) / 0.2e1 + cos(t155) / 0.2e1) * rSges(6,1) + (-t205 + pkin(4)) * t164) * m(6) * pkin(3) + t217;
t1 = [t166 ^ 2 * Icges(3,2) + t140 ^ 2 * Icges(4,2) + Icges(2,3) + (Icges(3,1) * t161 + 0.2e1 * Icges(3,4) * t166) * t161 + (Icges(4,1) * t138 + 0.2e1 * Icges(4,4) * t140) * t138 + (Icges(6,3) + Icges(5,2)) * t126 ^ 2 + ((t215 * Icges(6,1) + Icges(5,1) + (Icges(6,4) * t214 + Icges(6,2) * t158) * t158) * t124 + (-Icges(6,5) * t163 + Icges(6,6) * t158 + Icges(5,4)) * t193) * t124 + m(3) * (t90 ^ 2 + t91 ^ 2) + m(2) * (t107 ^ 2 + t110 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t86 ^ 2 + t87 ^ 2) + m(4) * (t88 ^ 2 + t89 ^ 2), t80, t81, t79; t80, (-t113 + t149 / 0.2e1 + t173 / 0.2e1 + t179 + t176 / 0.2e1 + t177 / 0.2e1 + t180) * m(6) + (-t116 + t207 / 0.2e1 + t218) * t151 + 0.2e1 * t217 + (m(4) + m(5)) * t180 + Icges(3,3) + (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) + t219, t78, t84; t81, t78, t181, t85; t79, t84, t85, t189;];
Mq = t1;
