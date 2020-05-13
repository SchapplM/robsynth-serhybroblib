% Calculate kinetic energy for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:23
% EndTime: 2020-04-12 19:33:26
% DurationCPUTime: 2.39s
% Computational Cost: add. (15516->233), mult. (21722->398), div. (1116->12), fcn. (6394->11), ass. (0->133)
t238 = pkin(4) ^ 2;
t237 = -2 * pkin(2);
t236 = (-pkin(3) - pkin(4));
t235 = (-pkin(3) + pkin(4));
t176 = sin(qJ(2));
t233 = pkin(1) * t176;
t178 = cos(qJ(2));
t232 = pkin(2) * t178;
t177 = sin(qJ(1));
t231 = Icges(2,4) * t177;
t230 = Icges(3,4) * t176;
t229 = Icges(3,4) * t178;
t185 = pkin(2) ^ 2;
t186 = pkin(1) ^ 2;
t218 = -0.2e1 * pkin(1) * t232 + t186;
t165 = t185 + t218;
t217 = pkin(3) ^ 2 - t238;
t146 = t165 + t217;
t169 = pkin(1) * t178 - pkin(2);
t144 = ((pkin(2) - t236) * (pkin(2) + t236)) + t218;
t145 = ((pkin(2) - t235) * (pkin(2) + t235)) + t218;
t187 = sqrt(-t144 * t145);
t220 = t176 * t187;
t124 = -pkin(1) * t220 - t146 * t169;
t127 = t146 * t233 - t169 * t187;
t160 = 0.1e1 / t165;
t184 = 0.1e1 / pkin(3);
t223 = t160 * t184;
t116 = qJ(2) + atan2(t127 * t223, t124 * t223);
t114 = sin(t116);
t228 = Icges(4,4) * t114;
t115 = cos(t116);
t227 = Icges(4,4) * t115;
t147 = t165 - t217;
t168 = pkin(1) - t232;
t125 = -pkin(2) * t220 + t147 * t168;
t226 = Icges(5,4) * t125;
t221 = t176 * t147;
t126 = pkin(2) * t221 + t168 * t187;
t225 = Icges(5,4) * t126;
t224 = t160 * t176 ^ 2;
t171 = V_base(6) + qJD(1);
t222 = t171 * t178;
t219 = t178 * t187;
t161 = 0.1e1 / t165 ^ 2;
t216 = t161 * t237;
t215 = V_base(5) * pkin(5) + V_base(1);
t123 = t126 ^ 2;
t181 = 0.1e1 / pkin(4);
t188 = t125 ^ 2;
t212 = ((t123 + t188) / t238 * t161) ^ (-0.1e1 / 0.2e1) * t160 * t181;
t211 = (-t144 - t145) * qJD(2) * pkin(2) * t233 / t187 * t160;
t167 = qJD(2) * t177 + V_base(4);
t208 = t211 / 0.2e1;
t207 = rSges(3,1) * t178 - rSges(3,2) * t176;
t206 = -rSges(4,1) * t115 + rSges(4,2) * t114;
t205 = -V_base(4) * pkin(5) + V_base(2);
t204 = Icges(3,1) * t178 - t230;
t203 = -Icges(4,1) * t115 + t228;
t202 = -Icges(3,2) * t176 + t229;
t201 = Icges(4,2) * t114 - t227;
t200 = Icges(3,5) * t178 - Icges(3,6) * t176;
t199 = -Icges(4,5) * t115 + Icges(4,6) * t114;
t179 = cos(qJ(1));
t192 = (-Icges(5,5) * t125 - Icges(5,6) * t126) * t212;
t122 = 0.1e1 / t188;
t88 = 0.2e1 * (-(t168 * t208 + (t185 * pkin(1) * t224 + ((t147 * t178 + t220) * t160 / 0.2e1 - t126 * t161 * t233) * pkin(2)) * qJD(2)) / t125 - (t176 * t208 + (-(-t219 + t221) * t160 / 0.2e1 + (t125 * t161 - t160 * t168) * t233) * qJD(2)) * pkin(2) * t126 * t122) * pkin(4) / (t122 * t123 + 0.1e1) * t165 * t181;
t86 = -t179 * t88 + V_base(5);
t87 = t177 * t88 + V_base(4);
t198 = t171 * (Icges(5,5) * t126 - Icges(5,6) * t125) * t212 + t86 * (-Icges(5,3) * t179 + t177 * t192) + t87 * (Icges(5,3) * t177 + t179 * t192);
t121 = 0.1e1 / t124 ^ 2;
t89 = ((-t169 * t211 + (0.2e1 * t186 * pkin(2) * t224 + ((t146 * t178 + t220) * t160 + t127 * t176 * t216) * pkin(1)) * qJD(2)) / t124 - (-t176 * t211 + (-t160 * t219 + ((t169 * t237 + t146) * t160 + t124 * t216) * t176) * qJD(2)) * pkin(1) * t127 * t121) * pkin(3) / (t121 * t127 ^ 2 + 0.1e1) * t165 * t184;
t84 = V_base(5) + (-qJD(2) - t89) * t179;
t85 = t177 * t89 + t167;
t197 = (-Icges(4,3) * t179 + t199 * t177) * t84 + (Icges(4,3) * t177 + t199 * t179) * t85 + (-Icges(4,5) * t114 - Icges(4,6) * t115) * t171;
t166 = -qJD(2) * t179 + V_base(5);
t196 = (-Icges(3,3) * t179 + t200 * t177) * t166 + (Icges(3,3) * t177 + t200 * t179) * t167 + (Icges(3,5) * t176 + Icges(3,6) * t178) * t171;
t195 = (-rSges(5,1) * t125 - rSges(5,2) * t126) * t212;
t194 = (-Icges(5,1) * t125 - t225) * t212;
t193 = (-Icges(5,2) * t126 - t226) * t212;
t104 = -Icges(4,6) * t179 + t201 * t177;
t105 = Icges(4,6) * t177 + t201 * t179;
t106 = -Icges(4,5) * t179 + t203 * t177;
t107 = Icges(4,5) * t177 + t203 * t179;
t111 = -Icges(4,2) * t115 - t228;
t112 = -Icges(4,1) * t114 - t227;
t191 = (t105 * t114 - t107 * t115) * t85 + (t104 * t114 - t106 * t115) * t84 + (t111 * t114 - t112 * t115) * t171;
t138 = -Icges(3,6) * t179 + t202 * t177;
t139 = Icges(3,6) * t177 + t202 * t179;
t140 = -Icges(3,5) * t179 + t204 * t177;
t141 = Icges(3,5) * t177 + t204 * t179;
t154 = Icges(3,2) * t178 + t230;
t157 = Icges(3,1) * t176 + t229;
t190 = (-t139 * t176 + t141 * t178) * t167 + (-t138 * t176 + t140 * t178) * t166 + (-t154 * t176 + t157 * t178) * t171;
t100 = (Icges(5,1) * t126 - t226) * t212;
t92 = -Icges(5,6) * t179 + t177 * t193;
t93 = Icges(5,6) * t177 + t179 * t193;
t94 = -Icges(5,5) * t179 + t177 * t194;
t95 = Icges(5,5) * t177 + t179 * t194;
t99 = (-Icges(5,2) * t125 + t225) * t212;
t189 = ((-t125 * t95 - t126 * t93) * t87 + (-t125 * t94 - t126 * t92) * t86 + (-t100 * t125 - t126 * t99) * t171) * t212;
t173 = Icges(2,4) * t179;
t164 = rSges(2,1) * t179 - rSges(2,2) * t177;
t163 = rSges(2,1) * t177 + rSges(2,2) * t179;
t162 = rSges(3,1) * t176 + rSges(3,2) * t178;
t159 = Icges(2,1) * t179 - t231;
t158 = Icges(2,1) * t177 + t173;
t156 = -Icges(2,2) * t177 + t173;
t155 = Icges(2,2) * t179 + t231;
t150 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t149 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t148 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t143 = rSges(3,3) * t177 + t207 * t179;
t142 = -rSges(3,3) * t179 + t207 * t177;
t135 = V_base(5) * rSges(2,3) - t163 * t171 + t215;
t134 = t164 * t171 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t133 = t163 * V_base(4) - t164 * V_base(5) + V_base(3);
t129 = -t142 * t171 + t162 * t166 + t215;
t128 = t143 * t171 - t162 * t167 + t205;
t120 = t142 * t167 - t143 * t166 + V_base(3);
t113 = -rSges(4,1) * t114 - rSges(4,2) * t115;
t109 = rSges(4,3) * t177 + t206 * t179;
t108 = -rSges(4,3) * t179 + t206 * t177;
t101 = (rSges(5,1) * t126 - rSges(5,2) * t125) * t212;
t97 = rSges(5,3) * t177 + t179 * t195;
t96 = -rSges(5,3) * t179 + t177 * t195;
t83 = -t108 * t171 + t113 * t84 + (t166 * t176 - t177 * t222) * pkin(2) + t215;
t82 = t109 * t171 - t113 * t85 + (-t167 * t176 + t179 * t222) * pkin(2) + t205;
t81 = t101 * t86 + (-pkin(1) * t177 - t96) * t171 + t215;
t80 = -t101 * t87 + (pkin(1) * t179 + t97) * t171 + t205;
t79 = t108 * t85 - t109 * t84 + V_base(3) + (-t166 * t179 + t167 * t177) * t232;
t78 = -t86 * t97 + t87 * t96 + V_base(3) + (t177 * V_base(4) - t179 * V_base(5)) * pkin(1);
t1 = m(1) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + t167 * (t196 * t177 + t190 * t179) / 0.2e1 + t166 * (t190 * t177 - t196 * t179) / 0.2e1 + m(4) * (t79 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t85 * (t197 * t177 + t191 * t179) / 0.2e1 + t84 * (t191 * t177 - t197 * t179) / 0.2e1 + m(5) * (t78 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + t87 * (t198 * t177 + t179 * t189) / 0.2e1 + t86 * (t177 * t189 - t198 * t179) / 0.2e1 + ((-t155 * t177 + t158 * t179 + Icges(1,4)) * V_base(5) + (-t156 * t177 + t159 * t179 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t155 * t179 + t158 * t177 + Icges(1,2)) * V_base(5) + (t156 * t179 + t159 * t177 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t139 * t178 + t141 * t176) * t167 + (t138 * t178 + t140 * t176) * t166 + (-t105 * t115 - t107 * t114) * t85 + (-t104 * t115 - t106 * t114) * t84 + ((-t125 * t93 + t126 * t95) * t87 + (-t125 * t92 + t126 * t94) * t86) * t212 + (Icges(2,3) + t154 * t178 + t157 * t176 - t111 * t115 - t112 * t114 + (t100 * t126 - t125 * t99) * t212) * t171) * t171 / 0.2e1 + t171 * V_base(4) * (Icges(2,5) * t179 - Icges(2,6) * t177) + V_base(5) * t171 * (Icges(2,5) * t177 + Icges(2,6) * t179) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
