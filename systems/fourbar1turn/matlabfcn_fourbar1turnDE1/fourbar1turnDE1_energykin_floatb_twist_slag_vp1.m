% Calculate kinetic energy for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:26
% EndTime: 2020-04-12 19:25:28
% DurationCPUTime: 2.85s
% Computational Cost: add. (22224->246), mult. (31914->418), div. (1740->13), fcn. (9254->10), ass. (0->135)
t254 = pkin(4) ^ 2;
t253 = pkin(3) ^ 2;
t203 = pkin(2) ^ 2;
t204 = pkin(1) ^ 2;
t195 = cos(qJ(2));
t246 = pkin(2) * t195;
t234 = -0.2e1 * pkin(1) * t246 + t204;
t182 = t203 + t234;
t233 = t253 - t254;
t163 = t182 + t233;
t186 = pkin(1) * t195 - pkin(2);
t193 = sin(qJ(2));
t250 = -pkin(3) - pkin(4);
t161 = (pkin(2) - t250) * (pkin(2) + t250) + t234;
t249 = -pkin(3) + pkin(4);
t162 = (pkin(2) - t249) * (pkin(2) + t249) + t234;
t205 = sqrt(-t161 * t162);
t237 = t193 * t205;
t141 = -pkin(1) * t237 - t163 * t186;
t247 = pkin(1) * t193;
t144 = t163 * t247 - t186 * t205;
t140 = t144 ^ 2;
t177 = 0.1e1 / t182;
t178 = 0.1e1 / t182 ^ 2;
t201 = 0.1e1 / pkin(3);
t206 = t141 ^ 2;
t227 = ((t140 + t206) / t253 * t178) ^ (-0.1e1 / 0.2e1) * t177 * t201;
t252 = (t141 * t193 + t144 * t195) * t227;
t251 = -0.2e1 * pkin(2);
t194 = sin(qJ(1));
t245 = Icges(2,4) * t194;
t244 = Icges(3,4) * t193;
t243 = Icges(3,4) * t195;
t164 = t182 - t233;
t185 = pkin(1) - t246;
t142 = -pkin(2) * t237 + t164 * t185;
t242 = Icges(5,4) * t142;
t238 = t193 * t164;
t143 = pkin(2) * t238 + t185 * t205;
t241 = Icges(5,4) * t143;
t240 = t144 * t193;
t239 = t177 * t193 ^ 2;
t188 = V_base(6) + qJD(1);
t236 = t195 * t188;
t235 = t195 * t205;
t232 = t178 * t251;
t231 = V_base(5) * pkin(5) + V_base(1);
t139 = t143 ^ 2;
t198 = 0.1e1 / pkin(4);
t207 = t142 ^ 2;
t228 = ((t139 + t207) / t254 * t178) ^ (-0.1e1 / 0.2e1) * t177 * t198;
t226 = (-t161 - t162) * qJD(2) * pkin(2) * t247 / t205 * t177;
t184 = qJD(2) * t194 + V_base(4);
t223 = t226 / 0.2e1;
t222 = rSges(3,1) * t195 - rSges(3,2) * t193;
t221 = -V_base(4) * pkin(5) + V_base(2);
t220 = Icges(3,1) * t195 - t244;
t219 = -Icges(3,2) * t193 + t243;
t218 = Icges(3,5) * t195 - Icges(3,6) * t193;
t196 = cos(qJ(1));
t211 = (-Icges(5,5) * t142 - Icges(5,6) * t143) * t228;
t138 = 0.1e1 / t207;
t100 = 0.2e1 * (-(t185 * t223 + (t203 * pkin(1) * t239 + ((t164 * t195 + t237) * t177 / 0.2e1 - t143 * t178 * t247) * pkin(2)) * qJD(2)) / t142 - (t193 * t223 + (-(-t235 + t238) * t177 / 0.2e1 + (t142 * t178 - t177 * t185) * t247) * qJD(2)) * pkin(2) * t143 * t138) * pkin(4) / (t138 * t139 + 0.1e1) * t182 * t198;
t98 = -t100 * t196 + V_base(5);
t99 = t100 * t194 + V_base(4);
t216 = (-Icges(5,3) * t196 + t194 * t211) * t98 + (Icges(5,3) * t194 + t196 * t211) * t99 + (Icges(5,5) * t143 - Icges(5,6) * t142) * t228 * t188;
t183 = -qJD(2) * t196 + V_base(5);
t215 = (-Icges(3,3) * t196 + t194 * t218) * t183 + (Icges(3,3) * t194 + t196 * t218) * t184 + (Icges(3,5) * t193 + Icges(3,6) * t195) * t188;
t214 = (-rSges(5,1) * t142 - rSges(5,2) * t143) * t228;
t213 = (-Icges(5,1) * t142 - t241) * t228;
t212 = (-Icges(5,2) * t143 - t242) * t228;
t127 = (-t141 * t195 + t240) * t227;
t155 = -Icges(3,6) * t196 + t194 * t219;
t156 = Icges(3,6) * t194 + t196 * t219;
t157 = -Icges(3,5) * t196 + t194 * t220;
t158 = Icges(3,5) * t194 + t196 * t220;
t171 = Icges(3,2) * t195 + t244;
t174 = Icges(3,1) * t193 + t243;
t209 = (-t156 * t193 + t158 * t195) * t184 + (-t155 * t193 + t157 * t195) * t183 + (-t171 * t193 + t174 * t195) * t188;
t116 = -Icges(5,6) * t196 + t194 * t212;
t117 = Icges(5,6) * t194 + t196 * t212;
t118 = -Icges(5,5) * t196 + t194 * t213;
t119 = Icges(5,5) * t194 + t196 * t213;
t129 = (-Icges(5,2) * t142 + t241) * t228;
t130 = (Icges(5,1) * t143 - t242) * t228;
t208 = ((-t117 * t143 - t119 * t142) * t99 + (-t116 * t143 - t118 * t142) * t98 + (-t129 * t143 - t130 * t142) * t188) * t228;
t190 = Icges(2,4) * t196;
t181 = rSges(2,1) * t196 - rSges(2,2) * t194;
t180 = rSges(2,1) * t194 + rSges(2,2) * t196;
t179 = rSges(3,1) * t193 + rSges(3,2) * t195;
t176 = Icges(2,1) * t196 - t245;
t175 = Icges(2,1) * t194 + t190;
t173 = -Icges(2,2) * t194 + t190;
t172 = Icges(2,2) * t196 + t245;
t167 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t166 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t165 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t160 = rSges(3,3) * t194 + t196 * t222;
t159 = -rSges(3,3) * t196 + t194 * t222;
t152 = V_base(5) * rSges(2,3) - t180 * t188 + t231;
t151 = t181 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t150 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t146 = -t159 * t188 + t179 * t183 + t231;
t145 = t160 * t188 - t179 * t184 + t221;
t137 = 0.1e1 / t206;
t136 = t159 * t184 - t160 * t183 + V_base(3);
t131 = (rSges(5,1) * t143 - rSges(5,2) * t142) * t228;
t125 = t196 * t252;
t124 = t196 * t127;
t123 = t194 * t252;
t122 = t194 * t127;
t121 = rSges(5,3) * t194 + t196 * t214;
t120 = -rSges(5,3) * t196 + t194 * t214;
t113 = -rSges(4,1) * t252 + rSges(4,2) * t127;
t112 = -Icges(4,1) * t252 + Icges(4,4) * t127;
t111 = -Icges(4,4) * t252 + Icges(4,2) * t127;
t110 = -Icges(4,5) * t252 + Icges(4,6) * t127;
t109 = rSges(4,1) * t124 + rSges(4,2) * t125 + rSges(4,3) * t194;
t108 = rSges(4,1) * t122 + rSges(4,2) * t123 - rSges(4,3) * t196;
t107 = Icges(4,1) * t124 + Icges(4,4) * t125 + Icges(4,5) * t194;
t106 = Icges(4,1) * t122 + Icges(4,4) * t123 - Icges(4,5) * t196;
t105 = Icges(4,4) * t124 + Icges(4,2) * t125 + Icges(4,6) * t194;
t104 = Icges(4,4) * t122 + Icges(4,2) * t123 - Icges(4,6) * t196;
t103 = Icges(4,5) * t124 + Icges(4,6) * t125 + Icges(4,3) * t194;
t102 = Icges(4,5) * t122 + Icges(4,6) * t123 - Icges(4,3) * t196;
t101 = ((-t186 * t226 + (0.2e1 * t204 * pkin(2) * t239 + ((t163 * t195 + t237) * t177 + t232 * t240) * pkin(1)) * qJD(2)) / t141 - (-t193 * t226 + (-t177 * t235 + ((t186 * t251 + t163) * t177 + t141 * t232) * t193) * qJD(2)) * pkin(1) * t144 * t137) * pkin(3) / (t137 * t140 + 0.1e1) * t182 * t201;
t97 = t101 * t194 + t184;
t96 = V_base(5) + (-qJD(2) - t101) * t196;
t95 = t131 * t98 + (-pkin(1) * t194 - t120) * t188 + t231;
t94 = -t131 * t99 + (pkin(1) * t196 + t121) * t188 + t221;
t93 = -t108 * t188 + t113 * t96 + (t193 * t183 - t194 * t236) * pkin(2) + t231;
t92 = t109 * t188 - t113 * t97 + (-t193 * t184 + t196 * t236) * pkin(2) + t221;
t91 = t120 * t99 - t121 * t98 + V_base(3) + (t194 * V_base(4) - V_base(5) * t196) * pkin(1);
t90 = t108 * t97 - t109 * t96 + V_base(3) + (-t196 * t183 + t194 * t184) * t246;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + t184 * (t194 * t215 + t196 * t209) / 0.2e1 + t183 * (t194 * t209 - t196 * t215) / 0.2e1 + m(4) * (t90 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t97 * ((t194 * t103 + t125 * t105 + t124 * t107) * t97 + (t102 * t194 + t104 * t125 + t106 * t124) * t96 + (t110 * t194 + t111 * t125 + t112 * t124) * t188) / 0.2e1 + t96 * ((-t103 * t196 + t105 * t123 + t107 * t122) * t97 + (-t196 * t102 + t123 * t104 + t122 * t106) * t96 + (-t110 * t196 + t111 * t123 + t112 * t122) * t188) / 0.2e1 + m(5) * (t91 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + t99 * (t194 * t216 + t196 * t208) / 0.2e1 + t98 * (t194 * t208 - t216 * t196) / 0.2e1 + ((-t172 * t194 + t175 * t196 + Icges(1,4)) * V_base(5) + (-t194 * t173 + t196 * t176 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t196 * t172 + t194 * t175 + Icges(1,2)) * V_base(5) + (t173 * t196 + t176 * t194 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t156 * t195 + t158 * t193) * t184 + (t155 * t195 + t157 * t193) * t183 + (t105 * t127 - t107 * t252) * t97 + (t104 * t127 - t106 * t252) * t96 + ((-t117 * t142 + t119 * t143) * t99 + (-t116 * t142 + t118 * t143) * t98) * t228 + (Icges(2,3) + t195 * t171 + t193 * t174 + t127 * t111 - t252 * t112 + (-t129 * t142 + t130 * t143) * t228) * t188) * t188 / 0.2e1 + V_base(4) * t188 * (Icges(2,5) * t196 - Icges(2,6) * t194) + V_base(5) * t188 * (Icges(2,5) * t194 + Icges(2,6) * t196) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
