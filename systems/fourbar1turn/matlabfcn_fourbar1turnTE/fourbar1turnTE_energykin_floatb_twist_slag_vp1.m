% Calculate kinetic energy for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnTE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:29
% EndTime: 2020-04-12 19:18:31
% DurationCPUTime: 2.54s
% Computational Cost: add. (12162->244), mult. (17172->421), div. (804->11), fcn. (5042->6), ass. (0->132)
t190 = pkin(2) ^ 2;
t191 = pkin(1) ^ 2;
t184 = cos(qJ(2));
t231 = pkin(2) * t184;
t218 = -0.2e1 * pkin(1) * t231 + t191;
t171 = t190 + t218;
t217 = pkin(3) ^ 2 - pkin(4) ^ 2;
t152 = t171 + t217;
t175 = pkin(1) * t184 - pkin(2);
t182 = sin(qJ(2));
t238 = -pkin(3) - pkin(4);
t150 = (pkin(2) - t238) * (pkin(2) + t238) + t218;
t237 = -pkin(3) + pkin(4);
t151 = (pkin(2) - t237) * (pkin(2) + t237) + t218;
t192 = sqrt(-t150 * t151);
t221 = t182 * t192;
t130 = -pkin(1) * t221 - t152 * t175;
t232 = pkin(1) * t182;
t133 = t152 * t232 - t175 * t192;
t166 = 0.1e1 / t171;
t189 = 0.1e1 / pkin(3);
t224 = t166 * t189;
t240 = (t184 * t133 / 0.2e1 + t182 * t130 / 0.2e1) * t224;
t239 = -0.2e1 * pkin(2);
t153 = t171 - t217;
t174 = pkin(1) - t231;
t131 = -pkin(2) * t221 + t153 * t174;
t236 = -t131 / 0.2e1;
t222 = t182 * t153;
t132 = pkin(2) * t222 + t174 * t192;
t235 = -t132 / 0.2e1;
t234 = t132 / 0.2e1;
t183 = sin(qJ(1));
t230 = Icges(2,4) * t183;
t229 = Icges(3,4) * t182;
t228 = Icges(3,4) * t184;
t227 = Icges(5,4) * t132;
t226 = t166 * t182 ^ 2;
t187 = 0.1e1 / pkin(4);
t225 = t166 * t187;
t223 = t182 * t133;
t177 = V_base(6) + qJD(1);
t220 = t184 * t177;
t219 = t184 * t192;
t167 = 0.1e1 / t171 ^ 2;
t216 = t167 * t239;
t215 = V_base(5) * pkin(5) + V_base(1);
t212 = (-t150 - t151) * qJD(2) * pkin(2) * t232 / t192 * t166;
t211 = Icges(5,4) * t236;
t173 = qJD(2) * t183 + V_base(4);
t208 = t212 / 0.2e1;
t207 = rSges(3,1) * t184 - rSges(3,2) * t182;
t206 = -V_base(4) * pkin(5) + V_base(2);
t205 = Icges(3,1) * t184 - t229;
t204 = -Icges(3,2) * t182 + t228;
t203 = Icges(3,5) * t184 - Icges(3,6) * t182;
t185 = cos(qJ(1));
t196 = (Icges(5,5) * t236 + Icges(5,6) * t235) * t225;
t129 = 0.1e1 / t131 ^ 2;
t93 = 0.2e1 * (-(t174 * t208 + (t190 * pkin(1) * t226 + ((t153 * t184 + t221) * t166 / 0.2e1 - t132 * t167 * t232) * pkin(2)) * qJD(2)) / t131 - (t182 * t208 + (-(-t219 + t222) * t166 / 0.2e1 + (t131 * t167 - t166 * t174) * t232) * qJD(2)) * pkin(2) * t132 * t129) * pkin(4) / (t129 * t132 ^ 2 + 0.1e1) * t171 * t187;
t91 = -t185 * t93 + V_base(5);
t92 = t183 * t93 + V_base(4);
t201 = (-Icges(5,3) * t185 + t183 * t196) * t91 + (Icges(5,3) * t183 + t185 * t196) * t92 + (Icges(5,5) * t234 + Icges(5,6) * t236) * t225 * t177;
t172 = -qJD(2) * t185 + V_base(5);
t200 = (-Icges(3,3) * t185 + t183 * t203) * t172 + (Icges(3,3) * t183 + t185 * t203) * t173 + (Icges(3,5) * t182 + Icges(3,6) * t184) * t177;
t199 = (rSges(5,1) * t236 + rSges(5,2) * t235) * t225;
t198 = (Icges(5,1) * t236 - t227 / 0.2e1) * t225;
t197 = (Icges(5,2) * t235 + t211) * t225;
t120 = (t223 / 0.2e1 - t184 * t130 / 0.2e1) * t224;
t144 = -Icges(3,6) * t185 + t183 * t204;
t145 = Icges(3,6) * t183 + t185 * t204;
t146 = -Icges(3,5) * t185 + t183 * t205;
t147 = Icges(3,5) * t183 + t185 * t205;
t160 = Icges(3,2) * t184 + t229;
t163 = Icges(3,1) * t182 + t228;
t194 = (-t145 * t182 + t147 * t184) * t173 + (-t144 * t182 + t146 * t184) * t172 + (-t160 * t182 + t163 * t184) * t177;
t109 = -Icges(5,6) * t185 + t183 * t197;
t110 = Icges(5,6) * t183 + t185 * t197;
t111 = -Icges(5,5) * t185 + t183 * t198;
t112 = Icges(5,5) * t183 + t185 * t198;
t122 = (t227 / 0.2e1 + Icges(5,2) * t236) * t225;
t123 = (Icges(5,1) * t234 + t211) * t225;
t193 = ((t110 * t235 + t112 * t236) * t92 + (t109 * t235 + t111 * t236) * t91 + (t122 * t235 + t123 * t236) * t177) * t225;
t179 = Icges(2,4) * t185;
t170 = rSges(2,1) * t185 - rSges(2,2) * t183;
t169 = rSges(2,1) * t183 + rSges(2,2) * t185;
t168 = rSges(3,1) * t182 + rSges(3,2) * t184;
t165 = Icges(2,1) * t185 - t230;
t164 = Icges(2,1) * t183 + t179;
t162 = -Icges(2,2) * t183 + t179;
t161 = Icges(2,2) * t185 + t230;
t156 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t155 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t154 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t149 = rSges(3,3) * t183 + t185 * t207;
t148 = -rSges(3,3) * t185 + t183 * t207;
t141 = V_base(5) * rSges(2,3) - t169 * t177 + t215;
t140 = t170 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t139 = t169 * V_base(4) - t170 * V_base(5) + V_base(3);
t135 = -t148 * t177 + t168 * t172 + t215;
t134 = t149 * t177 - t168 * t173 + t206;
t128 = 0.1e1 / t130 ^ 2;
t127 = t148 * t173 - t149 * t172 + V_base(3);
t124 = (rSges(5,1) * t234 + rSges(5,2) * t236) * t225;
t118 = t185 * t240;
t117 = t185 * t120;
t116 = t183 * t240;
t115 = t183 * t120;
t114 = t183 * rSges(5,3) + t185 * t199;
t113 = -t185 * rSges(5,3) + t183 * t199;
t106 = -rSges(4,1) * t240 + rSges(4,2) * t120;
t105 = -Icges(4,1) * t240 + Icges(4,4) * t120;
t104 = -Icges(4,4) * t240 + Icges(4,2) * t120;
t103 = -Icges(4,5) * t240 + Icges(4,6) * t120;
t102 = rSges(4,1) * t117 + rSges(4,2) * t118 + rSges(4,3) * t183;
t101 = rSges(4,1) * t115 + rSges(4,2) * t116 - rSges(4,3) * t185;
t100 = Icges(4,1) * t117 + Icges(4,4) * t118 + Icges(4,5) * t183;
t99 = Icges(4,1) * t115 + Icges(4,4) * t116 - Icges(4,5) * t185;
t98 = Icges(4,4) * t117 + Icges(4,2) * t118 + Icges(4,6) * t183;
t97 = Icges(4,4) * t115 + Icges(4,2) * t116 - Icges(4,6) * t185;
t96 = Icges(4,5) * t117 + Icges(4,6) * t118 + Icges(4,3) * t183;
t95 = Icges(4,5) * t115 + Icges(4,6) * t116 - Icges(4,3) * t185;
t94 = ((-t175 * t212 + (0.2e1 * t191 * pkin(2) * t226 + ((t152 * t184 + t221) * t166 + t216 * t223) * pkin(1)) * qJD(2)) / t130 - (-t182 * t212 + (-t166 * t219 + ((t175 * t239 + t152) * t166 + t130 * t216) * t182) * qJD(2)) * pkin(1) * t133 * t128) * pkin(3) / (t128 * t133 ^ 2 + 0.1e1) * t171 * t189;
t90 = t183 * t94 + t173;
t89 = V_base(5) + (-qJD(2) - t94) * t185;
t88 = t124 * t91 + (-pkin(1) * t183 - t113) * t177 + t215;
t87 = -t124 * t92 + (pkin(1) * t185 + t114) * t177 + t206;
t86 = -t101 * t177 + t106 * t89 + (t182 * t172 - t183 * t220) * pkin(2) + t215;
t85 = t102 * t177 - t106 * t90 + (-t182 * t173 + t185 * t220) * pkin(2) + t206;
t84 = t113 * t92 - t114 * t91 + V_base(3) + (t183 * V_base(4) - V_base(5) * t185) * pkin(1);
t83 = t101 * t90 - t102 * t89 + V_base(3) + (-t185 * t172 + t183 * t173) * t231;
t1 = m(1) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(2) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t173 * (t183 * t200 + t185 * t194) / 0.2e1 + t172 * (t183 * t194 - t185 * t200) / 0.2e1 + m(4) * (t83 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t90 * ((t117 * t100 + t118 * t98 + t183 * t96) * t90 + (t117 * t99 + t118 * t97 + t183 * t95) * t89 + (t103 * t183 + t104 * t118 + t105 * t117) * t177) / 0.2e1 + t89 * ((t100 * t115 + t116 * t98 - t185 * t96) * t90 + (t115 * t99 + t116 * t97 - t185 * t95) * t89 + (-t103 * t185 + t104 * t116 + t105 * t115) * t177) / 0.2e1 + m(5) * (t84 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t92 * (t183 * t201 + t185 * t193) / 0.2e1 + t91 * (t183 * t193 - t201 * t185) / 0.2e1 + ((-t161 * t183 + t164 * t185 + Icges(1,4)) * V_base(5) + (-t183 * t162 + t185 * t165 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t185 * t161 + t183 * t164 + Icges(1,2)) * V_base(5) + (t162 * t185 + t165 * t183 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t184 + t147 * t182) * t173 + (t144 * t184 + t146 * t182) * t172 + (-t100 * t240 + t120 * t98) * t90 + (t120 * t97 - t240 * t99) * t89 + ((t110 * t236 + t112 * t234) * t92 + (t109 * t236 + t111 * t234) * t91) * t225 + (Icges(2,3) + t184 * t160 + t182 * t163 + t104 * t120 - t105 * t240 + (t122 * t236 + t123 * t234) * t225) * t177) * t177 / 0.2e1 + V_base(4) * t177 * (Icges(2,5) * t185 - Icges(2,6) * t183) + V_base(5) * t177 * (Icges(2,5) * t183 + Icges(2,6) * t185) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
