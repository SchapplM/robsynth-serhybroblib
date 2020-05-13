% Calculate kinetic energy for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:44:52
% EndTime: 2020-05-07 09:44:56
% DurationCPUTime: 3.88s
% Computational Cost: add. (1420->220), mult. (1476->361), div. (166->8), fcn. (100->2), ass. (0->124)
t228 = Icges(3,4) + Icges(2,6);
t122 = qJ(1) + pkin(3);
t111 = 0.1e1 / t122;
t198 = t111 / pkin(1);
t164 = V_base(5) * t198;
t188 = -pkin(2) + t122;
t189 = -pkin(2) - t122;
t167 = (pkin(1) + t188) * (pkin(1) - t189) * (pkin(1) - t188) * (pkin(1) + t189);
t141 = sqrt(-t167);
t173 = -t198 / 0.2e1;
t163 = t141 * t173;
t135 = (pkin(3) ^ 2);
t130 = (qJ(1) ^ 2);
t207 = (pkin(3) * qJ(1));
t170 = t130 + 2 * t207;
t161 = -t135 - t170;
t136 = (pkin(2) ^ 2);
t138 = (pkin(1) ^ 2);
t192 = (t138 - t136);
t88 = -t161 + t192;
t205 = Icges(3,5) * t88;
t211 = -t141 / 0.2e1;
t166 = t88 * t173;
t81 = Icges(2,4) * t166;
t227 = (t205 / 0.2e1 + Icges(3,3) * t211) * t198 + Icges(2,1) * t163 + t81;
t201 = Icges(2,4) * t141;
t217 = -t88 / 0.2e1;
t78 = Icges(3,5) * t163;
t226 = Icges(3,3) * t166 + (Icges(2,1) * t217 + t201 / 0.2e1) * t198 + t78;
t210 = t141 / 0.2e1;
t216 = t88 / 0.2e1;
t225 = (Icges(2,5) * t217 + Icges(3,6) * t216 + t228 * t210) * t198;
t137 = t138 ^ 2;
t131 = rSges(4,2) ^ 2;
t132 = rSges(4,1) ^ 2;
t157 = (t132 / 0.4e1 - t131 / 0.4e1) * m(4) + Icges(4,2) / 0.4e1 - Icges(4,1) / 0.4e1;
t123 = pkin(2) + pkin(3);
t124 = pkin(2) - pkin(3);
t197 = t124 ^ 2 * t123 ^ 2;
t91 = (-t131 + t132) * m(4) + Icges(4,2) - Icges(4,1);
t223 = -t91 * t197 / 0.4e1 - t157 * t137;
t222 = -3 * t135;
t221 = t91 / 0.4e1;
t193 = (t136 + t138);
t175 = t135 - t193;
t220 = -t175 - t170;
t219 = pkin(1) * pkin(2);
t215 = (pkin(3) * rSges(3,3));
t214 = m(4) * rSges(4,3);
t156 = t138 + t161;
t213 = m(4) * (-t136 - t156);
t212 = t136 / 0.2e1;
t121 = t138 / 0.2e1;
t144 = t135 ^ 2;
t209 = t144 / 0.2e1;
t208 = pkin(1) * t122;
t206 = t138 * pkin(3);
t204 = t135 * rSges(3,3);
t203 = t91 * t135;
t115 = (V_base(4) ^ 2);
t202 = t175 * t115;
t200 = qJD(1) * (-t136 + t156);
t125 = pkin(1) + pkin(2);
t126 = pkin(1) - pkin(2);
t196 = t126 ^ 2 * t125 ^ 2;
t195 = t138 * t136;
t194 = (-t138 / 0.2e1 + t212) * rSges(3,3);
t114 = (V_base(5) ^ 2);
t190 = pkin(1) * V_base(3);
t186 = V_base(4) / 0.2e1;
t184 = V_base(6) / 0.2e1;
t183 = V_base(6) * pkin(1);
t182 = rSges(3,1) * t200;
t92 = (t131 + t132) * m(4) + Icges(4,3);
t180 = t92 * t195;
t179 = 0.2e1 * rSges(4,3) * t219;
t94 = -rSges(4,1) * t214 + Icges(4,5);
t178 = t94 * V_base(5);
t95 = -m(4) * rSges(4,1) * rSges(4,2) + Icges(4,4);
t177 = t95 * V_base(5);
t176 = t175 * V_base(5);
t174 = rSges(3,2) * t122;
t172 = t198 / 0.2e1;
t171 = -t135 / 0.2e1 + t212;
t169 = pkin(2) * t183;
t168 = 0.1e1 / t167;
t165 = V_base(4) * t198;
t162 = V_base(4) * V_base(5);
t93 = rSges(4,2) * t214 - Icges(4,6);
t160 = t93 * t169;
t159 = t94 * t169;
t158 = t122 * t167;
t155 = -rSges(4,1) * V_base(1) - rSges(4,2) * V_base(2);
t154 = V_base(5) * t160;
t153 = t220 * V_base(6);
t151 = t158 * t184;
t149 = t122 ^ 2;
t147 = t130 ^ 2;
t134 = pkin(3) * t135;
t129 = qJ(1) * t130;
t127 = rSges(3,3) / 0.2e1;
t120 = rSges(3,3) + qJ(1);
t113 = 2 * t215;
t98 = t222 + t193;
t87 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t86 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t85 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t83 = t95 * t176;
t82 = ((t132 / 0.2e1 + t131 / 0.2e1 + rSges(4,3) ^ 2) * m(4) + Icges(4,1) / 0.2e1 + Icges(4,2) / 0.2e1) * t136;
t77 = -t176 * t91 - 0.2e1 * t160;
t75 = -rSges(2,1) * t141 - (rSges(2,2) * t88);
t74 = -rSges(2,1) * t88 + rSges(2,2) * t141;
t73 = V_base(6) + t111 / t141 * t200;
t69 = (Icges(3,1) * t211 - t205 / 0.2e1) * t198;
t68 = Icges(3,1) * t172 * t88 + t78;
t67 = (-t201 / 0.2e1 + Icges(2,2) * t217) * t198;
t66 = Icges(2,2) * t141 * t172 + t81;
t59 = 0.2e1 * t122 * t190 + ((rSges(3,1) * t88) - t120 * t141) * V_base(4) + (rSges(3,1) * t141 + t129 + (2 * pkin(3) * t130) + ((t113 + t175) * qJ(1)) - (2 * t206) + ((t130 + t135 + t192) * rSges(3,3))) * V_base(5);
t58 = t173 * t73 * t75 + V_base(5) * rSges(2,3) + V_base(1);
t57 = t172 * t73 * t74 - V_base(4) * rSges(2,3) + t183 + V_base(2);
t56 = -V_base(5) * pkin(1) + V_base(3) + (t75 * t186 - V_base(5) * t74 / 0.2e1) * t198;
t55 = ((t129 + (0.5e1 / 0.2e1 * pkin(3) + t127) * t130 + ((2 * t135 + t215) * qJ(1)) + t134 / 0.2e1 + t204 / 0.2e1 + (t121 - t136 / 0.2e1) * pkin(3) + t194) * qJD(1) + ((t130 / 0.2e1 + t207 + t121 - t171) * V_base(6) * rSges(3,1) + (-t122 * V_base(1) + t174 * V_base(5)) * pkin(1)) * t122) * t141 + t182 * t216 + t120 * t151;
t54 = (-t182 / 0.2e1 + (t122 * V_base(2) + t174 * V_base(4)) * t208 + (-t129 / 0.2e1 - (pkin(3) + t127) * t130 + (t121 + t171 - t215) * qJ(1) + t206 - t204 / 0.2e1 + t194) * t122 * V_base(6)) * t141 + (t127 * t147 + ((t113 + 8 * t135 - t193) * t129) + ((3 * t135 * t130) + t209 - t196 / 0.2e1) * rSges(3,3) + ((7 * t130 - t193) * t134) + ((2 * rSges(3,3) * t134 + t193 * t222 + 3 * t144 + t147) * qJ(1)) + (0.9e1 / 0.2e1 * t147 - (3 * t193 * t130) + t209 + t196 / 0.2e1) * pkin(3)) * qJD(1) + rSges(3,1) * t151;
t1 = m(1) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6)) * t184 + m(2) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t225 * V_base(4) + (Icges(3,2) + Icges(2,3)) * t73) * t73 / 0.2e1 + (m(3) * (t59 ^ 2 / t149 / 0.4e1 + (-t54 ^ 2 - t55 ^ 2) / t149 ^ 2 * t168) / 0.2e1 + (-0.2e1 * (t95 * t202 - t77 * V_base(4) - t83 * V_base(5) + t170 * (t91 * t162 + (-t114 + t115) * t95) + 0.2e1 * (V_base(6) * t178 + ((V_base(4) * rSges(4,1) + V_base(5) * rSges(4,2)) * V_base(3) + t155 * V_base(6)) * m(4)) * t219) * t167 * t141 + 0.4e1 * (-pkin(2) * (-V_base(5) * rSges(4,1) + V_base(4) * rSges(4,2)) * t190 * t213 + ((-0.2e1 * t177 * t98 - t159) * V_base(4) + t154 + (t115 / 0.2e1 - t114 / 0.2e1) * t98 * t91) * t130 - 0.4e1 * (t202 * t221 + (-t83 + t159 / 0.2e1) * V_base(4) + t77 * V_base(5) / 0.4e1) * t207 + ((t82 + t203 / 0.2e1) * t138 + t223) * t115 + ((-(2 * t138 * t135) + t137 + t197) * t177 - t175 * t159) * V_base(4) + ((t82 - t203 / 0.2e1) * t138 - t223) * t114 + t175 * t154 + V_base(6) ^ 2 * t180 + (4 * pkin(3) * t129 + t147) * (t114 * t221 - t115 * t157 + t95 * t162) + (((-rSges(4,2) * t153 + t179 * V_base(5)) * V_base(1) - (-rSges(4,1) * t153 + t179 * V_base(4)) * V_base(2)) * t219 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * t195) * m(4)) * t167 + 0.8e1 * (-(-0.2e1 * t92 * t169 + (rSges(4,1) * V_base(2) - rSges(4,2) * V_base(1)) * t213 + t220 * (t93 * V_base(5) - t94 * V_base(4))) * pkin(2) * t208 * t141 + (m(4) * t155 + t93 * V_base(4) + t178) * t158 * t219 - 0.2e1 * t149 * t180 * qJD(1)) * qJD(1)) / pkin(2) ^ 2 * t168 / 0.8e1) / pkin(1) ^ 2 + (Icges(1,1) * V_base(4) + Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + t225 * t73 + (t210 * t66 + t211 * t69 + t226 * t217) * t165 + (t210 * t67 + t211 * t68 + t227 * t217) * t164) * t186 + (Icges(1,4) * V_base(4) + Icges(1,2) * V_base(5) + Icges(1,6) * V_base(6) + (t226 * t211 + t216 * t69 + t217 * t66) * t165 + (t227 * t211 + t216 * t68 + t217 * t67) * t164) * V_base(5) / 0.2e1 + t73 * (Icges(2,5) * t211 + Icges(3,6) * t210 + t228 * t217) * t164;
T = t1;
