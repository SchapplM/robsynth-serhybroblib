% Calculate kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:00:28
% EndTime: 2020-05-03 00:00:31
% DurationCPUTime: 2.95s
% Computational Cost: add. (1299->290), mult. (1401->444), div. (0->0), fcn. (1205->12), ass. (0->157)
t178 = V_base(6) + qJD(1);
t193 = cos(qJ(3));
t176 = pkin(3) * t193 + pkin(2);
t190 = sin(qJ(2));
t194 = cos(qJ(2));
t189 = sin(qJ(3));
t252 = pkin(3) * t189;
t260 = t178 * (-t176 * t194 + t190 * t252 - pkin(1));
t195 = cos(qJ(1));
t181 = qJD(2) * t195;
t171 = V_base(5) + t181;
t191 = sin(qJ(1));
t219 = -rSges(3,1) * t194 + rSges(3,2) * t190;
t205 = pkin(1) - t219;
t258 = t195 * rSges(3,3) + t205 * t191;
t257 = -t191 * rSges(3,3) + t205 * t195;
t254 = pkin(2) * t190;
t253 = pkin(2) * t194;
t249 = Icges(2,4) * t191;
t248 = Icges(3,4) * t190;
t247 = Icges(3,4) * t194;
t187 = qJ(2) + qJ(3);
t182 = sin(t187);
t246 = Icges(4,4) * t182;
t183 = cos(t187);
t245 = Icges(4,4) * t183;
t185 = qJ(4) + t187;
t174 = sin(t185);
t244 = Icges(5,4) * t174;
t175 = cos(t185);
t243 = Icges(5,4) * t175;
t241 = t174 * t191;
t240 = t174 * t195;
t239 = t176 * t190;
t188 = sin(qJ(5));
t238 = t188 * t195;
t237 = t189 * t194;
t236 = t191 * t188;
t192 = cos(qJ(5));
t235 = t191 * t192;
t234 = t192 * t195;
t180 = qJD(3) * t195;
t232 = qJD(5) * t174;
t231 = -qJD(2) - qJD(3);
t230 = V_base(5) * pkin(5) + V_base(1);
t177 = pkin(1) + t253;
t227 = t191 * V_base(4);
t226 = t194 * V_base(4);
t225 = V_base(5) * t195;
t224 = t191 * pkin(2) * t226 + pkin(1) * t227 + V_base(3);
t223 = t178 * t177;
t146 = t180 + t171;
t221 = pkin(4) * t175 + pkin(6) * t174;
t172 = -qJD(2) * t191 + V_base(4);
t220 = -qJD(2) - t225;
t218 = rSges(4,1) * t183 - rSges(4,2) * t182;
t217 = rSges(5,1) * t175 - rSges(5,2) * t174;
t167 = rSges(6,1) * t192 - rSges(6,2) * t188;
t216 = -V_base(4) * pkin(5) + V_base(2);
t215 = rSges(6,3) * t174 + t167 * t175;
t130 = qJD(4) * t195 + t146;
t214 = Icges(3,1) * t194 - t248;
t213 = Icges(4,1) * t183 - t246;
t212 = Icges(5,1) * t175 - t244;
t211 = -Icges(3,2) * t190 + t247;
t210 = -Icges(4,2) * t182 + t245;
t209 = -Icges(5,2) * t174 + t243;
t208 = Icges(3,5) * t194 - Icges(3,6) * t190;
t207 = Icges(4,5) * t183 - Icges(4,6) * t182;
t206 = Icges(5,5) * t175 - Icges(5,6) * t174;
t131 = V_base(4) + (-qJD(4) + t231) * t191;
t204 = (Icges(5,3) * t195 + t206 * t191) * t130 + (-Icges(5,3) * t191 + t206 * t195) * t131 + (-Icges(5,5) * t174 - Icges(5,6) * t175) * t178;
t147 = t231 * t191 + V_base(4);
t203 = (Icges(4,3) * t195 + t207 * t191) * t146 + (-Icges(4,3) * t191 + t207 * t195) * t147 + (-Icges(4,5) * t182 - Icges(4,6) * t183) * t178;
t202 = (Icges(3,3) * t195 + t208 * t191) * t171 + (-Icges(3,3) * t191 + t208 * t195) * t172 + (-Icges(3,5) * t190 - Icges(3,6) * t194) * t178;
t201 = -qJD(2) * t253 + pkin(3) * (-qJD(3) + t220 + t227) * t183 - t177 * t225 + t224;
t133 = pkin(3) * t237 + t239;
t149 = t190 * t193 + t237;
t200 = V_base(4) * t239 + t226 * t252 + (-pkin(3) * qJD(3) * t149 - qJD(2) * t133) * t191 + t216 - t260 * t195;
t102 = Icges(5,6) * t195 + t209 * t191;
t103 = -Icges(5,6) * t191 + t209 * t195;
t104 = Icges(5,5) * t195 + t212 * t191;
t105 = -Icges(5,5) * t191 + t212 * t195;
t135 = -Icges(5,2) * t175 - t244;
t136 = -Icges(5,1) * t174 - t243;
t199 = (-t103 * t174 + t105 * t175) * t131 + (-t102 * t174 + t104 * t175) * t130 + (-t135 * t174 + t136 * t175) * t178;
t110 = Icges(4,6) * t195 + t210 * t191;
t111 = -Icges(4,6) * t191 + t210 * t195;
t112 = Icges(4,5) * t195 + t213 * t191;
t113 = -Icges(4,5) * t191 + t213 * t195;
t143 = -Icges(4,2) * t183 - t246;
t144 = -Icges(4,1) * t182 - t245;
t198 = (-t111 * t182 + t113 * t183) * t147 + (-t110 * t182 + t112 * t183) * t146 + (-t143 * t182 + t144 * t183) * t178;
t122 = Icges(3,6) * t195 + t211 * t191;
t123 = -Icges(3,6) * t191 + t211 * t195;
t124 = Icges(3,5) * t195 + t214 * t191;
t125 = -Icges(3,5) * t191 + t214 * t195;
t158 = -Icges(3,2) * t194 - t248;
t161 = -Icges(3,1) * t190 - t247;
t197 = (-t123 * t190 + t125 * t194) * t172 + (-t122 * t190 + t124 * t194) * t171 + (-t158 * t190 + t161 * t194) * t178;
t196 = -V_base(5) * t239 + (-t149 * t180 - V_base(5) * t237) * pkin(3) - t133 * t181 + t230 + t260 * t191;
t184 = Icges(2,4) * t195;
t168 = rSges(2,1) * t195 - t191 * rSges(2,2);
t166 = rSges(6,1) * t188 + rSges(6,2) * t192;
t165 = t191 * rSges(2,1) + rSges(2,2) * t195;
t164 = -rSges(3,1) * t190 - rSges(3,2) * t194;
t163 = Icges(2,1) * t195 - t249;
t162 = Icges(2,1) * t191 + t184;
t160 = -Icges(2,2) * t191 + t184;
t159 = Icges(2,2) * t195 + t249;
t154 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t153 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t152 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t148 = qJD(5) * t175 + t178;
t145 = -rSges(4,1) * t182 - rSges(4,2) * t183;
t138 = -pkin(4) * t174 + pkin(6) * t175;
t137 = -rSges(5,1) * t174 - rSges(5,2) * t175;
t129 = t175 * t234 - t236;
t128 = -t175 * t238 - t235;
t127 = t175 * t235 + t238;
t126 = -t175 * t236 + t234;
t119 = t221 * t195;
t118 = t221 * t191;
t115 = -t191 * rSges(4,3) + t218 * t195;
t114 = rSges(4,3) * t195 + t218 * t191;
t107 = -t191 * rSges(5,3) + t217 * t195;
t106 = rSges(5,3) * t195 + t217 * t191;
t99 = V_base(5) * rSges(2,3) - t165 * t178 + t230;
t98 = t168 * t178 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t195 * t232 + t131;
t96 = t191 * t232 + t130;
t95 = t165 * V_base(4) - t168 * V_base(5) + V_base(3);
t94 = rSges(6,3) * t175 - t167 * t174;
t93 = Icges(6,5) * t175 + (-Icges(6,1) * t192 + Icges(6,4) * t188) * t174;
t92 = Icges(6,6) * t175 + (-Icges(6,4) * t192 + Icges(6,2) * t188) * t174;
t91 = Icges(6,3) * t175 + (-Icges(6,5) * t192 + Icges(6,6) * t188) * t174;
t88 = -t191 * t166 + t215 * t195;
t87 = t166 * t195 + t215 * t191;
t86 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t240;
t85 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t241;
t84 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t240;
t83 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t241;
t82 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t240;
t81 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t241;
t80 = t171 * t164 - t258 * t178 + t230;
t79 = -t172 * t164 + t257 * t178 + t216;
t78 = t219 * qJD(2) - t257 * V_base(5) + t258 * V_base(4) + V_base(3);
t77 = t178 * t115 - t147 * t145 + t172 * t254 + t223 * t195 + t216;
t76 = -t178 * t114 + t146 * t145 - t171 * t254 - t223 * t191 + t230;
t75 = -pkin(1) * t225 + t147 * t114 - t146 * t115 + t220 * t253 + t224;
t74 = t106 * t131 - t107 * t130 + t201;
t73 = t107 * t178 - t131 * t137 + t200;
t72 = -t178 * t106 + t130 * t137 + t196;
t71 = t118 * t131 - t119 * t130 + t87 * t97 - t88 * t96 + t201;
t70 = t119 * t178 - t131 * t138 + t148 * t88 - t94 * t97 + t200;
t69 = -t178 * t118 + t130 * t138 - t148 * t87 + t96 * t94 + t196;
t1 = m(1) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(2) * (t95 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t172 * (-t202 * t191 + t197 * t195) / 0.2e1 + t171 * (t197 * t191 + t202 * t195) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t147 * (-t203 * t191 + t198 * t195) / 0.2e1 + t146 * (t198 * t191 + t203 * t195) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t131 * (-t204 * t191 + t199 * t195) / 0.2e1 + t130 * (t199 * t191 + t204 * t195) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t97 * ((t128 * t84 + t129 * t86 + t82 * t240) * t97 + (t128 * t83 + t129 * t85 + t81 * t240) * t96 + (t128 * t92 + t129 * t93 + t91 * t240) * t148) / 0.2e1 + t96 * ((t126 * t84 + t127 * t86 + t82 * t241) * t97 + (t126 * t83 + t127 * t85 + t81 * t241) * t96 + (t126 * t92 + t127 * t93 + t91 * t241) * t148) / 0.2e1 + t148 * ((t148 * t91 + t81 * t96 + t82 * t97) * t175 + ((t188 * t84 - t192 * t86) * t97 + (t188 * t83 - t192 * t85) * t96 + (t188 * t92 - t192 * t93) * t148) * t174) / 0.2e1 + ((-t191 * t159 + t162 * t195 + Icges(1,4)) * V_base(5) + (-t191 * t160 + t195 * t163 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t195 * t159 + t191 * t162 + Icges(1,2)) * V_base(5) + (t160 * t195 + t191 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t123 * t194 - t125 * t190) * t172 + (-t122 * t194 - t124 * t190) * t171 + (-t111 * t183 - t113 * t182) * t147 + (-t110 * t183 - t112 * t182) * t146 + (-t103 * t175 - t105 * t174) * t131 + (-t102 * t175 - t104 * t174) * t130 + (-t175 * t135 - t174 * t136 - t183 * t143 - t182 * t144 - t194 * t158 - t190 * t161 + Icges(2,3)) * t178) * t178 / 0.2e1 + V_base(4) * t178 * (Icges(2,5) * t195 - Icges(2,6) * t191) + V_base(5) * t178 * (Icges(2,5) * t191 + Icges(2,6) * t195) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
