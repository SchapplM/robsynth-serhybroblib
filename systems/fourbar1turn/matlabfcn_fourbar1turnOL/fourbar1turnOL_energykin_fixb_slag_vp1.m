% Calculate kinetic energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnOL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:47
% EndTime: 2020-04-12 19:40:48
% DurationCPUTime: 0.87s
% Computational Cost: add. (342->139), mult. (633->257), div. (0->0), fcn. (534->8), ass. (0->91)
t190 = sin(qJ(1));
t231 = t190 ^ 2;
t193 = cos(qJ(1));
t230 = t193 ^ 2;
t189 = sin(qJ(2));
t227 = Icges(3,4) * t189;
t192 = cos(qJ(2));
t226 = Icges(3,4) * t192;
t187 = qJ(2) + qJ(3);
t185 = sin(t187);
t225 = Icges(4,4) * t185;
t186 = cos(t187);
t224 = Icges(4,4) * t186;
t188 = sin(qJ(4));
t223 = Icges(5,4) * t188;
t191 = cos(qJ(4));
t222 = Icges(5,4) * t191;
t221 = qJD(1) * t192;
t220 = qJD(2) * t190;
t219 = qJD(2) * t193;
t218 = qJD(4) * t190;
t217 = qJD(4) * t193;
t216 = t189 * qJD(2);
t215 = qJD(2) + qJD(3);
t214 = rSges(3,1) * t192 - rSges(3,2) * t189;
t213 = -rSges(4,1) * t186 + rSges(4,2) * t185;
t212 = rSges(5,1) * t191 - rSges(5,2) * t188;
t211 = Icges(3,1) * t192 - t227;
t210 = -Icges(4,1) * t186 + t225;
t209 = Icges(5,1) * t191 - t223;
t208 = -Icges(3,2) * t189 + t226;
t207 = Icges(4,2) * t185 - t224;
t206 = -Icges(5,2) * t188 + t222;
t205 = Icges(3,5) * t192 - Icges(3,6) * t189;
t204 = -Icges(4,5) * t186 + Icges(4,6) * t185;
t203 = Icges(5,5) * t191 - Icges(5,6) * t188;
t157 = -Icges(5,6) * t193 + t190 * t206;
t161 = -Icges(5,5) * t193 + t190 * t209;
t202 = t157 * t188 - t161 * t191;
t158 = Icges(5,6) * t190 + t193 * t206;
t162 = Icges(5,5) * t190 + t193 * t209;
t201 = -t158 * t188 + t162 * t191;
t159 = -Icges(3,6) * t193 + t190 * t208;
t163 = -Icges(3,5) * t193 + t190 * t211;
t200 = t159 * t189 - t163 * t192;
t160 = Icges(3,6) * t190 + t193 * t208;
t164 = Icges(3,5) * t190 + t193 * t211;
t199 = -t160 * t189 + t164 * t192;
t177 = Icges(5,2) * t191 + t223;
t179 = Icges(5,1) * t188 + t222;
t198 = -t177 * t188 + t179 * t191;
t178 = Icges(3,2) * t192 + t227;
t180 = Icges(3,1) * t189 + t226;
t197 = -t178 * t189 + t180 * t192;
t173 = t215 * t190;
t174 = t215 * t193;
t196 = (-Icges(4,5) * t185 - Icges(4,6) * t186) * qJD(1) - (-Icges(4,3) * t193 + t190 * t204) * t174 + (Icges(4,3) * t190 + t193 * t204) * t173;
t147 = -Icges(4,6) * t193 + t190 * t207;
t148 = Icges(4,6) * t190 + t193 * t207;
t149 = -Icges(4,5) * t193 + t190 * t210;
t150 = Icges(4,5) * t190 + t193 * t210;
t170 = -Icges(4,2) * t186 - t225;
t171 = -Icges(4,1) * t185 - t224;
t195 = (t148 * t185 - t150 * t186) * t173 - (t147 * t185 - t149 * t186) * t174 + (t170 * t185 - t171 * t186) * qJD(1);
t184 = rSges(2,1) * t193 - rSges(2,2) * t190;
t183 = rSges(2,1) * t190 + rSges(2,2) * t193;
t182 = rSges(3,1) * t189 + rSges(3,2) * t192;
t181 = rSges(5,1) * t188 + rSges(5,2) * t191;
t176 = Icges(3,5) * t189 + Icges(3,6) * t192;
t175 = Icges(5,5) * t188 + Icges(5,6) * t191;
t172 = -rSges(4,1) * t185 - rSges(4,2) * t186;
t168 = rSges(3,3) * t190 + t193 * t214;
t167 = rSges(5,3) * t190 + t193 * t212;
t166 = -rSges(3,3) * t193 + t190 * t214;
t165 = -rSges(5,3) * t193 + t190 * t212;
t156 = Icges(3,3) * t190 + t193 * t205;
t155 = -Icges(3,3) * t193 + t190 * t205;
t154 = Icges(5,3) * t190 + t193 * t203;
t153 = -Icges(5,3) * t193 + t190 * t203;
t152 = rSges(4,3) * t190 + t193 * t213;
t151 = -rSges(4,3) * t193 + t190 * t213;
t144 = -qJD(1) * t166 - t182 * t219;
t143 = qJD(1) * t168 - t182 * t220;
t142 = -t181 * t218 + (pkin(1) * t193 + t167) * qJD(1);
t141 = -t181 * t217 + (-pkin(1) * t190 - t165) * qJD(1);
t140 = (t166 * t190 + t168 * t193) * qJD(2);
t139 = (t165 * t190 + t167 * t193) * qJD(4);
t138 = -qJD(1) * t151 - t172 * t174 + (-t190 * t221 - t193 * t216) * pkin(2);
t137 = qJD(1) * t152 - t172 * t173 + (-t190 * t216 + t193 * t221) * pkin(2);
t136 = t151 * t173 + t152 * t174 + (t230 + t231) * t192 * qJD(2) * pkin(2);
t1 = m(3) * (t140 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + ((t190 * t176 + t193 * t197) * qJD(1) + (t231 * t156 + (t200 * t193 + (-t155 + t199) * t190) * t193) * qJD(2)) * t220 / 0.2e1 - ((-t193 * t176 + t190 * t197) * qJD(1) + (t230 * t155 + (t199 * t190 + (-t156 + t200) * t193) * t190) * qJD(2)) * t219 / 0.2e1 + m(4) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + t173 * (t196 * t190 + t195 * t193) / 0.2e1 - t174 * (t195 * t190 - t196 * t193) / 0.2e1 + m(5) * (t139 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + ((t190 * t175 + t193 * t198) * qJD(1) + (t231 * t154 + (t202 * t193 + (-t153 + t201) * t190) * t193) * qJD(4)) * t218 / 0.2e1 - ((-t193 * t175 + t190 * t198) * qJD(1) + (t230 * t153 + (t201 * t190 + (-t154 + t202) * t193) * t190) * qJD(4)) * t217 / 0.2e1 + (m(2) * (t183 ^ 2 + t184 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t160 * t192 + t164 * t189) * t190 - (t159 * t192 + t163 * t189) * t193) * qJD(2) + (-t148 * t186 - t150 * t185) * t173 - (-t147 * t186 - t149 * t185) * t174 + ((t158 * t191 + t162 * t188) * t190 - (t157 * t191 + t161 * t188) * t193) * qJD(4) + (-t186 * t170 - t185 * t171 + t191 * t177 + t192 * t178 + t188 * t179 + t189 * t180) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
