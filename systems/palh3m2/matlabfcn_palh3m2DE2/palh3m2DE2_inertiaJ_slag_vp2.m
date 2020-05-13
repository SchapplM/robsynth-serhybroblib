% Calculate joint inertia matrix for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:13
% EndTime: 2020-05-07 04:22:14
% DurationCPUTime: 1.60s
% Computational Cost: add. (574->229), mult. (600->240), div. (0->0), fcn. (240->90), ass. (0->101)
t222 = pkin(17) + pkin(18);
t210 = pkin(15) + t222;
t200 = (pkin(16) + t210);
t276 = pkin(3) * m(9);
t239 = cos(qJ(2));
t277 = 0.2e1 * t239;
t243 = m(5) + m(6);
t274 = pkin(8) * mrSges(6,2);
t273 = (pkin(10) * mrSges(6,3));
t236 = sin(pkin(15));
t272 = pkin(8) * t236;
t240 = cos(pkin(15));
t271 = pkin(8) * t240;
t237 = cos(qJ(4));
t270 = mrSges(6,1) * t237;
t233 = sin(qJ(4));
t269 = t233 * mrSges(6,2);
t215 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t207 = m(8) + m(9) + m(4) + t243;
t268 = t207 * pkin(1) ^ 2;
t234 = sin(qJ(3));
t235 = sin(qJ(2));
t267 = t234 * t235;
t266 = t234 * t239;
t221 = t243 * pkin(4) ^ 2;
t226 = qJ(2) + qJ(4);
t265 = qJ(4) - qJ(2);
t238 = cos(qJ(3));
t204 = -pkin(4) * t238 + pkin(1);
t162 = pkin(4) * t267 + t204 * t239;
t223 = pkin(15) + pkin(18);
t198 = pkin(4) * t243 + mrSges(4,1) + mrSges(9,1);
t263 = t198 * pkin(1) * t238;
t261 = t221 + Ifges(4,3) + Ifges(9,3);
t260 = -pkin(10) * m(6) - mrSges(6,3);
t216 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t259 = 2 * t200;
t258 = -qJ(2) + t210;
t257 = qJ(2) + t210;
t191 = -qJ(2) + t200;
t190 = qJ(2) + t200;
t189 = -qJ(4) + t200;
t188 = qJ(4) + t200;
t256 = mrSges(6,1) * t272 + t216 * t240;
t185 = -qJ(2) + t189;
t184 = qJ(2) + t189;
t183 = -qJ(2) + t188;
t182 = qJ(2) + t188;
t227 = qJ(2) + qJ(3);
t211 = sin(t227);
t212 = cos(t227);
t217 = qJ(3) + t226;
t218 = qJ(3) - t265;
t255 = -t212 * (Ifges(4,6) + Ifges(9,6)) + (-Ifges(4,5) - Ifges(9,5)) * t211 + (mrSges(5,3) * t211 + (-cos(t217) / 0.2e1 + cos(t218) / 0.2e1) * mrSges(6,1) + (sin(t217) + sin(t218)) * mrSges(6,2) / 0.2e1) * pkin(4);
t251 = pkin(8) ^ 2;
t250 = pkin(10) ^ 2;
t248 = 0.2e1 * qJ(2);
t247 = 0.2e1 * qJ(4);
t242 = pkin(8) * mrSges(6,1);
t232 = mrSges(4,2) + mrSges(9,2);
t229 = cos(pkin(16));
t228 = sin(pkin(16));
t225 = qJ(3) + t248;
t220 = pkin(8) * m(6) + mrSges(5,1);
t219 = -qJ(2) - pkin(15) + pkin(14);
t214 = -qJ(2) + t223;
t213 = qJ(2) + t223;
t209 = cos(t222);
t208 = sin(t222);
t206 = mrSges(5,2) + t260;
t205 = 0.2e1 * t227;
t203 = cos(t219);
t202 = sin(t219);
t201 = 0.2e1 * t223;
t197 = 0.2e1 * t219;
t196 = pkin(1) * t232 * t234;
t195 = qJ(3) + t257;
t194 = -qJ(3) + t258;
t193 = -qJ(4) + t259;
t192 = qJ(4) + t259;
t187 = qJ(3) + t190;
t186 = -qJ(3) + t191;
t181 = 2 * t200;
t179 = mrSges(6,1) * t233 + mrSges(6,2) * t237;
t174 = -qJ(3) + t185;
t173 = -qJ(3) + t183;
t172 = qJ(3) + t184;
t171 = qJ(3) + t182;
t170 = 0.2e1 * t189;
t169 = 0.2e1 * t188;
t168 = t238 * t239 - t267;
t167 = t235 * t238 + t266;
t166 = mrSges(6,2) * t272 + t215 * t240;
t165 = mrSges(6,1) * t271 - t216 * t236;
t164 = mrSges(6,2) * t271 - t215 * t236;
t161 = -pkin(4) * t266 + t204 * t235;
t160 = -t167 * t236 + t168 * t240;
t159 = t167 * t240 + t168 * t236;
t158 = -t161 * t236 + t162 * t240;
t157 = t161 * t240 + t162 * t236;
t1 = [((cos(t186) + cos(t187)) * t220 + (sin(t186) + sin(t187)) * t206 + (-sin(t171) / 0.2e1 + sin(t172) / 0.2e1 - sin(t173) / 0.2e1 + sin(t174) / 0.2e1) * mrSges(6,2) + (cos(t171) / 0.2e1 + cos(t172) / 0.2e1 + cos(t173) / 0.2e1 + cos(t174) / 0.2e1) * mrSges(6,1)) * pkin(4) + (t260 * pkin(8) - Ifges(5,4)) * sin(t181) + (-t216 - t274) * sin(t192) / 0.2e1 + (-t216 + t274) * sin(t193) / 0.2e1 + (t207 * pkin(12) * t277 + t232 * sin(t225) + (-cos(t190) - cos(t191)) * t220 + (-sin(t190) - sin(t191)) * t206 + (-cos(t225) - t238) * t198 + (-sin(t213) - sin(t214)) * mrSges(8,2) + (-cos(t213) - cos(t214)) * mrSges(8,1) + (-cos(t258) - cos(t257)) * t276 + (sin(t182) / 0.2e1 + sin(t183) / 0.2e1 - sin(t184) / 0.2e1 - sin(t185) / 0.2e1) * mrSges(6,2) + (-cos(t182) / 0.2e1 - cos(t183) / 0.2e1 - cos(t184) / 0.2e1 - cos(t185) / 0.2e1) * mrSges(6,1)) * pkin(1) + (-0.2e1 * cos(t210) * t276 + mrSges(3,1) * t277 - 0.2e1 * mrSges(8,1) * cos(t223) - 0.2e1 * mrSges(3,2) * t235 - 0.2e1 * mrSges(8,2) * sin(t223) - 0.2e1 * t220 * cos(t200) - 0.2e1 * t206 * sin(t200) - 0.2e1 * t198 * t212 + 0.2e1 * t232 * t211 + (m(3) + t207) * pkin(12) + (sin(t188) - sin(t189)) * mrSges(6,2) + (-cos(t188) - cos(t189)) * mrSges(6,1)) * pkin(12) + t268 / 0.2e1 + (m(7) * pkin(6) - 0.2e1 * mrSges(7,1) * t203 - 0.2e1 * mrSges(7,2) * t202) * pkin(6) + (-Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2) + t221) * cos(t205) / 0.2e1 + (-Ifges(3,1) + Ifges(3,2) + t268) * cos(t248) / 0.2e1 + (t242 - t215) * cos(t192) / 0.2e1 + (t242 + t215) * cos(t193) / 0.2e1 - pkin(8) * t269 + t273 + (0.2e1 * (-t250 + t251) * m(6) - (4 * t273) - (2 * Ifges(5,1)) - Ifges(6,1) + (2 * Ifges(5,2)) - Ifges(6,2) + (2 * Ifges(6,3))) * cos(t181) / 0.4e1 + t196 + Ifges(3,4) * sin(t248) + ((sin(t194) - sin(t195)) * mrSges(9,2) + (cos(t194) + cos(t195)) * mrSges(9,1) + (0.1e1 / 0.2e1 + cos(0.2e1 * t210) / 0.2e1) * t276) * pkin(3) + (Ifges(4,4) + Ifges(9,4)) * sin(t205) + t221 / 0.2e1 + (Ifges(7,2) - Ifges(7,1)) * cos(t197) / 0.2e1 + (Ifges(8,2) - Ifges(8,1)) * cos(t201) / 0.2e1 + (t250 + t251) * m(6) / 0.2e1 + pkin(8) * t270 - Ifges(7,4) * sin(t197) - Ifges(8,4) * sin(t201) + (cos(t169) / 0.8e1 + cos(t170) / 0.8e1 - cos(t247) / 0.4e1) * (Ifges(6,2) - Ifges(6,1)) + (sin(t169) / 0.4e1 - sin(t247) / 0.2e1 - sin(t170) / 0.4e1) * Ifges(6,4) + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.2e1 + Ifges(8,1) / 0.2e1 + Ifges(9,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.2e1 + Ifges(8,2) / 0.2e1 + Ifges(9,2) / 0.2e1 + Ifges(2,3) + Ifges(6,3) / 0.2e1; Ifges(7,6) * t203 - Ifges(7,5) * t202 + Ifges(3,5) * t235 + Ifges(3,6) * t239 + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t235 + (-sin(-t265) / 0.2e1 - sin(t226) / 0.2e1) * mrSges(6,2) + (-cos(-t265) / 0.2e1 + cos(t226) / 0.2e1) * mrSges(6,1)) * pkin(1) + t255; Ifges(3,3) + Ifges(7,3) + 0.2e1 * t196 + t261 - 0.2e1 * t263 + t268; t255; t196 + t261 - t263; t261; ((t165 * t229 - t256 * t228) * t237 + (-t164 * t229 + t166 * t228) * t233 + Ifges(6,3) * (-t228 * t236 + t229 * t240)) * t209 + ((-t228 * t165 - t256 * t229) * t237 + (t164 * t228 + t166 * t229) * t233 - Ifges(6,3) * (t228 * t240 + t229 * t236)) * t208 - (pkin(12) + t162) * (-t269 + t270); -((t157 * t229 + t158 * t228) * t209 + t208 * (-t157 * t228 + t158 * t229)) * t179; ((t159 * t229 + t160 * t228) * t209 + t208 * (-t159 * t228 + t160 * t229)) * t179 * pkin(4); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
