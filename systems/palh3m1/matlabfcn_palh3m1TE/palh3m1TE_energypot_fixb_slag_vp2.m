% Calculate potential energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1TE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:18:55
% EndTime: 2020-04-17 15:19:05
% DurationCPUTime: 5.24s
% Computational Cost: add. (78745->148), mult. (118441->209), div. (5616->6), fcn. (74977->24), ass. (0->99)
t201 = cos(qJ(2));
t274 = pkin(1) * t201 + pkin(13);
t273 = pkin(3) * m(9);
t272 = -m(5) - m(6);
t271 = -m(4) - m(9) - m(8);
t196 = sin(qJ(2));
t198 = sin(pkin(16));
t203 = cos(pkin(16));
t181 = t196 * t198 - t201 * t203;
t259 = pkin(5) * t181;
t252 = (-0.2e1 * t259 + pkin(1)) * pkin(1);
t244 = pkin(5) ^ 2 + t252;
t251 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t169 = t244 + t251;
t177 = pkin(1) * t181 - pkin(5);
t264 = -pkin(6) + pkin(2);
t265 = -pkin(6) - pkin(2);
t168 = sqrt(-((pkin(5) - t264) * (pkin(5) + t264) + t252) * ((pkin(5) - t265) * (pkin(5) + t265) + t252));
t183 = t196 * t203 + t198 * t201;
t253 = t183 * t168;
t237 = -pkin(1) * t253 - t169 * t177;
t238 = pkin(1) * t169 * t183 - t168 * t177;
t248 = cos(pkin(15)) / 0.2e1;
t170 = 0.1e1 / t244;
t254 = t170 / pkin(6);
t261 = sin(pkin(15));
t162 = (t237 * t248 + t238 * t261 / 0.2e1) * t254;
t163 = (t238 * t248 - t237 * t261 / 0.2e1) * t254;
t270 = pkin(7) * m(7) - m(3) * pkin(13) - mrSges(3,1) * t201 - mrSges(7,1) * t162 + mrSges(3,2) * t196 + mrSges(7,2) * t163 - mrSges(2,1);
t269 = m(6) * pkin(11) - mrSges(5,2) + mrSges(6,3);
t194 = sin(qJ(4));
t199 = cos(qJ(4));
t268 = -pkin(9) * m(6) - mrSges(6,1) * t199 + mrSges(6,2) * t194 - mrSges(5,1);
t267 = t194 * mrSges(6,1) + t199 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t266 = 0.1e1 / pkin(2);
t263 = -pkin(8) - pkin(10);
t262 = -pkin(8) + pkin(10);
t258 = cos(pkin(19));
t257 = sin(pkin(19));
t188 = t196 * pkin(1) + pkin(12);
t256 = cos(pkin(17));
t255 = sin(pkin(17));
t197 = sin(qJ(1));
t184 = t274 * t197;
t202 = cos(qJ(1));
t185 = t274 * t202;
t250 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t249 = pkin(18) + pkin(19);
t243 = pkin(1) - t259;
t242 = cos(t249);
t241 = sin(t249);
t240 = t244 - t251;
t234 = t266 * (-pkin(5) * t253 + t240 * t243);
t232 = -t234 / 0.2e1;
t233 = t266 * (pkin(5) * t183 * t240 + t168 * t243) / 0.2e1;
t160 = (t232 * t258 + t233 * t257) * t170;
t231 = t234 / 0.2e1;
t161 = (t231 * t257 + t233 * t258) * t170;
t159 = t160 * t201 - t161 * t196;
t158 = t160 * t196 + t161 * t201;
t195 = sin(qJ(3));
t200 = cos(qJ(3));
t239 = t195 * t201 + t196 * t200;
t180 = t195 * t196 - t200 * t201;
t230 = t195 * t231 + t200 * t233;
t229 = t195 * t233 + t200 * t232;
t228 = t170 * (-t229 * t242 - t230 * t241);
t227 = t170 * (t229 * t241 - t230 * t242);
t226 = pkin(4) * t228;
t225 = pkin(3) - t226;
t224 = -pkin(3) * t228 + pkin(4);
t223 = -0.2e1 * pkin(3) * t226 + pkin(4) ^ 2;
t222 = pkin(3) ^ 2 + t223;
t221 = 0.1e1 / t222;
t220 = 0.1e1 / pkin(8) * t221;
t219 = 0.1e1 / pkin(10) * t221;
t218 = t222 + t250;
t217 = t222 - t250;
t216 = sqrt(-((pkin(3) - t262) * (pkin(3) + t262) + t223) * ((pkin(3) - t263) * (pkin(3) + t263) + t223));
t215 = t216 * t227;
t214 = (-pkin(4) * t215 + t217 * t225) * t220;
t213 = (-pkin(3) * t215 + t218 * t224) * t219;
t212 = -(pkin(4) * t217 * t227 + t216 * t225) * t220 / 0.2e1;
t211 = (pkin(3) * t218 * t227 + t216 * t224) * t219 / 0.2e1;
t193 = cos(pkin(18));
t192 = sin(pkin(18));
t176 = t180 * t202;
t175 = t239 * t202;
t174 = t180 * t197;
t173 = t239 * t197;
t157 = t158 * t202;
t156 = t159 * t202;
t155 = t158 * t197;
t154 = t159 * t197;
t149 = -t193 * t214 / 0.2e1 + t192 * t212;
t148 = t192 * t214 / 0.2e1 + t193 * t212;
t147 = t256 * t211 + t255 * t213 / 0.2e1;
t146 = -t256 * t213 / 0.2e1 + t255 * t211;
t1 = (-mrSges(3,1) * t196 - mrSges(3,2) * t201 + mrSges(4,1) * t239 - mrSges(4,2) * t180 - mrSges(8,1) * t158 - mrSges(8,2) * t159 - mrSges(2,3) - mrSges(1,3) - m(7) * pkin(14) - t163 * mrSges(7,1) - t162 * mrSges(7,2) - (t158 * t193 - t159 * t192) * t273 - (t148 * t159 + t149 * t158) * mrSges(9,1) - (-t148 * t158 + t149 * t159) * mrSges(9,2) + t272 * (-pkin(4) * t239 + t188) + t271 * t188 + t269 * (t146 * t180 + t147 * t239) + t268 * (-t146 * t239 + t147 * t180) + (-m(7) - m(3) - m(2)) * pkin(12)) * g(3) + (-mrSges(4,1) * t174 - mrSges(4,2) * t173 - mrSges(8,1) * t154 + mrSges(8,2) * t155 - (-t148 * t155 + t149 * t154) * mrSges(9,1) - (-t148 * t154 - t149 * t155) * mrSges(9,2) - (t154 * t193 + t155 * t192) * t273 - mrSges(1,2) + t272 * (t174 * pkin(4) + t184) + t271 * t184 + t269 * (t146 * t173 - t147 * t174) + t270 * t197 + t268 * (t146 * t174 + t147 * t173) + t267 * t202) * g(2) + (mrSges(8,2) * t157 - mrSges(4,1) * t176 - mrSges(4,2) * t175 - (-t148 * t157 + t149 * t156) * mrSges(9,1) - (-t148 * t156 - t149 * t157) * mrSges(9,2) - (t156 * t193 + t157 * t192) * t273 - mrSges(8,1) * t156 - mrSges(1,1) + t272 * (t176 * pkin(4) + t185) + t271 * t185 + t269 * (t146 * t175 - t147 * t176) + t270 * t202 + t268 * (t146 * t176 + t147 * t175) - t267 * t197) * g(1);
U = t1;
