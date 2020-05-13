% Calculate matrix of centrifugal and coriolis load on the joints for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2DE2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:38
% EndTime: 2020-05-02 20:58:52
% DurationCPUTime: 4.00s
% Computational Cost: add. (13391->225), mult. (23667->335), div. (0->0), fcn. (32694->22), ass. (0->156)
t167 = sin(qJ(3));
t273 = pkin(1) * t167;
t213 = pkin(5) + t273;
t157 = pkin(22) + pkin(21);
t151 = sin(t157);
t152 = cos(t157);
t168 = sin(qJ(2));
t172 = cos(qJ(2));
t164 = cos(pkin(20));
t169 = sin(pkin(18));
t173 = cos(pkin(18));
t231 = sin(pkin(20));
t131 = t164 * t169 - t173 * t231;
t135 = t173 * t164 + t169 * t231;
t171 = cos(qJ(3));
t299 = t167 * t131 + t135 * t171;
t98 = t131 * t171 - t135 * t167;
t211 = -t168 * t299 + t98 * t172;
t319 = t168 * t98 + t299 * t172;
t306 = -t151 * t319 + t152 * t211;
t272 = pkin(1) * t171;
t56 = t151 * t211 + t319 * t152;
t55 = t56 * t272;
t35 = -t213 * t306 + t55;
t330 = t35 * t56;
t170 = cos(qJ(4));
t166 = sin(qJ(4));
t83 = -t131 * t152 + t135 * t151;
t243 = t166 * t83;
t84 = -t131 * t151 - t135 * t152;
t68 = mrSges(6,2) * t84 - mrSges(6,3) * t243;
t241 = t170 * t68;
t239 = t170 * t83;
t69 = -mrSges(6,1) * t84 - mrSges(6,3) * t239;
t244 = t166 * t69;
t268 = mrSges(5,3) * t84;
t181 = t241 / 0.2e1 - t244 / 0.2e1 + t268 / 0.2e1;
t331 = pkin(5) * t56;
t332 = t181 * t331;
t329 = t213 * t56;
t158 = t166 ^ 2;
t159 = t170 ^ 2;
t223 = t158 + t159;
t321 = t223 * t35;
t328 = t321 * t56;
t327 = (-t223 + 0.1e1) * t56 * pkin(5) ^ 2 * t306;
t192 = pkin(1) * t306;
t33 = t171 * t192 + t329;
t325 = t306 * t33;
t269 = mrSges(5,3) * t83;
t254 = mrSges(6,2) * t170;
t257 = mrSges(6,1) * t166;
t200 = t254 + t257;
t67 = t200 * t83;
t302 = -t67 - t269;
t323 = t302 * t306;
t136 = -t172 * t167 - t168 * t171;
t138 = -t167 * t168 + t171 * t172;
t160 = qJ(2) + qJ(3);
t154 = sin(t160);
t155 = cos(t160);
t165 = cos(pkin(19));
t232 = sin(pkin(19));
t132 = t165 * t171 - t167 * t232;
t225 = t132 * t172;
t130 = t167 * t165 + t171 * t232;
t226 = t130 * t168;
t197 = t225 - t226;
t96 = t130 * t172 + t132 * t168;
t320 = Ifges(11,5) * t155 + Ifges(4,5) * t138 + Ifges(9,5) * t197 - Ifges(11,6) * t154 + Ifges(4,6) * t136 - Ifges(9,6) * t96;
t318 = Ifges(4,4) * t136 - (Ifges(4,2) - Ifges(4,1)) * t138;
t316 = Ifges(11,4) * t154 - t155 * (Ifges(11,1) - Ifges(11,2));
t146 = t168 * pkin(1) - pkin(15);
t161 = sin(pkin(22));
t163 = cos(pkin(22));
t133 = t161 * t169 + t163 * t173;
t134 = t161 * t173 - t169 * t163;
t78 = (cos(pkin(21)) * t133 - t134 * sin(pkin(21))) * pkin(4) + t146;
t315 = -Ifges(11,4) * t155 ^ 2 - Ifges(4,4) * t138 ^ 2 - t78 * (mrSges(11,1) * t154 + mrSges(11,2) * t155) - t146 * (-mrSges(4,1) * t136 + mrSges(4,2) * t138);
t267 = mrSges(6,3) * t83;
t240 = t170 * t69;
t245 = t166 * t68;
t298 = -t240 / 0.2e1 - t245 / 0.2e1;
t308 = t298 - t267 * t223 / 0.2e1;
t270 = pkin(2) * t132;
t271 = pkin(2) * t130;
t179 = (mrSges(11,2) + mrSges(4,2)) * t273 - mrSges(10,1) * t270 - mrSges(10,2) * t271 - (mrSges(11,1) + mrSges(4,1)) * t272;
t305 = pkin(1) * t172;
t297 = t223 * t306;
t85 = -pkin(2) * t197 - pkin(15);
t293 = pkin(2) * (m(10) * t85 + mrSges(10,1) * t168 + mrSges(10,2) * t172) - pkin(15) * mrSges(9,1);
t115 = -pkin(5) * t138 + t146;
t66 = -pkin(9) * t84 - pkin(11) * t83 + t115;
t291 = m(6) * t223 * t66 + m(5) * t115 - mrSges(5,1) * t84 + mrSges(5,2) * t83 + t240 + t245;
t290 = 4 * qJD(2);
t289 = 4 * qJD(3);
t288 = m(5) / 0.2e1;
t286 = m(6) / 0.2e1;
t285 = m(6) / 0.4e1;
t24 = -t167 * t192 - t35 + t55;
t263 = t24 * t33;
t284 = m(5) * (t331 * t35 + t263);
t283 = m(6) * (t321 * t331 + t263);
t280 = pkin(5) * t306;
t278 = m(6) * t327;
t274 = sin(pkin(17));
t262 = t35 * t33;
t261 = Ifges(10,4) + Ifges(3,4);
t258 = pkin(5) * qJD(3);
t256 = mrSges(6,1) * t170;
t255 = mrSges(6,2) * t166;
t174 = cos(pkin(17));
t137 = t169 * t174 - t173 * t274;
t139 = t274 * t169 + t174 * t173;
t108 = t137 * t168 + t139 * t172;
t123 = -pkin(5) * t136 + t305;
t196 = t137 * t172 - t139 * t168;
t204 = -pkin(15) * mrSges(9,2) + Ifges(9,1) * t96;
t5 = t204 * t197 - (Ifges(9,2) * t197 - t293) * t96 + (pkin(15) * mrSges(3,2) - t85 * mrSges(10,2) + t261 * t168) * t168 + t291 * t123 + (t197 ^ 2 - t96 ^ 2) * Ifges(9,4) + (-pkin(15) * mrSges(3,1) + t85 * mrSges(10,1) - t261 * t172 + (-Ifges(10,1) - Ifges(3,1) + Ifges(3,2) + Ifges(10,2)) * t168 + (m(11) * t78 - mrSges(11,1) * t155 - mrSges(4,1) * t138 + mrSges(8,1) * t133 + mrSges(8,2) * t134 + (m(8) + m(4)) * t146) * pkin(1)) * t172 + (mrSges(11,2) * t305 - t316) * t154 + (-mrSges(4,2) * t305 - t318) * t136 + (pkin(14) * mrSges(7,2) + Ifges(7,4) * t196) * t196 + (pkin(14) * mrSges(7,1) + (Ifges(7,1) - Ifges(7,2)) * t196 - t108 * Ifges(7,4)) * t108 - t315;
t236 = t5 * qJD(1);
t180 = (t159 / 0.2e1 + t158 / 0.2e1) * t267 - t298;
t193 = t256 / 0.2e1 - t255 / 0.2e1;
t184 = t193 * t123;
t201 = -t255 + t256;
t190 = t201 * t83;
t187 = -t190 / 0.2e1;
t6 = t180 * t35 + t187 * t33 + t184;
t235 = t6 * qJD(1);
t34 = t272 * t306 + t329;
t229 = qJD(2) * t34;
t12 = -(-Ifges(9,4) * t96 + t293) * t96 + t316 * t154 + (-Ifges(9,4) * t197 + Ifges(9,2) * t96 - t204) * t197 + (t291 * pkin(5) + t318) * t136 + t315;
t228 = t12 * qJD(1);
t185 = t193 * t136;
t186 = t190 / 0.2e1;
t13 = (-t180 * t306 - t186 * t56 - t185) * pkin(5);
t227 = t13 * qJD(1);
t198 = t241 - t244;
t15 = t198 * t66 + ((-Ifges(6,4) * t239 + Ifges(6,6) * t84) * t170 + (Ifges(6,4) * t243 + Ifges(6,5) * t84 + (-Ifges(6,1) + Ifges(6,2)) * t239) * t166) * t83;
t224 = t15 * qJD(1);
t222 = mrSges(10,3) * t271;
t221 = mrSges(10,3) * t270;
t209 = t168 * t222;
t207 = t67 / 0.2e1 + t269 / 0.2e1;
t203 = t24 * t56 + t325;
t202 = t325 - t330;
t178 = t207 * t24 + t332;
t3 = (-t181 * t56 - t207 * t306) * pkin(5) + t178;
t199 = t3 * qJD(1);
t189 = qJD(4) * t201;
t188 = t198 + t268;
t14 = (-t187 * t56 - t185) * pkin(5) - t308 * t280;
t10 = (-t257 / 0.2e1 - t254 / 0.2e1 - t200 / 0.2e1) * t331;
t7 = t186 * t33 + t308 * t35 + t184;
t2 = t209 / 0.2e1 + (t226 / 0.2e1 - t225) * mrSges(10,3) * pkin(2) + t178 - pkin(5) * t323 / 0.2e1 + t332 + t320;
t1 = -t179 + t283 + t284;
t4 = [qJD(2) * t5 - qJD(3) * t12 + qJD(4) * t15, t2 * qJD(3) + t7 * qJD(4) + t188 * t229 + t236 + (Ifges(7,5) * t196 - Ifges(7,6) * t108 + t302 * t35 + (-Ifges(3,6) - Ifges(10,6) - t221) * t172 + (-Ifges(3,5) - Ifges(10,5) + t222) * t168 + ((t133 ^ 2 + t134 ^ 2) * t168 * mrSges(8,3) + (-t136 * t171 - t138 * t167) * mrSges(4,3) + (t154 * t171 - t155 * t167) * mrSges(11,3)) * pkin(1) + t320) * qJD(2), -t228 + t2 * qJD(2) + (-t172 * t221 + t209 + t320) * qJD(3) + t14 * qJD(4) + (t188 * t56 - t323) * t258, t224 + t7 * qJD(2) + t14 * qJD(3) + ((-mrSges(6,2) * t66 - Ifges(6,6) * t83) * t170 + (-mrSges(6,1) * t66 - Ifges(6,5) * t83) * t166) * qJD(4); qJD(3) * t3 - qJD(4) * t6 - t236, t1 * qJD(3) + ((t321 * t34 - t262) * t285 + (t34 * t35 - t262) * m(5) / 0.4e1) * t290, t1 * qJD(2) - t179 * qJD(3) - t278 * t289 / 0.4e1 + 0.2e1 * ((-t297 * t331 + t203 + t328) * t286 + (-t306 * t331 + t203 + t330) * t288) * t258 + t199, -t189 * t35 - t235; -qJD(2) * t3 - qJD(4) * t13 + t228, t278 * qJD(3) + (-t283 / 0.4e1 - t284 / 0.4e1) * t290 - t199 + (t179 + 0.2e1 * ((-t297 * t34 + t202 + t328) * t286 + (-t306 * t34 + t202 + t330) * t288) * pkin(5)) * qJD(2), t327 * t285 * t289 + qJD(2) * t278, t189 * t280 - t227; qJD(2) * t6 + qJD(3) * t13 - t224, t10 * qJD(3) - t200 * t229 + t235, -t200 * t56 * t258 + t10 * qJD(2) + t227, 0;];
Cq = t4;
