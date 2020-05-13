% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fivebar1OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:12
% EndTime: 2020-04-27 06:13:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (76->50), mult. (124->70), div. (0->0), fcn. (64->8), ass. (0->27)
t209 = pkin(2) * m(3);
t208 = pkin(3) * m(5);
t192 = sin(qJ(4));
t207 = t192 * mrSges(5,2);
t194 = sin(qJ(2));
t206 = t194 * mrSges(3,2);
t196 = cos(qJ(4));
t205 = t196 * mrSges(5,1);
t198 = cos(qJ(2));
t204 = t198 * mrSges(3,1);
t203 = t194 * mrSges(3,1) + t198 * mrSges(3,2);
t202 = -t204 + t206;
t201 = t205 - t207;
t200 = t192 * mrSges(5,1) + t196 * mrSges(5,2);
t199 = cos(qJ(1));
t197 = cos(qJ(3));
t195 = sin(qJ(1));
t193 = sin(qJ(3));
t187 = -mrSges(2,2) + t203;
t186 = mrSges(4,2) + t200;
t185 = -qJD(1) ^ 2 * pkin(2) - t199 * g(1) - t195 * g(2);
t184 = -qJD(3) ^ 2 * pkin(3) - t197 * g(1) - t193 * g(2);
t183 = qJDD(1) * pkin(2) + t195 * g(1) - t199 * g(2);
t182 = qJDD(3) * pkin(3) + t193 * g(1) - t197 * g(2);
t181 = mrSges(2,1) + t202 + t209;
t180 = mrSges(4,1) + t201 + t208;
t1 = [-(-t195 * t181 + t187 * t199) * g(1) - (t181 * t199 + t195 * t187) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + (t202 * qJDD(2) + (-0.2e1 * t204 + 0.2e1 * t206 + t209) * qJDD(1) + (0.2e1 * qJD(1) + qJD(2)) * t203 * qJD(2)) * pkin(2); Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (-t198 * t183 + t194 * t185) - mrSges(3,2) * (-t194 * t183 - t198 * t185); -(-t193 * t180 - t186 * t197) * g(1) - (t180 * t197 - t186 * t193) * g(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + Ifges(5,3) * qJDD(4) + (t201 * qJDD(4) + (0.2e1 * t205 - 0.2e1 * t207 + t208) * qJDD(3) + (-0.2e1 * qJD(3) - qJD(4)) * t200 * qJD(4)) * pkin(3); Ifges(5,3) * (qJDD(3) + qJDD(4)) + mrSges(5,1) * (t196 * t182 - t192 * t184) - mrSges(5,2) * (t192 * t182 + t196 * t184); 0;];
tauJ = t1;
