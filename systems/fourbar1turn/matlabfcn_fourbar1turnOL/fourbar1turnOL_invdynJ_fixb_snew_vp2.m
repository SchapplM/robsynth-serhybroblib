% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% fourbar1turnOL
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
%   pkin=[l1,l2,l3,l4,l5]';
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
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar1turnOL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:52
% EndTime: 2020-04-12 19:40:52
% DurationCPUTime: 0.33s
% Computational Cost: add. (467->128), mult. (982->183), div. (0->0), fcn. (611->8), ass. (0->49)
t220 = qJD(1) * qJD(2);
t219 = qJD(1) * qJD(4);
t210 = sin(qJ(2));
t218 = t210 * t220;
t211 = sin(qJ(1));
t215 = cos(qJ(1));
t202 = t211 * g(1) - t215 * g(2);
t203 = -t215 * g(1) - t211 * g(2);
t214 = cos(qJ(2));
t186 = -t210 * g(3) + t214 * t203;
t185 = -t214 * g(3) - t210 * t203;
t216 = qJD(1) ^ 2;
t179 = (t210 * t214 * t216 + qJDD(2)) * pkin(2) + t185;
t180 = (-t214 ^ 2 * t216 - qJD(2) ^ 2) * pkin(2) + t186;
t209 = sin(qJ(3));
t213 = cos(qJ(3));
t167 = -t213 * t179 + t209 * t180;
t168 = -t209 * t179 - t213 * t180;
t194 = (-t214 * t209 - t210 * t213) * qJD(1);
t199 = t210 * qJDD(1) + t214 * t220;
t201 = t214 * qJDD(1) - t218;
t172 = -t194 * qJD(3) + t209 * t199 - t213 * t201;
t193 = (t210 * t209 - t214 * t213) * qJD(1);
t173 = t193 * qJD(3) - t213 * t199 - t209 * t201;
t207 = qJD(2) + qJD(3);
t175 = Ifges(4,4) * t194 + Ifges(4,2) * t193 + Ifges(4,6) * t207;
t176 = Ifges(4,1) * t194 + Ifges(4,4) * t193 + Ifges(4,5) * t207;
t206 = qJDD(2) + qJDD(3);
t217 = mrSges(4,1) * t167 - mrSges(4,2) * t168 + Ifges(4,5) * t173 + Ifges(4,6) * t172 + Ifges(4,3) * t206 + t194 * t175 - t193 * t176;
t212 = cos(qJ(4));
t208 = sin(qJ(4));
t200 = t212 * qJDD(1) - t208 * t219;
t198 = t208 * qJDD(1) + t212 * t219;
t197 = -t216 * pkin(1) + t203;
t196 = -qJDD(1) * pkin(1) - t202;
t192 = Ifges(3,5) * qJD(2) + (t210 * Ifges(3,1) + t214 * Ifges(3,4)) * qJD(1);
t191 = Ifges(5,5) * qJD(4) + (t208 * Ifges(5,1) + t212 * Ifges(5,4)) * qJD(1);
t190 = Ifges(3,6) * qJD(2) + (t210 * Ifges(3,4) + t214 * Ifges(3,2)) * qJD(1);
t189 = Ifges(5,6) * qJD(4) + (t208 * Ifges(5,4) + t212 * Ifges(5,2)) * qJD(1);
t184 = -t208 * g(3) + t212 * t197;
t183 = -t212 * g(3) - t208 * t197;
t182 = t207 * mrSges(4,1) - t194 * mrSges(4,3);
t181 = -t207 * mrSges(4,2) + t193 * mrSges(4,3);
t178 = (-t201 + t218) * pkin(2) - t202;
t177 = -t193 * mrSges(4,1) + t194 * mrSges(4,2);
t174 = Ifges(4,5) * t194 + Ifges(4,6) * t193 + Ifges(4,3) * t207;
t165 = mrSges(4,2) * t178 - mrSges(4,3) * t167 + Ifges(4,1) * t173 + Ifges(4,4) * t172 + Ifges(4,5) * t206 + t193 * t174 - t207 * t175;
t164 = -mrSges(4,1) * t178 + mrSges(4,3) * t168 + Ifges(4,4) * t173 + Ifges(4,2) * t172 + Ifges(4,6) * t206 - t194 * t174 + t207 * t176;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t202 - mrSges(2,2) * t203 + t210 * (-mrSges(3,2) * t202 - mrSges(3,3) * t185 + Ifges(3,1) * t199 + Ifges(3,4) * t201 + Ifges(3,5) * qJDD(2) - qJD(2) * t190 + t209 * t164 - t213 * t165) + t214 * (Ifges(3,4) * t199 + Ifges(3,2) * t201 + Ifges(3,6) * qJDD(2) + qJD(2) * t192 + mrSges(3,1) * t202 + mrSges(3,3) * t186 - t209 * t165 - t213 * t164 - pkin(2) * (m(4) * t178 - t172 * mrSges(4,1) + t173 * mrSges(4,2) - t193 * t181 + t194 * t182)) + t208 * (mrSges(5,2) * t196 - mrSges(5,3) * t183 + Ifges(5,1) * t198 + Ifges(5,4) * t200 + Ifges(5,5) * qJDD(4) - qJD(4) * t189) + t212 * (-mrSges(5,1) * t196 + mrSges(5,3) * t184 + Ifges(5,4) * t198 + Ifges(5,2) * t200 + Ifges(5,6) * qJDD(4) + qJD(4) * t191) + (-m(5) * t196 + t200 * mrSges(5,1) - t198 * mrSges(5,2) + ((-mrSges(5,1) * t208 - mrSges(5,2) * t212) * qJD(4) + (t208 ^ 2 + t212 ^ 2) * qJD(1) * mrSges(5,3)) * qJD(1)) * pkin(1); Ifges(3,5) * t199 + Ifges(3,6) * t201 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t185 - mrSges(3,2) * t186 + pkin(2) * (-t209 * (m(4) * t168 - t206 * mrSges(4,2) + t172 * mrSges(4,3) + t193 * t177 - t207 * t182) - t213 * (m(4) * t167 + t206 * mrSges(4,1) - t173 * mrSges(4,3) - t194 * t177 + t207 * t181)) + (t210 * t190 - t214 * t192) * qJD(1) + t217; t217; mrSges(5,1) * t183 - mrSges(5,2) * t184 + Ifges(5,5) * t198 + Ifges(5,6) * t200 + Ifges(5,3) * qJDD(4) + (t208 * t189 - t212 * t191) * qJD(1); 0;];
tauJ = t1;
