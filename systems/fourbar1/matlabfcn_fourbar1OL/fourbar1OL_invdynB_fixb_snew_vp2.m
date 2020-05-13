% Calculate vector of inverse dynamics base forces with Newton-Euler for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = fourbar1OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1OL_invdynB_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:27
% EndTime: 2020-04-24 20:10:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (194->71), mult. (286->86), div. (0->0), fcn. (184->6), ass. (0->41)
t154 = sin(qJ(1));
t153 = sin(qJ(2));
t156 = cos(qJ(2));
t143 = t153 * mrSges(3,1) + t156 * mrSges(3,2);
t179 = -mrSges(2,2) + t143;
t181 = t154 * t179;
t177 = (2 * qJD(1) + qJD(2)) * qJD(2);
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t140 = t152 * mrSges(4,1) + mrSges(4,2) * t155;
t146 = mrSges(4,1) * t155 - t152 * mrSges(4,2);
t159 = qJD(3) ^ 2;
t180 = t146 * qJDD(3) - t140 * t159;
t157 = cos(qJ(1));
t134 = t179 * t157;
t172 = t156 * mrSges(3,1);
t173 = t153 * mrSges(3,2);
t144 = t172 - t173;
t176 = m(3) * pkin(2);
t136 = -mrSges(2,1) - t176 + t144;
t129 = t154 * t136 + t134;
t133 = t136 * t157;
t178 = mrSges(1,1) + t146 - t133;
t142 = t156 * Ifges(3,5) - t153 * Ifges(3,6);
t145 = Ifges(3,5) * t153 + Ifges(3,6) * t156;
t135 = -mrSges(3,3) * pkin(2) + Ifges(2,5) - t142;
t138 = -Ifges(2,6) + t145;
t164 = t154 * t135 - t138 * t157;
t131 = -t142 * t157 + t154 * t145;
t163 = t154 * t142 + t145 * t157;
t132 = t143 * t157 + t154 * t144;
t130 = t154 * t143 - t144 * t157;
t162 = mrSges(1,2) + t140;
t161 = qJD(1) ^ 2;
t151 = m(2) + m(3) + m(4) + m(1);
t149 = mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3);
t141 = t155 * Ifges(4,5) - t152 * Ifges(4,6);
t139 = -t152 * Ifges(4,5) - t155 * Ifges(4,6);
t128 = -t133 + t181;
t127 = t135 * t157 + t154 * t138;
t1 = [-t151 * g(1) + t129 * qJDD(1) + t132 * qJDD(2) - t140 * qJDD(3) - t128 * t161 - t177 * t130 - t146 * t159; -t151 * g(2) + t128 * qJDD(1) + t130 * qJDD(2) + t129 * t161 + t177 * t132 + t180; -t151 * g(3); t149 * g(2) - (t162 - t129) * g(3) + t127 * qJDD(1) - t164 * t161 + t131 * qJDD(2) + t141 * qJDD(3) + t139 * t159 + t177 * t163; -(-m(4) * pkin(1) - t178 - t181) * g(3) - t149 * g(1) + t164 * qJDD(1) + t127 * t161 - t163 * qJDD(2) - t139 * qJDD(3) + t141 * t159 + t177 * t131; -t178 * g(2) - (t134 - t162) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + Ifges(4,3) * qJDD(3) + (-t136 * g(1) - g(2) * t179) * t154 + (-m(4) * g(2) + t180) * pkin(1) + (-t144 * qJDD(2) + t177 * t143 + (-0.2e1 * t172 + 0.2e1 * t173 + t176) * qJDD(1)) * pkin(2);];
tauB = t1;
