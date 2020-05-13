% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% fourbar1IC
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
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1IC_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:50
% EndTime: 2020-04-24 20:15:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (52->34), mult. (74->48), div. (4->3), fcn. (41->8), ass. (0->18)
t21 = m(3) * pkin(2);
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t20 = t11 * mrSges(3,1) + t14 * mrSges(3,2);
t19 = t11 * mrSges(3,2);
t18 = t14 * mrSges(3,1);
t17 = -qJ(3) + qJ(1);
t16 = -t18 + t19;
t15 = cos(qJ(1));
t13 = cos(qJ(3));
t12 = sin(qJ(1));
t10 = sin(qJ(3));
t7 = sin(qJ(2) + t17);
t4 = -mrSges(2,2) + t20;
t3 = -qJD(1) ^ 2 * pkin(2) - t15 * g(1) - t12 * g(2);
t2 = t12 * g(1) - t15 * g(2) + qJDD(1) * pkin(2);
t1 = mrSges(2,1) + t16 + t21;
t5 = [-(-t12 * t1 + t4 * t15) * g(1) - (t1 * t15 + t12 * t4) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + (t16 * qJDD(2) + (-0.2e1 * t18 + 0.2e1 * t19 + t21) * qJDD(1)) * pkin(2) + ((-pkin(3) * t7 + pkin(2) * sin(t17)) / pkin(3) * (Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (t11 * t3 - t14 * t2) - mrSges(3,2) * (-t11 * t2 - t14 * t3)) + t11 * pkin(2) / pkin(4) * (-(-t10 * mrSges(4,1) - mrSges(4,2) * t13) * g(1) - (mrSges(4,1) * t13 - t10 * mrSges(4,2)) * g(2) + Ifges(4,3) * qJDD(3))) / t7 + (0.2e1 * qJD(1) + qJD(2)) * t20 * pkin(2) * qJD(2);];
tau = t5(:);
