% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% fourbarprisOL
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbarprisOL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:06
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.03s
% Computational Cost: add. (28->20), mult. (44->28), div. (0->0), fcn. (16->4), ass. (0->10)
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t69 = t66 * g(1) + t64 * g(2);
t68 = -t64 * g(1) + t66 * g(2);
t67 = qJD(1) ^ 2;
t65 = cos(qJ(3));
t63 = sin(qJ(3));
t60 = -t67 * qJ(2) + qJDD(2) + t69;
t59 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t68;
t1 = [mrSges(2,1) * t68 - mrSges(2,2) * t69 + mrSges(3,1) * t60 - mrSges(3,3) * t59 + qJ(2) * (-m(3) * t59 + t67 * mrSges(3,1)) + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1); m(3) * t60 + qJDD(1) * mrSges(3,1) - t67 * mrSges(3,3); Ifges(4,3) * qJDD(3) + mrSges(4,1) * (-t63 * g(1) + t65 * g(2)) - mrSges(4,2) * (t65 * g(1) + t63 * g(2)); 0;];
tauJ = t1;
