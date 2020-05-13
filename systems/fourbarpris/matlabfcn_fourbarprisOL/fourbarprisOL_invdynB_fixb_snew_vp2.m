% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = fourbarprisOL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_invdynB_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:06
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (179->73), mult. (279->85), div. (0->0), fcn. (88->4), ass. (0->28)
t100 = -m(2) - m(3);
t99 = (-mrSges(2,1) - mrSges(3,3));
t98 = -mrSges(2,2) + mrSges(3,1);
t97 = (Ifges(3,4) + Ifges(2,6));
t96 = Ifges(2,5) - Ifges(3,6);
t89 = sin(qJ(1));
t91 = cos(qJ(1));
t83 = t91 * g(1) + t89 * g(2);
t81 = -t89 * g(1) + t91 * g(2);
t78 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t81;
t93 = qJD(1) ^ 2;
t95 = -m(3) * t78 + (t93 * mrSges(3,1)) + qJDD(1) * mrSges(3,3);
t72 = m(2) * t81 + qJDD(1) * mrSges(2,1) - (t93 * mrSges(2,2)) + t95;
t79 = -t93 * qJ(2) + qJDD(2) + t83;
t73 = m(2) * t83 + m(3) * t79 + t98 * qJDD(1) + (t99 * t93);
t94 = -t91 * t72 - t89 * t73;
t92 = qJD(3) ^ 2;
t90 = cos(qJ(3));
t88 = sin(qJ(3));
t82 = t90 * g(1) + t88 * g(2);
t80 = -t88 * g(1) + t90 * g(2);
t77 = m(4) * t82 - t92 * mrSges(4,1) - qJDD(3) * mrSges(4,2);
t76 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t92 * mrSges(4,2);
t75 = -mrSges(4,2) * g(3) - mrSges(4,3) * t80 + Ifges(4,5) * qJDD(3) - t92 * Ifges(4,6);
t74 = mrSges(4,1) * g(3) + mrSges(4,3) * t82 + t92 * Ifges(4,5) + Ifges(4,6) * qJDD(3);
t71 = -mrSges(3,2) * t78 - mrSges(2,3) * t81 + t98 * g(3) + t96 * qJDD(1) - (t97 * t93);
t70 = -mrSges(3,2) * t79 + mrSges(2,3) * t83 + t96 * t93 + t97 * qJDD(1) + (m(3) * qJ(2) - t99) * g(3);
t1 = [-m(1) * g(1) + t89 * t72 - t91 * t73 + t88 * t76 - t90 * t77; -m(1) * g(2) - t90 * t76 - t88 * t77 + t94; (-m(1) - m(4) + t100) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t89 * t70 - t91 * t71 + t88 * t74 - t90 * t75; -mrSges(1,3) * g(1) - t91 * t70 - t89 * t71 - t90 * t74 - t88 * t75 + (-pkin(1) * t100 + mrSges(1,1)) * g(3); qJ(2) * t95 + pkin(1) * t94 + mrSges(3,1) * t79 - mrSges(3,3) * t78 + mrSges(2,1) * t81 - mrSges(2,2) * t83 + mrSges(4,1) * t80 - mrSges(4,2) * t82 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * qJDD(3) + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
