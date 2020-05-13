% Calculate vector of inverse dynamics joint torques for
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnOL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:52
% DurationCPUTime: 1.13s
% Computational Cost: add. (466->172), mult. (1156->267), div. (0->0), fcn. (728->10), ass. (0->79)
t48 = qJ(2) + qJ(3);
t44 = sin(t48);
t45 = cos(t48);
t105 = mrSges(4,1) * t44 + mrSges(4,2) * t45;
t70 = -mrSges(4,1) * t45 + t44 * mrSges(4,2);
t52 = sin(qJ(1));
t56 = cos(qJ(1));
t104 = g(1) * t56 + g(2) * t52;
t103 = m(5) * pkin(1);
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t54 = cos(qJ(3));
t55 = cos(qJ(2));
t65 = t50 * t55 + t51 * t54;
t25 = t65 * qJD(1);
t101 = -t25 / 0.2e1;
t49 = sin(qJ(4));
t100 = -t49 / 0.2e1;
t99 = -t51 / 0.2e1;
t53 = cos(qJ(4));
t98 = t53 / 0.2e1;
t97 = t55 / 0.2e1;
t96 = m(4) * pkin(2) ^ 2;
t73 = qJD(1) * qJD(2);
t33 = qJDD(1) * t55 - t51 * t73;
t93 = t33 * pkin(2);
t64 = t50 * t51 - t54 * t55;
t24 = t64 * qJD(1);
t89 = mrSges(4,3) * t24;
t88 = mrSges(4,3) * t25;
t87 = Ifges(3,4) * t51;
t86 = Ifges(4,4) * t25;
t85 = Ifges(5,4) * t49;
t84 = Ifges(5,4) * t53;
t83 = Ifges(5,5) * t53;
t82 = Ifges(5,2) * t53;
t81 = Ifges(5,6) * t49;
t80 = t105 * t52;
t79 = t105 * t56;
t78 = pkin(2) * qJD(2);
t77 = Ifges(3,5) * qJD(2);
t76 = Ifges(3,6) * qJD(2);
t75 = qJD(1) * t51;
t74 = qJD(1) * t55;
t72 = qJD(1) * qJD(4);
t71 = qJD(2) * qJD(3);
t69 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t68 = g(1) * t52 - g(2) * t56;
t67 = mrSges(3,1) * t51 + mrSges(3,2) * t55;
t66 = mrSges(5,1) * t49 + mrSges(5,2) * t53;
t41 = qJD(1) * t84;
t63 = (Ifges(5,6) * qJD(4) + (t82 + t85) * qJD(1)) * t100 + (qJD(1) * t49 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t41) * t98;
t62 = t65 * qJD(3);
t61 = t64 * qJD(3);
t60 = mrSges(2,1) + t103 + t70;
t47 = qJD(2) + qJD(3);
t11 = Ifges(4,2) * t24 + Ifges(4,6) * t47 - t86;
t19 = Ifges(4,4) * t24;
t12 = -Ifges(4,1) * t25 + Ifges(4,5) * t47 + t19;
t27 = (-qJDD(2) * t54 + t50 * t71) * pkin(2);
t28 = (-qJDD(2) * t50 - t54 * t71) * pkin(2);
t46 = qJDD(2) + qJDD(3);
t34 = qJDD(1) * t51 + t55 * t73;
t7 = qJD(1) * t61 - t33 * t50 - t34 * t54;
t8 = qJD(1) * t62 - t33 * t54 + t34 * t50;
t59 = -t28 * mrSges(4,2) + t11 * t101 + pkin(2) * (-mrSges(4,1) * t25 + mrSges(4,2) * t24) * t74 + t50 * t78 * t88 + t27 * mrSges(4,1) + t25 * (Ifges(4,1) * t24 + t86) / 0.2e1 + Ifges(4,3) * t46 + Ifges(4,6) * t8 + Ifges(4,5) * t7 - t47 * (Ifges(4,5) * t24 + Ifges(4,6) * t25) / 0.2e1 - (Ifges(4,2) * t25 + t12 + t19) * t24 / 0.2e1;
t58 = pkin(1) * t66 + (Ifges(5,1) * t53 - t85) * t100;
t42 = Ifges(3,4) * t74;
t35 = -mrSges(5,1) * t53 + mrSges(5,2) * t49;
t32 = qJDD(1) * t49 + t53 * t72;
t31 = qJDD(1) * t53 - t49 * t72;
t23 = Ifges(3,1) * t75 + t42 + t77;
t21 = t76 + (Ifges(3,2) * t55 + t87) * qJD(1);
t17 = mrSges(4,1) * t47 + t88;
t16 = -mrSges(4,2) * t47 + t89;
t15 = t65 * qJD(2) + t62;
t14 = t64 * qJD(2) + t61;
t13 = -mrSges(4,1) * t24 - mrSges(4,2) * t25;
t1 = [t47 * (Ifges(4,5) * t14 + Ifges(4,6) * t15) / 0.2e1 + t24 * (Ifges(4,4) * t14 + Ifges(4,2) * t15) / 0.2e1 + (Ifges(4,1) * t14 + Ifges(4,4) * t15) * t101 + t14 * t12 / 0.2e1 + t15 * t11 / 0.2e1 + Ifges(2,3) * qJDD(1) + (t31 * mrSges(5,1) - t32 * mrSges(5,2) + (-t35 + t103) * qJDD(1)) * pkin(1) + (t69 * t52 - t60 * t56) * g(2) + (t60 * t52 + t69 * t56) * g(1) + (t68 * mrSges(5,1) + Ifges(5,4) * t32 + Ifges(5,2) * t31 + Ifges(5,6) * qJDD(4)) * t53 + (-t68 * mrSges(5,2) + Ifges(5,1) * t32 + Ifges(5,4) * t31 + Ifges(5,5) * qJDD(4)) * t49 - (-mrSges(4,2) * t93 - t27 * mrSges(4,3) + Ifges(4,1) * t7 + Ifges(4,4) * t8 + Ifges(4,5) * t46) * t65 + (mrSges(4,1) * t93 + t28 * mrSges(4,3) + Ifges(4,4) * t7 + Ifges(4,2) * t8 + Ifges(4,6) * t46) * t64 + ((t83 / 0.2e1 - t81 / 0.2e1) * qJD(4) + ((-Ifges(5,2) * t49 + t84) * t98 - t58) * qJD(1) + t63) * qJD(4) + (-t68 * mrSges(3,2) + Ifges(3,1) * t34 + Ifges(3,4) * t33 + Ifges(3,5) * qJDD(2)) * t51 + (Ifges(3,4) * t34 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + t96) * t33 + t68 * mrSges(3,1) + (mrSges(4,1) * t8 - mrSges(4,2) * t7 - qJD(1) * (-mrSges(4,1) * t15 + mrSges(4,2) * t14) + t68 * m(4)) * pkin(2)) * t55 + (t21 * t99 + t23 * t97 + (Ifges(3,5) * t97 + Ifges(3,6) * t99) * qJD(2) + (t51 * (Ifges(3,1) * t55 - t87) / 0.2e1 + (Ifges(3,4) * t55 - Ifges(3,2) * t51) * t97 - t55 * t51 * t96) * qJD(1) + (t51 * t13 + (t14 * t54 - t15 * t50) * mrSges(4,3)) * pkin(2)) * qJD(2); ((t46 * mrSges(4,2) - t8 * mrSges(4,3) + qJD(3) * t17) * t50 + (-t46 * mrSges(4,1) - qJD(3) * t16 + (-qJD(2) * t24 + t7) * mrSges(4,3)) * t54 + (-g(3) * t55 + t104 * t51 - t27 * t54 - t28 * t50) * m(4)) * pkin(2) + ((-t77 / 0.2e1 - t23 / 0.2e1 - t42 / 0.2e1) * t55 + (t76 / 0.2e1 + t21 / 0.2e1 - pkin(2) * t13 + Ifges(3,4) * t75 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1 + t96) * t74) * t51) * qJD(1) + (t67 * t56 - t79) * g(1) + (t67 * t52 - t80) * g(2) + (-mrSges(3,1) * t55 + mrSges(3,2) * t51 - t70) * g(3) + t59 + Ifges(3,6) * t33 + Ifges(3,5) * t34 + Ifges(3,3) * qJDD(2); t59 - g(2) * t80 - g(1) * t79 - g(3) * t70 + (-t17 * t50 + (t16 - t89) * t54) * t78; Ifges(5,5) * t32 + Ifges(5,6) * t31 + Ifges(5,3) * qJDD(4) + g(3) * t35 + (-t53 * t41 / 0.2e1 - qJD(4) * (-t81 + t83) / 0.2e1 + (t49 * t82 / 0.2e1 + t58) * qJD(1) - t63) * qJD(1) + t104 * t66; 0;];
tau = t1;
