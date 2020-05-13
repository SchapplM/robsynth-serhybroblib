% Calculate vector of inverse dynamics joint torques with ic for
% fourbar1turnIC
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnIC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:30
% EndTime: 2020-05-07 11:33:32
% DurationCPUTime: 1.58s
% Computational Cost: add. (476->175), mult. (1161->276), div. (4->3), fcn. (733->12), ass. (0->84)
t53 = sin(qJ(2));
t118 = -t53 / 0.2e1;
t50 = qJ(2) + qJ(3);
t46 = sin(t50);
t47 = cos(t50);
t117 = mrSges(4,1) * t46 + mrSges(4,2) * t47;
t54 = sin(qJ(1));
t58 = cos(qJ(1));
t116 = g(1) * t58 + g(2) * t54;
t106 = m(4) * pkin(2) ^ 2;
t57 = cos(qJ(2));
t98 = Ifges(3,4) * t53;
t115 = (Ifges(3,1) * t57 - t98) * t118 + t57 * t53 * t106;
t74 = -t47 * mrSges(4,1) + t46 * mrSges(4,2);
t113 = m(4) * pkin(2);
t112 = m(5) * pkin(1);
t52 = sin(qJ(3));
t56 = cos(qJ(3));
t68 = t52 * t57 + t53 * t56;
t25 = t68 * qJD(1);
t110 = -t25 / 0.2e1;
t51 = sin(qJ(4));
t109 = -t51 / 0.2e1;
t55 = cos(qJ(4));
t108 = t55 / 0.2e1;
t107 = t57 / 0.2e1;
t80 = qJD(1) * qJD(2);
t33 = qJDD(1) * t57 - t53 * t80;
t103 = pkin(2) * t33;
t67 = t52 * t53 - t56 * t57;
t24 = t67 * qJD(1);
t99 = mrSges(4,3) * t24;
t97 = Ifges(5,4) * t51;
t96 = Ifges(5,4) * t55;
t95 = Ifges(3,5) * t57;
t94 = Ifges(5,5) * t55;
t93 = Ifges(5,2) * t55;
t92 = Ifges(3,6) * t53;
t91 = Ifges(5,6) * t51;
t49 = qJD(2) + qJD(3);
t90 = (mrSges(4,1) * t49 + mrSges(4,3) * t25) * t52;
t89 = t25 * Ifges(4,4);
t87 = t117 * t54;
t86 = t117 * t58;
t85 = pkin(2) * qJD(2);
t83 = t53 * qJD(1);
t82 = t57 * qJD(1);
t81 = -qJ(4) + qJ(2);
t79 = qJD(1) * qJD(4);
t78 = qJD(2) * qJD(3);
t77 = t53 * t113;
t76 = mrSges(4,3) * t85;
t72 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t71 = g(1) * t54 - g(2) * t58;
t69 = mrSges(5,1) * t51 + mrSges(5,2) * t55;
t43 = qJD(1) * t96;
t66 = (t51 * qJD(1) * Ifges(5,1) + Ifges(5,5) * qJD(4) + t43) * t108 + (Ifges(5,6) * qJD(4) + (t93 + t97) * qJD(1)) * t109;
t64 = t68 * qJD(3);
t63 = t67 * qJD(3);
t62 = mrSges(2,1) + t112 + t74;
t11 = t24 * Ifges(4,2) + t49 * Ifges(4,6) - t89;
t19 = Ifges(4,4) * t24;
t12 = -t25 * Ifges(4,1) + t49 * Ifges(4,5) + t19;
t27 = (-qJDD(2) * t56 + t52 * t78) * pkin(2);
t28 = (-qJDD(2) * t52 - t56 * t78) * pkin(2);
t48 = qJDD(2) + qJDD(3);
t34 = qJDD(1) * t53 + t57 * t80;
t7 = qJD(1) * t63 - t33 * t52 - t34 * t56;
t8 = qJD(1) * t64 - t33 * t56 + t34 * t52;
t61 = -t28 * mrSges(4,2) + t11 * t110 + pkin(2) * (-mrSges(4,1) * t25 + mrSges(4,2) * t24) * t82 + t52 * t25 * t76 + t27 * mrSges(4,1) + t25 * (Ifges(4,1) * t24 + t89) / 0.2e1 + Ifges(4,3) * t48 + Ifges(4,6) * t8 + Ifges(4,5) * t7 - t49 * (Ifges(4,5) * t24 + Ifges(4,6) * t25) / 0.2e1 - (Ifges(4,2) * t25 + t12 + t19) * t24 / 0.2e1;
t60 = pkin(1) * t69 + (Ifges(5,1) * t55 - t97) * t109;
t44 = Ifges(3,4) * t82;
t42 = sin(qJ(3) + t81);
t40 = 0.1e1 / t42;
t35 = -mrSges(5,1) * t55 + mrSges(5,2) * t51;
t32 = qJDD(1) * t51 + t55 * t79;
t31 = qJDD(1) * t55 - t51 * t79;
t23 = Ifges(3,1) * t83 + Ifges(3,5) * qJD(2) + t44;
t21 = Ifges(3,6) * qJD(2) + (t57 * Ifges(3,2) + t98) * qJD(1);
t16 = -mrSges(4,2) * t49 + t99;
t15 = qJD(2) * t68 + t64;
t14 = qJD(2) * t67 + t63;
t13 = -mrSges(4,1) * t24 - mrSges(4,2) * t25;
t1 = [t49 * (Ifges(4,5) * t14 + Ifges(4,6) * t15) / 0.2e1 + (Ifges(4,1) * t14 + Ifges(4,4) * t15) * t110 + t24 * (Ifges(4,4) * t14 + Ifges(4,2) * t15) / 0.2e1 + t14 * t12 / 0.2e1 + t15 * t11 / 0.2e1 + Ifges(2,3) * qJDD(1) + (mrSges(5,1) * t31 - mrSges(5,2) * t32 + (-t35 + t112) * qJDD(1)) * pkin(1) + (t54 * t72 - t58 * t62) * g(2) + (t54 * t62 + t58 * t72) * g(1) + (mrSges(5,1) * t71 + Ifges(5,4) * t32 + Ifges(5,2) * t31 + Ifges(5,6) * qJDD(4)) * t55 + (-mrSges(5,2) * t71 + Ifges(5,1) * t32 + Ifges(5,4) * t31 + Ifges(5,5) * qJDD(4)) * t51 - (-mrSges(4,2) * t103 - mrSges(4,3) * t27 + Ifges(4,1) * t7 + Ifges(4,4) * t8 + Ifges(4,5) * t48) * t68 + (mrSges(4,1) * t103 + mrSges(4,3) * t28 + Ifges(4,4) * t7 + Ifges(4,2) * t8 + Ifges(4,6) * t48) * t67 + ((t94 / 0.2e1 - t91 / 0.2e1) * qJD(4) + ((-Ifges(5,2) * t51 + t96) * t108 - t60) * qJD(1) + t66) * qJD(4) + (-mrSges(3,2) * t71 + Ifges(3,1) * t34 + Ifges(3,4) * t33 + Ifges(3,5) * qJDD(2)) * t53 + (Ifges(3,4) * t34 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + t106) * t33 + t71 * mrSges(3,1) + (mrSges(4,1) * t8 - mrSges(4,2) * t7 - qJD(1) * (-mrSges(4,1) * t15 + mrSges(4,2) * t14) + t71 * m(4)) * pkin(2)) * t57 + (t23 * t107 + t21 * t118 + (t95 / 0.2e1 - t92 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t57 - Ifges(3,2) * t53) * t107 - t115) * qJD(1) + (t53 * t13 + (t14 * t56 - t15 * t52) * mrSges(4,3)) * pkin(2)) * qJD(2); (-pkin(3) * t42 + pkin(2) * sin(t81)) / pkin(3) * t40 * (t61 + (-t90 + (t16 - t99) * t56) * t85 - g(3) * t74 - g(2) * t87 - g(1) * t86) + t61 + (-t27 * t56 - t28 * t52) * t113 - g(3) * (t57 * t113 + t74) - (-t92 + t95) * t80 / 0.2e1 - g(1) * (-t58 * t77 + t86) - g(2) * (-t54 * t77 + t87) - g(3) * (t57 * mrSges(3,1) - t53 * mrSges(3,2)) - t56 * pkin(2) * (mrSges(4,1) * t48 - mrSges(4,3) * t7) + Ifges(3,6) * t33 + Ifges(3,5) * t34 - t56 * t24 * t76 + Ifges(3,3) * qJDD(2) + (-t16 * t56 + t90) * pkin(2) * qJD(3) + (t21 / 0.2e1 - pkin(2) * t13) * t83 - (-Ifges(3,2) * t83 + t23 + t44) * t82 / 0.2e1 + t116 * (mrSges(3,1) * t53 + mrSges(3,2) * t57) + t115 * qJD(1) ^ 2 + (0.1e1 / pkin(4) * t40 * (Ifges(5,5) * t32 + Ifges(5,6) * t31 + Ifges(5,3) * qJDD(4) + g(3) * t35 + (-t55 * t43 / 0.2e1 - qJD(4) * (-t91 + t94) / 0.2e1 + (t51 * t93 / 0.2e1 + t60) * qJD(1) - t66) * qJD(1) + t116 * t69) + mrSges(4,2) * t48 - mrSges(4,3) * t8) * pkin(2) * t52;];
tau = t1(:);
