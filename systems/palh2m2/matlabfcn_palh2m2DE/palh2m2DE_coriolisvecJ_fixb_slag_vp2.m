% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m2DE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:30
% DurationCPUTime: 0.87s
% Computational Cost: add. (560->138), mult. (1454->247), div. (0->0), fcn. (701->6), ass. (0->84)
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t98 = pkin(4) * qJD(2) ^ 2;
t70 = t47 * t98;
t95 = pkin(5) * qJD(3) ^ 2;
t28 = -t46 * t95 - t70;
t99 = pkin(4) * t47;
t33 = -pkin(1) - t99;
t32 = -pkin(2) + t33;
t96 = pkin(5) * t46;
t25 = t32 - t96;
t21 = -pkin(3) + t25;
t10 = t21 * qJD(1);
t45 = cos(qJ(4));
t44 = sin(qJ(2));
t79 = pkin(4) * qJD(2);
t43 = sin(qJ(3));
t97 = pkin(5) * t43;
t26 = -qJD(3) * t97 - t44 * t79;
t42 = sin(qJ(4));
t84 = t26 * t42;
t4 = -t10 * t45 - t84;
t78 = qJD(1) * t26;
t2 = qJD(4) * t4 + t28 * t45 + t42 * t78;
t83 = t26 * t45;
t62 = t42 * t10 - t83;
t3 = t62 * qJD(4) - t28 * t42 + t45 * t78;
t65 = t4 * t45 - t42 * t62;
t110 = t65 * qJD(4) - t2 * t45 + t3 * t42;
t73 = t44 * qJD(1);
t108 = qJD(2) - qJD(3);
t107 = 2 * m(5);
t106 = -0.2e1 * pkin(1);
t105 = m(7) / 0.2e1;
t104 = -t43 / 0.2e1;
t103 = -t44 / 0.2e1;
t102 = t46 / 0.2e1;
t101 = t47 / 0.2e1;
t100 = pkin(4) * t44;
t94 = m(4) * t33;
t93 = m(6) * t25;
t88 = Ifges(3,4) * t44;
t87 = Ifges(3,4) * t47;
t86 = Ifges(5,4) * t43;
t85 = Ifges(3,2) * t47;
t82 = t43 * t46;
t63 = -mrSges(5,1) * t46 + mrSges(5,2) * t43;
t81 = t63 * t73;
t80 = qJD(1) / 0.2e1;
t77 = qJD(2) * t47;
t76 = qJD(3) * t46;
t75 = qJD(4) * t21;
t74 = t43 * qJD(1);
t72 = t46 * qJD(1);
t69 = Ifges(5,2) * t82;
t67 = Ifges(3,6) * t103;
t66 = Ifges(5,6) * t104;
t64 = t4 * t42 + t45 * t62;
t23 = t43 * t47 - t46 * t44;
t30 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t74;
t31 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t72;
t59 = t43 * t44 + t46 * t47;
t61 = -t23 * t30 + t31 * t59;
t60 = -t23 * t43 - t46 * t59;
t38 = qJD(1) * t87;
t57 = (Ifges(3,1) * t73 + Ifges(3,5) * qJD(2) + t38) * t101 + (Ifges(3,6) * qJD(2) + (t85 + t88) * qJD(1)) * t103;
t37 = Ifges(5,4) * t72;
t56 = (Ifges(5,1) * t74 + Ifges(5,5) * qJD(3) + t37) * t102 + (Ifges(5,6) * qJD(3) + (Ifges(5,2) * t46 + t86) * qJD(1)) * t104;
t55 = t32 * (mrSges(5,1) * t43 + mrSges(5,2) * t46);
t54 = qJD(4) * (mrSges(7,1) * t45 - mrSges(7,2) * t42);
t53 = t60 * qJD(3);
t41 = qJD(1) + qJD(4);
t52 = m(7) * t64 + (mrSges(7,1) * t42 + mrSges(7,2) * t45) * t41;
t50 = qJD(1) ^ 2;
t27 = -t43 * t95 - t44 * t98;
t16 = t27 * t99;
t14 = (-t42 * t73 - t45 * t77) * pkin(4);
t13 = (-t42 * t74 - t45 * t76) * pkin(5);
t12 = (t42 * t77 - t45 * t73) * pkin(4);
t11 = (t42 * t76 - t45 * t74) * pkin(5);
t7 = t108 * t59;
t6 = t108 * t23;
t1 = t3 * mrSges(7,1);
t5 = [m(7) * (t65 * t26 + (t64 * qJD(4) - t2 * t42 - t3 * t45) * t21) + t28 * mrSges(6,3) - t2 * mrSges(7,2) + t1 - mrSges(4,3) * t70 + ((Ifges(5,5) * t102 + t66) * qJD(3) + t56) * qJD(3) + ((t45 * t75 - t84) * mrSges(7,2) + (t42 * t75 + t83) * mrSges(7,1)) * t41 + ((Ifges(3,5) * t101 + t67) * qJD(2) + (t81 + (-t43 * t6 - t46 * t7 + t53) * mrSges(5,3)) * pkin(4) + t57) * qJD(2) + (0.2e1 * (mrSges(6,1) - t93) * t26 + ((mrSges(3,2) * t106 + 0.3e1 / 0.2e1 * t87) * t47 + (mrSges(3,1) * t106 - 0.3e1 / 0.2e1 * t88 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t47 + (t32 * t107 - (2 * mrSges(4,1)) + t63 + 0.2e1 * t94) * pkin(4)) * t44) * qJD(2) + (0.2e1 * t55 + 0.3e1 / 0.2e1 * Ifges(5,1) * t82 - 0.3e1 / 0.2e1 * t69 + (-0.3e1 / 0.2e1 * t43 ^ 2 + 0.3e1 / 0.2e1 * t46 ^ 2) * Ifges(5,4)) * qJD(3)) * qJD(1); (-mrSges(7,1) * t12 + mrSges(7,2) * t14) * t41 + (t6 * t30 - t7 * t31 + (t41 * t54 + (mrSges(4,1) + mrSges(6,1) - t94) * t50) * t44) * pkin(4) + m(6) * (-t28 * t100 + t16) + (-t47 * t38 / 0.2e1 + (t44 * t85 / 0.2e1 + pkin(1) * (mrSges(3,1) * t44 + mrSges(3,2) * t47) + (Ifges(3,1) * t47 - t88) * t103) * qJD(1) + (mrSges(5,3) * t53 - t81) * pkin(4) - t57) * qJD(1) + (Ifges(3,5) * t80 * t47 + qJD(1) * t67 + (t52 * t47 + t61 + (-t23 * t7 + t59 * t6) * t107 * pkin(4)) * pkin(4)) * qJD(2) + (t110 * t100 - t4 * t12 + t62 * t14 + t16) * m(7) + (-m(5) * t32 - t93) * t50 * t100; t50 * mrSges(6,1) * t97 + (-t11 * mrSges(7,1) + t13 * mrSges(7,2) + t54 * t97) * t41 + (mrSges(5,1) * t6 + mrSges(5,2) * t7 - t61) * t79 - m(7) * (t4 * t11 - t13 * t62) + 0.2e1 * (m(6) / 0.2e1 + t105) * t27 * t96 + 0.2e1 * (t110 * t105 + (-t28 / 0.2e1 - t25 * t50 / 0.2e1) * m(6)) * t97 + (-t46 * t37 / 0.2e1 + (-t55 + (Ifges(5,1) * t46 - t86) * t104 + t69 / 0.2e1) * qJD(1) - t60 * mrSges(5,3) * t79 - t56) * qJD(1) + (qJD(1) * t66 + (t52 * pkin(5) + Ifges(5,5) * t80) * t46) * qJD(3); -mrSges(7,1) * t41 * t62 + t1 + (t4 * t41 - t2) * mrSges(7,2);];
tauc = t5(:);
