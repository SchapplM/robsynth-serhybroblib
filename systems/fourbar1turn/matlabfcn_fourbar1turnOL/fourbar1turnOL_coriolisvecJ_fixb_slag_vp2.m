% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnOL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:48
% EndTime: 2020-04-12 19:40:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (307->97), mult. (925->177), div. (0->0), fcn. (567->6), ass. (0->50)
t32 = qJD(2) + qJD(3);
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t37 = cos(qJ(3));
t38 = cos(qJ(2));
t26 = t34 * t35 - t37 * t38;
t24 = t26 * qJD(1);
t65 = t24 / 0.2e1;
t44 = t34 * t38 + t37 * t35;
t25 = t44 * qJD(1);
t64 = -t25 / 0.2e1;
t33 = sin(qJ(4));
t63 = -t33 / 0.2e1;
t62 = -t35 / 0.2e1;
t36 = cos(qJ(4));
t61 = t36 / 0.2e1;
t60 = t38 / 0.2e1;
t59 = (m(4) * pkin(2) ^ 2);
t58 = Ifges(3,4) * t35;
t57 = Ifges(4,4) * t25;
t56 = Ifges(5,4) * t33;
t53 = t25 * mrSges(4,3);
t55 = (t32 * mrSges(4,1) + t53) * t34;
t54 = t24 * mrSges(4,3);
t52 = t33 * t36;
t51 = qJD(2) * pkin(2);
t50 = Ifges(3,5) * qJD(2);
t49 = Ifges(3,6) * qJD(2);
t48 = qJD(1) * t35;
t47 = qJD(1) * t38;
t46 = Ifges(5,2) * t52;
t45 = qJD(3) * t51;
t30 = Ifges(5,4) * t36 * qJD(1);
t43 = (Ifges(5,6) * qJD(4) + (Ifges(5,2) * t36 + t56) * qJD(1)) * t63 + (Ifges(5,1) * t33 * qJD(1) + Ifges(5,5) * qJD(4) + t30) * t61;
t42 = pkin(1) * (mrSges(5,1) * t33 + mrSges(5,2) * t36);
t41 = (Ifges(5,5) * t61 + Ifges(5,6) * t63) * qJD(4);
t14 = t32 * t26;
t11 = t14 * qJD(1);
t15 = t32 * t44;
t12 = t15 * qJD(1);
t19 = Ifges(4,4) * t24;
t7 = Ifges(4,2) * t24 + Ifges(4,6) * t32 - t57;
t8 = -Ifges(4,1) * t25 + Ifges(4,5) * t32 + t19;
t40 = t7 * t64 + Ifges(4,5) * t11 + t37 * mrSges(4,2) * t45 + t25 * (Ifges(4,1) * t24 + t57) / 0.2e1 - t32 * (Ifges(4,5) * t24 + Ifges(4,6) * t25) / 0.2e1 + pkin(2) * (-t25 * mrSges(4,1) + t24 * mrSges(4,2)) * t47 + Ifges(4,6) * t12 - (Ifges(4,2) * t25 + t19 + t8) * t24 / 0.2e1 + (mrSges(4,1) * t45 + t51 * t53) * t34;
t31 = Ifges(3,4) * t47;
t23 = Ifges(3,1) * t48 + t31 + t50;
t21 = t49 + (Ifges(3,2) * t38 + t58) * qJD(1);
t16 = -t32 * mrSges(4,2) + t54;
t13 = -t24 * mrSges(4,1) - t25 * mrSges(4,2);
t1 = [-t38 * pkin(2) * (-t12 * mrSges(4,1) + t11 * mrSges(4,2)) + t14 * t8 / 0.2e1 + t15 * t7 / 0.2e1 + t32 * (Ifges(4,5) * t14 + Ifges(4,6) * t15) / 0.2e1 + (t26 * t12 + t15 * t65) * Ifges(4,2) + (-t11 * t44 + t14 * t64) * Ifges(4,1) + (t41 + t43) * qJD(4) + (t21 * t62 + t23 * t60 + (Ifges(3,5) * t60 + Ifges(3,6) * t62) * qJD(2) + (t35 * t13 + (t14 * t37 - t15 * t34 + (-t26 * t37 + t34 * t44) * qJD(3)) * mrSges(4,3)) * pkin(2)) * qJD(2) + (t26 * t11 - t12 * t44 + t14 * t65 + t15 * t64) * Ifges(4,4) + ((-pkin(2) * (-t15 * mrSges(4,1) + t14 * mrSges(4,2)) + 0.3e1 / 0.2e1 * Ifges(3,4) * qJD(2) * t38) * t38 + (pkin(2) * (-t26 * mrSges(4,1) - mrSges(4,2) * t44) - 0.3e1 / 0.2e1 * t58 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - (2 * t59)) * t38) * qJD(2) * t35 + (0.3e1 / 0.2e1 * Ifges(5,1) * t52 - 0.3e1 / 0.2e1 * t46 - 0.2e1 * t42 + (-0.3e1 / 0.2e1 * t33 ^ 2 + 0.3e1 / 0.2e1 * t36 ^ 2) * Ifges(5,4)) * qJD(4)) * qJD(1); t40 + ((-t23 / 0.2e1 - t31 / 0.2e1 + t50 / 0.2e1) * t38 + (t21 / 0.2e1 - pkin(2) * t13 - t49 / 0.2e1 + Ifges(3,4) * t48 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1 + t59) * t47) * t35) * qJD(1) + ((-t16 * t37 + t55) * qJD(3) + (-t12 * t34 + (-qJD(2) * t24 + t11) * t37) * mrSges(4,3)) * pkin(2); (-t55 + (t16 - t54) * t37) * t51 + t40; (-t36 * t30 / 0.2e1 + (t42 + (Ifges(5,1) * t36 - t56) * t63 + t46 / 0.2e1) * qJD(1) + t41 - t43) * qJD(1); 0;];
tauc = t1(:);
