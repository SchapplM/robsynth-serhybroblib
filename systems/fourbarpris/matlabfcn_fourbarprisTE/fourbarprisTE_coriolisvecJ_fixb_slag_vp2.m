% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
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
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisTE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisTE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:21
% EndTime: 2020-05-07 09:01:24
% DurationCPUTime: 0.41s
% Computational Cost: add. (1758->53), mult. (1002->112), div. (155->12), fcn. (90->2), ass. (0->54)
t39 = (qJ(1) + pkin(3));
t57 = -pkin(2) - t39;
t28 = pkin(1) - t57;
t56 = -pkin(2) + t39;
t29 = pkin(1) - t56;
t30 = pkin(1) + t56;
t31 = pkin(1) + t57;
t63 = t30 * t31;
t54 = t29 * t63;
t51 = (t28 * t54);
t10 = ((-t51) ^ (-0.1e1 / 0.2e1));
t41 = (qJ(1) ^ 2);
t13 = pkin(1) ^ 2 - pkin(2) ^ 2 - t41 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t68 = t10 * t13;
t35 = (t39 ^ 2);
t36 = 1 / t39;
t76 = t13 ^ 2;
t78 = t36 / t35 * t76;
t74 = 0.1e1 / t31;
t75 = 1 / t29;
t77 = (t75 * t74);
t27 = 0.1e1 / t29 ^ 2;
t25 = 0.1e1 / t31 ^ 2;
t15 = (t39 * qJD(1));
t73 = -2 * t15;
t72 = m(3) * t41;
t71 = m(3) * qJ(1);
t70 = Ifges(2,3) * t76;
t69 = Ifges(4,3) * t35;
t65 = (t25 * t75);
t37 = 1 / (t39 ^ 2);
t52 = (1 / t51 * t68) / 0.2e1;
t47 = -t54 + (t63 + (t30 - t31) * t29) * t28;
t7 = t47 * qJD(1);
t50 = t7 * t52;
t40 = (qJD(1) ^ 2);
t55 = t40 * t68;
t60 = qJD(1) * t36;
t3 = -(t37 * t55) + ((t10 * t73) + t50) * t60;
t64 = t3 * qJ(1);
t59 = t13 * qJD(1);
t53 = t37 * t59;
t4 = -(t10 * t53) + (-(2 * t39 * t10) + t47 * t52) * t60;
t62 = t4 * qJ(1);
t61 = mrSges(3,3) * qJ(1);
t58 = 2 * t69;
t49 = (Ifges(3,2) + 2 * t61 + t72) * t77;
t48 = ((t37 * t70) / 0.2e1 + t58) * t40 * t77;
t21 = 0.1e1 / t30 ^ 2;
t20 = 0.1e1 / t30;
t17 = 0.1e1 / t28 ^ 2;
t2 = -(t36 * t10 * t59) - t62;
t1 = -(t36 * t55) - t64;
t5 = [t3 * mrSges(3,1) + t20 * t17 * t48 + (-mrSges(3,1) * t37 * t40 + (-t1 * t71 + Ifges(3,2) * t3 + (-t1 + t64) * mrSges(3,3)) * t36) * t68 + (-t4 * mrSges(3,1) + (-t72 / 0.2e1 - Ifges(3,2) / 0.2e1 - t61) * t7 * t37 * t25 * t21 * t27 * t17 * t76 + (mrSges(3,1) * t50 + (mrSges(3,1) * t73 + (t2 * t71 - Ifges(3,2) * t4 + (t2 - t62) * mrSges(3,3)) * t13) * t10) * t36) * qJD(1) + (t21 * t48 + ((2 * (2 * Ifges(2,3) * t77 + t49) * t15 * t53) + (-(2 * t65 * t69) + (t49 * t78) - (t27 * t58 - ((-4 * Ifges(4,3) * t39 + (-2 * t13 * t36 + t78) * Ifges(2,3)) * t75)) * t74 + (((-mrSges(3,3) - t71) * t76 * t77) + (-t65 / 0.2e1 - t74 * t27 / 0.2e1) * t70) * t37) * t40) * t20) / t28;];
tauc = t5(:);
