% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
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
% Cq [1x1]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1TE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:17
% EndTime: 2020-04-24 19:49:19
% DurationCPUTime: 0.29s
% Computational Cost: add. (1117->45), mult. (1565->81), div. (48->10), fcn. (392->4), ass. (0->45)
t64 = pkin(3) ^ 2;
t36 = pkin(1) ^ 2;
t30 = cos(qJ(1));
t52 = pkin(1) * t30;
t46 = -0.2e1 * pkin(2) * t52 + t36;
t58 = -pkin(3) - pkin(4);
t19 = (pkin(2) - t58) * (pkin(2) + t58) + t46;
t57 = -pkin(3) + pkin(4);
t20 = (pkin(2) - t57) * (pkin(2) + t57) + t46;
t48 = t19 * t20;
t37 = sqrt(-t48);
t13 = 0.1e1 / t37;
t63 = t13 / pkin(3);
t26 = -pkin(2) + t52;
t51 = pkin(2) * t26;
t11 = t37 * t51;
t35 = pkin(2) ^ 2;
t25 = t35 + t46;
t45 = -pkin(4) ^ 2 + t64;
t21 = t25 + t45;
t29 = sin(qJ(1));
t60 = pkin(1) * pkin(2);
t44 = t29 * t60;
t8 = -t21 * t44 + t11;
t5 = t8 ^ 2;
t22 = t25 - t45;
t9 = t22 * t44 + t11;
t62 = Ifges(3,3) * t9 ^ 2 + Ifges(4,3) * t5;
t10 = (-t19 - t20) * t44;
t50 = t13 * t10;
t4 = t50 * t51;
t43 = t29 ^ 2 * t35 * t36;
t47 = t29 * t37;
t2 = -0.2e1 * t43 + t4 + (-t21 * t30 - t47) * t60;
t59 = t2 * t8;
t7 = (pkin(1) * t47 + t21 * t26) * pkin(2);
t56 = mrSges(3,1) * t7;
t55 = mrSges(3,2) * t8;
t15 = 0.1e1 / t19;
t17 = 0.1e1 / t20;
t49 = t15 * t17;
t42 = m(3) / t64 / 0.4e1;
t3 = 0.2e1 * t43 + t4 + (t22 * t30 - t47) * t60;
t1 = (0.2e1 * t35 * t26 * t29 + (t30 * t37 + (-t21 + t50) * t29) * pkin(2)) * pkin(1);
t6 = [((t1 * t7 + t59) * t42 + (-Ifges(3,3) * t3 * t9 - Ifges(4,3) * t59) * t49 + (-0.2e1 * (-t55 / 0.4e1 + t56 / 0.4e1) * t9 * t10 / t48 + (t2 * t9 / 0.2e1 + t8 * t3 / 0.2e1) * mrSges(3,2) + (-t1 * t9 / 0.2e1 - t7 * t3 / 0.2e1) * mrSges(3,1)) * t63 + (t62 * (t15 / t20 ^ 2 + 0.1e1 / t19 ^ 2 * t17) + 0.2e1 * ((-t7 ^ 2 - t5) * t42 + t62 * t49 + (-t55 + t56) * t9 * t63) / t25) * t44) / t25 ^ 2 * qJD(1);];
Cq = t6;
