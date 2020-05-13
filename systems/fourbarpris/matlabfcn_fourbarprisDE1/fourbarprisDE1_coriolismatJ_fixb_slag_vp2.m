% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbarprisDE1
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
% Cq [1x1]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbarprisDE1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE1_coriolismatJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:52
% EndTime: 2020-05-07 09:09:55
% DurationCPUTime: 0.29s
% Computational Cost: add. (634->38), mult. (355->72), div. (153->12), fcn. (6->2), ass. (0->34)
t29 = (qJ(1) + pkin(3));
t43 = (-pkin(2) - t29);
t21 = (pkin(1) + t43);
t15 = 1 / (t21 ^ 2);
t42 = (-pkin(2) + t29);
t19 = (pkin(1) - t42);
t53 = 1 / t19;
t56 = t15 * t53;
t25 = (t29 ^ 2);
t26 = 1 / t29;
t30 = (qJ(1) ^ 2);
t4 = pkin(1) ^ 2 - pkin(2) ^ 2 - t30 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t54 = t4 ^ 2;
t55 = t26 / t25 * t54;
t52 = 1 / t21;
t48 = t52 * t53;
t17 = 1 / (t19 ^ 2);
t50 = m(3) * t30;
t49 = Ifges(2,3) * t54;
t47 = t4 * t26;
t46 = Ifges(4,3) * t25;
t20 = pkin(1) + t42;
t45 = t20 * t21;
t44 = qJ(1) * mrSges(3,3);
t41 = t53 * t46;
t40 = t19 * t45;
t18 = pkin(1) - t43;
t39 = t18 * t40;
t38 = t44 + t50 / 0.2e1 + Ifges(3,2) / 0.2e1;
t37 = (Ifges(3,2) + (2 * t44) + t50) * t48;
t27 = 0.1e1 / t29 ^ 2;
t36 = (2 * t52 * t41) + (t49 / 0.2e1 + t38 * t54) * t48 * t27;
t6 = 1 / t18;
t1 = [((((-t40 + (t45 + (t20 - t21) * t19) * t18) / t39 * t47) / 0.2e1 - (2 * t26 * t29) - t4 * t27) * ((-t39) ^ (-0.1e1 / 0.2e1)) * mrSges(3,1) + t36 * t6 / (t20 ^ 2) + (t36 / (t18 ^ 2) + (-(2 * t15 * t41) + t37 * t55 - ((2 * t17 * t46 - (-4 * Ifges(4,3) * t29 + (2 * t47 + t55) * Ifges(2,3)) * t53) * t52) + ((-t56 / 0.2e1 - (t52 * t17) / 0.2e1) * t49 + 0.2e1 * t4 * t29 * t37 + (-t38 * t56 + (((-m(3) * qJ(1) - mrSges(3,3)) * t53) - t38 * t17) * t52) * t54) * t27) * t6) / t20) * qJD(1);];
Cq = t1;
