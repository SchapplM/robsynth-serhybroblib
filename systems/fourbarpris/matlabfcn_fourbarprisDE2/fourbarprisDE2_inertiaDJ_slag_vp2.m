% Calculate time derivative of joint inertia matrix for
% fourbarprisDE2
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
% MqD [1x1]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisDE2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE2_inertiaDJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:03
% EndTime: 2020-05-07 09:45:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (634->39), mult. (398->80), div. (153->12), fcn. (6->2), ass. (0->40)
t61 = 2 * qJ(1);
t30 = qJ(1) + pkin(3);
t46 = -pkin(2) - t30;
t22 = pkin(1) + t46;
t57 = 1 / t22;
t45 = -pkin(2) + t30;
t20 = pkin(1) - t45;
t58 = 1 / t20;
t51 = t57 * t58;
t31 = qJ(1) ^ 2;
t4 = pkin(1) ^ 2 - pkin(2) ^ 2 - t31 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t60 = t4 ^ 2;
t19 = pkin(1) - t46;
t21 = pkin(1) + t45;
t49 = t21 * t22;
t43 = t20 * t49;
t42 = t19 * t43;
t1 = (-t42) ^ (-0.1e1 / 0.2e1);
t59 = -2 * t1;
t27 = 1 / t30;
t10 = 1 / t20 ^ 2;
t56 = 1 / t22 ^ 2;
t26 = t30 ^ 2;
t55 = 4 * t26;
t7 = 1 / t19;
t54 = t58 * t7;
t53 = Ifges(2,3) * t60;
t52 = t27 / t26 * t7;
t48 = -2 * t54;
t47 = Ifges(4,3) * t55;
t44 = Ifges(2,3) * t51;
t8 = 1 / t19 ^ 2;
t41 = -t10 * t7 + t58 * t8;
t40 = m(3) * t31 + mrSges(3,3) * t61 + Ifges(3,2);
t28 = 1 / t30 ^ 2;
t39 = (t28 * t53 + t47) * t7;
t38 = t40 * t51;
t37 = 2 * t38;
t11 = 1 / t21;
t2 = [(2 * (mrSges(3,1) * t27 * t59 + (t37 + 2 * t44) * t7 * t4 * t28 * t11) * t30 + (t28 * t59 + t27 * t1 * (-t43 + (t49 + (t21 - t22) * t20) * t19) / t42) * t4 * mrSges(3,1) + (t47 * t51 + (t38 + t44) * t28 * t60) * t7 / t21 ^ 2 + (-t58 * t56 * t39 - (t10 * t39 - ((-8 * t30 * t7 + t8 * t55) * Ifges(4,3) + (t28 * t8 + 2 * t52) * t53) * t58) * t57 + (t37 * t52 + (-t40 * t56 * t54 + (m(3) * (qJ(1) * t48 + t41 * t31) + t41 * Ifges(3,2) + (t41 * t61 + t48) * mrSges(3,3)) * t57) * t28) * t60) * t11) * qJD(1);];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t2(1);];
Mq = res;
