% Calculate time derivative of joint inertia matrix for
% fourbar1DE2
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
% MqD [1x1]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1DE2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1DE2_inertiaDJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:05:04
% EndTime: 2020-04-24 20:05:05
% DurationCPUTime: 0.32s
% Computational Cost: add. (1117->46), mult. (1620->86), div. (48->10), fcn. (392->4), ass. (0->47)
t64 = pkin(3) ^ 2;
t36 = pkin(1) ^ 2;
t30 = cos(qJ(1));
t53 = pkin(1) * t30;
t47 = -0.2e1 * pkin(2) * t53 + t36;
t59 = -pkin(3) - pkin(4);
t19 = (pkin(2) - t59) * (pkin(2) + t59) + t47;
t58 = -pkin(3) + pkin(4);
t20 = (pkin(2) - t58) * (pkin(2) + t58) + t47;
t49 = t19 * t20;
t37 = sqrt(-t49);
t13 = 0.1e1 / t37;
t63 = t13 / pkin(3);
t61 = pkin(1) * pkin(2);
t35 = pkin(2) ^ 2;
t25 = t35 + t47;
t46 = -pkin(4) ^ 2 + t64;
t21 = t25 + t46;
t29 = sin(qJ(1));
t45 = t29 * t61;
t42 = qJD(1) * t45;
t10 = (-t19 - t20) * t42;
t51 = t10 * t13;
t26 = -pkin(2) + t53;
t52 = pkin(2) * t26;
t4 = t51 * t52;
t43 = t29 ^ 2 * t35 * t36;
t48 = t29 * t37;
t2 = t4 + (-0.2e1 * t43 + (-t21 * t30 - t48) * t61) * qJD(1);
t11 = t37 * t52;
t8 = -t21 * t45 + t11;
t60 = t2 * t8;
t57 = m(3) / t64;
t7 = (pkin(1) * t48 + t21 * t26) * pkin(2);
t56 = mrSges(3,1) * t7;
t55 = mrSges(3,2) * t8;
t22 = t25 - t46;
t9 = t22 * t45 + t11;
t54 = Ifges(3,3) * t9 ^ 2;
t15 = 0.1e1 / t19;
t17 = 0.1e1 / t20;
t50 = t15 * t17;
t44 = 0.4e1 * t50;
t5 = t8 ^ 2;
t3 = t4 + (0.2e1 * t43 + (t22 * t30 - t48) * t61) * qJD(1);
t1 = (0.2e1 * t35 * t26 * qJD(1) * t29 + (t29 * t51 + (-t29 * t21 + t30 * t37) * qJD(1)) * pkin(2)) * pkin(1);
t6 = [((t1 * t7 + t60) * t57 / 0.2e1 + 0.2e1 * (-Ifges(3,3) * t3 * t9 - Ifges(4,3) * t60) * t50 + (-0.2e1 * (-t55 / 0.2e1 + t56 / 0.2e1) * t9 * t10 / t49 + (t2 * t9 + t3 * t8) * mrSges(3,2) + (-t1 * t9 - t3 * t7) * mrSges(3,1)) * t63 + ((-t7 ^ 2 * t57 + t44 * t54 + 0.4e1 * (-t55 + t56) * t9 * t63 + (Ifges(4,3) * t44 - t57) * t5) / t25 + 0.2e1 * (Ifges(4,3) * t5 + t54) * (t15 / t20 ^ 2 + 0.1e1 / t19 ^ 2 * t17)) * t42) / t25 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t6(1);];
Mq = res;
